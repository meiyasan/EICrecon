// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// TimesliceBuffer_service: process-wide, frame-indexed store of
// slow-detector time-aligned hits, used for cross-frame boundary recovery.
//
// JANA2 creates one JOmniFactory instance per in-flight timeslice
// (jana:max_inflight_timeslices, default nthreads) and assigns frames to
// instances arbitrarily, so a per-factory rolling buffer would see a
// nondeterministic subset of frames. This store is keyed by absolute frame
// (timeslice) number instead, so it stays independent of both thread count
// and processing order.
//
// Producer: TrkTimeAlignment_factory calls deposit() from the parallel map
//   stage. Deposits are order-free and never wait on other frames.
// Consumer: EventUnfolder calls fetch() from the sequential unfold stage.
//   fetch() waits, bounded by a timeout, for frames still in flight, and
//   returns std::nullopt for frames that will never arrive (before stream
//   start, after stream end). A frame is always deposited before its own
//   unfold runs, so a backward fetch waits for at most one in-flight
//   neighbor.

#pragma once

#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <map>
#include <mutex>
#include <iostream>
#include <optional>
#include <set>
#include <vector>

#include <JANA/Services/JServiceLocator.h>

#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>

// Effective provenance of a particle, as a (generatorStatus, time) pair.
// Geant4 secondaries carry generatorStatus 0, which erases the class band
// the link weights encode -- a delta electron behind a TOF hit would be
// labelled "no class" even though its primary is e.g. a ncdisq1 product.
// Climb the parent chain to the first ancestor with a nonzero status and
// label with that ancestor's status and production time instead: the
// ancestor's time is the collision time, so the pile-up slot also lands on
// the right collision. Every status-0 secondary in production frames traces
// to exactly one banded ancestor (verified on real data, 58,731/58,731);
// the depth cap is defensive. A chain with no labelled ancestor keeps the
// particle's own (0, time).
inline std::pair<int32_t, double> effectiveProvenance(const edm4hep::MCParticle& p) {
  edm4hep::MCParticle cur = p;
  for (int depth = 0; depth < 64; ++depth) {
    if (cur.getGeneratorStatus() != 0)
      return {cur.getGeneratorStatus(), double(cur.getTime())};
    const auto parents = cur.getParents();
    if (parents.empty() || !parents[0].isAvailable())
      break;
    cur = parents[0];
  }
  return {p.getGeneratorStatus(), double(p.getTime())};
}

class TimesliceBuffer_service : public JService {
public:
  // One buffered slow-detector hit, carrying everything the EventUnfolder
  // needs to rebuild that hit's truth in a child event.
  //
  // The RecHit alone is not enough. Its rawHit relation, and the
  // RawHitLink -> SimTrackerHit -> MCParticle chain behind it, point into
  // the *producing* frame's collections, which are destroyed as soon as
  // that frame is done. The unfolder therefore has to clone adjacent-frame
  // hits with cloneRelations=false, which used to leave their rawHit
  // dangling (index -1) and emit no link at all -- so every adjacent-frame
  // Si/B0 hit reached the child event with no decodable class label. That
  // was 20-30% of silicon hits in a low-background campaign and ~50-67%
  // in a high-background one (the ratio is just n_adjacent/(n_adjacent+1),
  // so it also moved with forward_frames being active or not).
  //
  // The fix is to resolve the truth at deposit time, while the producing
  // frame is still alive, and carry it here: an untracked clone of the raw
  // hit (so the child can own a real RawTrackerHit to point at and to hang
  // a link from), plus the MCParticle's generatorStatus and time. The time
  // is kept rather than a finished label because the pile-up slot is a
  // per-candidate quantity -- only the unfolder knows the candidate's
  // trigger_classes list, so it adds slot*1e6 itself.
  struct SlowHit {
    edm4eic::MutableTrackerHit    hit;      // time-corrected, untracked clone
    edm4eic::MutableRawTrackerHit raw;      // untracked clone; valid iff has_raw
    float status  = 0.0f;                   // effectiveProvenance() status: own, or
                                            // first labelled ancestor's for secondaries
    float ptime   = 0.0f;                   // matching time [ns], for the slot
    bool  has_raw = false;
    bool  has_truth = false;                // status/ptime resolved

    // Full truth chain for this hit, buffered only when
    // eventbuilder:truth:simhits is on (see m_buffer_truth in
    // TimeAlignment_factory). The weight-encoded label above is enough for
    // segmentation truth; these two are what make
    // link -> SimTrackerHit -> MCParticle resolve offline for a hit that came
    // from a neighbouring frame. Without them, cross-frame recovery
    // (backward_frames / forward_frames) produces hits that no amount of
    // simhits=1 can chain, because both the sim hit and the particle belong
    // to a frame that has already been destroyed.
    //
    // The clones are untracked and shared: one sim hit backs several raw
    // hits, and one particle backs many sim hits, so the producing frame
    // deduplicates before it deposits. The keys are the ORIGINAL
    // (collectionID, index) / index from that frame, so the consumer can
    // deduplicate again per child event.
    edm4hep::MutableSimTrackerHit sim;      // untracked clone; valid iff has_sim
    edm4hep::MutableMCParticle    particle; // untracked clone; valid iff has_particle
    uint64_t sim_key      = 0;
    uint32_t mcp_key      = 0;
    bool     has_sim      = false;
    bool     has_particle = false;

    // Frame that produced this hit. It MUST take part in any deduplication
    // of sim_key / mcp_key: podio indices restart in every frame, so a
    // backward and a forward neighbour both contain index 5, and keying on
    // the index alone makes the second hit reuse the first one's particle.
    // That is silent -- the chain still resolves, just to the wrong particle.
    uint64_t frame = 0;
  };

  // Corrected slow-detector hits of one frame: [slow_det_index][hit].
  // Hits are untracked podio clones (owned by no collection, refcounted).
  using FrameHits = std::vector<std::vector<SlowHit>>;

  // Frames retained behind the newest deposit. Bounds memory; must exceed the
  // largest usable backward/forward window plus the number of in-flight frames.
  static constexpr uint64_t kRetentionFrames = 32;

  void deposit(uint64_t frame_nr, FrameHits hits) {
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      // Frame numbering restarted (e.g. a new input stream): drop stale state.
      if (!m_frames.empty() && frame_nr + kRetentionFrames < m_frames.rbegin()->first) {
        m_frames.clear();
        m_floor       = 0;
        m_contig_next = 0;
        m_deposited.clear();
        // A restart is also a fact that resolves pending fetches: a waiter
        // asking for the OLD stream's next frame can now be told it is never
        // coming, instead of discovering that by timeout at every file
        // boundary of a multi-file job.
        ++m_generation;
      }
      m_frames[frame_nr] = std::move(hits);

      // Advance the contiguous watermark: every frame BELOW m_contig_next has
      // been deposited at some point. Deposits arrive out of order once more
      // than one timeslice is in flight, so the watermark tracks a set rather
      // than the newest number -- taking the maximum would declare a frame
      // absent while it was still on its way, which is the nondeterminism this
      // whole mechanism exists to remove.
      m_deposited.insert(frame_nr);
      while (m_deposited.erase(m_contig_next) != 0)
        ++m_contig_next;
      const uint64_t newest = m_frames.rbegin()->first;
      if (newest > kRetentionFrames) {
        const uint64_t cutoff = newest - kRetentionFrames;
        m_frames.erase(m_frames.begin(), m_frames.lower_bound(cutoff));
        if (m_floor < cutoff)
          m_floor = cutoff; // pruned frames can no longer be fetched
      }
    }
    m_cv.notify_all();
  }

  // Fetch one frame's deposit. The answer must not depend on WHEN it is asked:
  // a candidate's cross-frame recovery is physics, so the same input has to
  // give the same child event at any thread count and under any load.
  //
  // The decision is therefore made from facts about the stream, never from a
  // clock. A frame is resolved when one of these holds:
  //
  //   * it is in the buffer                      -> return it
  //   * frame_nr < m_floor                       -> evicted past the retention
  //                                                 window; gone for good
  //   * frame_nr < m_contig_next                 -> already deposited and no
  //                                                 longer here, so it is not
  //                                                 coming back
  //
  // and otherwise it has simply not been reached yet, so we wait. Crucially a
  // timeout no longer writes m_floor. It used to: one fetch that waited too
  // long marked every lower frame permanently absent, GLOBALLY, so a single
  // scheduling hiccup silently changed the recovered hits of every candidate
  // afterwards. That conflated "I stopped waiting" with "this frame is gone",
  // and it was measurable -- 109 of 391 child events differed between a
  // 1-thread and a 4-thread run of identical input.
  //
  // The timeout survives only to break the genuine tail case: the last frame
  // of a stream asks for a successor that will never be deposited, and nothing
  // in the buffer can prove a stream has ended. That costs one wait per
  // forward_frames at end of job and resolves no other case -- in steady state
  // a forward neighbour arrives (that is what jana:max_inflight_timeslices >
  // forward_frames guarantees) and a backward one is already below
  // m_contig_next.
  /// Why a fetch came back empty. The distinction is what keeps the result
  /// deterministic: Absent is a FACT about the stream (the frame was evicted,
  /// or the numbering restarted for a new input file) and is safe to act on;
  /// Timeout means the answer was still unknown when we gave up, and acting
  /// on it silently would make the recovered physics depend on machine load
  /// and thread count -- the caller must treat it as an error, except at the
  /// provable stream tail.
  enum class Why { Present, Absent, Timeout };

  // The wait itself is cheap insurance, not the mechanism: at >= 3 worker
  // threads the neighbouring frame's deposit was measured to ALWAYS lead the
  // fetch in steady state (the only misses at the old 2 s timeout were the
  // handful of fetches issued while one-time geometry/ONNX initialisation
  // still stalled the pipeline, plus multi-file boundaries -- both of which
  // now resolve on facts: the first deposit and the generation bump). The
  // timeout only bounds how long a pathologically stalled pipeline is given
  // before the caller turns the unknown into a hard error.
  std::optional<FrameHits> fetch(uint64_t frame_nr, Why& why,
                                 std::chrono::milliseconds timeout = std::chrono::seconds(30)) {
    std::unique_lock<std::mutex> lock(m_mutex);
    const uint64_t gen0 = m_generation;
    m_cv.wait_for(lock, timeout, [&] {
      return m_frames.count(frame_nr) != 0 || frame_nr < m_floor || frame_nr < m_contig_next ||
             m_generation != gen0;
    });
    auto it = m_frames.find(frame_nr);
    if (it != m_frames.end()) {
      why = Why::Present;
      return it->second; // shallow podio handle copies
    }
    if (frame_nr < m_floor || frame_nr < m_contig_next || m_generation != gen0) {
      why = Why::Absent; // evicted, already-passed, or a new input file began
      return std::nullopt;
    }
    why = Why::Timeout;
    return std::nullopt;
  }

  /// True while no frame numbered above frame_nr has ever been deposited --
  /// the signature of the stream tail, where a timed-out forward fetch is
  /// expected (the successor genuinely does not exist) rather than an error.
  bool nothingBeyond(uint64_t frame_nr) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_deposited.empty() && *m_deposited.rbegin() > frame_nr)
      return false;
    return m_contig_next <= frame_nr + 1;
  }

private:
  void acquire_services(JServiceLocator*) override {}

  mutable std::mutex            m_mutex;
  std::condition_variable       m_cv;
  std::map<uint64_t, FrameHits> m_frames; // frame number -> slow-det hits
  uint64_t                      m_floor = 0; // frames below this never arrive
  /// Every frame below this has been deposited. Together with m_floor it makes
  /// "this frame is not coming" a fact about the stream instead of a timeout.
  uint64_t                      m_contig_next = 0;
  /// Deposited frames at or above m_contig_next, i.e. the out-of-order window.
  /// Erased as the watermark passes them, so this stays as small as the number
  /// of timeslices in flight.
  std::set<uint64_t>            m_deposited;
  /// Bumped when the frame numbering restarts (a new input file). Lets a
  /// pending fetch for the old stream fail on the restart FACT, not a clock.
  uint64_t                      m_generation = 0;
};

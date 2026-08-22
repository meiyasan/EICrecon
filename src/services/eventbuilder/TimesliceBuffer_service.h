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
#include <optional>
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
        m_floor = 0;
      }
      m_frames[frame_nr] = std::move(hits);
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

  // Fetch one frame's deposit. Blocks until the frame arrives, it is known to
  // never arrive, or the timeout expires. A timeout marks the frame as
  // permanently absent so later fetches below it return immediately.
  std::optional<FrameHits> fetch(uint64_t frame_nr,
                                 std::chrono::milliseconds timeout = std::chrono::milliseconds(2000)) {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_cv.wait_for(lock, timeout, [&] {
      return m_frames.count(frame_nr) != 0 || frame_nr < m_floor;
    });
    auto it = m_frames.find(frame_nr);
    if (it != m_frames.end())
      return it->second; // shallow podio handle copies
    if (m_floor <= frame_nr)
      m_floor = frame_nr + 1; // timed out: treat as absent from now on
    return std::nullopt;
  }

private:
  void acquire_services(JServiceLocator*) override {}

  std::mutex                    m_mutex;
  std::condition_variable       m_cv;
  std::map<uint64_t, FrameHits> m_frames; // frame number -> slow-det hits
  uint64_t                      m_floor = 0; // frames below this never arrive
};

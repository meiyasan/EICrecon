// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// TimesliceBuffer_service — process-wide, frame-indexed store of slow-detector
// time-aligned hits, used for cross-frame boundary recovery.
//
// JANA2 creates one JOmniFactory instance per in-flight timeslice 
// (jana:max_inflight_timeslices, default nthreads) and assigns frames to
// instances arbitrarily, so a per-factory rolling buffer sees a nondeterministic
//  subset of frames. This store is keyed by absolute frame (timeslice) number instead, 
// which makes it independent of both thread count and processing order.
//
// Producers: TrkTimeAlignment_factory calls deposit() from the parallel map
//   stage. Deposits are order-free and never wait on other frames.
// Consumer: EventUnfolder calls fetch() from the sequential unfold stage.
//   fetch() waits (bounded by timeout) for frames still in flight and returns
//   std::nullopt for frames that will never arrive (before stream start,
//   after stream end). A frame is guaranteed to be deposited before its own
//   unfold runs, so backward fetches at most wait for an in-flight neighbour.

#pragma once

#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <map>
#include <mutex>
#include <optional>
#include <vector>

#include <JANA/Services/JServiceLocator.h>

#include <edm4eic/TrackerHitCollection.h>

class TimesliceBuffer_service : public JService {
public:
  // Corrected slow-detector hits of one frame: [slow_det_index][hit].
  // Hits are untracked podio clones (owned by no collection, refcounted).
  using FrameHits = std::vector<std::vector<edm4eic::MutableTrackerHit>>;

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

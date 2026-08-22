// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Jefferson Science Associates, LLC.

#include "ONNXRuntime_service.h"

std::shared_ptr<Ort::Session> ONNXRuntime_service::session(const std::string& path) {
  std::lock_guard<std::mutex> lock(m_mutex);

  auto it = m_sessions.find(path);
  if (it != m_sessions.end())
    return it->second;

  // JANA owns cross-event parallelism; keep ORT single-threaded per session.
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  opts.SetInterOpNumThreads(1);

  auto sess = std::make_shared<Ort::Session>(m_env, path.c_str(), opts);
  m_sessions.emplace(path, sess);
  return sess;
}

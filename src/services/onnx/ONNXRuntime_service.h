// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Jefferson Science Associates, LLC.

#pragma once

#include <JANA/JApplicationFwd.h>
#include <JANA/Services/JServiceLocator.h>
#include <onnxruntime_cxx_api.h>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

// Singleton JService that owns one Ort::Env per process and caches
// Ort::Session handles keyed by model path.
//
// ONNX Runtime's own recommendation is a single Ort::Env for the process;
// multiple Env instances each spin up their own thread pool and logging
// infrastructure. Sessions are thread-safe for concurrent Run() calls, so
// all JANA worker threads may share a single session for the same model.
//
// Usage in a factory Configure():
//   auto svc = GetApplication()->GetService<ONNXRuntime_service>();
//   m_session = svc->session("/path/to/model.onnx");
//   // m_session is a shared_ptr<Ort::Session>; keep it alive for Process().

class ONNXRuntime_service : public JService {
public:
  ONNXRuntime_service(JApplication* app) : m_app(app) {}
  ~ONNXRuntime_service() override = default;

  // Returns (or creates) a shared session for the given model path.
  // Thread-safe: protected by a mutex during session creation only.
  std::shared_ptr<Ort::Session> session(const std::string& model_path);

private:
  ONNXRuntime_service() = default;
  void acquire_services(JServiceLocator*) override {}

  [[maybe_unused]] JApplication*                             m_app{nullptr};
  Ort::Env                                                   m_env{ORT_LOGGING_LEVEL_WARNING, "eicrecon_onnx"};
  std::mutex                                                 m_mutex;
  std::unordered_map<std::string, std::shared_ptr<Ort::Session>> m_sessions;
};

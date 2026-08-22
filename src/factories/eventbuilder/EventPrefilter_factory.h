// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// EventPrefilter_factory
// ----------------------
// Optional ONNX-based, frame-level background rejection between
// EventBuilder_factory and EventUnfolder. Runs once per Timeslice.
//
// GNN (4 inputs) MultiClassEventGNN → probs [1,7]:
//    Input[0]  x          [N, 6]  float32 normalised hit features
//                                  [x/4000, y/4000, z/5000, t/50,
//                                   clip(eDep/10,0,1), det_id/6]
//    Input[1]  edge_idx_0 [2, E]  int64  k-NN edges layer 0 (4D spacetime)
//    Input[2]  edge_idx_1 [2, E]  int64  k-NN edges layer 1 (spatial approx)
//    Input[3]  edge_idx_2 [2, E]  int64  k-NN edges layer 2 (spatial approx)
//    Output[0] probs      [1, 7]  float32  softmax class probabilities
//    Frame passes if probs[BKG=0] < score_threshold.
//
// Known model:
//   vendor/eicrecon/share/models/prefilter-gold-coating-100epochs-2500hits.onnx
//
// Activation:
//   -Peventbuilder:onnx_model=/path/to/model.onnx
//   -Peventbuilder:score_threshold=0.5   (default)
//   -Peventbuilder:max_hits=2500          (0 = dynamic)
//   -Peventbuilder:k_neighbors=16
//   -Peventbuilder:time_weight=1.0
//
// Default (empty model path): zero overhead — pass-through.

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include <onnxruntime_cxx_api.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/TrackerHitCollection.h>

#include <extensions/jana/JOmniFactory.h>
#include "services/onnx/ONNXRuntime_service.h"
#include "services/onnx/ONNXRuntime_gnn.h"

struct EventPrefilter_factory : public JOmniFactory<EventPrefilter_factory> {

  std::string m_onnx_model_path;
  float       m_score_threshold_val{0.5f};
  int32_t     m_max_hits_val{0};
  int32_t     m_k_neighbors_val{16};
  float       m_time_weight_val{1.0f};

  PodioInput<edm4hep::EventHeader>              m_candidates_in{this, "EventCandidates"};
  VariadicPodioInput<edm4eic::TrackerHit, true> m_fast_hits_in{this, {}};
  PodioOutput<edm4hep::EventHeader>             m_candidates_out{this, "EventCandidatesFiltered"};

  std::shared_ptr<Ort::Session> m_session;
  int64_t m_n_features{4};

  // Input/output name storage (owned strings, pointers valid while m_session alive).
  std::vector<std::string>    m_in_name_strs;
  std::vector<std::string>    m_out_name_strs;
  std::vector<const char*>    m_in_names;
  std::vector<const char*>    m_out_names;

  ONNXRuntime_gnn m_gnn;

  void Configure() {
    auto* pm = GetApplication()->GetJParameterManager();
    pm->SetDefaultParameter("eventbuilder:onnx_model",     m_onnx_model_path,     "Path to frame-level ONNX prefilter model. Empty = pass-through.");
    pm->SetDefaultParameter("eventbuilder:score_threshold", m_score_threshold_val, "Maximum BKG probability to pass (GNN class-0 score).");
    pm->SetDefaultParameter("eventbuilder:max_hits",        m_max_hits_val,        "Fixed input length for padded models (e.g. 2500). 0 = use actual hit count.");
    pm->SetDefaultParameter("eventbuilder:k_neighbors",     m_k_neighbors_val,     "k for kNN graph construction (per DGCNN layer).");
    pm->SetDefaultParameter("eventbuilder:time_weight",     m_time_weight_val,     "Scale factor on the time dimension in layer-0 spacetime kNN.");

    if (m_onnx_model_path.empty()) return;

    m_session = GetApplication()->GetService<ONNXRuntime_service>()->session(m_onnx_model_path);

    Ort::AllocatorWithDefaultOptions alloc;
    const size_t n_inputs = m_session->GetInputCount();
    for (size_t i = 0; i < n_inputs; ++i)
      m_in_name_strs.push_back(m_session->GetInputNameAllocated(i, alloc).get());
    const size_t n_outputs = m_session->GetOutputCount();
    for (size_t i = 0; i < n_outputs; ++i)
      m_out_name_strs.push_back(m_session->GetOutputNameAllocated(i, alloc).get());

    for (const auto& s : m_in_name_strs)  m_in_names.push_back(s.c_str());
    for (const auto& s : m_out_name_strs) m_out_names.push_back(s.c_str());

    const auto in_shape = m_session->GetInputTypeInfo(0)
                              .GetTensorTypeAndShapeInfo().GetShape();
    if (in_shape.size() >= 2 && in_shape.back() > 0)
      m_n_features = in_shape.back();

    fprintf(stderr,
            "[EventPrefilter] model: %s  features/hit: %ld\n",
            m_onnx_model_path.c_str(), static_cast<long>(m_n_features));
  }

  void ChangeRun(int32_t) {}

  void Process(int64_t, uint64_t) {
    const auto* cands = m_candidates_in();
    if (!cands) return;

    m_candidates_out()->setSubsetCollection(true);

    if (!m_session) {
      for (const auto& h : *cands) m_candidates_out()->push_back(h);
      return;
    }

    const int64_t max = static_cast<int64_t>(m_max_hits_val);

    // det_id: TOF (0,1) → 2/6, MPGD (2..5) → 1/6
    std::vector<float> feat;
    size_t n_hits = 0;
    const auto& colls = m_fast_hits_in();
    for (size_t ci = 0; ci < colls.size(); ++ci) {
      const auto* coll = colls[ci];
      if (!coll) continue;
      const float det_id_norm = (ci <= 1) ? 2.f / 6.f : 1.f / 6.f;
      for (const auto& h : *coll) {
        if (max > 0 && static_cast<int64_t>(n_hits) >= max) goto done_collecting;
        const auto& pos = h.getPosition();
        const float e_dep = std::min(1.f, std::max(0.f, h.getEdep() / 10.f));
        feat.push_back(static_cast<float>(pos[0]) / 4000.f);
        feat.push_back(static_cast<float>(pos[1]) / 4000.f);
        feat.push_back(static_cast<float>(pos[2]) / 5000.f);
        feat.push_back(h.getTime()                / 50.f);
        feat.push_back(e_dep);
        feat.push_back(det_id_norm);
        ++n_hits;
      }
    }
    done_collecting:;

    const size_t target = (max > 0) ? static_cast<size_t>(max) : n_hits;
    feat.resize(target * static_cast<size_t>(m_n_features), 0.0f);
    if (target == 0) return; // no fast hits → treat as background

    Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    run_gnn(cands, feat, static_cast<int>(target), mem);
  }

private:
  void run_gnn(const edm4hep::EventHeaderCollection* cands,
               std::vector<float>& feat, int N,
               Ort::MemoryInfo& mem)
  {
    // Build k-NN edge indices for all 3 DGCNN layers.
    auto edges = m_gnn.build_dgcnn_edges(
        feat.data(), N, m_k_neighbors_val, m_time_weight_val,
        static_cast<int>(m_n_features));
    const int64_t E = static_cast<int64_t>(edges[0].row.size());

    // Flatten each EdgeIndex into a [2, E] int64 buffer (row-first).
    // Build in-place: reuse edges[i].row as the combined [row||col] buffer.
    std::array<std::vector<int64_t>, 3> ei_buf;
    for (int l = 0; l < 3; ++l) {
      ei_buf[l].reserve(2 * E);
      ei_buf[l].insert(ei_buf[l].end(), edges[l].row.begin(), edges[l].row.end());
      ei_buf[l].insert(ei_buf[l].end(), edges[l].col.begin(), edges[l].col.end());
    }

    const std::vector<int64_t> x_sh{N, m_n_features};
    const std::vector<int64_t> ei_sh{2, E};

    std::vector<Ort::Value> ins;
    ins.push_back(Ort::Value::CreateTensor<float>(
        mem, feat.data(), feat.size(), x_sh.data(), x_sh.size()));
    for (int l = 0; l < 3; ++l)
      ins.push_back(Ort::Value::CreateTensor<int64_t>(
          mem, ei_buf[l].data(), ei_buf[l].size(), ei_sh.data(), ei_sh.size()));

    auto out = m_session->Run(
        Ort::RunOptions{nullptr}, m_in_names.data(), ins.data(), ins.size(),
        m_out_names.data(), 1);

    // Output[0] = probs [1, 7]; class 0 = BKG.
    const float bkg_prob = *out[0].GetTensorMutableData<float>();
    if (bkg_prob < m_score_threshold_val)
      for (const auto& h : *cands) m_candidates_out()->push_back(h);
  }
};

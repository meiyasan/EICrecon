// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// EventPrefilter_factory
// ----------------------
// Optional ONNX-based, frame-level background rejection between
// EventBuilder_factory and EventUnfolder. Runs once per Timeslice.
//
// GNN, up to 6 inputs, MultiClassEventGNN, producing probs [1,N]:
//    Input[0]  x          [N, 6]  float32. Normalized hit/cluster features:
//                                  [x/4000, y/4000, z/5000, t/50,
//                                   clip(eDep-or-energy/10,0,1), det_id/6].
//                                  N spans all wired detectors together
//                                  (fast and slow tracking hits, plus calo
//                                  clusters) as one unified graph,
//                                  distinguished only by det_id; see
//                                  m_trk_hit_names/m_calo_cluster_names
//                                  below for the category assignment.
//    Input[1]  edge_idx_0 [2, E]  int64. k-NN edges, layer 0 (4D spacetime).
//    Input[2]  edge_idx_1 [2, E]  int64. k-NN edges, layer 1 (spatial approx).
//    Input[3]  edge_idx_2 [2, E]  int64. k-NN edges, layer 2 (spatial approx).
//    Input[4]  edge_idx_3 [2, E]  int64. k-NN edges, layer 3 (spatial approx).
//    Input[5]  cand_info  [C, 9] float32. Optional per-candidate event-level
//                                  summary (see build_cand_info() below).
//                                  Sent only if the loaded model declares a
//                                  6th input; older 5-input models are
//                                  unaffected.
//    Output[0] probs       [1, N]      float32. Softmax class probabilities.
//                                       Frame passes if probs[BKG=0] < score_threshold.
//    Output[?] lambda_k    [1, M]      float32. Optional, stored raw (see
//                                       eventbuilder:prefilter:store_aux_outputs
//                                       below). Not yet trained/meaningful.
//    Output[?] seg_logits  [n_hits, K] float32. Optional, per-hit, stored raw
//                                       (same parameter). Not yet trained/
//                                       meaningful. Rows follow the same
//                                       flattened hit order as Input[0] x,
//                                       truncated to the first n_hits real
//                                       rows (padding rows added by max_hits
//                                       are dropped).
//    probs is always read from Output[0]; lambda_k/seg_logits are looked up
//    by NAME (not position), so their exact slot in the model's output list
//    does not matter and their absence is not an error.
//
// No model is currently trained against the 6-input (hits+calo+cand_info)
// contract described here. This factory only wires the richer data
// through, so a future training pass has it available; nothing here tunes
// thresholds or weights for it. A 5-input-only (hits-only) model still
// works unchanged, since this file matches ins.size() to whatever the
// loaded session actually declares.
//
// Known model (5-input, hits-only, no calo/cand_info, 4 DGCNN layers):
//   vendor/eicrecon/share/models/prefilter-20270803.onnx
//
// Activation:
//   -Peventbuilder:prefilter:onnx_model=/path/to/model.onnx
//   -Peventbuilder:score_threshold=0.5   (default)
//   -Peventbuilder:max_hits=2500          (0 = dynamic)
//   -Peventbuilder:k_neighbors=16
//   -Peventbuilder:time_weight=1.0
//   -Peventbuilder:prefilter:store_aux_outputs=1  (0 = default: probs only)
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
#include <edm4eic/ClusterCollection.h>

#include <extensions/jana/JOmniFactory.h>
#include "services/onnx/ONNXRuntime_service.h"
#include "services/onnx/ONNXRuntime_gnn.h"

struct EventPrefilter_factory : public JOmniFactory<EventPrefilter_factory> {

  std::string m_onnx_model_path;
  float       m_score_threshold_val{0.5f};
  int32_t     m_max_hits_val{0};
  int32_t     m_k_neighbors_val{16};
  float       m_time_weight_val{1.0f};
  bool        m_veto_val{true};
  bool        m_store_aux_val{false};

  PodioInput<edm4hep::EventHeader>              m_candidates_in{this, "EventCandidates"};

  // Tracking hits: fast (TOF+MPGD) and slow (Si/B0), both gated by
  // TrkTimeCoincidence (two instances, tags `trkcoincidence` and
  // `coincidence`), so both are already cleaned to each candidate's t0,
  // not the frame-wide *TimeAlignRecHits mirror. Both are edm4eic::TrackerHit,
  // so one variadic input covers both. Names are hardcoded here, matching
  // EventBuilder_factory's own m_trk_in pattern: JOmniFactory distributes
  // one flat wiring-tag list evenly across variadic inputs, and with two
  // variadic inputs on this factory (this one plus m_calo_in below), an
  // uneven split would fail at wiring time.
  std::vector<std::string> m_trk_hit_names = {
      "TOFBarrelTrkCoincRecHits",          "TOFEndcapTrkCoincRecHits",
      "MPGDBarrelTrkCoincRecHits",         "OuterMPGDBarrelTrkCoincRecHits",
      "BackwardMPGDEndcapTrkCoincRecHits", "ForwardMPGDEndcapTrkCoincRecHits",
      "SiBarrelVertexTimeCoincRecHits",    "SiBarrelTrackerTimeCoincRecHits",
      "SiEndcapTrackerTimeCoincRecHits",   "B0TrackerTimeCoincRecHits"};
  VariadicPodioInput<edm4eic::TrackerHit, true, true> m_trk_hits_in{this, m_trk_hit_names};

  // Calo clusters, gated by CalTimeCoincidence (time, plus an optional
  // energy floor). A different EDM4hep type than tracking hits, so this
  // needs a separate variadic input. 10 entries, matching EventBuilder_
  // factory's own m_calo_collection_names, since JOmniFactory distributes
  // one flat wiring-tag list evenly across variadic inputs (10 here plus 10
  // in m_trk_hit_names above); an uneven split fails at wiring time.
  // EcalLumiSpec has no plugin yet; as an optional input it resolves to
  // nullptr until one is wired.
  std::vector<std::string> m_calo_cluster_names = {
      "B0ECalCoincClusters",           "EcalBarrelCoincClusters",
      "EcalEndcapNCoincClusters",      "EcalEndcapPCoincClusters",
      "HcalBarrelCoincClusters",       "HcalEndcapPInsertCoincClusters",
      "LFHCALCoincClusters",           "EcalFarForwardZDCCoincClusters",
      "HcalFarForwardZDCCoincClusters", "EcalLumiSpecCoincClusters"};
  VariadicPodioInput<edm4eic::Cluster, true, true> m_calo_in{this, m_calo_cluster_names};

  PodioOutput<edm4hep::EventHeader>             m_candidates_out{this, "EventCandidatesFiltered"};
  // Frame-level GNN class probabilities. Always computed and stored when a
  // model runs: one entry per frame, weights[0..6] is the 7-class softmax,
  // class 0 is background. This lets downstream code apply the prefilter
  // cut offline, without reloading the ONNX model or rebuilding the kNN
  // graph. EventUnfolder copies this onto each child PhysicsEvent. Empty in
  // pass-through mode.
  PodioOutput<edm4hep::EventHeader>             m_scores_out{this, "PrefilterScores"};

  // Raw model outputs beyond the classifier, stored only when
  // eventbuilder:prefilter:store_aux_outputs is on and the model actually
  // declares them (looked up by name; see run_gnn). Neither is trained to
  // be physically meaningful yet -- these exist for offline inspection of
  // the model as it develops, not for any current physics use.
  //   PrefilterLambdaK:   one entry per frame, weights[] = raw lambda_k row.
  //   PrefilterSegLogits: one entry per frame, weights[] = raw seg_logits
  //                       rows flattened in the same hit order as Input[0]
  //                       x (see the file header), truncated to real hits.
  // Both are empty collections when the parameter is off or the model has
  // no such output, exactly as PrefilterScores is empty in pass-through.
  PodioOutput<edm4hep::EventHeader>             m_lambda_k_out{this, "PrefilterLambdaK"};
  PodioOutput<edm4hep::EventHeader>             m_seg_logits_out{this, "PrefilterSegLogits"};

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
    pm->SetDefaultParameter("eventbuilder:prefilter:onnx_model",     m_onnx_model_path,     "Path to frame-level ONNX prefilter model. Empty = pass-through.");
    pm->SetDefaultParameter("eventbuilder:prefilter:score_threshold", m_score_threshold_val, "Maximum BKG probability to pass (GNN class-0 score).");
    pm->SetDefaultParameter("eventbuilder:prefilter:max_hits",        m_max_hits_val,        "Fixed input length for padded models (e.g. 2500). 0 = use actual hit count.");
    pm->SetDefaultParameter("eventbuilder:prefilter:k_neighbors",     m_k_neighbors_val,     "k for kNN graph construction (per DGCNN layer).");
    pm->SetDefaultParameter("eventbuilder:prefilter:time_weight",     m_time_weight_val,     "Scale factor on the time dimension in layer-0 spacetime kNN.");
    pm->SetDefaultParameter("eventbuilder:prefilter:veto",            m_veto_val,            "1 = drop frames with BKG prob >= score_threshold (online veto); 0 = score-only, keep every candidate and store the scores for a reversible OFFLINE cut (recommended for dataset production).");
    pm->SetDefaultParameter("eventbuilder:prefilter:store_aux_outputs", m_store_aux_val,      "1 = also store PrefilterLambdaK/PrefilterSegLogits, the model's raw lambda_k/seg_logits outputs, when it declares them. Neither is trained/meaningful yet; stored as-is for offline inspection. 0 (default) = only the classifier (probs, PrefilterScores) is stored.");

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

  void Process(int32_t, uint64_t) {
    const auto* cands = m_candidates_in();
    if (!cands) return;

    m_candidates_out()->setSubsetCollection(true);

    if (!m_session) {
      for (const auto& h : *cands) m_candidates_out()->push_back(h);
      return;
    }

    const int64_t max = static_cast<int64_t>(m_max_hits_val);

    // det_id categories, all divided by 6 to match the normalization
    // scheme. TOF=2/6 and MPGD=1/6 match the original 2-category scheme
    // that the one known trained model uses. New categories are appended
    // at previously unused values, so they do not disturb what that model
    // already learned:
    //   TOF=2/6 (unchanged)  MPGD=1/6 (unchanged)  Si=3/6  B0=4/6  Calo=5/6
    // No model has yet been trained against the 3 new categories; see the
    // file header.
    std::vector<float> feat;
    size_t n_hits = 0;
    bool cap_reached = false;
    const auto at_cap = [&] { return max > 0 && static_cast<int64_t>(n_hits) >= max; };

    const auto& trk_colls = m_trk_hits_in();
    for (size_t ci = 0; ci < trk_colls.size() && !cap_reached; ++ci) {
      const auto* coll = trk_colls[ci];
      if (!coll) continue;
      float det_id_norm;
      if (ci <= 1)      det_id_norm = 2.f / 6.f; // TOF (unchanged)
      else if (ci <= 5) det_id_norm = 1.f / 6.f; // MPGD (unchanged)
      else if (ci <= 8) det_id_norm = 3.f / 6.f; // Si
      else              det_id_norm = 4.f / 6.f; // B0
      for (const auto& h : *coll) {
        if (at_cap()) { cap_reached = true; break; }
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
    const auto& calo_colls = m_calo_in();
    for (size_t ci = 0; ci < calo_colls.size() && !cap_reached; ++ci) {
      const auto* coll = calo_colls[ci];
      if (!coll) continue;
      const float det_id_norm = 5.f / 6.f; // Calo (all detectors, one category)
      for (const auto& c : *coll) {
        if (at_cap()) { cap_reached = true; break; }
        const auto& pos = c.getPosition();
        const float e = std::min(1.f, std::max(0.f, c.getEnergy() / 10.f));
        feat.push_back(static_cast<float>(pos[0]) / 4000.f);
        feat.push_back(static_cast<float>(pos[1]) / 4000.f);
        feat.push_back(static_cast<float>(pos[2]) / 5000.f);
        feat.push_back(c.getTime()                / 50.f);
        feat.push_back(e);
        feat.push_back(det_id_norm);
        ++n_hits;
      }
    }

    const size_t target = (max > 0) ? static_cast<size_t>(max) : n_hits;
    feat.resize(target * static_cast<size_t>(m_n_features), 0.0f);
    if (target == 0) {
      // No hits to feed the GNN. Treated as background only when the
      // online veto is active. In score-only mode (veto off), every
      // candidate is still kept, just without a GNN score for this frame.
      if (!m_veto_val)
        for (const auto& h : *cands) m_candidates_out()->push_back(h);
      return;
    }

    Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    run_gnn(cands, feat, static_cast<int>(target), n_hits, mem);
  }

private:
  // Per-candidate event-level summary: weights[0/1/4/18/19/8/7/9/10], the
  // same fields already stored in EventBuilderInfo (t0, t0sigma, S, mu, k,
  // cmax_tof, cmax, E_calo, n_tracklets — cmax_tof stands in for the
  // retired cluster_size, same feature slot). One row per candidate in this
  // frame, sent as an optional 6th ONNX input (see run_gnn), so a future
  // model can use frame/candidate context that a pure hit graph cannot see
  // on its own.
  static std::vector<float> build_cand_info(const edm4hep::EventHeaderCollection* cands) {
    static constexpr int IDX[9] = {0, 1, 4, 18, 19, 8, 7, 9, 10};
    std::vector<float> info;
    info.reserve(cands->size() * 9);
    for (const auto& c : *cands) {
      const auto w = c.getWeights();
      for (int idx : IDX)
        info.push_back(static_cast<size_t>(idx) < w.size()
                            ? static_cast<float>(w[idx]) : 0.0f);
    }
    return info;
  }

  void run_gnn(const edm4hep::EventHeaderCollection* cands,
               std::vector<float>& feat, int N, size_t n_hits,
               Ort::MemoryInfo& mem)
  {
    // Build k-NN edge indices for all 4 DGCNN layers.
    auto edges = m_gnn.build_dgcnn_edges(
        feat.data(), N, m_k_neighbors_val, m_time_weight_val,
        static_cast<int>(m_n_features));
    const int64_t E = static_cast<int64_t>(edges[0].row.size());

    // Flatten each EdgeIndex into a [2, E] int64 buffer (row-first).
    // Build in-place: reuse edges[i].row as the combined [row||col] buffer.
    std::array<std::vector<int64_t>, 4> ei_buf;
    for (int l = 0; l < 4; ++l) {
      ei_buf[l].reserve(2 * E);
      ei_buf[l].insert(ei_buf[l].end(), edges[l].row.begin(), edges[l].row.end());
      ei_buf[l].insert(ei_buf[l].end(), edges[l].col.begin(), edges[l].col.end());
    }

    const std::vector<int64_t> x_sh{N, m_n_features};
    const std::vector<int64_t> ei_sh{2, E};

    std::vector<Ort::Value> ins;
    ins.push_back(Ort::Value::CreateTensor<float>(
        mem, feat.data(), feat.size(), x_sh.data(), x_sh.size()));
    for (int l = 0; l < 4; ++l)
      ins.push_back(Ort::Value::CreateTensor<int64_t>(
          mem, ei_buf[l].data(), ei_buf[l].size(), ei_sh.data(), ei_sh.size()));

    // Optional 6th input (cand_info): sent only if the loaded model
    // actually declares that many inputs, so a 5-input-only (hits-only)
    // model still runs unmodified. ins.size() must equal m_in_names.size()
    // for Run().
    std::vector<float> cand_info;
    std::vector<int64_t> ci_sh;
    if (m_in_names.size() >= 6) {
      cand_info = build_cand_info(cands);
      const int64_t C = static_cast<int64_t>(cands->size());
      ci_sh = {C, 9};
      if (C == 0) cand_info.resize(0); // empty candidate list — zero-row tensor
      ins.push_back(Ort::Value::CreateTensor<float>(
          mem, cand_info.data(), cand_info.size(), ci_sh.data(), ci_sh.size()));
    }

    // Fetch every output the model declares. Only Output[0] (probs) is
    // required; lambda_k/seg_logits, if present, are optional and located
    // by name below, not by position.
    auto out = m_session->Run(
        Ort::RunOptions{nullptr}, m_in_names.data(), ins.data(), ins.size(),
        m_out_names.data(), m_out_names.size());

    // Output[0] = probs [1, N]; class 0 = background. Every class is
    // stored per frame, so the score works as both a reversible offline cut
    // and a physics-class tag. The optional online veto is decided after.
    const float* probs = out[0].GetTensorMutableData<float>();
    const auto shp = out[0].GetTensorTypeAndShapeInfo().GetShape();
    const int64_t n_cls = shp.empty() ? 7 : shp.back();
    edm4hep::MutableEventHeader scores;
    for (int64_t c = 0; c < n_cls; ++c)
      scores.addToWeights(static_cast<double>(probs[c]));
    m_scores_out()->push_back(scores);

    if (m_store_aux_val) {
      for (size_t oi = 0; oi < m_out_name_strs.size(); ++oi) {
        if (m_out_name_strs[oi] == "lambda_k") {
          const float* lk = out[oi].GetTensorMutableData<float>();
          const auto lk_shp = out[oi].GetTensorTypeAndShapeInfo().GetShape();
          const int64_t n_lk = lk_shp.empty() ? 0 : lk_shp.back();
          edm4hep::MutableEventHeader hdr;
          for (int64_t c = 0; c < n_lk; ++c)
            hdr.addToWeights(static_cast<double>(lk[c]));
          m_lambda_k_out()->push_back(hdr);
        } else if (m_out_name_strs[oi] == "seg_logits") {
          const float* sl = out[oi].GetTensorMutableData<float>();
          const auto sl_shp = out[oi].GetTensorTypeAndShapeInfo().GetShape();
          // sl_shp = [rows, per_hit_dim]; rows includes any max_hits zero-
          // padding, so only the first n_hits rows correspond to real hits.
          const int64_t per_hit = (sl_shp.size() >= 2) ? sl_shp.back() : 0;
          const int64_t n_rows  = sl_shp.empty() ? 0 : sl_shp[0];
          const int64_t n_real  = std::min<int64_t>(static_cast<int64_t>(n_hits), n_rows);
          edm4hep::MutableEventHeader hdr;
          for (int64_t r = 0; r < n_real; ++r)
            for (int64_t c = 0; c < per_hit; ++c)
              hdr.addToWeights(static_cast<double>(sl[r * per_hit + c]));
          m_seg_logits_out()->push_back(hdr);
        }
      }
    }

    const float bkg_prob = probs[0];
    // Score-only mode (veto off): keep every candidate; the stored scores
    // let the veto be applied offline instead, reversibly.
    if (!m_veto_val || bkg_prob < m_score_threshold_val)
      for (const auto& h : *cands) m_candidates_out()->push_back(h);
  }
};

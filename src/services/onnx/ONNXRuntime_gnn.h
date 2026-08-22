// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025, Jefferson Science Associates, LLC.

// ONNXRuntime_gnn — k-NN graph construction for GNN-based event prefilter.
//
// Builds directed edge indices in the nanoflann / pyg-lib ≥ 0.6 convention:
//   row[i] = source (neighbor), col[i] = target (query)
//   edge_index shape: [2, N*k], layout: all sources followed by all targets
//
// Two compile-time KDTree backends (selected by CMake):
//   USE_NANOFLANN  — uses nanoflann::KDTreeSingleIndexAdaptor via a thin
//                    PointCloud adaptor.  CMake fetches nanoflann v1.6 via
//                    FetchContent if not found as a system package.
//   (default)      — built-in compact k-d tree; no external dependency.

#pragma once

#include <array>
#include <cstdint>
#include <vector>

// Edge index: row = source (neighbor), col = target (query), length E = N*k each.
struct EdgeIndex {
    std::vector<int64_t> row;
    std::vector<int64_t> col;
};

class ONNXRuntime_gnn {
public:
    // Build kNN edge index for N points using the first D dimensions of each row.
    // stride: floats per row in x (0 = D, i.e. tightly packed). Allows distance
    //   computation over a coordinate prefix of a wider feature matrix,
    //   e.g. kNN in (x,y,z,t) over rows of [x,y,z,t,eDep,det_id].
    // metric_weights: optional per-dim scale applied before distance computation
    //   (e.g. set weights[3] = time_weight to scale the time dimension for spacetime kNN).
    // Edge convention: source=neighbor, target=query (pyg-lib source_to_target flow).
    EdgeIndex build(const float* x, int N, int D, int k,
                    const float* metric_weights = nullptr,
                    int stride = 0) const;

    // Build the 3 edge sets required by a 3-layer DGCNN GNN:
    //   [0]: kNN in 4D spacetime (x,y,z,t) with optional time_weight on dim 3
    //   [1]: same spatial graph (approximation; exact would need layer-0 activations)
    //   [2]: same spatial graph (approximation; exact would need layer-1 activations)
    // x is a row-major [N, n_features] matrix; the first 4 columns must be (x,y,z,t).
    std::array<EdgeIndex, 3> build_dgcnn_edges(
        const float* x, int N, int k = 16, float time_weight = 1.0f,
        int n_features = 6) const;
};

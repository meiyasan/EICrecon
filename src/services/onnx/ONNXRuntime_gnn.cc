// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025, Jefferson Science Associates, LLC.

#include "ONNXRuntime_gnn.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>

// ---------------------------------------------------------------------------
// kNN graph construction — two compile-time paths:
//
//   USE_NANOFLANN (defined by CMake when nanoflann is found):
//     Uses nanoflann::KDTreeSingleIndexAdaptor via a lightweight PointCloud
//     adaptor.  When metric_weights != nullptr the coordinates are pre-scaled
//     into a temporary buffer so the standard L2 metric produces weighted
//     distances without a custom distance type.
//     Requires nanoflann >= 1.5  (SearchParameters API).
//
//   Fallback (no nanoflann):
//     Compact axis-aligned k-d tree written to the nanoflann point-cloud
//     adapter contract (kdtree_get_point_count / kdtree_get_pt /
//     kdtree_get_bbox).  Complexity O(N log N) build, O(k log N) per query.
//     Drop in the real nanoflann header by defining USE_NANOFLANN.
// ---------------------------------------------------------------------------

#ifdef USE_NANOFLANN

#include <nanoflann.hpp>

namespace {

// Thin adaptor wrapping a row-major float* array for nanoflann.
struct PointCloud {
    const float* pts;   // row-major, `stride` floats per row
    int          N;
    int          D;
    int          stride;

    size_t kdtree_get_point_count() const { return static_cast<size_t>(N); }
    float  kdtree_get_pt(size_t idx, int dim) const {
        return pts[idx * stride + dim];
    }
    template<class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

} // namespace

EdgeIndex ONNXRuntime_gnn::build(const float* x, int N, int D, int k,
                                  const float* metric_weights, int stride) const
{
    EdgeIndex ei;
    if (N <= 1 || k <= 0) return ei;
    if (stride <= 0) stride = D;
    const int eff_k = std::min(k, N - 1);
    ei.row.resize(static_cast<size_t>(N) * eff_k);
    ei.col.resize(static_cast<size_t>(N) * eff_k);

    // When metric_weights are provided, pre-scale a D-column copy of the
    // coordinates.  The L2 metric on the scaled space equals the weighted L2
    // on the original space, avoiding a custom distance type.
    std::vector<float> scaled;
    const float* search_pts    = x;
    int          search_stride = stride;
    if (metric_weights) {
        scaled.resize(static_cast<size_t>(N) * D);
        for (int i = 0; i < N; ++i)
            for (int d = 0; d < D; ++d)
                scaled[i * D + d] = x[i * stride + d] * metric_weights[d];
        search_pts    = scaled.data();
        search_stride = D;
    }

    using KDTreeT = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, PointCloud>,
        PointCloud, -1 /* dynamic dims */>;

    PointCloud cloud{search_pts, N, D, search_stride};
    KDTreeT    tree(D, cloud, {16 /*max_leaf*/});
    tree.buildIndex();

    const int        query_k = eff_k + 1; // +1: self may appear in results
    std::vector<size_t> idx(query_k);
    std::vector<float>  dist2(query_k);
    std::vector<float>  q(D);

    for (int i = 0; i < N; ++i) {
        // Build query point in (possibly scaled) search space.
        if (metric_weights)
            for (int d = 0; d < D; ++d) q[d] = x[i * stride + d] * metric_weights[d];
        else
            for (int d = 0; d < D; ++d) q[d] = x[i * stride + d];

        nanoflann::KNNResultSet<float> result_set(query_k);
        result_set.init(idx.data(), dist2.data());
        tree.findNeighbors(result_set, q.data(), nanoflann::SearchParameters{});

        int written = 0;
        for (size_t m = 0; m < result_set.size() && written < eff_k; ++m) {
            const int j = static_cast<int>(idx[m]);
            if (j == i) continue;
            ei.row[i * eff_k + written] = j;
            ei.col[i * eff_k + written] = i;
            ++written;
        }
        // Pad with last valid neighbor if fewer than eff_k unique neighbors found.
        const int last = std::max(0, written - 1);
        for (int w = written; w < eff_k; ++w) {
            ei.row[i * eff_k + w] = ei.row[i * eff_k + last];
            ei.col[i * eff_k + w] = i;
        }
    }
    return ei;
}

#else // ----------------------------------------------------------------
// Fallback: built-in compact k-d tree (no external dependency).
// ---------------------------------------------------------------------------

namespace {

using HeapEntry = std::pair<float, int>; // (distance², point-index)

struct KDNode {
    int   split_dim{-1};   // -1 = leaf
    float split_val{0.f};
    int   left{-1};
    int   right{-1};
    int   pt_begin{-1};    // leaf: range [pt_begin, pt_end) into m_indices
    int   pt_end{-1};
};

class KDTree {
public:
    KDTree(const float* pts, int N, int D, int stride,
           const float* weights = nullptr,
           int leaf_size = 16)
        : m_pts(pts), m_N(N), m_D(D), m_stride(stride),
          m_weights(weights), m_leaf_size(leaf_size)
    {
        m_indices.resize(N);
        std::iota(m_indices.begin(), m_indices.end(), 0);
        if (N > 0) build(0, N, 0);
    }

    void knn(const float* query, int k, std::vector<HeapEntry>& result) const {
        std::priority_queue<HeapEntry> heap;
        float worst = std::numeric_limits<float>::infinity();
        search(0, query, k, heap, worst);
        result.resize(heap.size());
        for (int i = static_cast<int>(heap.size()) - 1; i >= 0; --i) {
            result[i] = heap.top();
            heap.pop();
        }
    }

private:
    const float*        m_pts;
    int                 m_N, m_D, m_stride;
    const float*        m_weights;
    int                 m_leaf_size;
    std::vector<int>    m_indices;
    std::vector<KDNode> m_nodes;

    float dist2(const float* a, const float* b) const {
        float d = 0.f;
        for (int i = 0; i < m_D; ++i) {
            float dx = a[i] - b[i];
            if (m_weights) dx *= m_weights[i];
            d += dx * dx;
        }
        return d;
    }

    int build(int begin, int end, int /*depth*/) {
        const int node_idx = static_cast<int>(m_nodes.size());
        m_nodes.emplace_back();
        KDNode& node = m_nodes.back();

        if (end - begin <= m_leaf_size) {
            node.pt_begin = begin;
            node.pt_end   = end;
            return node_idx;
        }

        std::vector<float> lo(m_D,  std::numeric_limits<float>::infinity());
        std::vector<float> hi(m_D, -std::numeric_limits<float>::infinity());
        for (int i = begin; i < end; ++i) {
            const float* p = m_pts + m_indices[i] * m_stride;
            for (int d = 0; d < m_D; ++d) {
                lo[d] = std::min(lo[d], p[d]);
                hi[d] = std::max(hi[d], p[d]);
            }
        }
        int   best_dim = 0;
        float best_rng = -1.f;
        for (int d = 0; d < m_D; ++d) {
            float rng = hi[d] - lo[d];
            if (m_weights) rng *= m_weights[d];
            if (rng > best_rng) { best_rng = rng; best_dim = d; }
        }

        int mid = (begin + end) / 2;
        std::nth_element(m_indices.begin() + begin,
                         m_indices.begin() + mid,
                         m_indices.begin() + end,
                         [&](int a, int b) {
                             return m_pts[a * m_stride + best_dim]
                                  < m_pts[b * m_stride + best_dim];
                         });

        node.split_dim = best_dim;
        node.split_val = m_pts[m_indices[mid] * m_stride + best_dim];

        int left_idx  = build(begin, mid, best_dim + 1);
        int right_idx = build(mid,   end, best_dim + 1);
        m_nodes[node_idx].left  = left_idx;
        m_nodes[node_idx].right = right_idx;
        return node_idx;
    }

    void search(int node_idx, const float* query, int k,
                std::priority_queue<HeapEntry>& heap, float& worst) const
    {
        const KDNode& node = m_nodes[node_idx];

        if (node.split_dim < 0) {
            for (int i = node.pt_begin; i < node.pt_end; ++i) {
                float d = dist2(query, m_pts + m_indices[i] * m_stride);
                if (static_cast<int>(heap.size()) < k) {
                    heap.push({d, m_indices[i]});
                    if (static_cast<int>(heap.size()) == k)
                        worst = heap.top().first;
                } else if (d < worst) {
                    heap.pop();
                    heap.push({d, m_indices[i]});
                    worst = heap.top().first;
                }
            }
            return;
        }

        float diff = query[node.split_dim] - node.split_val;
        if (m_weights) diff *= m_weights[node.split_dim];
        const bool go_left = (diff <= 0.f);
        const int  near    = go_left ? node.left  : node.right;
        const int  far     = go_left ? node.right : node.left;

        search(near, query, k, heap, worst);
        if (static_cast<int>(heap.size()) < k || diff * diff < worst)
            search(far, query, k, heap, worst);
    }
};

} // namespace

EdgeIndex ONNXRuntime_gnn::build(const float* x, int N, int D, int k,
                                  const float* metric_weights, int stride) const
{
    EdgeIndex ei;
    if (N <= 1 || k <= 0) return ei;
    if (stride <= 0) stride = D;
    const int eff_k = std::min(k, N - 1);
    ei.row.resize(static_cast<size_t>(N) * eff_k);
    ei.col.resize(static_cast<size_t>(N) * eff_k);

    KDTree tree(x, N, D, stride, metric_weights);
    std::vector<HeapEntry> nbrs;
    nbrs.reserve(eff_k);

    for (int i = 0; i < N; ++i) {
        nbrs.clear();
        tree.knn(x + i * stride, eff_k + 1, nbrs);
        int written = 0;
        for (const auto& [d2, j] : nbrs) {
            if (j == i) continue;
            if (written >= eff_k) break;
            ei.row[i * eff_k + written] = j;
            ei.col[i * eff_k + written] = i;
            ++written;
        }
        const int last = std::max(0, written - 1);
        for (int w = written; w < eff_k; ++w) {
            ei.row[i * eff_k + w] = ei.row[i * eff_k + last];
            ei.col[i * eff_k + w] = i;
        }
    }
    return ei;
}

#endif // USE_NANOFLANN

// ---------------------------------------------------------------------------
// Shared: build the 3 edge sets required by a 3-layer DGCNN.
// ---------------------------------------------------------------------------

std::array<EdgeIndex, 3> ONNXRuntime_gnn::build_dgcnn_edges(
    const float* x, int N, int k, float time_weight, int n_features) const
{
    const float weights[4] = {1.f, 1.f, 1.f, time_weight};
    EdgeIndex e0 = build(x, N, 4, k, weights, n_features);
    return {e0, e0, e0};
}

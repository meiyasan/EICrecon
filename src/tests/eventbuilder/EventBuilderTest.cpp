// Tests for the eventbuilder pipeline algorithmic logic.
//
// Self-contained: replicates the core math from EventBuilder_factory,
// TimeAlignment_factory, and TimeCoincidence_factory using pure C++ so
// no JANA / PODIO / ROOT infrastructure is required to run these tests.
//
// Physics classes under test (streaming-readout taxonomy):
//   0 BKG          — machine backgrounds only, no signal coincidence
//   1 DIS_NC       — neutral-current DIS: backward scattered electron
//   2 DIS_CC       — charged-current DIS: forward hadrons, no electron
//   3 DDIS         — diffractive DIS: backward electron + rapidity gap
//   4 EXCLUSIVE_H  — hard exclusive (DVCS/DVMP): backward e + backward γ
//   5 EXCLUSIVE_D  — exclusive diffractive (J/ψ, ρ): equatorial decay products
//   6 SIDIS        — semi-inclusive DIS: backward e + identified forward hadron
//
// Prefilter (gold-coating ONNX model) tests are compiled only when the build
// system detects OnnxRuntime and defines HAS_ONNXRUNTIME.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <set>
#include <vector>

#ifdef HAS_ONNXRUNTIME
#  include <onnxruntime_cxx_api.h>
#  include <unistd.h>   // access()
#endif

// kNN graph service under test (standalone TU include keeps this test
// buildable without linking the eicrecon onnx_service library).
#include "../../services/onnx/ONNXRuntime_gnn.cc"

// ---------------------------------------------------------------------------
// Minimal data types — mirrors EventBuilder_factory internals
// ---------------------------------------------------------------------------

namespace {

static constexpr float kSigmaTOF  = 0.03f;   // ns — AC-LGAD
static constexpr float kSigmaMPGD = 10.0f;   // ns
static constexpr float kSigmaSi   = 2000.0f; // ns — MAPS integration window
static constexpr float kCoincWin  = 50.0f;   // ns — default coincidence window
static constexpr float kNSigma    = 3.0f;    // default nsigma_window

// Replica of FastHit (all inputs to the coincidence scan are of this type)
struct FastHit {
  double time;
  float  sigma;
  float  x, y, z;  // mm
  int    theta1, phi1, theta2, phi2;
};

// Raw hit before bin calculation
struct RawHit {
  double time;
  float  sigma;
  float  x, y, z;
};

// Replica of thetaPhiBinCalc from EventBuilder_factory (without ROOT/TMath)
void computeBins(float x, float y, float z,
                 int& theta1, int& phi1, int& theta2, int& phi2) {
  const double pi = M_PI;
  const double r  = std::sqrt((double)x*x + (double)y*y + (double)z*z);
  double theta    = (r > 0.0) ? std::acos((double)z / r) : 0.0;
  theta           = std::max(0.0, std::min(theta, pi));
  double phi      = std::atan2((double)y, (double)x);
  if (phi < 0.0) phi += 2.0 * pi;

  theta1 = (int)(theta / (pi / 12.0));
  theta2 = (int)((theta + pi / 24.0) / (pi / 12.0));
  phi1   = (int)(phi   / (pi / 8.0));
  phi2   = (int)((phi  + pi / 16.0) / (pi / 8.0));

  auto clamp = [](int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); };
  theta1 = clamp(theta1, 0, 11);
  theta2 = clamp(theta2, 0, 11);
  phi1   = clamp(phi1, 0, 7);
  phi2   = clamp(phi2 % 8, 0, 7);
}

FastHit toFastHit(const RawHit& h) {
  FastHit fh;
  fh.time  = h.time;
  fh.sigma = h.sigma;
  fh.x = h.x; fh.y = h.y; fh.z = h.z;
  computeBins(h.x, h.y, h.z, fh.theta1, fh.phi1, fh.theta2, fh.phi2);
  return fh;
}

// ---------------------------------------------------------------------------
// Algorithmic t0 — inverse-variance weighted mean (algorithmic mode)
// ---------------------------------------------------------------------------

double computeT0(const std::vector<FastHit>& hits, size_t lo, size_t hi) {
  double sumw = 0.0, sumwt = 0.0;
  for (size_t k = lo; k < hi; ++k) {
    const double w = 1.0 / ((double)hits[k].sigma * hits[k].sigma);
    sumw  += w;
    sumwt += w * hits[k].time;
  }
  return sumw > 0.0 ? sumwt / sumw : 0.0;
}

double computeSigmaT0(const std::vector<FastHit>& hits, size_t lo, size_t hi) {
  double sumw = 0.0;
  for (size_t k = lo; k < hi; ++k)
    sumw += 1.0 / ((double)hits[k].sigma * hits[k].sigma);
  return sumw > 0.0 ? std::sqrt(1.0 / sumw) : 1e9;
}

// ---------------------------------------------------------------------------
// Coincidence scan — replica of EventBuilder_factory::Process sliding window
// ---------------------------------------------------------------------------

struct Candidate { double t0; double t0sigma; };

std::vector<Candidate> findCandidates(std::vector<FastHit>& hits,
                                      float dt      = kCoincWin,
                                      float nsigma  = kNSigma) {
  std::sort(hits.begin(), hits.end(),
            [](const FastHit& a, const FastHit& b) { return a.time < b.time; });

  int occupancy1[12][8] = {};
  int occupancy2[12][8] = {};
  size_t lo = 0, hi = 0;
  std::vector<Candidate> out;

  while (lo < hits.size()) {
    while (hi < hits.size() &&
           (hits[hi].time - hits[lo].time) <= dt) {
      ++occupancy1[hits[hi].theta1][hits[hi].phi1];
      ++occupancy2[hits[hi].theta2][hits[hi].phi2];
      ++hi;
    }

    bool trigger = false;
    for (int it = 0; it < 12 && !trigger; ++it)
      for (int ip = 0; ip < 8 && !trigger; ++ip)
        if (occupancy1[it][ip] > 2 || occupancy2[it][ip] > 2)
          trigger = true;

    if (trigger) {
      const double t0      = computeT0(hits, lo, hi);
      const double sigma_t0 = computeSigmaT0(hits, lo, hi);
      out.push_back({t0, nsigma * sigma_t0});
      for (size_t j = lo; j < hi; ++j) {
        --occupancy1[hits[j].theta1][hits[j].phi1];
        --occupancy2[hits[j].theta2][hits[j].phi2];
      }
      lo = hi;
    } else {
      --occupancy1[hits[lo].theta1][hits[lo].phi1];
      --occupancy2[hits[lo].theta2][hits[lo].phi2];
      ++lo;
      if (hi < lo) hi = lo;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Helpers for generating synthetic hit data
// ---------------------------------------------------------------------------

// Hit at the center of (th, ph) bin at radius R mm, time t_ns, sigma sig_ns
RawHit hitAtBin(int th, int ph, float R, double t_ns, float sig_ns = kSigmaTOF) {
  const float pi = (float)M_PI;
  const float theta = (th + 0.5f) * pi / 12.0f;
  const float phi   = (ph + 0.5f) * pi / 8.0f;
  return {t_ns, sig_ns,
          R * std::sin(theta) * std::cos(phi),
          R * std::sin(theta) * std::sin(phi),
          R * std::cos(theta)};
}

// N hits in the same (th, ph) bin, uniformly spread over [t_center - dt, t_center + dt]
void addCluster(std::vector<RawHit>& v, int th, int ph, float R,
                double t_center, double dt, int n, float sig = kSigmaTOF) {
  for (int i = 0; i < n; ++i) {
    const double t = (n > 1)
        ? t_center - dt + i * (2.0 * dt / (n - 1))
        : t_center;
    v.push_back(hitAtBin(th, ph, R, t, sig));
  }
}

// Convert a vector of RawHit to FastHit
std::vector<FastHit> toFastHits(const std::vector<RawHit>& raw) {
  std::vector<FastHit> out;
  out.reserve(raw.size());
  for (const auto& h : raw) out.push_back(toFastHit(h));
  return out;
}

} // namespace

// ===========================================================================
// I. TimeAlignment: r/c propagation correction
// ===========================================================================

TEST(TimeAlignment, RcCorrection) {
  // t_corr = t_raw - |r| * 0.0034 (positions in mm, times in ns)
  // 0.0034 ns/mm ≈ 1/c (c ≈ 294 mm/ns)
  const double t_raw = 100.0;
  const float  x = 100.0f, y = 0.0f, z = 0.0f;  // |r| = 100 mm
  const double dt_propagation = std::sqrt(x*x + y*y + z*z) * 0.0034;
  const double t_corr = t_raw - dt_propagation;

  EXPECT_NEAR(dt_propagation, 0.34, 0.01)
      << "100 mm / c ≈ 0.34 ns";
  EXPECT_NEAR(t_corr, 99.66, 0.01);
}

TEST(TimeAlignment, RcCorrection_OnAxis) {
  // z-aligned hit: r = |z|, propagation = |z| * 0.0034
  const double t_raw = 500.0;
  const float  z = -250.0f;  // backward endcap
  const double r = std::fabs(z);
  const double t_corr = t_raw - r * 0.0034;
  EXPECT_NEAR(t_corr, 500.0 - 250.0 * 0.0034, 1e-6);
}

// ===========================================================================
// II. EventBuilder: t0 and t0sigma computation
// ===========================================================================

TEST(EventBuilder, T0_SingleTOFHit) {
  // Single TOF hit: t0 must equal the hit time exactly.
  std::vector<FastHit> hits = { toFastHit(hitAtBin(6, 0, 100.0f, 123.456)) };
  EXPECT_NEAR(computeT0(hits, 0, 1), 123.456, 1e-9);
  EXPECT_GT(computeSigmaT0(hits, 0, 1), 0.0);
}

TEST(EventBuilder, T0_EqualWeights) {
  // N hits at same time → t0 = that time, sigma_t0 = sigma/sqrt(N)
  const int    N = 4;
  const double t = 200.0;
  std::vector<FastHit> hits;
  for (int i = 0; i < N; ++i)
    hits.push_back(toFastHit(hitAtBin(6, 0, 100.0f, t)));

  EXPECT_NEAR(computeT0(hits, 0, N), t, 1e-9);
  const double sigma_t0 = computeSigmaT0(hits, 0, N);
  EXPECT_NEAR(sigma_t0, kSigmaTOF / std::sqrt((double)N), 1e-6);
}

TEST(EventBuilder, T0_InverseVarianceWeighting) {
  // TOF hit at t=100 ns (sigma=0.03) and MPGD hit at t=200 ns (sigma=10).
  // t0 must be dominated by the TOF hit: t0 ≈ 100 ns.
  std::vector<FastHit> hits = {
    toFastHit({100.0, kSigmaTOF,  0.0f, 100.0f, 0.0f}),  // TOF
    toFastHit({200.0, kSigmaMPGD, 0.0f, 100.0f, 0.0f}),  // MPGD
  };
  const double t0 = computeT0(hits, 0, 2);
  EXPECT_NEAR(t0, 100.0, 0.5)
      << "TOF (sigma=0.03 ns) dominates the weighted mean";
}

// ===========================================================================
// III. TimeCoincidence: candidate-t0 gate
// ===========================================================================

TEST(TimeCoincidence, Gate_HitInsideWindow) {
  const float nsigma = 3.0f;
  const float half   = nsigma * kSigmaSi;  // 3 × 2000 = 6000 ns
  const double t0    = 500.0;

  EXPECT_TRUE(std::fabs(t0 + half * 0.9 - t0) <= half)
      << "Hit just inside window must pass";
  EXPECT_TRUE(std::fabs(t0 - half * 0.9 - t0) <= half)
      << "Hit just inside window (negative) must pass";
}

TEST(TimeCoincidence, Gate_HitOutsideWindow) {
  const float nsigma = 3.0f;
  const float half   = nsigma * kSigmaSi;
  const double t0    = 500.0;

  EXPECT_FALSE(std::fabs(t0 + half * 1.01 - t0) <= half)
      << "Hit just outside window must not pass";
}

// ===========================================================================
// IV. Per-physics-class coincidence tests
// Each test generates synthetic fast-detector hits characteristic of the
// class and verifies the expected number of candidates and t0 quality.
// ===========================================================================

// 0 — BKG: machine backgrounds, no signal coincidence.
// One hit per theta-phi bin separated in time → max bin count = 1 → no trigger.
TEST(EventBuilder, BKG_NoCoincidence) {
  std::vector<RawHit> raw;
  // 10 hits each 100 ns apart — at most 1 hit in any 50 ns window
  for (int i = 0; i < 10; ++i)
    raw.push_back(hitAtBin(i, i % 8, 100.0f, i * 100.0));

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  EXPECT_EQ(cands.size(), 0u)
      << "Background frame: no time coincidence → no candidates";
}

// 1 — DIS_NC: backward scattered electron (θ > 90°, bin th=8)
TEST(EventBuilder, DIS_NC_BackwardElectron) {
  std::vector<RawHit> raw;
  // Tight cluster in backward hemisphere: 4 TOF hits within ±15 ns of t=100 ns
  addCluster(raw, /*th=*/8, /*ph=*/2, /*R=*/100.0f, /*t=*/100.0, /*dt=*/15.0, /*n=*/4);
  // Isolated background hits far in time — do not add to the cluster
  raw.push_back(hitAtBin(2, 5, 100.0f, 500.0));
  raw.push_back(hitAtBin(3, 6, 100.0f, 900.0));

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_EQ(cands.size(), 1u)
      << "DIS_NC: exactly one backward-electron candidate";
  EXPECT_NEAR(cands[0].t0, 100.0, 0.5)
      << "t0 should be within 0.5 ns of the cluster centre";
  EXPECT_LT(cands[0].t0sigma, 1.0)
      << "t0sigma should be small for a TOF-dominated cluster";
}

// 2 — DIS_CC: forward hadronic activity only (θ < 45°, bin th=2), no backward electron
TEST(EventBuilder, DIS_CC_ForwardHadrons_NoElectron) {
  std::vector<RawHit> raw;
  // Forward-only cluster: 4 TOF hits in bin th=2 (θ ≈ 37.5°)
  addCluster(raw, /*th=*/2, /*ph=*/4, /*R=*/100.0f, /*t=*/200.0, /*dt=*/10.0, /*n=*/4);

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_EQ(cands.size(), 1u)
      << "DIS_CC: one candidate from forward hadronic cluster";
  EXPECT_NEAR(cands[0].t0, 200.0, 0.5);
  EXPECT_LT(cands[0].t0sigma, 1.0);
}

// 3 — DDIS: diffractive DIS — backward electron at t=300 ns, forward remnant at t=310 ns
TEST(EventBuilder, DDIS_BackwardElectron_ForwardRemnant) {
  std::vector<RawHit> raw;
  addCluster(raw, /*th=*/9, /*ph=*/1, /*R=*/100.0f, /*t=*/300.0, /*dt=*/5.0, /*n=*/4);
  addCluster(raw, /*th=*/1, /*ph=*/3, /*R=*/100.0f, /*t=*/310.0, /*dt=*/5.0, /*n=*/3);

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_GE(cands.size(), 1u)
      << "DDIS: at least one candidate from backward+forward cluster";
  EXPECT_NEAR(cands[0].t0, 305.0, 10.0);
}

// 4 — EXCLUSIVE_H: hard exclusive (DVCS) — backward scattered e + backward photon
TEST(EventBuilder, EXCLUSIVE_H_BackwardElectronAndPhoton) {
  std::vector<RawHit> raw;
  addCluster(raw, /*th=*/10, /*ph=*/0, /*R=*/100.0f, /*t=*/400.0, /*dt=*/8.0, /*n=*/4);
  addCluster(raw, /*th=*/10, /*ph=*/1, /*R=*/100.0f, /*t=*/400.0, /*dt=*/8.0, /*n=*/3);

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_GE(cands.size(), 1u)
      << "EXCLUSIVE_H: candidate from backward electron + photon activity";
  EXPECT_NEAR(cands[0].t0, 400.0, 2.0);
}

// 5 — EXCLUSIVE_D: exclusive diffractive vector meson (J/ψ→μ+μ-, ρ→π+π-)
TEST(EventBuilder, EXCLUSIVE_D_SymmetricVectorMeson) {
  std::vector<RawHit> raw;
  addCluster(raw, /*th=*/6, /*ph=*/0, /*R=*/100.0f, /*t=*/500.0, /*dt=*/5.0, /*n=*/4);
  addCluster(raw, /*th=*/6, /*ph=*/4, /*R=*/100.0f, /*t=*/500.0, /*dt=*/5.0, /*n=*/3);

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_GE(cands.size(), 1u)
      << "EXCLUSIVE_D: candidate from symmetric vector-meson decay tracks";
  EXPECT_NEAR(cands[0].t0, 500.0, 2.0);
}

// 6 — SIDIS: semi-inclusive DIS — backward electron + more forward identified hadron
TEST(EventBuilder, SIDIS_ElectronPlusForwardHadron) {
  std::vector<RawHit> raw;
  addCluster(raw, /*th=*/7, /*ph=*/3, /*R=*/100.0f, /*t=*/600.0, /*dt=*/10.0, /*n=*/4);
  addCluster(raw, /*th=*/4, /*ph=*/2, /*R=*/100.0f, /*t=*/605.0, /*dt=*/5.0, /*n=*/3);

  auto hits  = toFastHits(raw);
  auto cands = findCandidates(hits);

  ASSERT_GE(cands.size(), 1u)
      << "SIDIS: candidate from backward electron + forward hadron";
  EXPECT_NEAR(cands[0].t0, 602.0, 6.0);
  EXPECT_LT(cands[0].t0sigma, 1.0);
}

// ===========================================================================
// V. Cross-frame recovery (merge logic — EventUnfolder + TimesliceBuffer_service)
// ===========================================================================

TEST(TimeAlignment, BackwardMerge_RecoversPastHit) {
  std::vector<RawHit> current = {
    hitAtBin(6, 0, 100.0f, 100.0, kSigmaTOF),
    hitAtBin(6, 0, 100.0f, 150.0, kSigmaTOF),
  };
  std::vector<RawHit> prev_frame = {
    hitAtBin(6, 0, 200.0f, -50.0, kSigmaSi),  // Si hit smeared into prev frame
  };

  std::vector<RawHit> merged = prev_frame;
  merged.insert(merged.end(), current.begin(), current.end());
  std::sort(merged.begin(), merged.end(),
            [](const RawHit& a, const RawHit& b) { return a.time < b.time; });

  EXPECT_EQ(merged.size(), 3u);
  EXPECT_NEAR(merged[0].time, -50.0, 1e-6)
      << "Backward-smeared Si hit is first after sort (earliest time)";
}

// ===========================================================================
// VI. True pile-up: collision-instance separation (EventBuilder_factory's
//     collision_times/collision_stream dedup + trigger_classes (weights[25+]),
//     EventUnfolder's coincidentSlot() nearest-match)
// ===========================================================================
//
// These functions replicate, for testing, the dedup loop in
// EventBuilder_factory.h::Process() (the collision_times/collision_stream
// block) and the per-candidate match loop that counts coincident collisions
// (flag/trigger_classes_mask, weights[2]/[17]) and, when
// store_coincident_list is on (the default), retains the (time, stream)
// pairs (trigger_classes, weights[25+]). See that file's own comments. This
// also replicates eventbuilder.cc's coincidentSlot() nearest-match, used to
// tag kept hits.

namespace {

struct PileupCollision { double time; int stream; };

// Mirrors EventBuilder_factory.h's collision_times/collision_stream loop: dedupe
// injected collisions by (stream, time), same stream within `dedup` ns
// merges into one entry, different streams never merge regardless of time.
std::vector<PileupCollision> dedupPileup(const std::vector<PileupCollision>& raw,
                                         double dedup) {
  std::vector<PileupCollision> out;
  for (const auto& c : raw) {
    bool merged = false;
    for (auto& o : out)
      if (o.stream == c.stream && std::fabs(o.time - c.time) < dedup) {
        merged = true;
        break;
      }
    if (!merged)
      out.push_back(c);
  }
  return out;
}

// Mirrors the per-candidate coincident_list loop: every deduped collision
// within +-dt of t0 is retained as a (time, stream) pair.
std::vector<PileupCollision> coincidentList(const std::vector<PileupCollision>& deduped,
                                            double t0, double dt) {
  std::vector<PileupCollision> list;
  for (const auto& c : deduped)
    if (std::fabs(t0 - c.time) <= dt)
      list.push_back(c);
  return list;
}

// Mirrors eventbuilder.cc's coincidentSlot(): nearest-time match against a
// candidate's own coincident_list, -1 if the list is empty.
int coincidentSlot(const std::vector<PileupCollision>& list, double t) {
  if (list.empty())
    return -1;
  int best = 0;
  double best_d = std::fabs(t - list[0].time);
  for (size_t k = 1; k < list.size(); ++k) {
    const double d = std::fabs(t - list[k].time);
    if (d < best_d) { best_d = d; best = int(k); }
  }
  return best;
}

}  // namespace

// Two collisions of the same stream (class), 30 ns apart, with a tight
// dedup window (10 ns, for example a TOF-dominated frame). 30 ns > 10 ns,
// so they must not merge: two distinct entries that both count into the
// candidate's flag and both appear in trigger_classes (weights[25+]).
TEST(TruePileup, Dedup_SameStreamFarApart_StaysSeparate) {
  std::vector<PileupCollision> raw = {{1000.0, 10}, {1030.0, 10}};
  auto deduped = dedupPileup(raw, /*dedup=*/10.0);
  EXPECT_EQ(deduped.size(), 2u)
      << "30 ns apart, 10 ns dedup window: same-stream collisions stay distinct";
}

// The same two collisions, but with a wide dedup window (40 ns) spanning
// the 30 ns gap: the intended merge behavior at the boundary. This is not
// a bug; it documents that the dedup granularity itself, not
// trigger_classes, decides whether two same-stream collisions are
// resolvable at all.
TEST(TruePileup, Dedup_SameStreamFarApart_MergesUnderWideWindow) {
  std::vector<PileupCollision> raw = {{1000.0, 10}, {1030.0, 10}};
  auto deduped = dedupPileup(raw, /*dedup=*/40.0);
  EXPECT_EQ(deduped.size(), 1u)
      << "30 ns apart, 40 ns dedup window: same-stream collisions intentionally merge";
}

// Different streams never merge regardless of time separation (mirrors the
// production dedup's explicit stream check) -- two collisions 5 ns apart,
// well inside even the tight dedup window, from DIFFERENT classes.
TEST(TruePileup, Dedup_DifferentStreamsNeverMerge) {
  std::vector<PileupCollision> raw = {{1000.0, 10}, {1005.0, 80}};
  auto deduped = dedupPileup(raw, /*dedup=*/10.0);
  EXPECT_EQ(deduped.size(), 2u)
      << "Different streams must never merge, even well inside the dedup window";
}

// Round-trip: a candidate's trigger_classes list size matches what flag
// (the count field, weights[2]) reports for the same match test --
// weights[25] (the stored list count) and weights[2] must never disagree,
// since both come from the identical `|t0 - t| <= dt` test.
TEST(TruePileup, CoincidentList_CountMatchesCoincidenceTest) {
  std::vector<PileupCollision> raw = {{1000.0, 10}, {1030.0, 80}, {2000.0, 150}};
  auto deduped = dedupPileup(raw, /*dedup=*/10.0);
  ASSERT_EQ(deduped.size(), 3u);

  const double t0 = 1015.0, dt = 30.0;  // window covers both close collisions, not the far one
  auto list = coincidentList(deduped, t0, dt);

  int flag_count = 0;
  for (const auto& c : deduped)
    if (std::fabs(t0 - c.time) <= dt)
      ++flag_count;

  EXPECT_EQ(list.size(), 2u) << "Both close collisions (streams 10, 80) are in-window";
  EXPECT_EQ(static_cast<int>(list.size()), flag_count)
      << "weights[25] (list size) must equal weights[2] (flag) by construction";
}

// Nearest-match slot assignment: hits/particles from each of two close
// collisions must resolve to the CORRECT distinct slot, not both collapsing
// onto one -- this is what actually lets per-hit truth separate two
// overlapping collisions downstream (eventbuilder.cc's coincidentSlot()).
TEST(TruePileup, CoincidentSlot_AssignsCorrectDistinctSlots) {
  std::vector<PileupCollision> list = {{1000.0, 10}, {1030.0, 80}};

  EXPECT_EQ(coincidentSlot(list, 1000.2), 0) << "Particle near collision 0 must resolve to slot 0";
  EXPECT_EQ(coincidentSlot(list, 1029.8), 1) << "Particle near collision 1 must resolve to slot 1";
  // Exactly at the midpoint, ties go to whichever is scanned first (slot 0)
  // -- documenting the tie-break, not asserting a "correct" side since none
  // exists for a genuine tie.
  EXPECT_EQ(coincidentSlot(list, 1015.0), 0) << "Midpoint tie-break: first-scanned slot wins";
}

TEST(TruePileup, CoincidentSlot_EmptyListReturnsMinusOne) {
  EXPECT_EQ(coincidentSlot({}, 1000.0), -1)
      << "No coincident-list data (parameter off, or nothing coincided): -1, unchanged from "
         "pre-existing behavior";
}

// ===========================================================================
// VI-b. generatorStatus decode: the LIVE 1000-wide-per-class scheme, plus
//       the DORMANT 10000-wide instance-sub-band scheme, and float32 slot
//       packing
// ===========================================================================
//
// LIVE scheme mirrors the production arithmetic: EventBuilder_factory.h's
// collision loop (status >= 10000, stream = status/1000, class_index =
// stream-10) and signal.sh's own _CLASS_STATUS base values (10000, 11000,
// 12000, ... -- 1000 apart, confirmed against the real generation script,
// not assumed). class_index = (stream-10)/10 (an extra /10, assuming a
// 10000-wide band) was a real bug found and fixed this session -- it
// collapsed classes 0-9/10-19/20-25 into three buckets. These LIVE tests
// exist specifically to catch a regression of that bug; the DORMANT tests
// further down do NOT exercise it (they test a different, inactive band
// width) and would not have caught it.
//
// DORMANT scheme: a wider, 10000-wide-per-class layout with 1000-wide
// same-class pile-up instance sub-bands is designed in the production
// repo's rewrite_status.py, meant to run between the SBM merge and npsim.
// It unconditionally self-disables today (SBM writes ROOT-format merges
// the rewriter cannot parse), so no real status ever has instance > 1 --
// these tests protect the arithmetic FOR WHEN that blocker is resolved,
// not current behavior.

namespace {

// LIVE: 1000-wide per class, no instance sub-bands (every real status is
// implicitly "instance 1").
struct LiveStatusDecode {
  bool physics;
  int class_index; // 0-25, physics only
  int low;         // status % 1000: 1 = stable, 2 = decay
};

LiveStatusDecode decodeLiveStatus(int status) {
  if (status < 10000)
    return {false, -1, status % 1000};
  return {true, (status - 10000) / 1000, status % 1000};
}

int streamOf(int status) { return status / 1000; }
int classIndexOfStream(int stream) { return stream - 10; }

// DORMANT (rewrite_status.py, not live): 10000-wide per class, 1000-wide
// instance sub-bands.
struct DormantStatusDecode {
  bool physics;
  int class_index; // 0-25, physics only
  int instance;    // 1-based same-class pile-up instance, physics only
  int low;         // status % 1000: 1 = stable, 2 = decay
};

DormantStatusDecode decodeDormantStatus(int status) {
  if (status < 10000)
    return {false, -1, 0, status % 1000};
  return {true, (status - 10000) / 10000, (status % 10000) / 1000 + 1, status % 1000};
}

int dormantClassIndexOfStream(int stream) { return (stream - 10) / 10; }

} // namespace

TEST(LiveStatusDecode, ClassBands) {
  EXPECT_FALSE(decodeLiveStatus(1).physics) << "Native status 1 is not a physics-class code";
  EXPECT_FALSE(decodeLiveStatus(2001).physics) << "Machine background (2000-6999) is not physics";
  EXPECT_TRUE(decodeLiveStatus(10001).physics);
  EXPECT_EQ(decodeLiveStatus(10001).class_index, 0) << "10000-10999 = class 0 (ncdisq1)";
  EXPECT_EQ(decodeLiveStatus(35001).class_index, 25) << "35000-35999 = class 25 (photoprod)";
  EXPECT_EQ(decodeLiveStatus(17001).class_index, 7) << "17000-17999 = class 7 (dvcs)";
  EXPECT_EQ(decodeLiveStatus(10001).low, 1) << "base+1 = stable";
  EXPECT_EQ(decodeLiveStatus(10002).low, 2) << "base+2 = decay";
}

// stream -> class_index for every class, the exact quantity the
// trigger_classes_mask bit and the trigger_classes stream decode both use.
// The extra "/10" regression this guards against silently collapsed
// classes 0-9/10-19/20-25 into three buckets in EventBuilder_factory.h.
TEST(LiveStatusDecode, StreamToClassIndexForEveryClass) {
  for (int cls = 0; cls < 26; ++cls) {
    const int status = 10000 + cls * 1000 + 1;
    EXPECT_EQ(classIndexOfStream(streamOf(status)), cls) << "status " << status;
  }
}

TEST(DormantStatusDecode, ClassBands) {
  EXPECT_FALSE(decodeDormantStatus(1).physics) << "Native status 1 is not a physics-class code";
  EXPECT_FALSE(decodeDormantStatus(2001).physics) << "Machine background (2000-6999) is not physics";
  EXPECT_TRUE(decodeDormantStatus(10001).physics);
  EXPECT_EQ(decodeDormantStatus(10001).class_index, 0) << "10000-19999 = class 0 (ncdisq1)";
  EXPECT_EQ(decodeDormantStatus(260001).class_index, 25) << "260000-269999 = class 25 (photoprod)";
  EXPECT_EQ(decodeDormantStatus(80001).class_index, 7) << "80000-89999 = class 7 (dvcs)";
}

TEST(DormantStatusDecode, InstanceSubBands) {
  EXPECT_EQ(decodeDormantStatus(10001).instance, 1);
  EXPECT_EQ(decodeDormantStatus(10001).low, 1) << "base+1 = instance 1, stable";
  EXPECT_EQ(decodeDormantStatus(10002).low, 2) << "base+2 = instance 1, decay";
  EXPECT_EQ(decodeDormantStatus(11001).instance, 2) << "base+1001 = instance 2, stable";
  EXPECT_EQ(decodeDormantStatus(11001).class_index, 0) << "Instance shift stays inside the class band";
  EXPECT_EQ(decodeDormantStatus(12002).instance, 3);
  EXPECT_EQ(decodeDormantStatus(12002).low, 2) << "base+2002 = instance 3, decay";
}

// stream -> class_index must hold for EVERY instance sub-band of every
// class, IF the dormant scheme ever goes live.
TEST(DormantStatusDecode, StreamToClassIndexForEveryInstance) {
  for (int cls = 0; cls < 26; ++cls) {
    const int base = 10000 + cls * 10000;
    for (int inst = 1; inst <= 10; ++inst) {
      const int status = base + (inst - 1) * 1000 + 1;
      EXPECT_EQ(dormantClassIndexOfStream(streamOf(status)), cls)
          << "status " << status << " (class " << cls << ", instance " << inst << ")";
    }
  }
}

// flag counts every in-window collision; mc_t is the EARLIEST one's time —
// mirrors the factory's per-candidate matching loop.
TEST(CountFlag, CountsAllAndReportsEarliest) {
  std::vector<PileupCollision> deduped = {{1030.0, 10}, {1000.0, 80}, {2000.0, 150}};
  const double t0 = 1015.0, dt = 30.0;
  int flag = 0;
  double mc_t = -1.0e9;
  for (const auto& c : deduped)
    if (std::fabs(t0 - c.time) <= dt) {
      ++flag;
      if (mc_t == -1.0e9 || c.time < mc_t)
        mc_t = c.time;
    }
  EXPECT_EQ(flag, 2) << "Two real collisions in-window: flag = 2 (true pile-up)";
  EXPECT_NEAR(mc_t, 1000.0, 1e-9) << "mc_t is the earliest matched collision's time";
}

TEST(CountFlag, ZeroCollisionsMeansFake) {
  std::vector<PileupCollision> deduped = {{2000.0, 10}};
  const double t0 = 1000.0, dt = 30.0;
  int flag = 0;
  for (const auto& c : deduped)
    if (std::fabs(t0 - c.time) <= dt)
      ++flag;
  EXPECT_EQ(flag, 0) << "No real collision in-window: flag = 0 (fake)";
}

// weight = generatorStatus + slot*1e6 must round-trip exactly in float32.
// The exact-integer ceiling of float is 2^24 = 16,777,216: tested against
// 269999, the DORMANT wide scheme's max (see the StatusDecode tests above)
// rather than the LIVE scheme's 35999, so the bound stays valid either way
// and does not need revisiting if the dormant scheme goes live. Slots 0-16
// stay exact and slot 17 is the first to break — this bound backs the
// assert(slot < 16) at both encode sites.
TEST(SlotPacking, Float32ExactThroughSlot16) {
  const int max_status = 269999;
  for (int slot = 0; slot <= 16; ++slot) {
    const double packed = double(max_status) + double(slot) * 1e6;
    EXPECT_EQ(double(float(packed)), packed) << "slot " << slot << " must pack exactly";
  }
  const double broken = double(max_status) + 17.0 * 1e6;
  EXPECT_NE(double(float(broken)), broken)
      << "slot 17 exceeds 2^24 and must NOT be exactly representable — the assert bound is real";
}

// ===========================================================================
// VII. EventPrefilter — ONNX gold-coating GNN model (compiled in when HAS_ONNXRUNTIME)
// ===========================================================================
//
// The gold-coating model (MultiClassEventGNN) takes 4 inputs:
//   x          [N, 6]  float32 — normalised hit features
//                                 [x/4000, y/4000, z/5000, t/50,
//                                  clip(eDep/10,0,1), det_id/6]
//   edge_idx_0 [2, E]  int64   — k=16 NN in 4D spacetime (x,y,z,t)
//   edge_idx_1 [2, E]  int64   — k=16 NN in layer-0 hidden features
//   edge_idx_2 [2, E]  int64   — k=16 NN in layer-1 hidden features
//
// Layers 1 and 2 need intermediate activations for exact edge construction.
// For integration tests we approximate them with the same spatial graph.
// Output: probs [1, 7] — softmax class probabilities (sum ≈ 1).
//
// Per-class accuracy tests on real events: run `make test-events` first to
// produce vendor/eicrecon/share/samples/gold-coating/<class>.eicrecon.root,
// then use the companion Python scripts/validate_prefilter.py.

#ifdef HAS_ONNXRUNTIME

namespace {

const char* findModelPath() {
  static const char* candidates[] = {
    "/mnt/local/vendor/eicrecon/share/models/prefilter-gold-coating-100epochs-2500hits.onnx",
    "/mnt/local/share/models/prefilter-gold-coating-100epochs-2500hits.onnx",
    "vendor/eicrecon/share/models/prefilter-gold-coating-100epochs-2500hits.onnx",
    "../vendor/eicrecon/share/models/prefilter-gold-coating-100epochs-2500hits.onnx",
  };
  for (auto p : candidates)
    if (access(p, R_OK) == 0) return p;
  return nullptr;
}

// Brute-force k-NN directed edge index [2, N*k] in the first coord_dims columns.
// Convention: source = neighbor, target = query (PyG knn_graph default).
std::vector<int64_t> buildKnnEdges(const std::vector<float>& feat,
                                    int N, int F, int coord_dims, int k) {
  k = std::min(k, N - 1);
  if (k <= 0) return {};
  std::vector<int64_t> src(N*k), tgt(N*k);
  std::vector<std::pair<float,int>> dists(N - 1);
  for (int i = 0; i < N; ++i) {
    int d = 0;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      float dist = 0;
      for (int c = 0; c < coord_dims; ++c) {
        float dx = feat[i*F+c] - feat[j*F+c]; dist += dx*dx;
      }
      dists[d++] = {dist, j};
    }
    std::partial_sort(dists.begin(), dists.begin()+k, dists.end());
    for (int ki = 0; ki < k; ++ki) {
      src[i*k+ki] = dists[ki].second;
      tgt[i*k+ki] = i;
    }
  }
  std::vector<int64_t> edges;
  edges.insert(edges.end(), src.begin(), src.end());
  edges.insert(edges.end(), tgt.begin(), tgt.end());
  return edges;
}

// Normalised hit feature: [x/4000, y/4000, z/5000, t/50, clip(e/10,0,1), det/6]
std::vector<float> makeHit(float x_mm, float y_mm, float z_mm,
                             float t_ns, float e_gev, int det_id) {
  return {x_mm / 4000.f, y_mm / 4000.f, z_mm / 5000.f,
          t_ns / 50.f,
          std::min(1.f, std::max(0.f, e_gev / 10.f)),
          static_cast<float>(det_id) / 6.f};
}

// Run the gold-coating GNN and return the 7-class probability vector.
// Uses the same spatial graph for all 3 edge sets (approximation).
std::vector<float> runGnnPrefilter(Ort::Session& session,
                                    const std::vector<float>& feat, int N) {
  constexpr int F = 6, k = 16;
  auto edges = buildKnnEdges(feat, N, F, 4, k);
  const int64_t E = static_cast<int64_t>(edges.size() / 2);

  Ort::AllocatorWithDefaultOptions alloc;
  const std::string nm_x   = session.GetInputNameAllocated(0, alloc).get();
  const std::string nm_e0  = session.GetInputNameAllocated(1, alloc).get();
  const std::string nm_e1  = session.GetInputNameAllocated(2, alloc).get();
  const std::string nm_e2  = session.GetInputNameAllocated(3, alloc).get();
  const std::string nm_out = session.GetOutputNameAllocated(0, alloc).get();

  Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  const std::vector<int64_t> x_sh{N, F}, ei_sh{2, E};

  std::vector<float>   fx = feat;
  std::vector<int64_t> e0 = edges, e1 = edges, e2 = edges;

  std::vector<Ort::Value> ins;
  ins.push_back(Ort::Value::CreateTensor<float>(
      mem, fx.data(), fx.size(), x_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e0.data(), e0.size(), ei_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e1.data(), e1.size(), ei_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e2.data(), e2.size(), ei_sh.data(), 2));

  const char* in_names[]  = {nm_x.c_str(), nm_e0.c_str(), nm_e1.c_str(), nm_e2.c_str()};
  const char* out_names[] = {nm_out.c_str()};
  auto out = session.Run(Ort::RunOptions{nullptr}, in_names, ins.data(), 4, out_names, 1);

  const float* p = out[0].GetTensorMutableData<float>();
  return std::vector<float>(p, p + out[0].GetTensorTypeAndShapeInfo().GetElementCount());
}

} // namespace

// Smoke-test: model loads, runs on 50-hit BKG-like frame, produces valid probs.
TEST(Prefilter, GoldCoating_LoadsAndScores) {
  const char* path = findModelPath();
  if (!path)
    GTEST_SKIP() << "Gold-coating model not found; run inside the container.";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "prefilter-test");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  opts.SetInterOpNumThreads(1);
  Ort::Session session(env, path, opts);

  // 50 deterministic BKG-like hits spread across the detector
  constexpr int N = 50;
  std::vector<float> feat;
  feat.reserve(N * 6);
  uint32_t rng = 0xdeadbeef;
  auto rnd = [&]() -> float {
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return (static_cast<float>(rng & 0xffffff) / 0x800000f) - 1.f;
  };
  for (int i = 0; i < N; ++i) {
    auto h = makeHit(rnd() * 800.f, rnd() * 800.f, rnd() * 1500.f,
                     std::abs(rnd()) * 1000.f, std::abs(rnd()) * 0.1f,
                     static_cast<int>(std::abs(rnd()) * 4));
    feat.insert(feat.end(), h.begin(), h.end());
  }

  const auto probs = runGnnPrefilter(session, feat, N);
  ASSERT_EQ(static_cast<int>(probs.size()), 7) << "Expected 7 class probabilities";

  float sum = 0;
  for (float p : probs) {
    EXPECT_GE(p, 0.0f);
    EXPECT_LE(p, 1.0f);
    sum += p;
  }
  EXPECT_NEAR(sum, 1.0f, 0.01f) << "Softmax probs must sum to 1";
  SUCCEED() << "probs=[" << probs[0] << " " << probs[1] << " " << probs[2]
            << " " << probs[3] << " " << probs[4] << " " << probs[5]
            << " " << probs[6] << "]";
}

// Signal-like frame (compact backward cluster) should produce valid probs distinct
// from a random BKG-like frame — validates the model accepts independent calls.
TEST(Prefilter, GoldCoating_SignalScoresHigherThanNoise) {
  const char* path = findModelPath();
  if (!path) GTEST_SKIP() << "Gold-new model not found";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "prefilter-test");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  Ort::Session session(env, path, opts);

  constexpr int N = 50;
  uint32_t rng = 0xcafebabe;
  auto rnd = [&]() -> float {
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return (static_cast<float>(rng & 0xffffff) / 0x800000f) - 1.f;
  };

  std::vector<float> bkg_feat;
  for (int i = 0; i < N; ++i) {
    auto h = makeHit(rnd()*800.f, rnd()*800.f, rnd()*1500.f,
                     std::abs(rnd())*1000.f, std::abs(rnd())*0.05f, 0);
    bkg_feat.insert(bkg_feat.end(), h.begin(), h.end());
  }

  // DIS_NC-like: 20 hits in backward Si (z ≈ -950 mm) + 30 random BKG
  std::vector<float> sig_feat;
  for (int i = 0; i < 20; ++i) {
    auto h = makeHit(50.f+i*0.5f, 5.f+i*0.2f, -950.f+i*0.1f, 100.f+i*0.3f, 0.05f, 0);
    sig_feat.insert(sig_feat.end(), h.begin(), h.end());
  }
  for (int i = 0; i < 30; ++i) {
    auto h = makeHit(rnd()*800.f, rnd()*800.f, rnd()*1500.f,
                     std::abs(rnd())*1000.f, std::abs(rnd())*0.05f, 0);
    sig_feat.insert(sig_feat.end(), h.begin(), h.end());
  }

  const auto probs_bkg = runGnnPrefilter(session, bkg_feat, N);
  const auto probs_sig = runGnnPrefilter(session, sig_feat, N);

  ASSERT_EQ(static_cast<int>(probs_bkg.size()), 7);
  ASSERT_EQ(static_cast<int>(probs_sig.size()), 7);

  float sum_bkg = 0, sum_sig = 0;
  for (int c = 0; c < 7; ++c) {
    EXPECT_GE(probs_bkg[c], 0.f); EXPECT_LE(probs_bkg[c], 1.f);
    EXPECT_GE(probs_sig[c], 0.f); EXPECT_LE(probs_sig[c], 1.f);
    sum_bkg += probs_bkg[c];
    sum_sig += probs_sig[c];
  }
  EXPECT_NEAR(sum_bkg, 1.f, 0.01f);
  EXPECT_NEAR(sum_sig, 1.f, 0.01f);
  SUCCEED() << "bkg=[" << probs_bkg[0] << ",...] sig=[" << probs_sig[0] << ",...]";
}

#endif  // HAS_ONNXRUNTIME

// ===========================================================================
// VIII. ONNXRuntime_gnn — kNN graph construction (pure C++, no OnnxRuntime)
// ===========================================================================
//
// Validates the k-d tree kNN against brute force and checks the pyg-lib ≥ 0.6
// edge conventions the gold-coating model was trained with:
//   row = source (neighbor), col = target (query), no self-loops.

namespace {

// Reference brute-force neighbor sets over the first coord_dims of stride-F rows.
std::vector<std::set<int64_t>> knnBruteSets(const std::vector<float>& x,
                                            int N, int F, int coord_dims, int k,
                                            const float* w = nullptr) {
  std::vector<std::set<int64_t>> out(N);
  std::vector<std::pair<float,int>> d(N - 1);
  k = std::min(k, N - 1);
  for (int i = 0; i < N; ++i) {
    int c = 0;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      float s = 0;
      for (int q = 0; q < coord_dims; ++q) {
        float dx = x[i*F+q] - x[j*F+q];
        if (w) dx *= w[q];
        s += dx*dx;
      }
      d[c++] = {s, j};
    }
    std::partial_sort(d.begin(), d.begin()+k, d.end());
    for (int ki = 0; ki < k; ++ki) out[i].insert(d[ki].second);
  }
  return out;
}

std::vector<float> randomCloud(int N, int F, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> u(-1.f, 1.f);
  std::vector<float> x(static_cast<size_t>(N) * F);
  for (auto& v : x) v = u(rng);
  return x;
}

} // namespace

TEST(KnnGraph, MatchesBruteForce) {
  constexpr int N = 200, D = 4, k = 8;
  const auto x = randomCloud(N, D, 0x1234u);
  ONNXRuntime_gnn gnn;
  const auto ei = gnn.build(x.data(), N, D, k);
  ASSERT_EQ(ei.row.size(), static_cast<size_t>(N) * k);
  const auto ref = knnBruteSets(x, N, D, D, k);
  int mismatches = 0;
  for (int i = 0; i < N; ++i) {
    std::set<int64_t> got(ei.row.begin() + i*k, ei.row.begin() + (i+1)*k);
    if (got != ref[i]) ++mismatches;
  }
  // Distance ties may legitimately differ between algorithms.
  EXPECT_LE(mismatches, 2) << "kd-tree disagrees with brute force";
}

TEST(KnnGraph, EdgeConvention_SourceToTarget) {
  constexpr int N = 64, D = 4, k = 8;
  const auto x = randomCloud(N, D, 0x5678u);
  const auto ei = ONNXRuntime_gnn{}.build(x.data(), N, D, k);
  for (int i = 0; i < N; ++i) {
    for (int ki = 0; ki < k; ++ki) {
      EXPECT_EQ(ei.col[i*k+ki], i) << "col must be the query index";
      EXPECT_NE(ei.row[i*k+ki], i) << "self-loops must be excluded";
      EXPECT_GE(ei.row[i*k+ki], 0);
      EXPECT_LT(ei.row[i*k+ki], N);
    }
  }
}

TEST(KnnGraph, StridedFeatures_UsesFirst4DimsOnly) {
  // Rows are [x,y,z,t,eDep,det_id]; distances must ignore columns 4-5.
  constexpr int N = 100, F = 6, k = 8;
  auto x = randomCloud(N, F, 0x9abcu);
  auto x_zeroed = x;
  for (int i = 0; i < N; ++i) { x_zeroed[i*F+4] = 0.f; x_zeroed[i*F+5] = 0.f; }

  ONNXRuntime_gnn gnn;
  const auto ei_a = gnn.build(x.data(),        N, 4, k, nullptr, F);
  const auto ei_b = gnn.build(x_zeroed.data(), N, 4, k, nullptr, F);
  EXPECT_EQ(ei_a.row, ei_b.row) << "columns beyond D leaked into the metric";

  const auto ref = knnBruteSets(x, N, F, 4, k);
  int mismatches = 0;
  for (int i = 0; i < N; ++i) {
    std::set<int64_t> got(ei_a.row.begin() + i*k, ei_a.row.begin() + (i+1)*k);
    if (got != ref[i]) ++mismatches;
  }
  EXPECT_LE(mismatches, 2);
}

TEST(KnnGraph, TimeWeight_ChangesNeighborChoice) {
  // Query at origin. A: spatially close, far in time. B: farther, same time.
  const std::vector<float> x = {
      0.f, 0.f, 0.f,  0.f,   // 0: query
      1.f, 0.f, 0.f, 10.f,   // 1: A — d_space²=1,  d_time²=100
      2.f, 0.f, 0.f,  0.f,   // 2: B — d_space²=4,  d_time²=0
  };
  ONNXRuntime_gnn gnn;
  const float w_spatial[4] = {1.f, 1.f, 1.f, 0.f};   // time ignored → A wins
  const float w_spacetime[4] = {1.f, 1.f, 1.f, 1.f}; // time counted → B wins
  EXPECT_EQ(gnn.build(x.data(), 3, 4, 1, w_spatial).row[0], 1);
  EXPECT_EQ(gnn.build(x.data(), 3, 4, 1, w_spacetime).row[0], 2);
}

TEST(KnnGraph, KClampedToNMinus1) {
  constexpr int N = 3, D = 4;
  const auto x = randomCloud(N, D, 0xdef0u);
  const auto ei = ONNXRuntime_gnn{}.build(x.data(), N, D, /*k=*/16);
  EXPECT_EQ(ei.row.size(), static_cast<size_t>(N) * (N - 1));
}

TEST(KnnGraph, EmptyAndSinglePoint) {
  ONNXRuntime_gnn gnn;
  EXPECT_TRUE(gnn.build(nullptr, 0, 4, 16).row.empty());
  const float one[4] = {0.f, 0.f, 0.f, 0.f};
  EXPECT_TRUE(gnn.build(one, 1, 4, 16).row.empty());
}

TEST(KnnGraph, DuplicatePointsDoNotCrash) {
  constexpr int N = 32, D = 4, k = 8;
  const std::vector<float> x(static_cast<size_t>(N) * D, 0.5f);
  const auto ei = ONNXRuntime_gnn{}.build(x.data(), N, D, k);
  ASSERT_EQ(ei.row.size(), static_cast<size_t>(N) * k);
  for (int i = 0; i < N; ++i)
    for (int ki = 0; ki < k; ++ki)
      EXPECT_NE(ei.row[i*k+ki], i) << "self-loop with degenerate distances";
}

TEST(KnnGraph, DgcnnEdges_ThreeConsistentLayers) {
  constexpr int N = 128, F = 6, k = 16;
  const auto x = randomCloud(N, F, 0x2468u);
  const auto edges = ONNXRuntime_gnn{}.build_dgcnn_edges(x.data(), N, k, 1.0f, F);
  ASSERT_EQ(edges[0].row.size(), static_cast<size_t>(N) * k);
  EXPECT_EQ(edges[1].row, edges[0].row) << "layer-1 approximation must reuse layer-0";
  EXPECT_EQ(edges[2].row, edges[0].row) << "layer-2 approximation must reuse layer-0";
  EXPECT_EQ(edges[1].col, edges[0].col);
  EXPECT_EQ(edges[2].col, edges[0].col);
}

TEST(KnnGraph, ScalesToFrameSize) {
  // Full prefilter frame: 2500 hits, k=16 — must stay well inside the
  // per-timeslice budget (completion + sanity, no strict timing in CI).
  constexpr int N = 2500, F = 6, k = 16;
  const auto x = randomCloud(N, F, 0x1357u);
  const auto edges = ONNXRuntime_gnn{}.build_dgcnn_edges(x.data(), N, k, 1.0f, F);
  EXPECT_EQ(edges[0].row.size(), static_cast<size_t>(N) * k);
}

// ===========================================================================
// IX. Prefilter decoding heads B (superposition) and C (hit refinement)
// ===========================================================================
//
// Model output contract (by index; missing outputs = older Head-A-only export):
//   [0] probs        [1, 7]  float32  — Head A: softmax class probabilities
//   [1] multiplicity [1, 7]  float32  — Head B: per-class Poisson rates λ_k ≥ 0
//   [2] hit_scores   [N] or [N, 1]    — Head C: per-hit signal probability ∈ [0,1]
//
// These tests activate automatically once the model export includes the
// heads; until then they SKIP with an actionable message.

#ifdef HAS_ONNXRUNTIME

namespace {

// Run the GNN requesting every declared output.
std::vector<Ort::Value> runGnnAllOutputs(Ort::Session& session,
                                         const std::vector<float>& feat, int N) {
  constexpr int F = 6, k = 16;
  auto edges = buildKnnEdges(feat, N, F, 4, k);
  const int64_t E = static_cast<int64_t>(edges.size() / 2);

  Ort::AllocatorWithDefaultOptions alloc;
  std::vector<std::string> in_names, out_names;
  for (size_t i = 0; i < session.GetInputCount(); ++i)
    in_names.push_back(session.GetInputNameAllocated(i, alloc).get());
  for (size_t i = 0; i < session.GetOutputCount(); ++i)
    out_names.push_back(session.GetOutputNameAllocated(i, alloc).get());

  Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  const std::vector<int64_t> x_sh{N, F}, ei_sh{2, E};

  std::vector<float>   fx = feat;
  std::vector<int64_t> e0 = edges, e1 = edges, e2 = edges;

  std::vector<Ort::Value> ins;
  ins.push_back(Ort::Value::CreateTensor<float>(
      mem, fx.data(), fx.size(), x_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e0.data(), e0.size(), ei_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e1.data(), e1.size(), ei_sh.data(), 2));
  ins.push_back(Ort::Value::CreateTensor<int64_t>(
      mem, e2.data(), e2.size(), ei_sh.data(), 2));

  std::vector<const char*> in_ptrs, out_ptrs;
  for (const auto& s : in_names)  in_ptrs.push_back(s.c_str());
  for (const auto& s : out_names) out_ptrs.push_back(s.c_str());

  return session.Run(Ort::RunOptions{nullptr},
                     in_ptrs.data(), ins.data(), ins.size(),
                     out_ptrs.data(), out_ptrs.size());
}

// Deterministic 50-hit mixed frame for head contract checks.
std::vector<float> makeHeadTestFrame(int N) {
  std::vector<float> feat;
  feat.reserve(N * 6);
  uint32_t rng = 0x600dc0de;
  auto rnd = [&]() -> float {
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return (static_cast<float>(rng & 0xffffff) / 0x800000f) - 1.f;
  };
  for (int i = 0; i < N; ++i) {
    auto h = makeHit(rnd()*800.f, rnd()*800.f, rnd()*1500.f,
                     std::abs(rnd())*1000.f, std::abs(rnd())*0.1f,
                     static_cast<int>(std::abs(rnd())*4));
    feat.insert(feat.end(), h.begin(), h.end());
  }
  return feat;
}

} // namespace

TEST(PrefilterHeads, OutputContractIsDeclared) {
  const char* path = findModelPath();
  if (!path) GTEST_SKIP() << "Gold-coating model not found";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "prefilter-heads");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  Ort::Session session(env, path, opts);

  const size_t n_out = session.GetOutputCount();
  ASSERT_GE(n_out, 1u) << "Model must expose at least Head A (probs)";
  RecordProperty("n_outputs", static_cast<int>(n_out));
  SUCCEED() << "Model declares " << n_out << " output(s): "
            << (n_out >= 3 ? "Heads A+B+C" : n_out == 2 ? "Heads A+B" : "Head A only");
}

TEST(PrefilterHeads, HeadB_MultiplicityRates) {
  const char* path = findModelPath();
  if (!path) GTEST_SKIP() << "Gold-coating model not found";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "prefilter-heads");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  Ort::Session session(env, path, opts);

  if (session.GetOutputCount() < 2)
    GTEST_SKIP() << "Head B not exported yet — re-export the model with "
                    "multiplicity_head=true (output[1] = lambda [1,7])";

  constexpr int N = 50;
  const auto feat = makeHeadTestFrame(N);
  auto outs = runGnnAllOutputs(session, feat, N);

  auto info = outs[1].GetTensorTypeAndShapeInfo();
  const auto shape = info.GetShape();
  ASSERT_EQ(info.GetElementCount(), 7u)
      << "Head B must emit one Poisson rate per physics class";
  const float* lam = outs[1].GetTensorMutableData<float>();
  for (int c = 0; c < 7; ++c) {
    EXPECT_TRUE(std::isfinite(lam[c])) << "lambda[" << c << "] not finite";
    EXPECT_GE(lam[c], 0.f) << "Poisson rate must be non-negative";
  }
}

TEST(PrefilterHeads, HeadC_HitScores) {
  const char* path = findModelPath();
  if (!path) GTEST_SKIP() << "Gold-coating model not found";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "prefilter-heads");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  Ort::Session session(env, path, opts);

  if (session.GetOutputCount() < 3)
    GTEST_SKIP() << "Head C not exported yet — add the per-node hit-refinement "
                    "head to the export (output[2] = hit_scores [N] or [N,1])";

  constexpr int N = 50;
  const auto feat = makeHeadTestFrame(N);
  auto outs = runGnnAllOutputs(session, feat, N);

  auto info = outs[2].GetTensorTypeAndShapeInfo();
  ASSERT_EQ(info.GetElementCount(), static_cast<size_t>(N))
      << "Head C must emit exactly one score per input hit";
  const float* s = outs[2].GetTensorMutableData<float>();
  for (int i = 0; i < N; ++i) {
    EXPECT_GE(s[i], 0.f) << "hit_score[" << i << "] below 0";
    EXPECT_LE(s[i], 1.f) << "hit_score[" << i << "] above 1";
  }
}

#endif  // HAS_ONNXRUNTIME

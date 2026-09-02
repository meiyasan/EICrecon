// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// TriggerBenchmark_service
// ------------------------
// Process-wide accumulator behind the eventbuilder_trigger benchmark. Both
// taps (FrameBenchmark_processor at Timeslice level, EventBenchmark_processor
// at PhysicsEvent level) fold their counts in here; this service owns the
// table rendering, the periodic progress report, and the CSV sidecar.
//
// Why a service and not one processor: the two levels see different things.
// Only the Timeslice tap sees frames that produced ZERO candidates, and those
// frames carry injected collisions the trigger missed entirely -- exactly the
// events that belong in an efficiency denominator. Only the PhysicsEvent tap
// can pull the per-candidate ACTS/calo reconstruction. Neither level can
// reach the other, so they meet here.
//
// Ordering: both taps call add*() from their ProcessSequential callback with
// ordering enabled, so JANA serializes them and no lock is needed for
// correctness of the counts themselves. The mutex guards only the case of the
// two taps (different levels, different arrows) touching the store at once.

#pragma once

#include <chrono>
#include <cstdint>
#include <map>
#include <mutex>
#include <set>
#include <string>

#include <utility>
#include <tuple>
#include <vector>

#include <JANA/JApplicationFwd.h>
#include <JANA/JService.h>
#include <spdlog/fwd.h>

#include "PdfReport.h"
#include "ResolutionHists.h"
#include "factories/eventbuilder/TruthClassLabels.h"

namespace eicrecon::eb {

/// Per-physics-class counters. Efficiency is per class; PURITY IS NOT --
/// fakes carry no class (mask 0, empty trigger_classes tail), so there is no
/// per-class fake population to divide by. Purity lives on the totals only.
struct ClassRow {
  std::uint64_t injected     = 0; ///< collisions injected into the stream
  std::uint64_t found        = 0; ///< distinct injected collisions a real candidate recovered
  std::uint64_t candidates   = 0; ///< candidates whose trigger_classes_mask carries this class
  std::uint64_t real_cands   = 0; ///< real candidates carrying it (real_cands - found = duplicates)
  std::uint64_t trk_matched  = 0; ///< STAGE 3: majority-matched in-acceptance charged MC
  std::uint64_t trk_expected = 0;
  std::uint64_t cal_matched  = 0; ///< STAGE 3: matched neutral truth groups
  std::uint64_t cal_expected = 0;
  std::uint64_t labeled_hits = 0;
  /// STAGE 1 per class: of the REAL candidates carrying this class, how many
  /// passed each fast-primitive term. x_ok / real_cands is that stage's OWN
  /// acceptance -- deliberately conditional on a candidate existing, because
  /// the primitives are per-candidate quantities: a collision the clustering
  /// never recovered has no n_tracklets and no E_calo, so it cannot be said to
  /// have failed the calo term. This is what a threshold is tuned against.
  /// It does mean the row is NOT a left-to-right chain; the composed figure is
  /// the "chain" footer line, from found_trig below.
  std::uint64_t trk_ok = 0, cal_ok = 0, trig_ok = 0;
  /// Distinct injected collisions recovered by a candidate that ALSO passed
  /// each term. A collision counts for a term if ANY candidate that recovered
  /// it passed that term (see addFrame): the terms are OR-ed across duplicate
  /// claimants, not read off whichever candidate happened to be seen first.
  ///
  /// found_x / x_ok is the per-class purity: the shortfall is duplicate
  /// candidates claiming the same collision, which IS class-attributable
  /// (unlike a fake, which carries no class at all). found_trig / injected is
  /// the trigger column: TIME && (TRACK || CALO), end to end.
  std::uint64_t found_trk = 0, found_cal = 0, found_trig = 0;
  /// Per-class purity, counted in RECONSTRUCTED OBJECTS on both sides (the
  /// efficiency columns count MC particles -- mixing the two would divide
  /// tracks by particles). *_any counts objects whose majority owner is this
  /// class at all; *_good the subset whose owner is also credited to the
  /// candidate and in acceptance. The shortfall is contamination within the
  /// class: secondaries, out-of-acceptance, or another collision of the same
  /// class inside the gate.
  /// Wall time spent in factories for candidates of this class, and how many
  /// candidates that covers. Only filled when profiling is on.
  double        proc_seconds = 0.0;
  std::uint64_t proc_cands   = 0;
  std::uint64_t trk_reco_any = 0, trk_reco_good = 0;
  std::uint64_t cal_reco_any = 0, cal_reco_good = 0; ///< tracker RawHitLink hits carrying this class band
};

/// Stream-wide counters that are not per class.
struct Totals {
  std::uint64_t frames     = 0;
  std::uint64_t candidates = 0;
  std::uint64_t real       = 0; ///< weights[FLAG] > 0
  std::uint64_t fake       = 0; ///< weights[FLAG] == 0

  // STAGE 1 -- the fast primitives the trigger definition actually specifies
  // (n_tracklets / E_calo thresholded), NOT ACTS or CaloIsland output.
  std::uint64_t trk_ok = 0, trk_fake = 0;
  std::uint64_t cal_ok = 0, cal_fake = 0;
  std::uint64_t trig_ok = 0, trig_fake = 0; ///< TRK || CAL

  // STAGE 2 -- GNN prefilter, one decision per frame.
  std::uint64_t gnn_frames = 0, gnn_passed = 0;
  /// Frames with a GNN score and ZERO real candidates: the false-accept
  /// baseline that stands in for the retired dedicated bkg class.
  std::uint64_t gnn_fake_frames = 0, gnn_fake_accepted = 0;

  // STAGE 3 -- ACTS tracks and calo clusters, aggregate.
  // *_matched / *_expected is EFFICIENCY (truth found).
  // *_reco_matched / *_reco is PURITY (reconstructed objects that are real):
  // the denominator counts every reconstructed object, the numerator only
  // those carrying a majority match, so the shortfall is the ghost rate.
  std::uint64_t trk_matched = 0, trk_expected = 0;
  std::uint64_t cal_matched = 0, cal_expected = 0;
  std::uint64_t trk_reco = 0, trk_reco_matched = 0;
  std::uint64_t cal_reco = 0, cal_reco_matched = 0;

  // Collision recovery.
  std::uint64_t collisions_injected = 0;
  std::uint64_t collisions_found    = 0;
  /// Stream-wide counterparts of ClassRow::found_x: distinct collisions
  /// recovered by a candidate that also passed each term. found_trig over
  /// collisions_injected is the end-to-end TIME && (TRACK || CALO) figure the
  /// "chain" footer line reports -- the one number the stage columns, being
  /// per-stage, cannot be multiplied out to give.
  std::uint64_t found_trk = 0, found_cal = 0, found_trig = 0;

  /// true once any candidate's stored weights[N_EXPECTED] was used as the
  /// tracking denominator, so the table can say the number is the file's own
  /// acceptance rather than the benchmark's recomputation.
  bool used_stored_nexp = false;
  std::uint64_t stored_nexp = 0;
  /// Child PhysicsEvents the STAGE 3 tap actually saw. The frame tap sees
  /// every frame and every candidate, so `candidates` is deterministic; this
  /// counts what reached the other level. The two should agree -- one child
  /// per candidate -- and a shortfall means the run ended before the child
  /// stream drained, which is exactly what makes BLIND and the STAGE 3
  /// columns wobble between otherwise identical runs.
  std::uint64_t children = 0;
  /// The trigger thresholds the file was produced with (weights[12]/[14]).
  /// A calo threshold of 0 means the calo term is DISABLED: cal_pass is then
  /// true for every candidate, so the calo and trigger columns collapse onto
  /// the time column and report the same measurement three times.
  double trk_threshold = -1.0;
  double cal_threshold = -1.0;
  bool saw_stage3 = false; ///< false => STAGE 3 columns render "--", not 0
  bool saw_gnn    = false;
};

class TriggerBenchmark_service : public JService {
public:
  explicit TriggerBenchmark_service(JApplication* app) : m_app(app) {}

  void acquire_services(JServiceLocator* locator) override;

  // -- configuration (read once, in acquire_services) ------------------------
  int    reportEvery() const { return m_report_every; }
  bool   stage3()      const { return m_stage3; }
  double frameNs()     const { return m_frame_ns; }
  double caloDr()      const { return m_calo_dr; }
  /// N-sigma multiplier the trigger folded into weights[T0SIGMA].
  double nsigmaWindow() const { return m_nsigma_window; }

  /// Residual histograms, shared: the tap fills them, the report plots them.
  ResolutionHists& res() { return m_res; }
  /// Records EVERY source read, not just the first: a run over many files
  /// otherwise reports one of them as if it were the whole input.
  void addInputFile(const std::string& f);

  /// Each tap registers at construction; the final table is printed when the
  /// last of them finishes, so Finish() ordering between taps does not matter.
  void registerTap();
  void tapFinished();

  // -- accumulation ----------------------------------------------------------
  /// One recovered collision, with the stage flags of the candidate that
  /// recovered it, so distinct-vs-duplicate can be resolved per stage.
  struct FoundKey {
    std::uint64_t frame = 0;
    int           cls   = -1;
    double        time  = 0.0;
    bool          trk = false, cal = false;
  };

  void addFrame(const Totals& delta, const std::map<int, ClassRow>& per_class,
                const std::vector<FoundKey>& found_keys);
  void addEvent(const Totals& delta, const std::map<int, ClassRow>& per_class);

  /// Frames that produced no child event and so never reached the tap. Their
  /// injected collisions cannot be counted; reporting the count keeps the
  /// efficiency figure honest about what it could not see.
  void addBlindFrames(std::uint64_t n);

  /// Wall-clock timing of the run, so the table can report throughput per
  /// frame and per candidate alongside the physics numbers. Started on the
  /// first event seen and advanced on every one after it, so it measures the
  /// processing window rather than including JANA start-up and geometry load.
  void tickClock();

  /// Per-factory wall time, accumulated from JANA's OWN call graph
  /// (JCallGraphRecorder timestamps every factory call). Off unless
  /// eventbuilder:benchmark:profile=1, because recording the call stack costs
  /// throughput -- JANA's own parameter description says so.
  bool profiling() const { return m_profile; }
  void addFactoryTime(const std::string& factory, double seconds);
  /// Slowest factories first: {factory, total seconds, call count}.
  std::vector<std::tuple<std::string, double, std::uint64_t>> factoryProfile() const;
  /// Same, but per physics class: {class, total seconds, candidates}.
  std::vector<std::tuple<std::string, double, std::uint64_t>> classProfile() const;

  /// Distinct parent frames the PhysicsEvent tap reached, reported once at
  /// its Finish(). Blind frames are the frame tap's count minus this.
  void addChildFrames(std::uint64_t n);

  /// Standalone mode = replaying a written *.eicrecon.root, where no Timeslice
  /// parent exists. Denominators degrade (zero-candidate frames are invisible);
  /// the header says so rather than quietly reporting a different number.
  void markStandalone() { m_standalone = true; }

  /// Called by the frame tap; emits a progress table every reportEvery frames.
  void maybeReport();

  void report(bool final_report);

  /// The statistical table as CELLS rather than pre-formatted lines, so the
  /// PDF can align and colour it instead of dumping monospace text.
  /// Returns {header, rows}; an empty row is a spacer.
  std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>> tableCells() const;

  /// Footer facts (FAR, recovery, blind frames) as ready-to-print lines.
  std::vector<std::string> tableFooter() const;

private:
  std::string renderTable(bool final_report) const;
  void        writeCsv() const;

  JApplication*                   m_app{};
  std::shared_ptr<spdlog::logger> m_log;

  mutable std::mutex      m_mutex;
  Totals                  m_totals;
  std::map<int, ClassRow> m_classes;
  /// Bits in m_found's value: the stage terms a collision has already been
  /// credited with, so a second candidate claiming it cannot double-count.
  static constexpr unsigned kCreditedTrk = 1U, kCreditedCal = 2U, kCreditedTrig = 4U;
  /// (frame, class_index, quantised time) of every recovered collision ->
  /// the stage terms already credited to it, as a bitmask of kCredited*. A
  /// collision seen by two overlapping candidates counts once, but each term
  /// counts if ANY of those candidates passed it, so the value has to be
  /// carried rather than decided at first insert.
  std::map<std::tuple<std::uint64_t, int, long>, unsigned> m_found;

  int  m_taps_registered = 0;
  int  m_taps_finished   = 0;
  bool m_standalone      = false;
  std::uint64_t m_last_report_frames = 0;
  std::uint64_t m_blind_frames      = 0;
  bool                                  m_clock_started = false;
  std::chrono::steady_clock::time_point m_t_first{};
  std::chrono::steady_clock::time_point m_t_last{};
  std::uint64_t m_child_frames      = 0;
  bool          m_have_child_frames = false;

  int         m_report_every = 250;
  bool        m_stage3       = true;
  double      m_frame_ns     = 2000.0;
  double      m_calo_dr      = 0.05;
  double      m_nsigma_window = 3.0;
  std::string m_csv_path;
  bool        m_profile = false;
  std::map<std::string, std::pair<double, std::uint64_t>> m_factory_time;
  std::string m_pdf_path;
  std::vector<std::string> m_input_files;
  ResolutionHists m_res;
};

} // namespace eicrecon::eb

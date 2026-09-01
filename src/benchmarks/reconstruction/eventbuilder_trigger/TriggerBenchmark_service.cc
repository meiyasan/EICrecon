// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "TriggerBenchmark_service.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <JANA/JApplication.h>
#include <JANA/Services/JParameterManager.h>

#include "services/log/Log_service.h"

namespace eicrecon::eb {

namespace {

/// Percentage cell. Renders "--" for an empty denominator rather than 0.0,
/// so "nothing to find" never reads as "found nothing".
std::string pct(std::uint64_t num, std::uint64_t den, int width = 7) {
  std::ostringstream os;
  if (den == 0)
    os << std::setw(width) << "--";
  else
    os << std::setw(width) << std::fixed << std::setprecision(1)
       << (100.0 * static_cast<double>(num) / static_cast<double>(den));
  return os.str();
}

std::string num(std::uint64_t v, int width = 7) {
  std::ostringstream os;
  os << std::setw(width) << v;
  return os.str();
}

/// "eff/pur" in one fixed-width cell. Either half renders "--" on an empty
/// denominator, so "nothing to find" and "nothing to be wrong about" stay
/// visually distinct from a genuine 0.0%.
std::string pair_cell(std::uint64_t en, std::uint64_t ed, std::uint64_t pn, std::uint64_t pd,
                      int width = 13) {
  auto half = [](std::uint64_t n, std::uint64_t d) {
    std::ostringstream h;
    if (d == 0)
      h << "--";
    else
      h << std::fixed << std::setprecision(1) << (100.0 * static_cast<double>(n) /
                                                  static_cast<double>(d));
    return h.str();
  };
  std::ostringstream os;
  os << std::setw(width) << (half(en, ed) + "/" + half(pn, pd));
  return os.str();
}

std::string dash_pair(int width = 13) {
  std::ostringstream os;
  os << std::setw(width) << "--/--";
  return os.str();
}

std::string dash(int width = 7) {
  std::ostringstream os;
  os << std::setw(width) << "--";
  return os.str();
}

} // namespace

void TriggerBenchmark_service::acquire_services(JServiceLocator* /*locator*/) {
  m_log = m_app->GetService<Log_service>()->logger("EventBuilderBenchmark");

  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:report_every", m_report_every,
      "Emit an interim benchmark table every N frames (timeslices). 0 = final table only.");
  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:stage3", m_stage3,
      "Score ACTS tracks and calo clusters per class (STAGE 3). Off = STAGE 1 + GNN only.");
  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:frame_ns", m_frame_ns,
      "Time-frame length in ns, used to turn fakes/frame into a false-alarm rate.");
  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:calo_dr", m_calo_dr,
      "eta-phi radius within which neutral truth particles merge into one calo target.");
  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:csv", m_csv_path,
      "Write the final table as CSV to this path. Empty = no CSV.");
  m_app->SetDefaultParameter(
      "eventbuilder:benchmark:pdf", m_pdf_path,
      "Write a multi-page PDF report (title, software provenance, table, plots) to this "
      "path. Empty = no PDF.");

  // weights[T0SIGMA] is NOT a 1-sigma uncertainty: EventBuilder_factory stores
  //   t0sigma = nsigma_window * sqrt(1 / sum(1/sigma_i^2))
  // i.e. the N-sigma HALF-WIDTH of the confidence interval, the same
  // nsigma_window that sets the coincidence windows. Dividing a residual by it
  // therefore under-states the pull by exactly that factor. Recover the true
  // 1-sigma by dividing it back out; the parameter is read from whichever name
  // the factory registered so this cannot silently drift from the trigger.
  for (const char* name : {"eventbuilder:nsigma_window", "trigger:nsigma_window",
                           "eventbuilder:trigger:nsigma_window"}) {
    try {
      const auto v = m_app->GetParameterValue<double>(name);
      if (v > 0.0) {
        m_nsigma_window = v;
        break;
      }
    } catch (...) {
      // not this name; fall through to the documented default of 3
    }
  }

  // A PDF needs the residual plots even when no -Phistsfile was requested, so
  // book them in memory in that case.
  m_res.init(m_app, !m_pdf_path.empty());
}

void TriggerBenchmark_service::addInputFile(const std::string& f) {
  std::lock_guard<std::mutex> lock(m_mutex);
  if (std::find(m_input_files.begin(), m_input_files.end(), f) == m_input_files.end())
    m_input_files.push_back(f);
}

void TriggerBenchmark_service::registerTap() {
  std::lock_guard<std::mutex> lock(m_mutex);
  ++m_taps_registered;
}

void TriggerBenchmark_service::tapFinished() {
  bool last = false;
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    last = (++m_taps_finished >= m_taps_registered);
  }
  if (last)
    report(true);
}

void TriggerBenchmark_service::addFrame(const Totals& d,
                                        const std::map<int, ClassRow>& per_class,
                                        const std::vector<FoundKey>& found_keys) {
  std::lock_guard<std::mutex> lock(m_mutex);

  m_totals.frames += d.frames;
  m_totals.candidates += d.candidates;
  m_totals.real += d.real;
  m_totals.fake += d.fake;
  m_totals.trk_ok += d.trk_ok;
  m_totals.trk_fake += d.trk_fake;
  m_totals.cal_ok += d.cal_ok;
  m_totals.cal_fake += d.cal_fake;
  m_totals.trig_ok += d.trig_ok;
  m_totals.trig_fake += d.trig_fake;
  m_totals.gnn_frames += d.gnn_frames;
  m_totals.gnn_passed += d.gnn_passed;
  m_totals.gnn_fake_frames += d.gnn_fake_frames;
  m_totals.gnn_fake_accepted += d.gnn_fake_accepted;
  m_totals.collisions_injected += d.collisions_injected;
  m_totals.saw_gnn = m_totals.saw_gnn || d.saw_gnn;
  if (d.cal_threshold >= 0.0) {
    m_totals.trk_threshold = d.trk_threshold;
    m_totals.cal_threshold = d.cal_threshold;
  }

  for (const auto& [ci, row] : per_class) {
    auto& dst = m_classes[ci];
    dst.injected += row.injected;
    dst.candidates += row.candidates;
    dst.real_cands += row.real_cands;
    dst.trk_ok += row.trk_ok;
    dst.cal_ok += row.cal_ok;
    dst.trig_ok += row.trig_ok;
  }

  // Deduplicate recovered collisions: overlapping candidate windows can each
  // claim the same injected collision, and it must count once.
  //
  // The stage terms are OR-ed across those claimants instead of being read off
  // whichever candidate was seen first. A collision IS track-triggered if any
  // candidate that recovered it cleared the track threshold, and first-wins
  // would have charged it to the first claimant's flags -- an arbitrary
  // choice, and a biased one wherever duplicates are common (the time purity
  // here runs near 76%, so roughly one collision in four has a second
  // claimant whose flags were being discarded).
  for (const auto& fk : found_keys) {
    const auto   ci = fk.cls;
    const double t  = fk.time;
    const auto   frame = fk.frame;
    // Quantise to 1 ps. The trigger_classes time and the MCParticle time come
    // from the same source, so an exact-equality key would work in principle;
    // quantising protects against float32 round-trips through the weights.
    const long key = static_cast<long>(t * 1000.0);

    const auto [slot, first_claim] = m_found.emplace(std::make_tuple(frame, ci, key), 0U);
    if (first_claim) {
      ++m_classes[ci].found;
      ++m_totals.collisions_found;
    }
    unsigned& credited = slot->second;

    if (fk.trk && (credited & kCreditedTrk) == 0U) {
      credited |= kCreditedTrk;
      ++m_classes[ci].found_trk;
      ++m_totals.found_trk;
    }
    if (fk.cal && (credited & kCreditedCal) == 0U) {
      credited |= kCreditedCal;
      ++m_classes[ci].found_cal;
      ++m_totals.found_cal;
    }
    if ((fk.trk || fk.cal) && (credited & kCreditedTrig) == 0U) {
      credited |= kCreditedTrig;
      ++m_classes[ci].found_trig;
      ++m_totals.found_trig;
    }
  }
}

void TriggerBenchmark_service::addEvent(const Totals& d,
                                        const std::map<int, ClassRow>& per_class) {
  std::lock_guard<std::mutex> lock(m_mutex);
  m_totals.trk_matched += d.trk_matched;
  m_totals.trk_expected += d.trk_expected;
  m_totals.cal_matched += d.cal_matched;
  m_totals.cal_expected += d.cal_expected;
  m_totals.trk_reco += d.trk_reco;
  m_totals.trk_reco_matched += d.trk_reco_matched;
  m_totals.cal_reco += d.cal_reco;
  m_totals.cal_reco_matched += d.cal_reco_matched;
  m_totals.children += d.children;
  m_totals.saw_stage3 = m_totals.saw_stage3 || d.saw_stage3;
  m_totals.used_stored_nexp = m_totals.used_stored_nexp || d.used_stored_nexp;
  m_totals.stored_nexp += d.stored_nexp;

  for (const auto& [ci, row] : per_class) {
    auto& dst = m_classes[ci];
    dst.trk_matched += row.trk_matched;
    dst.trk_expected += row.trk_expected;
    dst.cal_matched += row.cal_matched;
    dst.cal_expected += row.cal_expected;
    dst.trk_reco_any += row.trk_reco_any;
    dst.trk_reco_good += row.trk_reco_good;
    dst.cal_reco_any += row.cal_reco_any;
    dst.cal_reco_good += row.cal_reco_good;
    dst.labeled_hits += row.labeled_hits;
  }
}

void TriggerBenchmark_service::addChildFrames(std::uint64_t n) {
  std::lock_guard<std::mutex> lock(m_mutex);
  m_child_frames += n;
  m_have_child_frames = true;
  // Frames the STAGE 3 tap never saw produced no child event, so their
  // injected collisions are missing from that tap's denominators. Counting
  // them by difference is order-independent; inferring them from gaps in the
  // frame numbering was not (see EventBenchmark_processor.cc).
  m_blind_frames = (m_totals.frames > m_child_frames) ? m_totals.frames - m_child_frames : 0;
}

void TriggerBenchmark_service::addBlindFrames(std::uint64_t n) {
  std::lock_guard<std::mutex> lock(m_mutex);
  m_blind_frames += n;
}

void TriggerBenchmark_service::maybeReport() {
  if (m_report_every <= 0)
    return;
  bool due = false;
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_totals.frames - m_last_report_frames >= static_cast<std::uint64_t>(m_report_every)) {
      m_last_report_frames = m_totals.frames;
      due                  = true;
    }
  }
  if (due)
    report(false);
}

void TriggerBenchmark_service::report(bool final_report) {
  const std::string table = renderTable(final_report);
  if (m_log)
    m_log->info("\n{}", table);
  if (final_report && !m_csv_path.empty())
    writeCsv();
  if (final_report && !m_pdf_path.empty()) {
    ReportContext ctx;
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      ctx.input_files = m_input_files;
    }
    // The full invocation, straight from the kernel: reporting a hand-picked
    // subset of parameters (as this page used to) hides exactly the flag that
    // turns out to matter when a number looks wrong.
    {
      std::ifstream cl("/proc/self/cmdline", std::ios::binary);
      std::string   raw((std::istreambuf_iterator<char>(cl)), std::istreambuf_iterator<char>());
      for (auto& ch : raw)
        if (ch == '\0')
          ch = ' ';
      while (!raw.empty() && raw.back() == ' ')
        raw.pop_back();
      ctx.command_line = raw;
    }
    ctx.mode       = m_standalone ? "standalone" : "inline";
    ctx.nsigma_window = m_nsigma_window;
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      ctx.frames     = m_totals.frames;
      ctx.candidates = m_totals.candidates;
      ctx.real       = m_totals.real;
      ctx.fake       = m_totals.fake;
      if (m_totals.collisions_injected > 0)
        ctx.overall_eff = static_cast<double>(m_totals.collisions_found) /
                          static_cast<double>(m_totals.collisions_injected);
      if (m_totals.candidates > 0)
        ctx.overall_pur =
            static_cast<double>(m_totals.real) / static_cast<double>(m_totals.candidates);
      for (int ci = 0; ci < kNumClasses; ++ci) {
        auto it = m_classes.find(ci);
        if (it == m_classes.end() || it->second.injected == 0)
          continue; // a class with nothing injected has no efficiency to show
        ClassEff ce;
        ce.name     = std::string(className(ci));
        ce.injected = it->second.injected;
        ce.found    = it->second.found;
        ce.eff      = static_cast<double>(ce.found) / static_cast<double>(ce.injected);
        ctx.classes.push_back(ce);
      }
    }
    ctx.table  = tableCells();
    ctx.footer = tableFooter();
    if (writePdfReport(m_pdf_path, table, m_res.hists(), ctx) && m_log)
      m_log->info("benchmark: wrote {}", m_pdf_path);
  }
}

std::string TriggerBenchmark_service::renderTable(bool final_report) const {
  std::lock_guard<std::mutex> lock(m_mutex);
  std::ostringstream o;

  const char* mode = m_standalone ? "standalone (reduced denominators -- "
                                    "zero-candidate frames invisible)"
                                  : "inline";
  o << (final_report ? "== eventbuilder trigger benchmark (final) "
                     : "== eventbuilder trigger benchmark (progress) ")
    << std::string(30, '=') << "\n"
    << "   frames " << m_totals.frames << " | candidates " << m_totals.candidates << " | real "
    << m_totals.real << " | fake " << m_totals.fake << " | children " << m_totals.children
    << " | mode: " << mode << "\n\n";

  o << "   " << std::left << std::setw(14) << "class" << std::right << std::setw(6) << "inj"
    << std::setw(7) << "found" << std::setw(13) << "time" << std::setw(13) << "track"
    << std::setw(13) << "calo" << std::setw(13) << "trigger" << std::setw(13) << "gnn"
    << std::setw(13) << "acts-trk" << std::setw(13) << "calo-island" << "\n";
  o << "   " << std::left << std::setw(14) << "" << std::right << std::setw(6) << "" << std::setw(7)
    << "" << std::setw(13) << "eff/pur" << std::setw(13) << "eff/pur" << std::setw(13) << "eff/pur"
    << std::setw(13) << "eff/pur" << std::setw(13) << "eff/pur" << std::setw(13) << "eff/pur"
    << std::setw(13) << "eff/pur" << "\n";
  o << "   " << std::string(118, '-') << "\n";

  for (int ci = 0; ci < kNumClasses; ++ci) {
    auto it = m_classes.find(ci);
    if (it == m_classes.end() || (it->second.injected == 0 && it->second.candidates == 0))
      continue; // class absent from this stream; a 0/0 row is noise
    const ClassRow& r = it->second;
    // Same cells tableCells() hands the PDF: the two renderers read the same
    // counters through the same expressions, so the text table and the report
    // cannot drift apart. Every efficiency half divides by `injected`; the
    // purity halves are the duplicate shortfall. `gnn` has no per-class form --
    // its decision is per frame, not per candidate.
    o << "   " << std::left << std::setw(14) << className(ci) << std::right << num(r.injected, 6)
      << num(r.found, 7) << pair_cell(r.found, r.injected, r.found, r.real_cands)
      << pair_cell(r.trk_ok, r.real_cands, r.found_trk, r.trk_ok)
      << pair_cell(r.cal_ok, r.real_cands, r.found_cal, r.cal_ok)
      << pair_cell(r.found_trig, r.injected, r.found_trig, r.trig_ok) << dash_pair()
      << (m_totals.saw_stage3
              ? pair_cell(r.trk_matched, r.trk_expected, r.trk_reco_good, r.trk_reco_any)
              : dash_pair())
      << (m_totals.saw_stage3
              ? pair_cell(r.cal_matched, r.cal_expected, r.cal_reco_good, r.cal_reco_any)
              : dash_pair())
      << "\n";
  }

  o << "   " << std::string(118, '-') << "\n";
  o << "   " << std::left << std::setw(14) << "TOTAL" << std::right
    << num(m_totals.collisions_injected, 6) << num(m_totals.collisions_found, 7)
    // time: collisions recovered / injected, and real / all candidates
    << pair_cell(m_totals.collisions_found, m_totals.collisions_injected, m_totals.real,
                 m_totals.candidates)
    // track/calo: each primitive's OWN acceptance over the candidates it
    // judged -- independent of the time column and of each other. trigger is
    // the one COMPOSED column: TIME && (TRACK || CALO), counted in distinct
    // injected collisions, so trigger <= time always holds and the row reads
    // as the trigger chain end to end.
    << pair_cell(m_totals.trk_ok, m_totals.real, m_totals.trk_ok,
                 m_totals.trk_ok + m_totals.trk_fake)
    << pair_cell(m_totals.cal_ok, m_totals.real, m_totals.cal_ok,
                 m_totals.cal_ok + m_totals.cal_fake)
    << pair_cell(m_totals.found_trig, m_totals.collisions_injected, m_totals.trig_ok,
                 m_totals.trig_ok + m_totals.trig_fake)
    << (m_totals.saw_gnn ? pair_cell(m_totals.gnn_passed, m_totals.gnn_frames, m_totals.gnn_passed,
                                     m_totals.gnn_passed + m_totals.gnn_fake_accepted)
                         : dash_pair())
    << (m_totals.saw_stage3 ? pair_cell(m_totals.trk_matched, m_totals.trk_expected,
                                        m_totals.trk_reco_matched, m_totals.trk_reco)
                            : dash_pair())
    << (m_totals.saw_stage3 ? pair_cell(m_totals.cal_matched, m_totals.cal_expected,
                                        m_totals.cal_reco_matched, m_totals.cal_reco)
                            : dash_pair())
    << "\n";

  // How often each stage ACCEPTS a fake. Not an efficiency: the denominator is
  // the fake population, so lower is better and there is no purity half.
  o << "   " << std::left << std::setw(14) << "FAKE-ACCEPT" << std::right << dash(6) << dash(7)
    << dash_pair() << pair_cell(m_totals.trk_fake, m_totals.fake, 0, 0)
    << pair_cell(m_totals.cal_fake, m_totals.fake, 0, 0)
    << pair_cell(m_totals.trig_fake, m_totals.fake, 0, 0)
    << (m_totals.saw_gnn ? pair_cell(m_totals.gnn_fake_accepted, m_totals.gnn_fake_frames, 0, 0)
                         : dash_pair())
    << dash_pair() << dash_pair() << "\n\n";

  // Per-class purity is absent by construction, not by omission: a fake
  // candidate has trigger_classes_mask == 0 and an empty trigger_classes
  // tail, and a ghost track/cluster matches no MC particle, so neither has a
  // class to be charged against.

  if (m_totals.frames > 0) {
    const double per_frame =
        static_cast<double>(m_totals.fake) / static_cast<double>(m_totals.frames);
    const double hz = (m_frame_ns > 0.0) ? per_frame / (m_frame_ns * 1e-9) : 0.0;
    o << "   FAR     " << m_totals.fake << " fakes / " << m_totals.frames << " frames = "
      << std::fixed << std::setprecision(3) << per_frame << "/frame = " << std::setprecision(1)
      << hz * 1e-3 << " kHz  (frame = " << m_frame_ns << " ns)\n";
  }

  const std::uint64_t missed = (m_totals.collisions_injected > m_totals.collisions_found)
                                   ? m_totals.collisions_injected - m_totals.collisions_found
                                   : 0;
  o << "   recovery " << m_totals.collisions_found << " / " << m_totals.collisions_injected
    << " collisions found (" << missed << " missed)\n";

  if (m_totals.used_stored_nexp)
    o << "   in-acceptance charged   " << m_totals.trk_expected << " recomputed   |   "
      << m_totals.stored_nexp << " stored (N_EXPECTED)\n";
  if (m_blind_frames > 0)
    o << "   BLIND   " << m_blind_frames
      << " frame(s) had no candidate; their collisions are not in the denominator (%eff is "
         "an over-estimate)\n";
  if (!m_totals.saw_stage3)
    o << "   note    STAGE 3 columns '--': no ACTS/calo associations seen\n";

  return o.str();
}

std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
TriggerBenchmark_service::tableCells() const {
  std::lock_guard<std::mutex> lock(m_mutex);
  auto cell = [](std::uint64_t n, std::uint64_t d, std::uint64_t pn, std::uint64_t pd) {
    auto half = [](std::uint64_t a, std::uint64_t b) {
      if (b == 0)
        return std::string("--");
      std::ostringstream o;
      o << std::fixed << std::setprecision(1) << (100.0 * double(a) / double(b));
      return o.str();
    };
    return half(n, d) + "/" + half(pn, pd);
  };

  // With min_cal_energy = 0 the calo term accepts everything, so the calo and
  // trigger columns are algebraically identical to time -- one measurement
  // printed three times. Say "off" instead of repeating it.
  const bool cal_off = (m_totals.cal_threshold == 0.0);
  std::vector<std::string> header = {"class", "inj",  "found", "time",    "track",
                                     "calo",  "trigger", "gnn", "acts-trk", "calo-island"};
  std::vector<std::vector<std::string>> rows;

  for (int ci = 0; ci < kNumClasses; ++ci) {
    auto it = m_classes.find(ci);
    if (it == m_classes.end() || (it->second.injected == 0 && it->second.candidates == 0))
      continue;
    const ClassRow& r = it->second;
    rows.push_back({std::string(className(ci)), std::to_string(r.injected),
                    std::to_string(r.found), cell(r.found, r.injected, r.found, r.real_cands),
                    cell(r.trk_ok, r.real_cands, r.found_trk, r.trk_ok),
                    cell(r.cal_ok, r.real_cands, r.found_cal, r.cal_ok),
                    cell(r.found_trig, r.injected, r.found_trig, r.trig_ok), "--/--",
                    m_totals.saw_stage3
                        ? cell(r.trk_matched, r.trk_expected, r.trk_reco_good, r.trk_reco_any)
                        : "--/--",
                    m_totals.saw_stage3
                        ? cell(r.cal_matched, r.cal_expected, r.cal_reco_good, r.cal_reco_any)
                        : "--/--"});
  }
  rows.push_back({});
  rows.push_back(
      {"TOTAL", std::to_string(m_totals.collisions_injected),
       std::to_string(m_totals.collisions_found),
       cell(m_totals.collisions_found, m_totals.collisions_injected, m_totals.real,
            m_totals.candidates),
       cell(m_totals.trk_ok, m_totals.real, m_totals.trk_ok, m_totals.trk_ok + m_totals.trk_fake),
       cell(m_totals.cal_ok, m_totals.real, m_totals.cal_ok, m_totals.cal_ok + m_totals.cal_fake),
       cell(m_totals.found_trig, m_totals.collisions_injected, m_totals.trig_ok,
            m_totals.trig_ok + m_totals.trig_fake),
       m_totals.saw_gnn ? cell(m_totals.gnn_passed, m_totals.gnn_frames, m_totals.gnn_passed,
                               m_totals.gnn_passed + m_totals.gnn_fake_accepted)
                        : "--/--",
       m_totals.saw_stage3 ? cell(m_totals.trk_matched, m_totals.trk_expected,
                                  m_totals.trk_reco_matched, m_totals.trk_reco)
                           : "--/--",
       m_totals.saw_stage3 ? cell(m_totals.cal_matched, m_totals.cal_expected,
                                  m_totals.cal_reco_matched, m_totals.cal_reco)
                           : "--/--"});
  rows.push_back({"FAKE-ACCEPT", "--", "--", "--/--",
                  cell(m_totals.trk_fake, m_totals.fake, 0, 0),
                  cal_off ? "off" : cell(m_totals.cal_fake, m_totals.fake, 0, 0),
                  cal_off ? "off" : cell(m_totals.trig_fake, m_totals.fake, 0, 0),
                  m_totals.saw_gnn ? cell(m_totals.gnn_fake_accepted, m_totals.gnn_fake_frames, 0, 0)
                                   : "--/--",
                  "--/--", "--/--"});
  return {header, rows};
}

std::vector<std::string> TriggerBenchmark_service::tableFooter() const {
  std::lock_guard<std::mutex> lock(m_mutex);
  std::vector<std::string> f;
  std::ostringstream o;
  if (m_totals.frames > 0) {
    const double per_frame = double(m_totals.fake) / double(m_totals.frames);
    const double hz        = (m_frame_ns > 0.0) ? per_frame / (m_frame_ns * 1e-9) : 0.0;
    o << "FAR   " << m_totals.fake << " fakes / " << m_totals.frames << " frames = " << std::fixed
      << std::setprecision(3) << per_frame << "/frame = " << std::setprecision(1) << hz * 1e-3
      << " kHz";
    f.push_back(o.str());
  }
  std::ostringstream r;
  r << "recovery   " << m_totals.collisions_found << " / " << m_totals.collisions_injected
    << " collisions found";
  f.push_back(r.str());
  // BLIND and the recomputed-vs-stored in-acceptance count belong in the
  // REPORT, not only the console: they qualify the very numbers printed above
  // them -- BLIND says the efficiency on this page is an over-estimate, and
  // the denominator comparison says by roughly how much the ACTS column could
  // move. A reader of the PDF alone must see both.
  if (m_totals.cal_threshold >= 0.0) {
    std::ostringstream t;
    t << "thresholds   min_tracklets = " << static_cast<long>(m_totals.trk_threshold)
      << ",  min_cal_energy = " << std::fixed << std::setprecision(3) << m_totals.cal_threshold
      << " GeV";
    if (m_totals.cal_threshold == 0.0)
      t << "  -> calo term DISABLED (accepts every candidate); calo and trigger "
           "columns would duplicate time";
    f.push_back(t.str());
  }
  if (m_totals.used_stored_nexp) {
    std::ostringstream n;
    n << "in-acceptance charged   " << m_totals.trk_expected << " recomputed  |  "
      << m_totals.stored_nexp << " stored (N_EXPECTED)";
    f.push_back(n.str());
  }
  if (m_blind_frames > 0) {
    std::ostringstream b;
    b << "BLIND   " << m_blind_frames
      << " frame(s) had no candidate; their collisions are not in the denominator, so "
         "efficiency is an over-estimate";
    f.push_back(b.str());
  }
  return f;
}

void TriggerBenchmark_service::writeCsv() const {
  std::lock_guard<std::mutex> lock(m_mutex);
  std::ofstream f(m_csv_path);
  if (!f) {
    if (m_log)
      m_log->error("benchmark: cannot write CSV to '{}'", m_csv_path);
    return;
  }
  f << "class,injected,found,candidates,real_cands,trk_matched,trk_expected,"
       "cal_matched,cal_expected,labeled_hits,trk_reco_matched,trk_reco,"
       "cal_reco_matched,cal_reco\n";
  for (int ci = 0; ci < kNumClasses; ++ci) {
    auto it = m_classes.find(ci);
    if (it == m_classes.end())
      continue;
    const ClassRow& r = it->second;
    f << className(ci) << ',' << r.injected << ',' << r.found << ',' << r.candidates << ','
      << r.real_cands << ',' << r.trk_matched << ',' << r.trk_expected << ',' << r.cal_matched
      << ',' << r.cal_expected << ',' << r.labeled_hits << ",,,,\n";
  }
  f << "__TOTAL__," << m_totals.collisions_injected << ',' << m_totals.collisions_found << ','
    << m_totals.candidates << ',' << m_totals.real << ',' << m_totals.trk_matched << ','
    << m_totals.trk_expected << ',' << m_totals.cal_matched << ',' << m_totals.cal_expected
    << ",0," << m_totals.trk_reco_matched << ',' << m_totals.trk_reco << ','
    << m_totals.cal_reco_matched << ',' << m_totals.cal_reco << '\n';
  if (m_log)
    m_log->info("benchmark: wrote {}", m_csv_path);
}

} // namespace eicrecon::eb

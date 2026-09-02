// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// PdfReport -- the benchmark's multi-page PDF.
//
//   1  title page: project, generation time
//   2  software provenance: versions and the run's own configuration
//   3+ the statistical table, verbatim, paginated
//   4+ residual plots, grouped by detector, four to a page
//
// Written once at Finish(), from the same numbers the terminal table shows.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class TH2D;

namespace eicrecon::eb {

/// One bar of the per-class efficiency chart.
struct ClassEff {
  std::string   name;
  double        eff = 0.0; ///< 0..1
  std::uint64_t found = 0;
  std::uint64_t injected = 0;
};

struct ReportContext {
  std::vector<std::string> input_files; ///< every source actually read
  std::string              command_line; ///< the full eicrecon invocation
  std::string detector_config;
  std::string mode;      ///< "inline" / "standalone"
  std::uint64_t frames = 0;
  std::uint64_t candidates = 0;
  std::uint64_t real = 0;
  std::uint64_t fake = 0;
  double        overall_eff = -1.0; ///< headline number; <0 = unavailable
  double        overall_pur = -1.0;
  double        nsigma_window = 3.0;
  std::vector<ClassEff> classes;
  std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>> table;
  std::vector<std::string> footer;
  /// {factory, total seconds, calls}, slowest first. Empty unless profiling.
  std::vector<std::tuple<std::string, double, std::uint64_t>> profile;
  /// Same shape, per physics class: {class, total seconds, candidates}.
  std::vector<std::tuple<std::string, double, std::uint64_t>> class_profile;
};

/// Writes the report. Returns false if the file could not be opened.
bool writePdfReport(const std::string& path, const std::string& table,
                    const std::map<std::string, TH2D*>& hists, const ReportContext& ctx);

} // namespace eicrecon::eb

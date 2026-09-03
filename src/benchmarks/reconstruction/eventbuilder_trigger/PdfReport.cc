// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "PdfReport.h"

#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <sstream>
#include <vector>

#include <JANA/JVersion.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <TH2D.h>
#include <TLine.h>
#include <TBox.h>
#include <TPad.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <edm4eic/EDM4eicVersion.h>
#include <podio/podioVersion.h>

#include "ResolutionHists.h"

namespace eicrecon::eb {

namespace {

// Portrait, in the proportions of the octofit reports (399 x 567 pt).
constexpr int kCanvasW = 850;
constexpr int kCanvasH = 1208;
// Font ids taken from octofit's Report::Style (tools/xfitter-report.cc:302,
// Report/ExperimentComparisonPage.cc:574): 132 = Times-Roman, 32 = the bold
// face its headings use, 82 = Courier for table cells.
constexpr int kFont     = 132;
constexpr int kFontBold = 32;
constexpr int kMono     = 82;

/// Draw pre-formatted lines top-down. Returns nothing; callers paginate.
void drawLines(const std::vector<std::string>& lines, double x0, double y0, double dy,
               double size, int font = kMono) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(font);
  t.SetTextSize(size);
  double y = y0;
  for (const auto& l : lines) {
    // TLatex treats '#' and '^' as markup; the table has neither, but a
    // detector name could. SetTextFont(82) still parses them, so draw with
    // TLatex only after neutralising the one character we actually emit.
    std::string safe = l;
    for (auto& c : safe)
      if (c == '#')
        c = '+';
    t.DrawLatex(x0, y, safe.c_str());
    y -= dy;
  }
}


/// Gaussian sigma of the CORE of a distribution: fit within +-2 RMS of the
/// peak, not the full range. For the time residuals this matters -- the tails
/// are the coincidence gate, flat out to +-dt, and a full-range fit would
/// report the gate width instead of the detector.
/// Fit and leave the function attached to the histogram so it can be drawn.
/// Returns nullptr if the fit did not produce a usable core. coreSigma() below
/// is the measurement-only wrapper.
TF1* fitCore(TH1* h, double& sigma, double& err, double* mu_out = nullptr) {
  if (h == nullptr || h->GetEntries() < 200)
    return nullptr;
  const double axis_lo   = h->GetXaxis()->GetXmin();
  const double axis_hi   = h->GetXaxis()->GetXmax();
  const double half_span = 0.5 * (axis_hi - axis_lo);
  const double bin_w     = h->GetBinWidth(1);
  const double rms       = h->GetRMS();
  const double peak      = h->GetBinContent(h->GetMaximumBin());
  const double peak_x    = h->GetBinCenter(h->GetMaximumBin());
  if (!(rms > 0) || !(peak > 0))
    return nullptr;

  // A SINGLE Gaussian, iterated over +-2.5 sigma of the peak -- which is what
  // the comment above always described, and what the code did not do.
  //
  // It used to fit gaus(0)+gaus(3) across the whole axis and report whichever
  // component came out NARROWER, guarded only by that component holding 5% of
  // the peak height. On a genuinely double-humped distribution that is fine.
  // On a clean single Gaussian the two components are degenerate: the fit
  // parks a narrow spike on the top few bins, lets the broad one absorb the
  // rest, and reports the spike. That is how the TOF sensor residual -- a pure
  // 25 ps digitisation smear, about as Gaussian as it gets -- came back as
  // 9.78 ps, well under the width actually injected into it.
  //
  // Restricting the range instead of adding a second component gets the same
  // tail rejection honestly: the coincidence-gate tails sit far outside
  // 2.5 sigma and are simply not fitted.
  auto* f = new TF1("core", "gaus", axis_lo, axis_hi);
  double mu = peak_x;
  double sg = std::max(rms, 2.0 * bin_w);
  for (int pass = 0; pass < 4; ++pass) {
    const double lo = std::max(mu - 2.5 * sg, axis_lo);
    const double hi = std::min(mu + 2.5 * sg, axis_hi);
    if (!(hi > lo)) {
      delete f;
      return nullptr;
    }
    f->SetParameters(peak, mu, sg);
    f->SetParLimits(0, 0.0, 10.0 * peak);
    f->SetParLimits(1, lo, hi);
    f->SetParLimits(2, 0.5 * bin_w, 4.0 * half_span);
    if (h->Fit(f, "QNR", "", lo, hi) != 0) {
      delete f;
      return nullptr;
    }
    mu                 = f->GetParameter(1);
    const double new_s = std::fabs(f->GetParameter(2));
    if (!(new_s > 0)) {
      delete f;
      return nullptr;
    }
    const bool settled = std::fabs(new_s - sg) < 0.01 * sg;
    sg                 = new_s;
    if (settled)
      break; // converged: another pass would move it by <1%
  }
  sigma = sg;
  err   = f->GetParError(2);
  // Unresolved (narrower than a bin) or so wide it is the axis, not a core.
  if (sigma <= bin_w || sigma >= 0.5 * half_span) {
    delete f;
    return nullptr;
  }
  if (mu_out != nullptr)
    *mu_out = mu;
  f->SetRange(mu - 3.0 * sg, mu + 3.0 * sg); // draw over what was fitted
  return f;
}

bool coreSigma(TH1* h, double& sigma, double& err) {
  TF1* f = fitCore(h, sigma, err);
  if (f == nullptr)
    return false;
  delete f;
  return true;
}

/// Two Gaussians sharing one mean: a narrow core (the resolution) and a wide
/// tail. This is the traditional two-scale model for a residual -- a single
/// Gaussian is dragged wide by the tail and reports neither. Returns the
/// fitted TF1 with core/tail widths and the tail-only amplitude, so the caller
/// can draw the TAIL alone and leave the core as visible excess above it.
/// sigma_core <= sigma_tail by construction.
TF1* fitDoubleGaussian(TH1* h, double& core, double& tail, double& mu, double& tail_amp) {
  // 100 was too strict: a sparse panel then showed no fit at all, which
  // reads as "the fit failed" rather than "there are 63 entries".
  if (h == nullptr || h->GetEntries() < 30)
    return nullptr;
  const double lo = h->GetXaxis()->GetXmin();
  const double hi = h->GetXaxis()->GetXmax();
  const double w  = 0.5 * (hi - lo);
  const double bw = h->GetBinWidth(1);
  const double s0 = h->GetRMS();
  const double a0 = h->GetMaximum();
  const double m0 = h->GetBinCenter(h->GetMaximumBin());
  if (!(s0 > 0) || !(a0 > 0))
    return nullptr;

  // Shared mean (par 4): the two components describe one population, so
  // letting the means float independently fits two unrelated peaks instead.
  auto* f = new TF1((std::string(h->GetName()) + "_dg").c_str(),
                    "[0]*exp(-0.5*((x-[4])/[1])^2)+[2]*exp(-0.5*((x-[4])/[3])^2)", lo, hi);
  f->SetParameters(0.7 * a0, 0.5 * s0, 0.3 * a0, 3.0 * s0, m0);
  f->SetParLimits(0, 0.0, 10.0 * a0);
  f->SetParLimits(1, bw, w);
  f->SetParLimits(2, 0.0, 10.0 * a0);
  f->SetParLimits(3, bw, w);
  f->SetParLimits(4, m0 - s0, m0 + s0);
  if (h->Fit(f, "QNR") != 0) {
    // Two components need enough entries to separate. Fall back to a single
    // Gaussian rather than reporting nothing: core is then the whole width
    // and tail is reported as absent.
    delete f;
    auto* g = new TF1((std::string(h->GetName()) + "_sg").c_str(), "gaus", lo, hi);
    g->SetParameters(a0, m0, s0);
    if (h->Fit(g, "QNR") != 0 || !(g->GetParameter(2) > 0)) {
      delete g;
      return nullptr;
    }
    core     = std::fabs(g->GetParameter(2));
    tail     = 0.0;
    mu       = g->GetParameter(1);
    tail_amp = 0.0;
    return g;
  }
  double a1 = f->GetParameter(0), s1 = std::fabs(f->GetParameter(1));
  double a2 = f->GetParameter(2), s2 = std::fabs(f->GetParameter(3));
  mu        = f->GetParameter(4);
  if (!(s1 > 0) || !(s2 > 0)) {
    delete f;
    return nullptr;
  }
  if (s1 > s2) { // convention: first component is the core
    std::swap(a1, a2);
    std::swap(s1, s2);
  }
  core     = s1;
  tail     = s2;
  tail_amp = a2;
  return f;
}

/// Background under the resolution peak: a wide Gaussian fitted to the
/// SIDEBANDS only, with the core region masked out. Fitting core+background
/// together and drawing the sum tells you the fit converged; fitting the
/// background alone and drawing only that leaves the core standing above it
/// as visible excess, which is what a resolution plot is actually asked.
/// Returns nullptr when the sidebands carry too little to constrain a shape.
TF1* fitBackground(TH1* h, double mu, double core_sigma) {
  if (h == nullptr || !(core_sigma > 0))
    return nullptr;
  // Mask +-3 sigma around the core; ROOT has no fit-range exclusion, so mask
  // by zeroing a clone rather than fitting a piecewise range.
  auto* side = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_side").c_str()));
  side->SetDirectory(nullptr);
  double kept = 0.0;
  for (int b = 1; b <= side->GetNbinsX(); ++b) {
    if (std::fabs(side->GetBinCenter(b) - mu) < 3.0 * core_sigma)
      side->SetBinContent(b, 0.0);
    else
      kept += side->GetBinContent(b);
  }
  if (kept < 50.0) { // no meaningful background to describe
    delete side;
    return nullptr;
  }
  auto* bkg = new TF1((std::string(h->GetName()) + "_bkg").c_str(), "gaus",
                      h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  bkg->SetParameters(side->GetMaximum(), mu, 4.0 * core_sigma);
  bkg->SetParLimits(0, 0.0, 10.0 * std::max(side->GetMaximum(), 1.0));
  bkg->SetParLimits(2, core_sigma, 20.0 * core_sigma);
  const int rc = side->Fit(bkg, "QNR");
  delete side;
  if (rc != 0 || !(bkg->GetParameter(2) > 0)) {
    delete bkg;
    return nullptr;
  }
  return bkg;
}

/// 68% containment. The position residuals are |rec - sim| and the calo
/// angular residual is dR: both are positive-definite and not Gaussian, so a
/// sigma is meaningless and the quantile is the honest summary.
bool quantile68(TH1* h, double& q) {
  if (h == nullptr || h->GetEntries() < 20)
    return false;
  double prob = 0.68, out = 0.0;
  h->GetQuantiles(1, &out, &prob);
  q = out;
  return true;
}

std::string fmt(double v, int prec, bool ok) {
  if (!ok)
    return "     --";
  std::ostringstream os;
  os << std::fixed << std::setprecision(prec) << v;
  return os.str();
}

/// The unit out of an axis title's trailing "[...]"; empty if dimensionless.
/// Reading it off the histogram is what lets one column hold two units without
/// lying: the tracker `position` cells are microns and the calo ones radians,
/// because those axes say so.
std::string axisUnit(const TH1* h) {
  if (h == nullptr)
    return {};
  const std::string at = h->GetXaxis()->GetTitle();
  const auto        ob = at.rfind('[');
  const auto        cb = at.rfind(']');
  if (ob == std::string::npos || cb == std::string::npos || cb <= ob + 1)
    return {};
  return at.substr(ob + 1, cb - ob - 1);
}

/// "value unit", with the precision matched to the magnitude so 265 microns
/// does not print as 265.0000 and 0.445 rad does not print as 0.
std::string fmtUnit(double v, bool ok, const std::string& unit) {
  if (!ok)
    return "--";
  const double a    = std::fabs(v);
  const int    prec = (a >= 100.0) ? 0 : (a >= 10.0 ? 1 : (a >= 1.0 ? 2 : 3));
  std::ostringstream os;
  os << std::fixed << std::setprecision(prec) << v;
  if (!unit.empty())
    os << " " << unit;
  return os.str();
}

/// A fraction as a percentage. The energy residual is (E_rec - E_MC)/E_MC, so
/// "-85.8%" says the reconstructed energy is 86% below truth -- far easier to
/// read at a glance than "-0.858".
std::string fmtPercent(double v, bool ok) {
  if (!ok)
    return "--";
  std::ostringstream os;
  os << std::fixed << std::setprecision(1) << 100.0 * v << "%";
  return os.str();
}


/// One table cell: text, plus the colour that carries its meaning.
struct Cell {
  std::string text;
  Color_t     color = kBlack;
};

/// Column geometry: NDC x, and TLatex alignment (11 = left, 31 = right).
struct Col {
  double x;
  short  align;
};

/// Draw a table the way the octofit reports do: header row, a rule spanning
/// the columns, then rows with per-cell colour. Drawing cell by cell rather
/// than one monospace string per line is what makes both the alignment and
/// the colouring possible.
/// octofit ParameterTablePage metrics: rows start at yTop, are baseDY apart,
/// and the body is drawn at 0.024. When there are more rows than fit, the row
/// height shrinks and the text shrinks with it by the same factor -- rather
/// than rows being silently dropped off the bottom.
constexpr double kTableYTop  = 0.818;
/// Table top on the resolution page only, which needs the space above for the
/// legend explaining what each column measures.
constexpr double kResTableYTop = 0.780;
constexpr double kTableYBot  = 0.05;
constexpr double kTableBaseDY = 0.032;
constexpr double kTableTextSize = 0.024;

/// Largest text size at which `ncols` columns of at most `maxchars` monospace
/// characters still fit inside `width` (NDC). octofit's 0.024 body size is
/// sized for its 8 narrow columns; a 10-column eff/pur table needs less, and
/// hard-coding 0.024 there just overlaps the columns.
double fitTextSize(int ncols, int maxchars, double width, double cap = kTableTextSize) {
  if (ncols <= 0 || maxchars <= 0)
    return cap;
  const double pitch = width / ncols;
  return std::min(cap, pitch / (0.62 * maxchars));
}

void drawTable(const std::vector<Col>& cols, const std::vector<std::string>& header,
               const std::vector<std::vector<Cell>>& rows, double y0, double dy, double size,
               double rule_x0, double rule_x1) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(kMono);
  t.SetTextSize(size);

  double y = y0;
  t.SetTextColor(kBlack);
  for (std::size_t i = 0; i < header.size() && i < cols.size(); ++i) {
    t.SetTextAlign(cols[i].align);
    t.DrawLatex(cols[i].x, y, header[i].c_str());
  }
  y -= 0.45 * dy;
  TLine rule;
  rule.SetNDC();
  rule.DrawLineNDC(rule_x0, y, rule_x1, y);
  y -= 0.75 * dy;

  for (const auto& row : rows) {
    if (row.empty()) { // blank spacer row
      y -= 0.5 * dy;
      continue;
    }
    for (std::size_t i = 0; i < row.size() && i < cols.size(); ++i) {
      t.SetTextAlign(cols[i].align);
      t.SetTextColor(row[i].color);
      t.DrawLatex(cols[i].x, y, row[i].text.c_str());
    }
    y -= dy;
  }
  t.SetTextColor(kBlack);
}

/// Bold serif centred page title with an optional subtitle and a rule under
/// it -- the header every c6 page carries.
double pageTitle(const std::string& title, const std::string& subtitle = "",
                 double y = 0.90) {
  TLatex t;
  t.SetNDC();
  t.SetTextAlign(22);
  t.SetTextFont(kFontBold);
  t.SetTextSize(0.044);
  t.DrawLatex(0.5, y, title.c_str());
  double next = y - 0.040;
  if (!subtitle.empty()) {
    t.SetTextFont(132); // Times-Roman
    t.SetTextSize(0.022);
    t.DrawLatex(0.5, next, subtitle.c_str());
    next -= 0.030;
  }
  TLine rule;
  rule.SetNDC();
  rule.DrawLineNDC(0.13, next, 0.87, next);
  return next - 0.045;
}

/// Longest monospace string that fits between x and the right margin.
int monoFit(double x, double size, double right = 0.94) {
  return std::max(8, static_cast<int>((right - x) / (0.62 * size)));
}

/// One colour per detector family, so a reader scanning the residual pages
/// can tell at a glance which subsystem a panel belongs to without reading the
/// title. Matched longest-prefix-first: OuterMPGDBarrel must not fall through
/// to the generic MPGD entry ahead of its own.
Color_t detectorColor(const std::string& name) {
  static const std::vector<std::pair<const char*, Color_t>> kFamilies = {
      {"TOF", kAzure + 1},        {"MPGD", kOrange + 7},
      {"SiBarrelVertex", kViolet + 1}, {"SiBarrel", kViolet - 5},
      {"SiEndcap", kViolet - 5},  {"B0Tracker", kTeal + 2},
      {"B0ECal", kTeal - 1},      {"EcalFarForward", kMagenta - 4},
      {"HcalFarForward", kMagenta - 7}, {"LFHCAL", kSpring - 6},
      {"Ecal", kRed - 4},         {"Hcal", kBlue - 4},
  };
  for (const auto& [prefix, col] : kFamilies)
    if (name.find(prefix) != std::string::npos)
      return col;
  return kGray + 2;
}

/// Split a path into at most `w`-character lines, breaking after a '/' so each
/// line is a run of whole directory names. Counting characters only works in
/// the MONO font -- the serif value column is proportional, which is why a
/// 48-character DETECTOR_PATH still ran off the right edge of the page.
std::vector<std::string> wrapPath(const std::string& path, int w) {
  std::vector<std::string> out;
  std::size_t              i = 0;
  while (i < path.size()) {
    std::size_t take = std::min(static_cast<std::size_t>(w), path.size() - i);
    if (i + take < path.size()) {
      // Back up to just after the last '/' that fits, unless that would leave
      // a stub -- a spack hash is one 60-character name with no break in it.
      const std::size_t cut = path.rfind('/', i + take);
      if (cut != std::string::npos && cut > i && (cut - i) > static_cast<std::size_t>(w) / 3)
        take = cut - i + 1;
    }
    out.push_back(path.substr(i, take));
    i += take;
  }
  return out;
}

/// Two-column label/value block in serif, as on c6's reference-set page.
double drawKeyValueEnd(const std::vector<std::pair<std::string, std::string>>& kv, double y0,
                       double dy = 0.034, double size = 0.026, double x_label = 0.20,
                       double x_value = 0.46);

void drawKeyValue(const std::vector<std::pair<std::string, std::string>>& kv, double y0,
                  double dy = 0.034, double size = 0.026, double x_label = 0.20,
                  double x_value = 0.46) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(132);
  t.SetTextSize(size);
  double y = y0;
  for (const auto& [k, v] : kv) {
    if (k.empty() && v.empty()) {
      y -= 0.5 * dy;
      continue;
    }
    t.SetTextAlign(11);
    t.DrawLatex(x_label, y, k.c_str());
    t.DrawLatex(x_value, y, v.c_str());
    y -= dy;
  }
}

/// Project the class axis away and rebin toward ~15 entries per bin. BOTH the
/// summary table and the plots go through this: fitting a raw 2000-bin
/// projection in one place and a rebinned one in the other made the same
/// quantity read 39.4 ps on its plot and "--" in the table.
/// drawKeyValue, returning the y it finished at so a caller can place
/// something underneath without recounting the rows.
double drawKeyValueEnd(const std::vector<std::pair<std::string, std::string>>& kv, double y0,
                       double dy, double size, double x_label, double x_value) {
  drawKeyValue(kv, y0, dy, size, x_label, x_value);
  double y = y0;
  for (const auto& [k, v] : kv)
    y -= (k.empty() && v.empty()) ? 0.5 * dy : dy;
  return y - 0.5 * dy;
}

/// Clip the x axis to the central 99% of the entries, so a peak that occupies
/// a percent of a deliberately generous axis still fills its pad.
void zoomToData(TH1* px) {
  if (px == nullptr)
    return;
  double q[2] = {0, 0}, pr[2] = {0.005, 0.995};
  px->GetQuantiles(2, q, pr);
  if (!(q[1] > q[0]))
    return;
  const double pad = 0.10 * (q[1] - q[0]);
  double lo = q[0] - pad, hi = q[1] + pad;
  if (px->GetXaxis()->GetXmin() >= 0.0)
    lo = std::max(0.0, lo);
  px->GetXaxis()->SetRangeUser(std::max(lo, px->GetXaxis()->GetXmin()),
                               std::min(hi, px->GetXaxis()->GetXmax()));
}

TH1D* projection(TH2D* h, const std::string& suffix, int ybin_lo = 1,
                 int ybin_hi = kNumClasses) {
  if (h == nullptr)
    return nullptr;
  auto* px =
      h->ProjectionX((std::string(h->GetName()) + suffix).c_str(), ybin_lo, ybin_hi);
  if (px == nullptr)
    return nullptr;

  // Rebin for the window the reader will actually SEE, not the whole axis.
  // Ranges are booked generously (a TOF time axis is +-2 ns for a ~20 ps
  // core), so the data occupies a percent or so of it. Targeting entries per
  // bin across the full axis then merged a 63-entry histogram down to about
  // four fat bars. Find the occupied window first, count how many bins fall
  // inside it, and choose the bin WIDTH there by Freedman-Diaconis,
  //
  //     w = 2 * IQR / N^(1/3)
  //
  // the standard rule for exactly this question. It reads the width off the
  // data's own interquartile spread, so a narrow core gets fine bins and a
  // broad one coarse bins, and it backs off as N^(-1/3) so a 60-entry panel is
  // not cut into forty one-count spikes. No fixed target can do both: 40 bins
  // shattered the sparse panels, and few enough bins to fix those threw away
  // the resolution on the 8000-entry ones.
  double qv[4] = {0, 0, 0, 0}, pr[4] = {0.005, 0.25, 0.75, 0.995};
  px->GetQuantiles(4, qv, pr);
  const double lo = qv[0], iqr = qv[2] - qv[1], hi = qv[3];
  if (hi > lo) {
    const double bw   = px->GetBinWidth(1);
    const double neff = std::max(1.0, px->GetEffectiveEntries());
    // A degenerate IQR -- the whole core inside one bin -- offers no scale to
    // read, so fall back to the visible span.
    const double w_fd = (iqr > 0.0) ? 2.0 * iqr / std::cbrt(neff) : (hi - lo) / 20.0;
    const int target = std::clamp(static_cast<int>((hi - lo) / std::max(w_fd, 1e-12)), 10, 60);
    const int in_view = std::max(1, static_cast<int>((hi - lo) / bw));
    int       factor  = std::max(1, in_view / target);
    // Rebin() needs a divisor of the bin count.
    while (factor > 1 && px->GetNbinsX() % factor != 0)
      --factor;
    if (factor > 1)
      px->Rebin(factor);
  }
  return px;
}

std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream is(s);
  std::string         line;
  while (std::getline(is, line))
    out.push_back(line);
  return out;
}

std::string nowString() {
  const std::time_t t = std::time(nullptr);
  char buf[64];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S %Z", std::localtime(&t));
  return buf;
}

std::string envOr(const char* k, const char* fallback) {
  const char* v = std::getenv(k);
  return (v != nullptr && *v != '\0') ? v : fallback;
}

} // namespace

bool writePdfReport(const std::string& path, const std::string& table,
                    const std::map<std::string, TH2D*>& hists, const ReportContext& ctx) {
  const bool batch_before = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);
  // octofit Report::Style::Apply(), same values.
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFont(kFont, "T");
  gStyle->SetTitleSize(0.05, "T");
  gStyle->SetLabelFont(kFont, "XYZ");
  gStyle->SetTitleFont(kFont, "XYZ");
  gStyle->SetTextFont(kFont);
  gStyle->SetStatFont(kFont);
  gStyle->SetLegendFont(kFont);
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.3, "Y");
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetEndErrorSize(3);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetTickLength(0.025, "XYZ");
  gStyle->SetGridStyle(3);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetPalette(kBird);

  TCanvas c("eb_report", "report", kCanvasW, kCanvasH);

  // ---- page 1: title -------------------------------------------------------
  c.Clear();
  {
    TLatex t;
    t.SetNDC();
    t.SetTextAlign(22);
    t.SetTextFont(kFontBold);
    t.SetTextSize(0.055);
    t.DrawLatex(0.5, 0.80, "Event Builder Trigger Benchmark");
    t.SetTextFont(kFont);
    t.SetTextSize(0.034);
    // "inline"/"standalone" alone means nothing to a reader of the PDF.
    t.DrawLatex(0.5, 0.745, "Software trigger performance");
    if (ctx.overall_eff >= 0.0) {
      std::ostringstream os;
      os << std::fixed << std::setprecision(1) << "Efficiency = " << 100.0 * ctx.overall_eff
         << "%   Purity = " << 100.0 * ctx.overall_pur << "%";
      t.SetTextSize(0.030);
      t.DrawLatex(0.5, 0.705, os.str().c_str());
    }
    TLine rule;
    rule.SetNDC();
    rule.DrawLineNDC(0.13, 0.675, 0.87, 0.675);

    std::vector<std::pair<std::string, std::string>> kv;
    kv.emplace_back("Frames:", std::to_string(ctx.frames));
    {
      std::ostringstream o;
      o << ctx.candidates << "   (real " << ctx.real << ", fake " << ctx.fake << ")";
      kv.emplace_back("Candidates:", o.str());
    }
    kv.emplace_back("Classes seen:", std::to_string(ctx.classes.size()));
    kv.emplace_back("", "");
    // Every input, not just the first: a run over many files used to report
    // one of them as if it were the whole input.
    {
      std::ostringstream n;
      n << ctx.input_files.size() << " file" << (ctx.input_files.size() == 1 ? "" : "s");
      kv.emplace_back("Input:", n.str());
    }
    const std::size_t kShow = 3;
    for (std::size_t f = 0; f < ctx.input_files.size() && f < kShow; ++f) {
      const std::string& in = ctx.input_files[f];
      // Show the tail: the leading path is common to all of them, the part
      // that identifies the file is at the end.
      // Wrap at the width that actually fits between the value column and the
      // right margin; a fixed character count overflowed the page.
      const int    wrap  = monoFit(0.46, 0.026);
      const std::string shown =
          in.size() > static_cast<std::size_t>(3 * wrap)
              ? "..." + in.substr(in.size() - static_cast<std::size_t>(3 * wrap))
              : in;
      for (std::size_t i = 0; i < shown.size(); i += wrap)
        kv.emplace_back("", shown.substr(i, wrap));
    }
    if (ctx.input_files.size() > kShow) {
      std::ostringstream m;
      m << "(+" << ctx.input_files.size() - kShow << " more)";
      kv.emplace_back("", m.str());
    }

    drawKeyValue(kv, 0.615);

    t.SetTextFont(kFont);
    t.SetTextSize(0.020);
    t.DrawLatex(0.5, 0.10,
                ("Generated by EICrecon eventbuilder_trigger: " + nowString()).c_str());
  }
  c.Print((path + "(").c_str(), "pdf");

  // ---- page 2: software provenance ----------------------------------------
  c.Clear();
  {
    const double y0 = pageTitle("Software and configuration",
                                "versions this report was produced with");
    std::vector<std::pair<std::string, std::string>> kv;
    kv.emplace_back("ROOT:", gROOT->GetVersion());
    kv.emplace_back("JANA2:", JVersion::GetVersion());
    {
      std::ostringstream os;
      os << podio::version::build_version;
      kv.emplace_back("podio:", os.str());
    }
    {
      std::ostringstream os;
      os << EDM4EIC_VERSION_MAJOR << "." << EDM4EIC_VERSION_MINOR << "."
         << EDM4EIC_VERSION_PATCH;
      kv.emplace_back("EDM4eic:", os.str());
    }
    std::string compiler_str =
#if defined(__clang__)
        "clang " __clang_version__
#elif defined(__GNUC__)
        "gcc " __VERSION__
#else
        "unknown"
#endif
        ;
    if (compiler_str.size() > 42)
      compiler_str = compiler_str.substr(0, 42) + "...";
    kv.emplace_back("Compiler:", compiler_str);
    kv.emplace_back("Built:", std::string(__DATE__) + " " + __TIME__);
    kv.emplace_back("", "");
    // The actual compact FILE, not just the config name: "Config: epic" says
    // nothing about which geometry was loaded.
    std::string geometry, geometry_label;
    {
      const std::string dpath = envOr("DETECTOR_PATH", "");
      const std::string dcfg  = ctx.detector_config.empty() ? envOr("DETECTOR_CONFIG", "")
                                                            : ctx.detector_config;
      std::string xml;
      if (!dpath.empty() && !dcfg.empty()) {
        const std::string cand = dpath + "/" + dcfg + ".xml";
        if (::access(cand.c_str(), R_OK) == 0)
          xml = cand;
      }
      // DETECTOR_PATH is the directory of this same file, so printing both
      // spent four lines saying one thing. Show the resolved file when there
      // is one, the search path only when there is not.
      geometry = xml.empty() ? dpath : xml;
      if (xml.empty() && dpath.empty())
        geometry = dcfg.empty() ? "(unset)" : dcfg + " (file not found)";
      geometry_label = xml.empty() ? "DETECTOR_PATH:" : "Geometry:";
    }
    kv.emplace_back("", "");

    double after = drawKeyValueEnd(kv, y0);

    if (!geometry.empty()) {
      TLatex t;
      t.SetNDC();
      t.SetTextFont(kFont);
      t.SetTextSize(0.026);
      t.SetTextAlign(11);
      t.DrawLatex(0.20, after, geometry_label.c_str());
      TLatex m;
      m.SetNDC();
      m.SetTextFont(kMono);
      m.SetTextSize(0.016);
      m.SetTextAlign(11);
      double y = after;
      for (const auto& line : wrapPath(geometry, monoFit(0.40, 0.016))) {
        m.DrawLatex(0.40, y, line.c_str());
        y -= 0.021;
      }
      after = std::min(after - 0.030, y - 0.012);
    }

    // The whole invocation. Listing selected parameters (this page used to
    // show only nsigma_window) hides the one flag that matters when a number
    // looks wrong.
    if (!ctx.command_line.empty()) {
      TLatex t;
      t.SetNDC();
      t.SetTextFont(kFont);
      t.SetTextSize(0.026);
      t.SetTextAlign(11);
      t.DrawLatex(0.20, after, "Command:");
      TLatex m;
      m.SetNDC();
      m.SetTextFont(kMono);
      m.SetTextSize(0.016);
      m.SetTextAlign(11);
      // Break on argument boundaries, not mid-token: a command chopped every
      // 78 characters splits flags and paths and is unreadable. Continuations
      // are indented so the argument list reads as one block.
      const int         wrap = monoFit(0.22, 0.016);
      std::vector<std::string> lines;
      {
        std::istringstream is(ctx.command_line);
        std::string        tok, cur;
        while (is >> tok) {
          if (!cur.empty() && cur.size() + 1 + tok.size() > static_cast<std::size_t>(wrap)) {
            lines.push_back(cur);
            cur.clear();
          }
          if (tok.size() > static_cast<std::size_t>(wrap)) { // a very long path
            if (!cur.empty()) {
              lines.push_back(cur);
              cur.clear();
            }
            for (std::size_t i = 0; i < tok.size(); i += wrap)
              lines.push_back(tok.substr(i, wrap));
            continue;
          }
          cur += (cur.empty() ? "" : " ") + tok;
        }
        if (!cur.empty())
          lines.push_back(cur);
      }
      double y = after - 0.030;
      for (std::size_t i = 0; i < lines.size() && y > 0.06; ++i) {
        m.DrawLatex(i == 0 ? 0.20 : 0.22, y, lines[i].c_str());
        y -= 0.021;
      }
    }
  }
  c.Print(path.c_str(), "pdf");

  // ---- per-class efficiency, as a bar chart -------------------------------
  // Traffic-light coded and annotated with the raw counts, so a low bar can be
  // told apart from a bar that is low because almost nothing was injected.
  if (!ctx.classes.empty()) {
    c.Clear();
    c.SetLeftMargin(0.30);
    c.SetRightMargin(0.06);
    const int    n   = static_cast<int>(ctx.classes.size());
    const double x0  = 0.0;
    const double x1  = 145.0; // headroom for the annotation past 100%
    // Annotations sit in one column just right of the 100% reference line
    // rather than chasing each bar's end. Ragged labels make the values hard
    // to compare down the page, and a short bar pushes its own label into the
    // middle of the plot where it reads as data.
    const double x_ann = 101.5;
    TH2D frame_h("eff_frame", ";efficiency [%];", 100, x0, x1, n, -0.5, n - 0.5);
    frame_h.SetStats(0);
    for (int i = 0; i < n; ++i)
      frame_h.GetYaxis()->SetBinLabel(i + 1, ctx.classes[i].name.c_str());
    frame_h.GetYaxis()->SetLabelSize(std::min(0.022, 0.65 / std::max(n, 1)));
    frame_h.GetYaxis()->SetLabelFont(kFont);
    frame_h.GetXaxis()->SetLabelSize(0.020);
    frame_h.GetXaxis()->SetLabelFont(kFont);
    frame_h.GetXaxis()->SetTitleSize(0.026);
    frame_h.GetXaxis()->SetTitleFont(kFont);
    frame_h.GetXaxis()->SetTitleOffset(1.3);
    c.SetGridx(); // dotted verticals, as octofit draws behind its chi2 bars
    frame_h.Draw();

    TBox   box;
    TLatex ann;
    ann.SetTextFont(kFont);
    ann.SetTextSize(std::min(0.020, 0.60 / std::max(n, 1)));
    ann.SetTextAlign(12);
    for (int i = 0; i < n; ++i) {
      const auto&  ce  = ctx.classes[i];
      const double pct = 100.0 * ce.eff;
      // green >= 90%, orange >= 60%, red below: the same three-band reading as
      // the octofit chi2/N bars.
      box.SetFillColor(pct >= 90.0 ? kGreen + 1 : (pct >= 60.0 ? kOrange + 7 : kRed + 1));
      box.DrawBox(0.0, i - 0.28, pct, i + 0.28);
      std::ostringstream os;
      os << std::fixed << std::setprecision(1) << pct << "  (" << ce.found << "/" << ce.injected
         << ")";
      ann.DrawLatex(x_ann, i, os.str().c_str());
    }
    TLine ref;
    ref.SetLineStyle(2);
    ref.DrawLine(100.0, -0.5, 100.0, n - 0.5);
    c.Print(path.c_str(), "pdf");
    c.SetGridx(0);
    c.SetLeftMargin(0.10);
    c.SetRightMargin(0.10);
  }

  // ---- fitted resolution summary, one row per subsystem ------------------
  c.Clear();
  {
    TLatex h;
    h.SetNDC();
    h.SetTextAlign(22);
    h.SetTextFont(132);
    h.SetTextSize(0.044);
    h.DrawLatex(0.5, 0.930, "Fitted resolution by subsystem");
    h.SetTextFont(kFont);
    // Spell the columns out. "time sigma" and "space 68%" assumed the reader
    // already knew which width, of what, in what unit -- and the position
    // column silently changes both quantity and unit between trackers and
    // calorimeters, which no header can carry on its own.
    h.SetTextSize(0.017);
    h.DrawLatex(0.5, 0.893, "residuals measured against MC truth: the sim hit, not the collision");
    h.SetTextColor(static_cast<Color_t>(kGray + 3));
    h.SetTextSize(0.015);
    h.DrawLatex(0.5, 0.868,
                "time = sensor resolution: core width of (t_{rec} + |r|/c - t_{sim}), "
                "TOF in ps and gaseous/silicon in ns");
    h.DrawLatex(0.5, 0.848,
                "position = 68% of hits land within this of truth: trackers in #mum, "
                "calorimeters in rad (the cluster-to-particle opening angle)");
    h.DrawLatex(0.5, 0.828,
                "energy = Gaussian width of (E_{rec} - E_{MC}) / E_{MC},  and E bias its mean "
                "(negative = reconstructed energy low)");
    h.DrawLatex(0.5, 0.810, "-- = the fit did not converge, usually for want of entries");
    h.SetTextColor(kBlack);

    const std::vector<Col> cols = {{0.09, 11}, {0.42, 31}, {0.56, 31},
                                   {0.70, 31}, {0.82, 31}, {0.93, 31}};
    const std::vector<std::string> header = {"Subsystem", "entries", "time",
                                             "position",  "energy",  "E bias"};
    std::vector<std::vector<Cell>> rows;

    auto add = [&](const std::string& det, bool calo) {
      auto grab = [&](const std::string& key) -> TH2D* {
        auto it = hists.find(key);
        return (it != hists.end()) ? it->second : nullptr;
      };
      TH2D* ht = grab("time_" + det);
      TH2D* hs = grab("space_" + det);
      TH2D* he = calo ? grab("energy_" + det) : nullptr;
      if (ht == nullptr && hs == nullptr && he == nullptr)
        return;

      double ts = 0, te = 0, sq = 0, es = 0, ee = 0, emean = 0;
      bool   ok_t = false, ok_s = false, ok_e = false;
      long   n    = 0;
      if (ht != nullptr && ht->GetEntries() > 0) {
        ok_t = coreSigma(projection(ht, "_t_fit"), ts, te);
        n    = static_cast<long>(ht->GetEntries());
      }
      if (hs != nullptr && hs->GetEntries() > 0)
        ok_s = quantile68(projection(hs, "_s_fit"), sq);
      if (he != nullptr && he->GetEntries() > 0) {
        auto* px = projection(he, "_e_fit");
        ok_e     = coreSigma(px, es, ee);
        emean    = px->GetMean();
      }
      // Green where a fit converged, grey where it did not: a blank cell then
      // reads as "not resolved" rather than as a missing detector.
      const auto ok_c   = static_cast<Color_t>(kGreen + 2);
      const auto none_c = static_cast<Color_t>(kGray + 2);
      // Thin statistics are the usual reason a fit is blank, so flag the count
      // itself rather than leaving the reader to guess.
      const auto n_c = static_cast<Color_t>((n > 0 && n < 500) ? kOrange + 7 : kBlack);
      // Units come off each histogram's own axis, so a row cannot claim a unit
      // its histogram does not use.
      rows.push_back({{det},
                      {std::to_string(n), n_c},
                      {fmtUnit(ts, ok_t, axisUnit(ht)), ok_t ? ok_c : none_c},
                      {fmtUnit(sq, ok_s, axisUnit(hs)), ok_s ? ok_c : none_c},
                      {fmtPercent(es, ok_e), ok_e ? ok_c : none_c},
                      {fmtPercent(emean, ok_e), ok_e ? ok_c : none_c}});
    };

    for (const auto& d : ResolutionHists::trackerNames())
      add(d, false);
    rows.push_back({});
    for (const auto& d : ResolutionHists::caloNames())
      add(d, true);

    {
      // Shrink to fit, octofit-style: dY = min(baseDY, available / rows).
      // Lower than kTableYTop: this page carries a five-line legend above the
      // table, and the header row has to clear it.
      const double avail  = kResTableYTop - kTableYBot;
      const double dY = std::min(kTableBaseDY, avail / std::max<double>(rows.size() + 12, 1));
      const double sz = std::min(kTableTextSize * (dY / kTableBaseDY),
                                 fitTextSize(static_cast<int>(cols.size()), 12, 0.88));
      drawTable(cols, header, rows, kResTableYTop, dY, sz, 0.075, 0.945);
    }

    // ---- candidate t0: measured spread against the quoted uncertainty ------
    auto grab = [&](const std::string& k) -> TH2D* {
      auto it = hists.find(k);
      return (it != hists.end()) ? it->second : nullptr;
    };
    TH2D* hr = grab("t0_residual");
    if (hr != nullptr && hr->GetEntries() > 0) {
      double rs = 0, re = 0, ps = 0, pe = 0, mean_sig = 0;
      const bool ok_r = coreSigma(projection(hr, "_t0r_fit"), rs, re);
      if (TH2D* hs2 = grab("t0_sigma"); hs2 != nullptr && hs2->GetEntries() > 0)
        mean_sig = projection(hs2, "_t0s_fit")->GetMean();
      bool ok_p = false;
      if (TH2D* hp = grab("t0_pull"); hp != nullptr && hp->GetEntries() > 0)
        ok_p = coreSigma(projection(hp, "_t0p_fit"), ps, pe);

      // Anchor to the row height actually used above; with octofit's taller
      // rows a stale constant put this block on top of the table.
      const double avail_r = kResTableYTop - kTableYBot;
      const double dY_r    = std::min(kTableBaseDY, avail_r / std::max<double>(rows.size() + 12, 1));
      const double ybase   = kResTableYTop - dY_r * static_cast<double>(rows.size() + 3);
      h.SetTextFont(132);
      h.SetTextAlign(11);
      h.SetTextSize(0.020);
      h.DrawLatex(0.09, ybase, "Candidate t0");

      const std::vector<Col>         c2  = {{0.11, 11}, {0.62, 31}, {0.72, 11}};
      const std::vector<std::string> hd2 = {"", "", ""};
      std::vector<std::vector<Cell>> r2;
      // rs / mean_sig are already picoseconds: the t0 histograms are booked
      // in ps, so no conversion here (there used to be a x1000).
      r2.push_back({{"measured   sigma(t0 - t_MC)"},
                    {fmt(rs, 1, ok_r), static_cast<Color_t>(ok_r ? kGreen + 2 : kGray + 2)},
                    {"ps"}});
      r2.push_back({{"reported   mean delta_t0 / nsigma"},
                    {fmt(mean_sig, 1, mean_sig > 0), kBlack},
                    {"ps"}});
      // A calibrated uncertainty gives a unit-width pull. Colour says which way
      // it is wrong without the reader having to remember the convention.
      const auto pull_c = static_cast<Color_t>(
          !ok_p ? kGray + 2 : (ps > 1.25 ? kRed + 1 : (ps < 0.8 ? kOrange + 7 : kGreen + 2)));
      r2.push_back({{"pull       sigma[(t0-t_MC)/delta_t0]"},
                    {fmt(ps, 3, ok_p), pull_c},
                    {ok_p ? (ps > 1.25 ? "underestimated" : (ps < 0.8 ? "conservative" : "calibrated"))
                          : ""}});
      {
        std::ostringstream note;
        note << "weights[T0SIGMA] is a " << std::fixed << std::setprecision(1)
             << ctx.nsigma_window << "-sigma half-width; divided out above";
        r2.push_back({});
        r2.push_back({{note.str(), static_cast<Color_t>(kGray + 2)}});
      }
      drawTable(c2, hd2, r2, ybase - 0.030, dY_r, kTableTextSize * (dY_r / kTableBaseDY),
                0.075, 0.945);
    }
  }
  c.Print(path.c_str(), "pdf");

  // ---- the per-class statistics table --------------------------------------
  // Rendered as cells, not as dumped monospace lines: only that way do the
  // columns align and the numbers carry colour.
  {
    auto header = ctx.table.first;
    auto rows   = ctx.table.second;
    // Drop columns with nothing in them (gnn when no model ran). They are pure
    // width, and width is why the last real column fell off the page.
    for (std::size_t col = header.size(); col-- > 3;) {
      bool any = false;
      for (const auto& r : rows)
        if (col < r.size() && r[col].rfind("--", 0) != 0 && r[col] != "n/a" && !r[col].empty()) {
          any = true;
          break;
        }
      if (any)
        continue;
      header.erase(header.begin() + static_cast<long>(col));
      for (auto& r : rows)
        if (col < r.size())
          r.erase(r.begin() + static_cast<long>(col));
    }
    // Column x positions: class name left-aligned, everything else right.
    std::vector<Col> cols;
    cols.push_back({0.050, 11}); // class name, left
    const std::size_t nstage = header.size() - 4; // class, inj, found, fake
    cols.push_back({0.165, 31});
    cols.push_back({0.225, 31});
    cols.push_back({0.285, 31});
    for (std::size_t i = 0; i < nstage; ++i)
      cols.push_back({0.325 + (0.965 - 0.325) * double(i + 1) / double(nstage), 31});

    constexpr std::size_t kPerPage = 34;
    std::size_t           page     = 0;
    for (std::size_t i = 0; i < rows.size(); i += kPerPage, ++page) {
      c.Clear();
      TLatex h;
      h.SetNDC();
      h.SetTextAlign(22);
      h.SetTextFont(132);
      h.SetTextSize(0.044);
      h.DrawLatex(0.5, 0.930, page == 0 ? "Per-class trigger performance"
                                        : "Per-class trigger performance (cont.)");
      h.SetTextFont(kFont);
      h.SetTextSize(0.028);
      h.DrawLatex(0.5, 0.888,
                  "each column is a stage: efficiency / purity, in percent");

      // Second header line naming what the pair means, so a bare number is
      // never ambiguous about which stage it belongs to.
      std::vector<std::vector<Cell>> body;
      {
        std::vector<Cell> sub;
        // Counts are counts: labelling them eff/pur too said inj and found
        // were percentages.
        sub.push_back({""});
        sub.push_back({"count", static_cast<Color_t>(kGray + 2)});
        sub.push_back({"x/inj", static_cast<Color_t>(kGray + 2)});
        sub.push_back({"count", static_cast<Color_t>(kGray + 2)});
        for (std::size_t k = 0; k < header.size() - 4; ++k)
          sub.push_back({"eff/pur", static_cast<Color_t>(kGray + 2)});
        body.push_back(sub);
      }
      for (std::size_t r = i; r < rows.size() && r < i + kPerPage; ++r) {
        if (rows[r].empty()) {
          body.push_back({});
          continue;
        }
        std::vector<Cell> cells;
        for (std::size_t k = 0; k < rows[r].size(); ++k) {
          Color_t col = kBlack;
          if (k >= 3) {
            // Colour the efficiency half: green >= 90, orange >= 60, red below,
            // grey when the stage has nothing to say for this class.
            const std::string& t = rows[r][k];
            if (t.rfind("--", 0) == 0)
              col = static_cast<Color_t>(kGray + 2);
            else {
              const double v = std::atof(t.c_str());
              col            = static_cast<Color_t>(v >= 90.0 ? kGreen + 2
                                                    : (v >= 60.0 ? kOrange + 7 : kRed + 1));
            }
          }
          cells.push_back({rows[r][k], col});
        }
        body.push_back(cells);
      }
      const double avail = kTableYTop - kTableYBot;
      // Budget for what is drawn AFTER the table too -- the FAR / recovery /
      // in-acceptance / BLIND footer and its note all use this same dY. Counting
      // only the rows pushed the last footer lines below yBot and off the page
      // (octofit's ParameterTablePage carries a comment about hitting exactly
      // this with its trailing Minuit commands).
      const double footerUnits = ctx.footer.empty() ? 0.0 : 1.5 + 0.9 * ctx.footer.size();
      const double dY =
          std::min(kTableBaseDY, avail / std::max<double>(body.size() + 4 + footerUnits, 1));
      // "100.0/100.0" is the widest cell; 10 columns of it will not fit at
      // octofit's body size, so take the smaller of the row- and width-limited
      // sizes instead of overlapping the columns.
      const double sz = std::min(kTableTextSize * (dY / kTableBaseDY),
                                 fitTextSize(static_cast<int>(cols.size()), 11, 0.93));
      drawTable(cols, header, body, kTableYTop, dY, sz, 0.045, 0.975);

      // Footer facts once, on the last page of the table.
      // Set as a table, not as prose. As free text the values never lined up
      // under each other and the long lines had to be ellipsised, which cut
      // "dt = 42.4 ns" off exactly the line that needed it.
      if (i + kPerPage >= rows.size() && !ctx.footer.empty()) {
        // Anchored to the rows actually drawn: a fixed constant here put the
        // footer on top of the last few class rows once the rows grew.
        const double y0 = kTableYTop - dY * static_cast<double>(body.size() + 2.5);
        const std::vector<Col>         fcols = {{0.050, 11}, {0.400, 31}, {0.425, 11}};
        const std::vector<std::string> fhdr  = {"", "", ""};
        std::vector<std::vector<Cell>> frows;
        for (const auto& [k, v, det] : ctx.footer)
          frows.push_back({{k}, {v}, {det, static_cast<Color_t>(kGray + 2)}});
        drawTable(fcols, fhdr, frows, y0, dY * 0.85, std::min(0.017, sz * 1.1), 0.045, 0.975);
      }
      c.Print(path.c_str(), "pdf");
    }
  }

  // ---- per-factory timing, when profiling is on -------------------------
  if (!ctx.profile.empty()) {
    c.Clear();
    TLatex h;
    h.SetNDC();
    h.SetTextAlign(22);
    h.SetTextFont(kFontBold);
    h.SetTextSize(0.044);
    h.DrawLatex(0.5, 0.930, "Where the time goes");
    h.SetTextFont(kFont);
    h.SetTextSize(0.028);
    h.DrawLatex(0.5, 0.888, "per-factory self time -- callees subtracted -- from JANA's call graph");

    double total = 0.0;
    for (const auto& [n, sec, calls] : ctx.profile)
      total += sec;

    const std::vector<Col> cols = {{0.07, 11}, {0.62, 31}, {0.74, 31}, {0.86, 31}, {0.95, 31}};
    const std::vector<std::string> hdr = {"Factory", "total s", "share", "calls", "ms/call"};
    std::vector<std::vector<Cell>> rows;
    std::size_t shown = 0;
    for (const auto& [name, sec, calls] : ctx.profile) {
      if (shown++ >= 28)
        break;
      const double share = (total > 0.0) ? 100.0 * sec / total : 0.0;
      std::ostringstream a, b, cc, d;
      a << std::fixed << std::setprecision(3) << sec;
      b << std::fixed << std::setprecision(1) << share << "%";
      cc << calls;
      d << std::fixed << std::setprecision(3) << (calls ? 1000.0 * sec / calls : 0.0);
      // Red where one factory dominates: that is the thing to look at.
      const auto col = static_cast<Color_t>(share >= 20.0 ? kRed + 1
                                            : (share >= 5.0 ? kOrange + 7 : kBlack));
      rows.push_back({{name}, {a.str(), col}, {b.str(), col}, {cc.str()}, {d.str()}});
    }
    const double avail = kTableYTop - kTableYBot;
    const double dY    = std::min(kTableBaseDY, avail / std::max<double>(rows.size() + 4, 1));
    drawTable(cols, hdr, rows, kTableYTop, dY,
              std::min(kTableTextSize * (dY / kTableBaseDY),
                       fitTextSize(static_cast<int>(cols.size()), 14, 0.90)),
              0.06, 0.96);

    double yAfter = kTableYTop - dY * (rows.size() + 3.0);

    {
      TLatex n;
      n.SetNDC();
      n.SetTextFont(kMono);
      n.SetTextSize(0.016);
      n.SetTextColor(static_cast<Color_t>(kGray + 2));
      n.DrawLatex(0.07, std::max(0.04, yAfter),
                  "RECORD_CALL_STACK was on: absolute throughput here is not representative");
    }
    c.Print(path.c_str(), "pdf");

    // Cost per physics class, as its own table: which classes are expensive
    // is a different question from which factories are, and merging them into
    // one table would make both unreadable. It also needs a page of its own --
    // the factory list is long enough to leave no room, and the two collided.
    if (!ctx.class_profile.empty()) {
      c.Clear();
      yAfter = kTableYTop;
      TLatex ch;
      ch.SetNDC();
      ch.SetTextAlign(22);
      ch.SetTextFont(kFontBold);
      ch.SetTextSize(0.032);
      ch.DrawLatex(0.5, 0.960, "Cost per physics class");
      ch.SetTextFont(kFont);
      ch.SetTextSize(0.020);
      ch.DrawLatex(0.5, 0.932, "frame wall time charged to the classes the frame carried");

      const std::vector<Col> ccols = {{0.09, 11}, {0.46, 31}, {0.62, 31}, {0.78, 31}};
      const std::vector<std::string> chdr = {"Class", "frames", "total s", "ms/frame"};
      std::vector<std::vector<Cell>> crows;
      double cmax = 0.0;
      for (const auto& [n2, sec, cands] : ctx.class_profile)
        cmax = std::max(cmax, cands ? 1000.0 * sec / cands : 0.0);
      for (const auto& [name, sec, cands] : ctx.class_profile) {
        const double per = cands ? 1000.0 * sec / cands : 0.0;
        std::ostringstream a, b, cc;
        a << cands;
        b << std::fixed << std::setprecision(2) << sec;
        cc << std::fixed << std::setprecision(1) << per;
        // Flag the classes costing most per candidate.
        const auto col = static_cast<Color_t>(
            (cmax > 0.0 && per >= 0.75 * cmax) ? kRed + 1
                                               : ((cmax > 0.0 && per >= 0.5 * cmax) ? kOrange + 7
                                                                                    : kBlack));
        crows.push_back({{name}, {a.str()}, {b.str()}, {cc.str(), col}});
      }
      const double cdY = std::min(kTableBaseDY * 0.8,
                                  std::max(0.012, (yAfter - 0.10) / (crows.size() + 3.0)));
      drawTable(ccols, chdr, crows, yAfter - 0.030, cdY,
                std::min(kTableTextSize * (cdY / kTableBaseDY),
                         fitTextSize(static_cast<int>(ccols.size()), 14, 0.80)),
                0.06, 0.90);
      yAfter -= 0.030 + cdY * (crows.size() + 3.0);
      TLatex n;
      n.SetNDC();
      n.SetTextFont(kMono);
      n.SetTextSize(0.016);
      n.SetTextColor(static_cast<Color_t>(kGray + 2));
      n.DrawLatex(0.07, std::max(0.04, yAfter),
                  "a frame's time is charged in full to every class it carried, so the "
                  "column sums past the wall time");
      c.Print(path.c_str(), "pdf");
    }
  }

  // ---- residual grid: one row per detector, one column per quantity ------
  // The panel title names the detector and the quantity, so it has to be
  // readable; the style default is smaller than the axis numbers beneath it.
  const float saved_title_size = gStyle->GetTitleFontSize();
  gStyle->SetTitleFontSize(0.105);
  // Trackers and calo get their own pages: they do not share quantities (dx/dy
  // vs dR/energy), and interleaving them made a reader hunt for a detector.
  {
    struct Panel { std::string key; };
    // `ybin` selects the class-axis row: the all-real aggregate, or the fake
    // one. Each detector contributes two Rows, drawn one above the other.
    struct Row { std::string name; std::vector<std::string> keys; int ybin; };

    // Below this a panel is a few counts smeared over the whole gate: no fit is
    // meaningful and the bars do not even render at this pad size.
    constexpr double kMinEntries = 20;
    auto drawGrid = [&](const std::string& title, const std::vector<std::string>& colnames,
                        const std::vector<Row>& all_rows) {
      // Drop detectors with nothing drawable rather than spending a sixth of
      // the page on three "no data" boxes in a row.
      const auto binEntries = [&](const std::string& k, int ybin) {
        auto i = hists.find(k);
        if (i == hists.end() || i->second == nullptr)
          return 0.0;
        return i->second->Integral(1, i->second->GetNbinsX(), ybin, ybin);
      };
      std::vector<Row> rows;
      for (const auto& r : all_rows) {
        // Count in THIS row's class-axis bin. Judging a row by the whole
        // histogram kept detectors alive whose every panel then blanked --
        // a page of "no data" boxes under a real title.
        const bool any = std::any_of(r.keys.begin(), r.keys.end(), [&](const std::string& k) {
          return binEntries(k, r.ybin) >= kMinEntries;
        });
        if (any)
          rows.push_back(r);
      }
      if (rows.empty())
        return;

      // ROOT shrinks a pad title that does not fit its box, so "TOFBarrel dt"
      // came out large and "BackwardMPGDEndcap dt" small on the same page.
      // Size every title from the LONGEST label instead, so a page reads as
      // one set of plots rather than a dozen different ones.
      std::size_t longest = 1;
      for (const auto& r : rows)
        for (const auto& q : colnames)
          longest = std::max(longest, r.name.size() + 1 + q.size());
      gStyle->SetTitleFontSize(
          std::clamp(0.105 * 17.0 / static_cast<double>(longest), 0.050, 0.105));
      // ONE DETECTOR PER ROW. Packing two detectors side by side fitted each
      // system onto a single page, but at six panels across nothing could be
      // read off them: the point of these pages is the width of a peak, and a
      // pad too small to show the peak's shape does not carry it. Four rows a
      // page keeps the pads square and large; a system spans as many pages as
      // it needs.
      constexpr std::size_t kMaxRows = 4;
      const std::size_t     ncol     = colnames.size();
      const std::size_t     blocks   = 1;
      const std::size_t     ntot     = ncol;
      for (std::size_t r0 = 0; r0 < rows.size(); r0 += kMaxRows) {
        const std::size_t nrow = std::min(kMaxRows, rows.size() - r0);
        c.Clear();
        TLatex pt;
        pt.SetNDC();
        pt.SetTextAlign(22);
        pt.SetTextFont(kFontBold);
        pt.SetTextSize(0.032);
        pt.DrawLatex(0.5, 0.960, title.c_str());
        pt.SetTextFont(kFont);
        pt.SetTextSize(0.020);
        pt.DrawLatex(0.5, 0.932, "curve = fitted tail (wide) component only");

        const double x0 = 0.070, x1 = 0.980, yTop = 0.900, yBot = 0.030;
        const double pw = (x1 - x0) / ntot;
        // Divide by the page capacity, not by how many rows this page got:
        // a last page with two rows was stretching them to half a page each.
        const double ph = (yTop - yBot) / kMaxRows;
        std::vector<TPad*> pads;
        for (std::size_t r = 0; r < nrow; ++r) {
          for (std::size_t k = 0; k < ntot; ++k) {
            const double px0 = x0 + pw * k;
            const double py1 = yTop - ph * r;
            auto* pad = new TPad(("g" + std::to_string(r0 + r) + "_" + std::to_string(k)).c_str(),
                                 "", px0, py1 - ph, px0 + pw, py1);
            pad->SetLeftMargin(0.17);
            pad->SetBottomMargin(0.28); // the centred axis title sits here
            pad->SetRightMargin(0.04);
            pad->SetTopMargin(0.15);
            pad->Draw();
            pads.push_back(pad);
          }
        }
        std::size_t ip = 0;
        for (std::size_t r = 0; r < nrow; ++r) {
          for (std::size_t k = 0; k < ntot; ++k, ++ip) {
            pads[ip]->cd();
            const std::size_t q  = k % ncol;                 // quantity within a block
            const std::size_t di = r0 + (k / ncol) * nrow + r; // detector this pad shows
            if (di >= rows.size())
              continue;                                      // ragged last block
            // Say WHY a panel is blank rather than leaving dead space.
            auto blank = [&] {
              TLatex e;
              e.SetNDC();
              e.SetTextAlign(22);
              e.SetTextFont(kFont);
              e.SetTextSize(0.13);
              e.SetTextColor(static_cast<Color_t>(kGray + 1));
              e.DrawLatex(0.5, 0.55, "no data");
              e.SetTextSize(0.10);
              e.DrawLatex(0.5, 0.35, (rows[di].name + " " + colnames[q]).c_str());
            };
            auto it = hists.find(rows[di].keys[q]);
            if (it == hists.end() || it->second == nullptr ||
                binEntries(rows[di].keys[q], rows[di].ybin) < kMinEntries) {
              blank();
              continue;
            }
            auto* px = projection(it->second, "_g" + std::to_string(rows[di].ybin),
                                  rows[di].ybin, rows[di].ybin);
            // The parent can hold entries this slice does not: a shared-hit TH2
            // counts every class row, so a per-detector projection of it came
            // out empty and drew a bare +-2 ns frame with no bars in it.
            // A peak of one or two counts is a flat smear, not a distribution:
            // the shared-hit time residual is uniform across the whole
            // coincidence gate by construction and rendered as a bare frame.
            if (px == nullptr || px->GetEntries() < kMinEntries || px->GetMaximum() < 3) {
              blank();
              continue;
            }
            px->SetTitle((rows[di].name + " " + colnames[q]).c_str());
            px->SetLineWidth(2);
            const Color_t dcol = detectorColor(rows[di].name);
            px->SetLineColor(dcol);
            px->SetFillColorAlpha(dcol, 0.35);
            // Default ~10 divisions collide at this pad size.
            px->GetXaxis()->SetNdivisions(505);
            px->GetYaxis()->SetNdivisions(505);
            // Times Roman (132) everywhere, matching the body text; ROOT's
            // default 42 is Helvetica and read as a different document.
            px->GetXaxis()->SetLabelFont(kFont);
            px->GetXaxis()->SetTitleFont(kFont);
            px->GetYaxis()->SetLabelFont(kFont);
            px->GetYaxis()->SetTitleFont(kFont);
            px->GetXaxis()->SetLabelSize(0.070);
            px->GetXaxis()->SetTitleSize(0.085);
            px->GetXaxis()->SetTitleOffset(1.20);
            px->GetXaxis()->CenterTitle(true);
            px->GetYaxis()->SetLabelSize(0.065);
            px->GetYaxis()->CenterTitle(true);
            // NOT px->SetTitleSize(): TH1::SetTitleSize defaults to axis "X",
            // so it silently resized the axis title -- pushing it off the pad
            // -- and left the panel title at the style default. The pad title
            // is a gStyle property, set once around the grid below.
            zoomToData(px);
            px->Draw("hist");

            // A residual booked over [0, x] is a distance, not a signed
            // residual: it has no Gaussian core about zero, and fitting one
            // gave an upward-arcing "tail" through the bulk. Quote the 68%
            // containment radius instead, the standard for such a quantity.
            if (it->second->GetXaxis()->GetXmin() >= 0.0) {
              double r68 = 0, pr68[1] = {0.68};
              px->GetQuantiles(1, &r68, pr68);
              TLine q68;
              q68.SetLineColor(kBlue + 2);
              q68.SetLineStyle(2);
              // The marker line carries it; the number itself belongs on the
              // fitted-resolution page, with its unit and its uncertainty.
              q68.DrawLine(r68, 0, r68, px->GetMaximum() * 1.05);
              continue;
            }

            double fs = 0, ts = 0, mu = 0, ta = 0;
            if (TF1* f = fitDoubleGaussian(px, fs, ts, mu, ta); f != nullptr) {
              if (ts > 0.0 && ta > 0.0) {
                auto* bg = new TF1((std::string(px->GetName()) + "_tl").c_str(), "gaus",
                                   px->GetXaxis()->GetXmin(), px->GetXaxis()->GetXmax());
                bg->SetParameters(ta, mu, ts);
                bg->SetLineColor(kBlack);
                bg->SetLineWidth(2);
                bg->SetNpx(400);
                bg->Draw("same");
              }
              TLine ln;
              ln.SetLineColor(kBlack);
              ln.DrawLine(0, 0, 0, px->GetMaximum() * 1.05);
              ln.SetLineColor(kBlue + 2);
              ln.SetLineStyle(2);
              for (int sg = -1; sg <= 1; sg += 2) {
                const double xv = mu + sg * fs;
                if (xv > px->GetXaxis()->GetXmin() && xv < px->GetXaxis()->GetXmax())
                  ln.DrawLine(xv, 0, xv, px->GetMaximum() * 1.05);
              }
              // No printed width here. "core 8.00 / tail >range" was two
              // unitless numbers and a qualifier about the fit range, sitting
              // where the title should be; the widths belong on the fitted
              // resolution page, which carries their units. The drawn curve
              // and the +-sigma markers say the same thing graphically.
            }
          }
        }
        c.Print(path.c_str(), "pdf");
      }
    };

    // Signal+background first, the fake row directly beneath it, so the two
    // are read against each other for one detector before the eye moves on.
    // ROOT bins are 1-based, hence the +1.
    constexpr int kBinReal = ResolutionHists::kAllReal + 1;
    constexpr int kBinFake = ResolutionHists::kFake + 1;

    std::vector<Row> trk;
    for (const auto& d : ResolutionHists::trackerNames()) {
      trk.push_back({d, {"time_" + d, "dx_" + d, "dy_" + d}, kBinReal});
      trk.push_back({d + " (fake)", {"time_" + d, "dx_" + d, "dy_" + d}, kBinFake});
    }
    drawGrid("Tracker residuals", {"#Deltat", "#Deltax", "#Deltay"}, trk);

    std::vector<Row> cal;
    for (const auto& d : ResolutionHists::caloNames()) {
      cal.push_back({d, {"time_" + d, "space_" + d, "energy_" + d}, kBinReal});
      cal.push_back({d + " (fake)", {"time_" + d, "space_" + d, "energy_" + d}, kBinFake});
    }
    drawGrid("Calorimeter residuals", {"#Deltat", "#DeltaR", "#DeltaE/E"}, cal);
    gStyle->SetTitleFontSize(saved_title_size);
  }

  c.Clear();
  c.Print((path + ")").c_str(), "pdf"); // close the document
  gROOT->SetBatch(batch_before);
  return true;
}

} // namespace eicrecon::eb

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "ResolutionHists.h"

#include <JANA/JApplication.h>
#include <JANA/Services/JGlobalRootLock.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TH2D.h>

#include "factories/eventbuilder/TruthClassLabels.h"
#include "services/rootfile/RootFile_service.h"

namespace eicrecon::eb {

const std::vector<std::string>& ResolutionHists::trackerNames() {
  static const std::vector<std::string> n = {
      "TOFBarrel",      "TOFEndcap",          "MPGDBarrel",         "OuterMPGDBarrel",
      "BackwardMPGDEndcap", "ForwardMPGDEndcap", "SiBarrelVertex",  "SiBarrelTracker",
      "SiEndcapTracker", "B0Tracker",
      // The TOF charge-sharing chain, measured separately: these are the pads
      // ACTS actually consumes (SiliconChargeSharing -> LGADHitClustering ->
      // Measurement2D), whereas indices 0/1 above are the plain one-hit-per-
      // deposit chain that carries the GNN's per-hit truth.
      "TOFBarrelShared", "TOFEndcapShared"};
  return n;
}

const std::vector<std::string>& ResolutionHists::caloNames() {
  static const std::vector<std::string> n = {
      "B0ECal",       "EcalBarrel",         "EcalEndcapN", "EcalEndcapP", "HcalBarrel",
      "HcalEndcapPInsert", "LFHCAL",        "EcalFarForwardZDC", "HcalFarForwardZDC"};
  return n;
}

void ResolutionHists::init(JApplication* app, bool want_anyway) {
  // Only book when a histogram file was actually requested. RootFile_service
  // creates the file on demand, so asking for the directory unconditionally
  // would produce an empty ROOT file on every run.
  std::string hists_file;
  app->SetDefaultParameter("histsfile", hists_file);
  if (hists_file.empty()) {
    // No ROOT file asked for, but a PDF still needs the plots: book them
    // detached from any directory so nothing is written to disk.
    m_in_memory = want_anyway;
    return;
  }

  auto rootfile_svc = app->GetService<RootFile_service>();
  auto lock         = app->GetService<JGlobalRootLock>();
  lock->acquire_write_lock();
  TDirectory* top = rootfile_svc->GetHistFile();
  if (top != nullptr)
    m_dir = top->mkdir("benchmark_eventbuilder");
  lock->release_lock();
}

TH2D* ResolutionHists::get(const std::string& key, const std::string& title,
                           const std::string& xlabel, int nbins, double lo, double hi) {
  auto it = m_hists.find(key);
  if (it != m_hists.end())
    return it->second;
  TDirectory* prev = gDirectory;
  if (m_dir != nullptr)
    m_dir->cd();
  // y axis is the class index: one bin per physics class, labelled, so the
  // projection for a single class is a plain ProjectionX on a named bin.
  auto* h = new TH2D(key.c_str(), (title + ";" + xlabel + ";class").c_str(), nbins, lo, hi,
                     kNumClasses, -0.5, kNumClasses - 0.5);
  if (m_dir == nullptr)
    h->SetDirectory(nullptr); // in-memory only; never written
  for (int c = 0; c < kNumClasses; ++c)
    h->GetYaxis()->SetBinLabel(c + 1, std::string(className(c)).c_str());
  if (prev != nullptr)
    prev->cd();
  m_hists.emplace(key, h);
  return h;
}

void ResolutionHists::fillTrackerTime(std::size_t det, int cls, double dt_ns) {
  if (!enabled() || det >= trackerNames().size() || cls < 0)
    return;
  // Binning has to match the sensor, not the gate. AC-LGAD TOF resolves ~30 ps;
  // at the 0.5 ns/bin a common +-50 ns axis would give, its entire peak falls
  // inside ONE bin and any fitted sigma just measures the coincidence window.
  // TOF therefore gets a narrow, finely binned axis; the ~10 ns MPGDs and the
  // slow silicon keep a wide one.
  // Bin width tracks the expected resolution, not the axis: ~2.5 ps for the
  // AC-LGAD TOF (30 ps core spans ~12 bins) and 0.15 ns for the ~10 ns gaseous
  // and silicon detectors (~65 bins across the core). Finer than that just
  // spreads limited statistics over empty bins.
  const bool   tof   = (det < 2 || det >= 10);
  // TOF is plotted in PICOSECONDS, the unit the device is specified in: an
  // AC-LGAD resolves ~30 ps, which on a nanosecond axis prints as "0.0300"
  // and reads as noise. The slow gaseous and silicon detectors stay in ns,
  // where their ~10 ns spread is the natural unit. The axis title carries the
  // unit and every consumer (the fit label, the summary table) lifts it from
  // there, so the two families cannot be silently mixed up.
  //
  // TOF gets +-5000 ps, not +-1000: referenced to the TRUE collision time
  // rather than to t0, any residual per-layer offset the r/c alignment did
  // not remove shifts the peak off zero, and a tight window would push the
  // whole distribution into overflow. 10 ps bins still resolve a 30 ps core.
  //
  // NOTE: this residual carries the particle's TIME OF FLIGHT (~2-6 ns in the
  // barrel), so for TOF it is not a sensor resolution and the core fit finds
  // nothing to fit -- the summary table reads "--" for both TOF rows. That is
  // honest, not a binning failure. Referencing the SIM HIT time instead was
  // tried and does NOT fix it: measured on TOFBarrel, t_rec - t_sim has a
  // -2.95 ns median with a ~1 ns spread (866 ps over the central 90%), i.e. a
  // systematic offset far larger than the 25 ps the reconstruction itself
  // quotes as timeError -- and only 469 of 3886 TOF hits carry a sim link, so
  // it also throws away 88% of the statistics. Getting a ps-scale number out
  // of this needs the offset understood first; see the README.
  const double range = tof ? 2000.0 : 45.0; // ps for TOF, ns for the rest
  const int    nbins = tof ? 2000 : 600;    // 2 ps  |  0.15 ns per bin
  get("time_" + trackerNames()[det], trackerNames()[det] + " hit time residual",
      tof ? "t_{rec} - t_{sim} [ps]" : "t_{rec} - t_{sim} [ns]", nbins, -range, range)
      ->Fill(tof ? dt_ns * 1000.0 : dt_ns, cls);
}

void ResolutionHists::fillTrackerSpace(std::size_t det, int cls, double dr_mm) {
  if (!enabled() || det >= trackerNames().size() || cls < 0)
    return;
  // MICRONS, not mm: every one of these residuals is sub-millimetre, and a mm
  // axis prints the interesting ones as a string of leading zeros.
  // MAPS silicon and B0 resolve to tens of microns; TOF strips and the gaseous
  // MPGDs to a few hundred. A common 0-10 mm axis would put every silicon hit
  // in the first bin or two.
  const bool   fine  = (det >= 6 && det < 10); // Si + B0; TOF/MPGD are coarser
  const double range = fine ? 500.0 : 6000.0;  // um
  const int    nbins = fine ? 250 : 300;       // 2 um  |  20 um per bin
  get("space_" + trackerNames()[det], trackerNames()[det] + " in-plane position residual",
      "in-plane |rec - sim| [#mum]", nbins, 0.0, range)
      ->Fill(dr_mm * 1000.0, cls);
}

void ResolutionHists::fillTrackerSpaceRadial(std::size_t det, int cls, double dr_mm) {
  if (!enabled() || det >= trackerNames().size() || cls < 0)
    return;
  const bool   fine  = (det >= 6 && det < 10);
  const double range = fine ? 500.0 : 6000.0; // um, as for the in-plane residual
  const int    nbins = fine ? 500 : 600;      // 2 um  |  20 um per bin
  get("radial_" + trackerNames()[det], trackerNames()[det] + " radial position residual",
      "r_{rec} - r_{sim} [#mum]", nbins, -range, range)
      ->Fill(dr_mm * 1000.0, cls);
}

void ResolutionHists::fillTrackerX(std::size_t det, int cls, double dx_mm) {
  if (!enabled() || det >= trackerNames().size() || cls < 0)
    return;
  const bool   fine  = (det >= 6 && det < 10);
  const double range = fine ? 0.5 : 6.0;
  get("dx_" + trackerNames()[det], trackerNames()[det] + " #Deltax",
      "x_{rec} - x_{sim} [mm]", fine ? 500 : 600, -range, range)
      ->Fill(dx_mm, cls);
}

void ResolutionHists::fillTrackerY(std::size_t det, int cls, double dy_mm) {
  if (!enabled() || det >= trackerNames().size() || cls < 0)
    return;
  const bool   fine  = (det >= 6 && det < 10);
  const double range = fine ? 0.5 : 6.0;
  get("dy_" + trackerNames()[det], trackerNames()[det] + " #Deltay",
      "y_{rec} - y_{sim} [mm]", fine ? 500 : 600, -range, range)
      ->Fill(dy_mm, cls);
}

void ResolutionHists::fillCaloTime(std::size_t sys, int cls, double dt_ns) {
  if (!enabled() || sys >= caloNames().size() || cls < 0)
    return;
  get("time_" + caloNames()[sys], caloNames()[sys] + " cluster time residual",
      "t_{cluster} - t_{MC} [ns]", 800, -40.0, 40.0)
      ->Fill(dt_ns, cls);
}

void ResolutionHists::fillCaloAngle(std::size_t sys, int cls, double dr) {
  if (!enabled() || sys >= caloNames().size() || cls < 0)
    return;
  get("space_" + caloNames()[sys], caloNames()[sys] + " cluster angular residual",
      "#DeltaR(cluster, MC) [rad]", 250, 0.0, 0.5)
      ->Fill(dr, cls);
}

void ResolutionHists::fillCaloEnergy(std::size_t sys, int cls, double frac) {
  if (!enabled() || sys >= caloNames().size() || cls < 0)
    return;
  get("energy_" + caloNames()[sys], caloNames()[sys] + " cluster energy residual",
      "(E_{rec} - E_{MC}) / E_{MC}", 200, -1.0, 1.0)
      ->Fill(frac, cls);
}

void ResolutionHists::fillT0Residual(int cls, double dt_ns) {
  if (!enabled() || cls < 0)
    return;
  // +-400 ps at 2 ps/bin. t0 is TOF-dominated, so its accuracy is a ps-scale
  // quantity even though the window around it is tens of ns -- and it is
  // plotted in ps for the same reason the TOF hit times are. A +-2 ns axis
  // spread the entries so thin (<1 per bin) that the core fit could not
  // converge even though the pull, on the same data, could.
  get("t0_residual", "candidate t0 accuracy", "t_{0} - t_{MC} [ps]", 400, -400.0, 400.0)
      ->Fill(dt_ns * 1000.0, cls);
}

void ResolutionHists::fillT0Sigma(int cls, double sigma_ns) {
  if (!enabled() || cls < 0)
    return;
  get("t0_sigma", "reported t0 uncertainty", "#delta t_{0} [ps]", 1000, 0.0, 1000.0)
      ->Fill(sigma_ns * 1000.0, cls);
}

void ResolutionHists::fillT0Pull(int cls, double pull) {
  if (!enabled() || cls < 0)
    return;
  get("t0_pull", "candidate t0 pull", "(t_{0} - t_{MC}) / #delta t_{0}", 400, -10.0, 10.0)
      ->Fill(pull, cls);
}

} // namespace eicrecon::eb

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// ResolutionHists -- per-detector, per-class resolution histograms.
//
// Three residuals, each booked once per detector as a TH2 with the physics
// class on the y axis, so a 21-detector x 26-class breakdown is 21 objects
// rather than 546:
//
//   time    t_hit - t0            [ns]   every tracker detector and calo system
//   space   |rec - sim| position  [mm]   trackers (needs truth:simhits=1)
//           dR(cluster, MC)       [rad]  calo, angular instead of metric
//   energy  (E_rec - E_mc)/E_mc   [1]    calo only
//
// The class axis is the per-HIT label where one exists (the RawHitLink weight
// carries generatorStatus, so a hit is attributed to the collision that made
// it, not to whatever mix the candidate was credited with) and the matched MC
// particle's own status on the calo side.
//
// Everything here is inert unless a histogram file was requested: booking goes
// through RootFile_service, and fill() is a no-op when no directory exists.

#pragma once

#include <map>
#include <string>
#include <vector>

#include <JANA/JApplicationFwd.h>

class TDirectory;
class TH2D;

namespace eicrecon::eb {

class ResolutionHists {
public:
  /// Wired detector names, index-aligned with the collection lists the
  /// eventbuilder emits into each child.
  static const std::vector<std::string>& trackerNames();
  static const std::vector<std::string>& caloNames();

  /// `want_anyway` books histograms even with no -Phistsfile, holding them in
  /// memory so a PDF report can still plot them.
  void init(JApplication* app, bool want_anyway = false);

  /// Booked histograms, keyed "<kind>_<detector>".
  const std::map<std::string, TH2D*>& hists() const { return m_hists; }
  bool enabled() const { return m_dir != nullptr || m_in_memory; }

  /// Hit time against the TRUE COLLISION time. Carries the particle's time of
  /// flight, so for TOF it is not a sensor resolution -- see the booking in
  /// ResolutionHists.cc for why referencing the sim hit time does not fix it.
  void fillTrackerTime(std::size_t det, int cls, double dt_ns);
  /// In-plane residual: the component on the measuring surface,
  /// sqrt((r*dphi)^2 + dz^2). This is the reconstruction resolution.
  void fillTrackerSpace(std::size_t det, int cls, double dr_mm);
  /// Radial residual: r_rec - r_sim. A barrel segmentation pins rec hits to a
  /// FIXED cylinder, so a sub-volume displaced from its segmentation radius
  /// shows up here as a fixed offset -- geometry, not resolution. Kept
  /// separate so it cannot contaminate the number above.
  void fillTrackerSpaceRadial(std::size_t det, int cls, double dr_mm);
  /// Cartesian residuals, the quantities res_check plots: one column each.
  void fillTrackerX(std::size_t det, int cls, double dx_mm);
  void fillTrackerY(std::size_t det, int cls, double dy_mm);
  void fillCaloTime(std::size_t sys, int cls, double dt_ns);
  void fillCaloAngle(std::size_t sys, int cls, double dr);
  void fillCaloEnergy(std::size_t sys, int cls, double frac);

  /// Candidate-level timing: the trigger's own t0 against MC truth, the
  /// uncertainty it reports for that t0, and the pull of one against the
  /// other. A pull width of 1 means the reported uncertainty is calibrated.
  void fillT0Residual(int cls, double dt_ns);
  void fillT0Sigma(int cls, double sigma_ns);
  void fillT0Pull(int cls, double pull);

private:
  /// Book on first use: a detector absent from the run never gets an object.
  TH2D* get(const std::string& key, const std::string& title, const std::string& xlabel, int nbins,
            double lo, double hi);

  TDirectory*                  m_dir{};
  bool                         m_in_memory = false;
  std::map<std::string, TH2D*> m_hists;
};

} // namespace eicrecon::eb

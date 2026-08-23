# `eventbuilder` — software event building for the ePIC streaming readout

> **Start here:** `src/factories/eventbuilder/README.md` is the compact,
> authoritative structural guide (file map, status bands, pile-up semantics,
> weights tables, invariants). This README holds the deep derivations and
> measurements behind it.

The ePIC detector has **no hardware trigger**. The readout electronics
stream all data continuously. The stream is sliced into fixed **time
frames** of 2 microseconds. These frames are also called **timeslices**.
This document uses the two words, frame and timeslice, interchangeably.
JANA2 also uses the two words interchangeably. A frame is not a physics
event. A frame contains whatever happened during those 2 us. Usually a
frame contains a lot of machine background, such as synchrotron photons
and beam-gas interactions. Sometimes a frame contains one or more genuine
electron-ion collisions.

This plugin finds the collisions inside the frames, in software. The
plugin sends each collision to the standard reconstruction as its own
event. The plugin turns one input frame into zero, one, or several
**PhysicsEvents**. Each PhysicsEvent carries only the hits that belong to
it.

## Quick start

```sh
# standard eventbuilder run
eicrecon -Peventbuilder=true -Pplugins=eventbuilder \
         -Ppodio:output_collections="EventHeader,TimesliceHeader,EventBuilderInfo,MCParticles,..." \
         input_sim.edm4hep.root
```

The flag `-Peventbuilder=true` makes the file reader deliver each entry as
a 2 us frame (`JEventLevel::Timeslice`). Without this flag, the file
reader delivers each entry as a single event. Without this flag, the
plugin is loaded but stays idle. Without this flag, reconstruction behaves
like stock eicrecon.

Narrow `podio:output_collections`. If you do not narrow it, the writer
probes the full set of standard factories on unfolded events. This aborts
the process for factories whose inputs exist only at the Timeslice level.
See the Notes section.

The CKF tracking chain itself works correctly on unfolded children. The
unfolder projects the gated RecHits, the raw-hit associations, and the TOF
Measurement2D collections (`TOFBarrelClusterHits`/`TOFEndcapClusterHits`).
Request, for example, `CentralCKFTracks,CentralCKFTrackAssociations`. The
chain then resolves entirely from projected data.

## Vocabulary

| Term | Meaning |
|---|---|
| **frame / timeslice** | one 2 us slab of streamed data; the unit the file reader delivers |
| **fast detectors** | TOF (sigma ~ 30 ps) and MPGD (sigma ~ 10 ns) — good enough timing to *find* events |
| **slow detectors** | the silicon trackers and B0 (MAPS sensors) — too slow to find events, but essential hits once an event is found |
| **candidate** | a group of time-coincident fast hits that looks like a collision; what the trigger produces |
| **t0** | the candidate's collision time, estimated from its fast hits |
| **t0sigma** | the uncertainty on t0 (typically ~ 25 ps, because TOF dominates) |
| **gate** | the per-detector time window `t0 +/- 3 x sigma_hit` used to decide which hits belong to the candidate |
| **trigger** | the decision that a hit group is a candidate (this plugin's `EventBuilder_factory`, tag `trigger`) |
| **real / fake** | in simulation, a candidate either matches an injected MC collision (real) or is a background accident (fake) |
| **frame background** | the density of fast hits in the frame unrelated to any collision; sets how often fakes happen |
| **S, S_TOF** | *significance* — how far a candidate's hit count sticks out above the frame background, standardized by the Poisson mean and standard deviation (a Pearson residual). It reads as an actual Gaussian sigma only asymptotically, at large `mu` — see "Stored significance" for why it's named `S`, not `SNR`. `S` (unsubscripted — the primary trigger variable) is computed over the COMBINED TOF+MPGD hit population; `S_TOF` is computed over a narrower TOF-ONLY sub-cluster drawn from that same combined population — a subset, not a disjoint sibling group. There is no `S_Si`: the slow detectors (Si trackers, B0) never enter the hit-counting/coincidence step at all, they're gated by `t0` only after a candidate already exists |
| **false alarm rate (FAR)** | fake candidates per frame that survive a given cut; the realized, measured counterpart of `p_fake` (which only targets a rate). Reported at the applied trigger threshold (cut = 0) and, offline, as a function of a tighter `S`/`S_TOF` cut — see "Stored significance" |

## The detectors — what feeds the event builder

ePIC is a layered 4pi detector around the interaction point. Silicon
vertexing sits closest to the beam. Gaseous tracking surrounds the
silicon vertexing. Time-of-flight detectors follow, for particle ID.
Cherenkov counters and calorimeters follow the time-of-flight detectors.
Far-forward and far-backward instrumentation runs along the beamlines.
For event building, one number matters for each subsystem: how precisely
the subsystem timestamps a hit. This single number decides each
detector's role.

The event builder uses 10 tracker subsystems and 4 EM calorimeters today.
The list below orders them from fastest to slowest:

| idx | subsystem | technology | physics role | sigma_t | eventbuilder role |
|---|---|---|---|---|---|
| 1 | TOFEndcap | AC-LGAD silicon | forward hadron PID by time of flight | ~25 ps | **TOF trigger group** — sets t0; feeds both `S` and `S_TOF` |
| 0 | TOFBarrel | AC-LGAD silicon | central hadron PID by time of flight | ~30 ps | **TOF trigger group** — sets t0; feeds both `S` and `S_TOF` |
| 2 | MPGDBarrel | Micromegas (gas) | tracking point + timing between Si layers | ~10 ns | **combined trigger group** — feeds `S` only, not `S_TOF`\* |
| 3 | OuterMPGDBarrel | uRWELL (gas) | outer tracking point before TOF/calo | ~10 ns | **combined trigger group** — feeds `S` only, not `S_TOF`\* |
| 4 | BackwardMPGDEndcap | uRWELL (gas) | electron-side disk tracking | ~10 ns | **combined trigger group** — feeds `S` only, not `S_TOF`\* |
| 5 | ForwardMPGDEndcap | uRWELL (gas) | hadron-side disk tracking | ~10 ns | **combined trigger group** — feeds `S` only, not `S_TOF`\* |
| 6 | SiBarrelVertex | MAPS silicon | vertex position (impact parameter) | 10 ns sim / ~2 us real | slow collector, cross-frame buffered |
| 7 | SiBarrelTracker | MAPS silicon | momentum measurement (sagitta) | 10 ns sim / ~2 us real | slow collector, cross-frame buffered |
| 8 | SiEndcapTracker | MAPS silicon | forward/backward momentum | 10 ns sim / ~2 us real | slow collector, cross-frame buffered |
| 9 | B0Tracker | AC-LGAD in B0 dipole | very-forward proton/photon tagging | 8 ns sim | slow collector, cross-frame buffered |
| — | B0ECal, EcalBarrel, EcalEndcapN/P | crystals / imaging+ScFi / W-SciGlass | electron & photon energy | ~ns | clusters time-aligned and carried through |

\* MPGD has no test of its own. MPGD only ever contributes to the
combined, unsubscripted, `S`, alongside TOF. See "Stored significance"
below.

Remember this pattern: the fastest detectors define WHEN a collision
happened. The slowest detectors define WHAT happened. AC-LGAD TOF
measures time in tens of picoseconds. AC-LGAD TOF is designed for pi/K/p
separation by flight time. AC-LGAD TOF anchors t0 almost single-handedly,
through inverse-variance weighting. The MPGDs use gas amplification, with
a resolution of about 10 ns. The MPGDs are too coarse to fix t0. The
MPGDs are plentiful enough to make the coincidence and topology trigger
robust. The MAPS silicon detector actually measures vertex position and
momentum. The MAPS silicon detector measures most of the physics. The
MAPS silicon detector integrates over about 2 us, because of its rolling
shutter. Therefore, the MAPS silicon detector cannot vote on WHEN an
event happened. The MAPS silicon detector is a pure hit collector. The
event builder gates the MAPS silicon detector after t0 exists. The
TimesliceBuffer recovers MAPS hits across frame boundaries. This README
flags an honesty gap here. The current simulation smears MAPS times with
a 10 ns Gaussian. The current simulation does not model the 2 us shutter.
Therefore, the slow gates in simulation are much tighter than they can be
in real data.

The event builder does not yet use these detectors. The list below gives
the reasons:

- **Cherenkov PID** (DIRC, DRICH, pfRICH): these are single-photon
  detectors with good intrinsic timing. Their reconstruction is heavy.
  Their hits are ring images, not space points. Therefore, they have no
  trigger value. A pfRICH registration exists, in
  `src/detectors/PFRICH/PFRICH.cc`. This registration aborts on the
  craterlake geometry. See the Notes section. DIRC and DRICH would need
  new frame-level registrations of their own.
- **HcalEndcapN**: this handles energy for backward-going jets and
  neutrons. The barrel calorimeter (BHCAL) and forward calorimeter
  (FHCAL) already have cluster carry-through. See "Missing subsystems"
  below. HcalEndcapN is the one HCal region still without a frame-level
  plugin.
- **Far-forward/backward** (Roman Pots, Off-Momentum, LOW-Q2 tagger,
  LumiSpec): these detectors handle tagging and luminosity measurement,
  along the beamlines, several metres from the IP. They operate in a
  different timing regime. The plugin deliberately keeps them outside
  central event finding. ZDC sits in the same far-forward region. ZDC
  already has cluster carry-through. See "Missing subsystems" below.

So the honest count is this: **14 subsystems participate today**. This
total is 6 trigger detectors, plus 4 slow trackers, plus 4 EM
calorimeters. Full ePIC has about 25 subsystems. When you append a new
detector, existing indices stay stable. Any detector without genuine
trigger timing joins as a collector, after index 9. Such a detector
never joins the fast or mid groups.

## How an event gets built

The input is one frame with all its hits. Six steps build an event. Each
step is a factory in this plugin:

1. **Digitize at frame level** (`src/detectors/<DET>/<DET>.cc`, each
   detector's own file). Each detector file re-registers its standard
   ePIC digitization and hit reconstruction factories, in place, to also
   run at `JEventLevel::Timeslice`. JANA binds each factory to one event
   level. The stock registrations run per event, so the plugin cannot
   reuse them directly. This is the only reason these Timeslice-level
   registrations exist. Their outputs carry a `Digi` suffix, for example
   `SiBarrelVertexRecHitDigi`. This suffix means the outputs never
   collide with the standard per-event names. This step sets the
   per-detector time resolutions. This step is the ONLY place that sets
   them. Every later stage reads the resolution
   off each hit, through `getTimeError()`. This gives a single source of
   truth: TOF 30/25 ps, MPGD 10 ns, Si 10 ns, B0 8 ns.

2. **Align hit times** (`TrkTimeAlignment_factory`, tag `align`). This
   step corrects each hit's time by its straight-line flight time from the
   origin (r/c, with c = 299.79 mm/ns). After correction, all times refer
   to the collision moment, not the arrival moment. The output is
   `<Detector>TimeAlignRecHits`. The name states the template type,
   `edm4eic::TrackerHit`, not the tracking stage. This step runs long
   before any tracking. Its sibling factory, `CalTimeAlignment_factory`
   (tag `CalTimeAlignment`, template `edm4eic::Cluster`), applies the
   same r/c correction to calorimeter clusters. For the slow detectors,
   this step also deposits the aligned hits into a process-wide buffer.
   The buffer is keyed by frame number, through `TimesliceBuffer_service`.
   This deposit supports cross-frame recovery later, in step 6.

   **Origin = the NOMINAL IP (0,0,0), not the true per-event vertex.**
   This r/c flight-time correction and the topology angles theta/phi, in
   step 3, are both computed from the fixed global origin. The event
   builder runs BEFORE any tracking. Only downstream tracking and
   vertexing can find the true vertex. So at frame level, (0,0,0) is all
   there is. The trigger is deliberately built to tolerate this
   approximation. The luminous region spreads the real vertex by about a
   centimeter along the beam, and about a few microns transverse. This
   spread displaces a hit's apparent theta by only a few degrees. This
   displacement is less than one topology cell. The coarse grid and the
   half-cell stagger absorb it. The spread also shifts flight time by
   about tens of ps per cm. This shift is negligible against the 50 ns /
   1 ns coincidence windows. So event FINDING is robust to the vertex
   spread. The t0 computed here is an estimate for gating. The standard
   reconstruction recovers the precise vertex and event time. This
   reconstruction runs on each unfolded child event.

3. **Find candidates — the trigger** (`EventBuilder_factory`, tag
   `trigger`). The next section describes this step in detail. The
   output is `EventCandidates`. Each entry holds one candidate, with t0,
   t0sigma, its significance, and, in simulation, whether it matches a
   true collision.

4. **Optional GNN prefilter** (`EventPrefilter_factory`, tag
   `prefilter`). If you configure an ONNX model, this step rejects
   frames whose hit pattern looks like pure background. This rejection
   happens before any event is unfolded. This step is off by default.
   The output is `EventCandidatesFiltered`.

5. **Gate the slow hits** (`TrkTimeCoincidence`, tag `coincidence` — the
   same `TimeCoincidence_factory<edm4eic::TrackerHit>` type is wired a
   second time, tag `trkcoincidence`, to gate the fast TOF/MPGD hits;
   `CalTimeCoincidence` = `TimeCoincidence_factory<edm4eic::Cluster>` does
   the same for calo clusters, tag `calcoincidence`).
   For each candidate, this step keeps only the slow-detector hits
   inside `t0 +/- 3 x sigma_hit`. This is about +/-30 ns for Si and
   +/-24 ns for B0. Each hit uses its own resolution, read from the hit.
   The output is `<Detector>TimeCoincRecHits`, a subset of the aligned
   collections.

6. **Unfold** (`EventUnfolder`). This step creates one PhysicsEvent for
   each surviving candidate. Each PhysicsEvent contains: fast hits gated
   at `t0 +/- 3 x sigma_hit` (tight, about +/-90 ps for TOF); slow hits
   from step 5, plus the neighbouring frames' deposits (see Cross-frame
   recovery); calorimeter clusters; MC particles within `mc_time_window`
   of t0; and the three header collections described under "What each
   PhysicsEvent contains". The standard reconstruction then runs on each
   child event.

## The trigger, in detail

The trigger asks a simple physical question: did several fast hits
happen at the same time AND in the same direction? Background hits
spread over the whole frame and the whole detector. Collision products
bunch together in both time and angle.

**Same time — the sliding window.** The plugin sorts all fast hits of
the frame by aligned time. Two indices into that sorted list, `lo` and
`hi`, bound the current candidate window `[lo, hi)`. `hi` moves right as
far as it can go, while every hit up to it stays within `dt` of
`hits[lo]`. The value `dt` is a per-frame *derived* quantity, not a fixed
constant. See "Stored significance" below for the formula. `dt` is
typically about 40-50 ns. Each step of the scan does one of two things:

- **below threshold: slide.** The scan drops `hits[lo]` (`++lo`) and
  re-grows `hi` from the next hit. This is the "advances by ONE HIT"
  behavior. There is no fixed time step and no binning. The anchor just
  moves to wherever the next real hit is. So the scan tests every
  distinct hit grouping. This makes the scan as fine as the data itself.
- **at or above threshold: consume and jump.** The scan emits one
  candidate over *all* of `[lo, hi)`. The scan then jumps straight past
  it, `lo = hi`. The scan never re-examines those hits as the start of a
  second candidate. This is what keeps one physical cluster from
  producing a run of near-duplicate candidates, as the scan continues.

Here is a worked example, at threshold 3, with a 50 ns window, and hits
at 10/12/15/60/62/65/68/200 ns. Set cell binning aside for a moment. See
the eic-shell repo's `marco.md` §5.2 for the full step-by-step trace and
an ASCII table. The scan fires once on `{10,12,15,60}`: four hits, all
within 50 ns of the 10 ns anchor. The scan consumes them. The scan
restarts at 62 ns. The scan fires again on `{62,65,68}`. The lone hit at
200 ns never becomes a candidate. Note that `60` belongs to the *first*
cluster, within 50 ns of `10`, even though it is also within 50 ns of
`62`. Once consumed, a hit cannot also join the second candidate.

Why not a fixed 50 ns grid, or a 25 ns overlapping grid, instead of
anchoring on hits? A fixed grid can split a genuine cluster across a
boundary it had no way of anticipating. For example, four hits at
40/45/55/60 ns form a 20 ns-wide cluster, comfortably inside any
reasonable window. Under a naive `[0,50), [50,100), ...` scheme, these
hits land as 2+2 across the bin edge at 50 ns. This split misses a
threshold-4 coincidence that is obviously there. Halving the bin width
shrinks the odds of this problem, but cannot eliminate it. The
hit-anchored scan has no edges to straddle. The window's membership can
only change *at* a hit's own time. So testing an anchor at every hit is
not an approximation of a continuous slide — it *is* one, at no extra
cost. This is also what keeps the whole pass `O(N)`: each hit enters the
window via `hi` and leaves via `lo` exactly once over the full scan.

Here is why 50 ns specifically: the window must contain one event's
genuine spread of fast-hit times. The slowest trigger detector sets this
spread. Two MPGD hits from the same collision differ by
sqrt(2) x 10 ns = 14 ns, one sigma. So a 3-sigma catch needs about 42 ns.
The extra margin, up to 50 ns, allows for slow (beta < 1) particles
arriving late. Below about 40 ns, the window starts cutting genuine
coincidences in half.

**Same direction — the staggered angular grids.** The plugin bins every
fast hit by its direction from the origin. Polar angle theta uses 12
bins of EQUAL SOLID ANGLE. These bins are uniform in `cos(theta)`, so
each row spans an equal `dcos(theta)` = 2/12. The rows are NOT equal
15-degree width. See "What S/S_TOF are" for why this matters. Azimuth
phi uses 8 bins of 45 degrees each. Together these form a 12 x 8 =
96-cell equal-solid-angle map. Hits from one collision product, or a jet
of them, fall into the same cell. Background photons scatter across all
96 cells. One grid has a problem: a genuine cluster sitting exactly on a
bin edge can split between two cells. This split can cause the cluster
to miss the count threshold. So the plugin overlays a SECOND, identical
grid. This second grid is shifted by half a cell in both angles: half a
`cos(theta)` step in theta, 22.5 degrees in phi. Whatever straddles an
edge of grid 1 sits squarely inside a cell of grid 2. This is the
"12x8 + 12x8 staggered" test. The granularity is a runtime choice:
`topology:theta_bins` and `topology:phi_bins`, default 12 x 8. The
half-cell stagger follows automatically from these values. One cell
should be comparable in size to a plausible cluster's angular spread.
Finer cells cut accidentals harder. Finer cells also risk splitting
broad clusters. The default value has not been tuned systematically, but
it can be tuned, since it is now a configurable knob rather than a
hardcoded value.

**Why equal solid angle, and why it affects ACCEPTANCE (not just the
fake rate).** A cell spanning `[theta1,theta2] x [phi1,phi2]` subtends
solid angle `dOmega = dphi (cos theta1 - cos theta2)`. With equal-WIDTH
theta bins, this solid angle is tiny at the poles, because
`sin theta -> 0` there. This solid angle is maximal at the equator. For
12 equal-`dtheta` rows, the equatorial rows cover about 7.6 times more
sky than the polar rows. This single geometric fact biases the trigger
in two ways.

1. *Background normalization.* The flat `mu` (see "What S/S_TOF are")
   assigns every cell an equal 1/ncells share. But equatorial cells
   genuinely collect about 7.6 times more background. So ordinary
   equatorial fluctuations are scored against a too-small `mu`, and they
   look significant. Fakes pile up at the equator. The measured
   theta-row occupancy non-uniformity is about 32 times, under
   equal-width binning.
2. *Signal acceptance.* The trigger fires only if a cluster's hits
   CONCENTRATE into one cell, to reach `topology:threshold`. Whether a
   few-hit cluster lands in one cell, or splits across a boundary,
   depends on the cell shape at that theta. So equal-width binning makes
   the efficiency theta-dependent. It drops low-multiplicity events that
   fall in badly-shaped cells.

Binning uniformly in `cos(theta)` gives every cell the SAME solid angle,
`dOmega = dphi . dcos(theta)`. This binning makes the flat `mu`
unbiased, for an isotropic background. The occupancy non-uniformity
drops from about 32 times to about 7 times. The residual 7 times comes
from real detector geometry. This binning also makes the acceptance
theta-uniform. Here is the direct evidence: on the gold-coating scan,
DIS NC efficiency recovers from 96% to 100%, at the default
`topology:threshold=3`, purely from the binning change. A
low-multiplicity event, which the mis-sized cells were splitting, now
concentrates and fires. This binning change fixes the geometric
theta-bias only. It does NOT turn the low-`mu` Pearson-residual score
into a Gaussian sigma. The residual 7-times per-cell non-uniformity
would need a per-cell `mu_j` estimate. See "What S/S_TOF are".

**The decision.** Within the current window, the plugin counts hits per
cell in both grids. If any cell reaches `topology:threshold` hits
(default 2), the window triggers. So `topology:threshold` is a
GEOMETRICAL requirement: at least N time-coincident hits must point into
the same equal-solid-angle cell. At the default grid, this is one of 96
cells, each spanning `4pi/96` sr. When the window triggers, all hits of
the window form the candidate. The plugin computes t0. The scan then
jumps past the consumed cluster. The plugin maintains the two grids'
counters incrementally. This keeps the whole scan O(N).

**t0 and t0sigma.** The candidate's time, t0, is the inverse-variance
weighted mean of its fast hit times:
`t0 = sum(t_i/sigma_i^2) / sum(1/sigma_i^2)`. Each hit's sigma_i is its
own `getTimeError()`. TOF's 30 ps enters the formula as 1/sigma^2. So
TOF dominates completely whenever it is present. The measured t0sigma is
about 22-25 ps. Silicon never enters this calculation. Silicon cannot
time-resolve events. Silicon only contributes hits after t0 already
exists.

**Real or fake (simulation only).** A candidate is labelled real, flag
1, if its t0 lies within the coincidence window of an injected MC
signal-collision time. A real candidate records that collision's true
time. Otherwise, the plugin labels the candidate fake, flag 2. The
plugin identifies signal collisions generator-agnostically. Signal-stream
particles keep their raw generator status. The merger offsets background
streams by +1000×stream. So the plugin reads collision times off
status-1, final-state, particles. The plugin deduplicates these times
per collision. This method works for every generator: pythia, rapgap,
lAger, and DVCS. The previous method used the pythia-only beam-remnant
status 61, which does not work for every generator. The plugin stores
this label per event. This label is the ground truth for trigger-purity
studies and GNN training. Storing it costs nothing extra.

## Statistical foundation — the trigger as a hypothesis test

Everything above is, formally, a per-cell hypothesis test. This section
states that test compactly. The detailed prose on what the stored
numbers do and do not mean is in "Stored significance" and "What
S/S_TOF are" below.

**Null and alternative.** Consider an angular cell with occupancy `C`:
the number of hits in that cell within the window. The test is:

    H0 : C ~ Poisson(mu)          background only
    H1 : C ~ Poisson(mu + s), s>0  a collision adds a local excess

`mu` is the expected accidental count per (cell × window). The plugin
does NOT assume `mu` or take it from simulation. The plugin estimates
`mu` from the frame itself: `mu = (N · dt/dT) / ncells`. Here, N is the
frame's fast-hit count, dt is the coincidence window, and dT is the
frame's hit-time span. This dT value is the same self-describing value
stored per candidate in `weights[18]`. A frame is about 99.9% background. So this is a
consistent, self-calibrating estimate of the background rate. The
signal in the frame inflates `mu` only at the `O(N_sig/N_bkg) ~ 0.1%`
level.

**What level S lives at, and how a quoted number is built.** S is a
PER-CANDIDATE scalar, one value per built PhysicsEvent, stored in
`weights[4]`. The plugin computes S once at trigger time, from hit
COUNTS: `S = (cmax − mu)/sqrt(mu)`. Here, `cmax` is the hit count in the
busiest `(theta,phi)` cell, stored in `weights[7]`. S is NOT a per-hit
quantity. Worked example: with `cmax=7` and `mu=0.18`,
`S = (7−0.18)/sqrt(0.18) ≈ 16`. A range quoted per class, for example
"7.3–18.3σ true", is a two-level summary. First, the plugin takes the
median S over the class's true candidates. Then it takes the min-max of
those medians across classes. The chain is: hits → (cmax, mu) → one S
per event → median per class → range across classes.

**The test is optimal for what it uses.** The Poisson family has a
monotone likelihood ratio in `C`. So, by the Karlin-Rubin theorem, the
count threshold `C >= k` is the UNIFORMLY MOST POWERFUL test of H0 vs
H1, among all tests that use only the count. No cleverer scalar on the
same count does better. Improving on it requires using the hit PATTERN.
This is exactly what the GNN prefilter does: it is a learned surrogate
for the full likelihood ratio. The counting trigger is the optimal cheap
stage. The GNN is the way past its ceiling.

**`p_fake` is a per-trial size, not a rate.** The trigger reports the
MAXIMUM count over all cells and window positions. This is about
`M ~ N × 2 grids × ncells ~ 1e5` trials per frame. So the per-cell tail
`p_fake` is NOT the per-frame fake probability. Under H0, the per-frame
fake expectation is the look-elsewhere product:

    E[fakes / frame]  ≈  M × p_fake .

This is why a tiny per-trial value, `p_fake = 1e-5 … 1e-8`, still leaves
a manageable number of accidental candidates per frame. The realized
number is the measured FAR. `p_fake` is only its per-trial target. The
single most common misreading is to interpret `p_fake` as a per-frame or
per-candidate probability.

**The Poisson assumption must be validated, because it may fail.** H0
assumes independent hits. But one background interaction makes SEVERAL
correlated hits. So the count is compound Poisson, and it is generally
OVER-DISPERSED. The Fano factor is
`F = Var[C]/E[C] = E[M] + Var[M]/E[M] >= 1`, for a random cluster size
`M`. The trigger lives in the far tail of the distribution. So
over-dispersion would make the threshold optimistic: it would produce
more fakes than `p_fake` implies. The check compares the measured
upper-tail probability against the Poisson prediction. This check is a
Pearson χ² goodness-of-fit test, on the occupancy distribution:
`χ² = Σ (O_n−E_n)²/E_n`, over count bins. This test is a ONE-TIME model
audit. It is NOT a per-candidate quantity. Note that `χ² = Σ
residuals²`. So this test shares the name "Pearson" with the per-cell
residual `(C−mu)/sqrt(mu)`. But it is a different object: one value per
sample, not one value per candidate. If the check fails, the negative
binomial distribution is the natural over-dispersed replacement. See the
per-cell background estimate in the outlook.

## Adaptive-threshold mode — holding the fake rate constant

With the default fixed threshold, the same rule of "3 hits in a cell"
means different things in different frames. In a quiet vacuum frame, 3
hits is about a 5 sigma fluctuation. In a dense gold-coating frame, 3
hits is only about a 1.7 sigma fluctuation. So a fixed count is too
tight for quiet beams and too loose for dense beams, at the same time.
Setting `adaptive_threshold=1` replaces the fixed count with a
statistics test. The plugin computes this test per frame:

```
mu = (N x dt / dT) / ncells   # N = frame's fast-hit count, dt = derived combined-group window (see below), dT = frame hit-time span; ncells = theta_bins x phi_bins
k  = smallest k such that P(X >= k | Poisson(mu)) < p_fake      # floored at topology:threshold
```

You choose `p_fake`. This is the probability that pure background fires
one (cell x window) trial. At the historical 12x8 = 96-cell grid, a
frame performs roughly n_hits x 2 grids x 96 cells, about 1e5 trials.
The current production default is a 1x1 = 1 cell grid instead
(`theta_bins=phi_bins=1`). See `EB_TOPOLOGY_THETA_BINS` and
`EB_TOPOLOGY_PHI_BINS` in the production Makefile (external repo). So `ncells` and the trial
count scale down accordingly, in production.
Setting `p_fake = 1e-5` targets about one fake per frame, at the
historical grid. This target stays true whether the frame is
gold-coating, vacuum, or vanilla, because the plugin measures mu from
the frame itself. Across all six signal classes, at p_fake 1e-5, the
rate is 12 fakes/frame, with only 66-100% efficiency by class (DVCS is
worst, at 66%). At p_fake 1e-8, the rate is 2.5 fakes/frame, with 29-94%
efficiency (DVCS at 29%). Adaptive
mode discards low-multiplicity signal FIRST as it tightens. The two
events lost at a fixed threshold of 4 are the same low-multiplicity
population. But adaptive mode's per-frame adaptivity makes this effect
worse in dense frames. Adaptive mode is the right knob for a
bandwidth-limited ONLINE trigger, holding a fixed fake-rate budget.
Adaptive mode is the wrong knob for training-dataset production. In
production, a lost real event is unrecoverable, while extra labelled
fakes are useful. So adaptive mode is off (`adaptive_threshold=0`) in
the production baseline. See the Makefile `EB_TOPOLOGY_THRESHOLD` and
`EB_TRIGGER_PFAKE` docs, and the "Configuration reference" section
below, for the full parameter sweep.

**Two resolution groups.** MPGD forces the combined window `dt`. This
wastes TOF's timing, which is 300 times better. A TOF pair inside
`dt_tof` is far more significant than the same pair spread over `dt`.
Therefore, in adaptive mode, the trigger fires on EITHER of the
following:

- **combined group** (unsubscripted): all fast hits within `dt`. The
  plugin Poisson-tests this group against mu, as described above. This
  group is floored at `topology:threshold`.
- **TOF group**: TOF hits only, within `dt_tof`. The plugin
  Poisson-tests this group against its own, much smaller, mu_TOF. This
  group is floored at `topology:threshold_tof` (2).

Neither `dt` nor `dt_tof` is a fixed constant or a separate parameter.
The plugin derives both, once per frame, from the single
`nsigma_window` confidence level and each pool's own measured hit
resolution: `dt = nsigma_window * sqrt(2) * max(sigma over that pool's
hits this frame)`. The plugin stores this value per candidate, in
the per-frame EventBuilderFrameInfo record (indices 3/4), described below, rather than assuming a constant value
across a file. The `sqrt(2)` factor is the quadrature width of the
worst-case *pair*: two hits both at the pool's own worst resolution.
This is provably the widest pairwise gap the pool can produce, since
mixing in any better-resolution hit only narrows it. At the default
`nsigma_window=3`, this typically lands near `dt≈42 ns`, MPGD-dominated,
and `dt_tof≈0.13 ns`, TOF-only. See `EventBuilder_factory.h`'s
derivation comment for the full reasoning.

Both groups use the same `p_fake`. Silicon has no trigger group, by
construction. Per-detector significances are not used: a
per-detector-per-cell mu is too small to estimate stably per frame, and
the trigger counts combined cell occupancy anyway.

`topology:threshold` and S are complementary, not alternatives. The
threshold is the PHYSICAL floor. It sets a minimum cluster size for a
coincidence to mean anything. It protects the quiet-frame limit, where
mu approaches 0 and any bare pair would otherwise look "significant".
The parameter p_fake supplies the STATISTICAL requirement, which scales
with background. `k`, the applied threshold, is always the maximum of
the two values.

**Stored significance.** In BOTH modes, every candidate records how far
it stuck out, at two nested scopes. These are NOT two disjoint detector
groups:

- `S = (cmax - mu)/sqrt(mu)`: the significance of `cmax`, the
  coincidence count over the FULL flattened hit population, TOF plus
  MPGD together. See `EventBuilder_factory::Process`, "Flatten the fast
  (TOF+MPGD) hits". This is the primary trigger variable.
  `topology:threshold` mode fires on `cmax >= k` alone.
- `S_TOF = (cmax_TOF - mu_TOF)/sqrt(mu_TOF)`: the significance of
  `cmax_TOF`, a narrower count restricted to the TOF-only SUBSET of that
  same population, in its own, usually shorter, derived window `dt_tof`.
  `S_TOF` can independently fire the trigger only in adaptive mode,
  through `cmax_TOF >= k_TOF`. Otherwise, `S_TOF` is an extra
  offline-cuttable feature.

### How cmax_TOF is computed, and what it costs

`cmax` (combined channel) is cheap by construction: the two-pointer scan's
own `lo`/`hi` bookkeeping already maintains `occupancy1`/`occupancy2` per
cell incrementally (one increment on entry via `hi`, one decrement on exit
via `lo`), so reading `cmax` each step is just an O(ncells) argmax over
already-current counts — bounded purely by grid size, never by how dense
the background is.

`cmax_TOF` is harder, because `dt_tof < dt`: a genuine TOF-only cluster can
be tighter than the combined window reveals, so `cmax_TOF` is a max over a
family of `dt_tof`-wide *sub*-windows nested inside `[lo, hi)`, not a plain
count. The naive approach — re-filter the TOF-only hits out of `[lo, hi)`
and rescan them from scratch on every outer-scan step — costs O(N x n_tof)
per frame, where `n_tof` is the number of TOF hits active in a typical
window. Unlike the combined channel's O(ncells) argmax, this scales with
*background density*, not grid size.

Production uses a per-cell tracker (`CellTracker` in
`EventBuilder_factory.h`) that maintains each cell's own best TOF-window
count as hits enter or leave, turning the query into an O(ncells) read
instead of an O(n_tof) rescan. Its machinery generalizes to a future
second timing-sensitive detector group, should one ever need its own
nested coincidence test the way TOF does today.

**Caveat for pile-up.** Whenever background hits concentrate into few
cells, which is exactly what happens at the coarsest grid
(`theta_bins=phi_bins=1`, the grid adaptive-threshold mode uses by
default, since it needs one global occupancy estimate rather than
per-cell ones), `n_tof` grows, and grows further as pile-up raises the
background hit rate. This is not currently on the production hot path,
since the trigger normally runs on the finer grid, but is worth
remembering if adaptive/1x1 configurations become routine under heavier
pile-up. Note that the combined channel has no nested-window problem to
solve: every hit in `[lo, hi)` already satisfies `dt` by construction, so
its O(ncells) argmax cost is already density-independent and would not
benefit from `CellTracker`'s machinery.

### What S/S_TOF are, statistically — read before using them

There are three separate statements, often conflated. Keeping them apart
is the whole point:

**1. The trigger's decision is exact and correct.** The online
accept/reject decision is an integer count test, `cmax >= k`. In
adaptive mode, `k` comes from the EXACT Poisson survival function,
`poissonMinK`: the smallest `k` with `P(X >= k | mu) < p_fake`. This
calculation uses no approximation. So nothing below affects whether a
candidate fires. It affects only the interpretation of the stored score.

**2. The stored score is a Pearson residual, not a Gaussian sigma.** The
formula `S = (c - mu)/sqrt(mu)` standardizes the observed count by a
Poisson's own mean, `mu`, and standard deviation, `sqrt(mu)`. This is a
*Pearson residual*. It is exact as a standardization of the first two
moments. But its value equals a number of Gaussian sigmas ONLY when
Poisson(`mu`) is itself approximately Gaussian. This holds when `mu` is
of order a few or larger, say `mu >= 5`. The Poisson distribution is not
"inaccurate" here. The Gaussian APPROXIMATION to it is inaccurate. At
`mu << 1`, the Poisson distribution is sharply right-skewed and
discrete. It looks nothing like a bell curve. So converting the residual
to a p-value, through the normal distribution, overstates the rarity.

  Here is a concrete example. With `mu = 0.01` and `c = 2`, the formula
  gives `S = 19.9`. But the exact tail is `P(X >= 2 | 0.01) ~ 5e-5`, a
  rarity of about 4 sigma, not 20. At the low `mu < 1` typical of every
  candidate, the stored score over-reads the exact one-sided
  Poisson-tail sigma, by a median factor of about 2 times, and up to
  about 9 times. For example, a "6.8 sigma" fake is really about 3.2
  sigma. The code also returns a 99.0 sentinel value when `mu == 0`.
  This sentinel is a symptom of the same low-mu pathology.

  The consequence is this: treat `S` and `S_TOF` as monotone RANKING
  variables, for offline cuts and as GNN input features. For those uses,
  only the ordering matters, and the ordering is fine. Where you need an
  actual probability at low `mu`, compute the Poisson survival
  `P(X >= c | mu)` directly. Both `cmax` (weights[7]) and `mu`
  (weights[18]) are stored for exactly this purpose. The external
  diagnostics tooling renders the Pearson-vs-exact-tail comparison,
  and the P(X=n)/P(X>=k) curves, per dataset.

**3. The `mu` estimate assumes an equal per-cell background share.**
`mu` is the frame-wide fast-hit rate, divided equally over all
`theta_bins x phi_bins` cells. The cells are POSITION-angle bins:
`theta = acos(z/|r|)` and `phi = atan2(y,x)`, of the HIT POSITION
relative to the origin (0,0,0), the nominal IP. This is where hits
landed, not their emission directions. For the equal-share `mu` to be
right, the cells must receive equal background on average. There are
two reasons they might not, with the status of each:

  - *Solid angle (fixed in code).* Theta is now binned in `cos(theta)`
    over [-1,1]. So each cell subtends an equal solid angle,
    `dOmega = dcos(theta) dphi`. The previous behaviour, equal-WIDTH
    theta bins, did not have this property. Polar cells subtend far
    less solid angle than equatorial ones. So an isotropic background
    piled up about 7.6 times more in the central theta rows. The flat
    `mu` under-estimated the local rate there, inflating those cells'
    S. The cos(theta) binning removes this purely geometric bias.
  - *Detector geometry and beam structure (not fully modelled).*
    Different theta rings still cross different detector layers, such
    as barrel versus endcap, with different radii and channel
    densities. Synchrotron background is also not azimuthally uniform.
    Under the pre-fix binning, per-theta-row background occupancy
    spanned about 32 times, well beyond the 7.6-times solid-angle part.
    A frame-wide scalar `mu` cannot capture this. A per-cell `mu_j`, a
    rolling or historical background estimate, is the proper fix. This
    fix is future work. The winning-cell record, weights[20..22],
    exists to measure the residual non-uniformity empirically.

Despite the parallel `S`/`S_TOF` naming, `S_TOF`'s hit pool is a subset
of `S`'s: TOF hits ⊂ TOF+MPGD hits. Cutting on both together is not
double-counting independent detectors. It is a coarse cut, S, plus a
finer cut on the timing-cleanest slice of the same cluster, S_TOF. There
is no third `S_Si`. Si/B0 hits are gated by the already-fixed `t0` (see
"gate" above). Si/B0 hits never enter `hits` or the coincidence count.
So there is no cluster for them to "stick out" of.

One measurement of this separation, on gold-coating dis_nc data with
p_fake=1e-8, used truth labels and theta binning both since superseded
(see "Real or fake" and "Why equal solid angle" above), so treat it as
indicative only. The median S was 32 for real candidates, versus 19 for
fake candidates. The median S_TOF was 117 for
real candidates, versus 28 for fake candidates. Separation varies by
class and background variant. The current production baseline's
separation, per class and variant, is in the external dataset
diagnostics report. You can choose the operating point OFFLINE. You can produce
candidates loosely, then cut on S or S_TOF later. This cut is
reversible. Both values are also free input features for the GNN.

**False alarm rate, made concrete.** Be precise about what `p_fake` is.
It is a target tail probability per (cell x window-position) TRIAL. The
name `p_local_tail` would describe it equally well. It is NOT
fakes-per-frame. It is NOT a probability per candidate. It is NOT a rate
in Hz. A frame performs roughly n_hits x 2 grids x ncells, about 1e5
such trials. So the frame-level fake count is about the effective trial
count, times p_fake. With 2 us frames, the wall-clock rate is
fakes/frame times 500 kHz. FAR, fakes per frame at a given cut, is the
measured, physical counterpart of the targeted `p_fake`. The calibration
from p_fake to FAR is an empirical curve, per background condition. It
is not an identity. `p_fake` is a target, set before production. FAR is
what actually came out of one file, at whatever cut is applied. Every
surviving candidate carries its own `S` and `S_TOF`. So a tighter
offline cut can be swept from 0, the applied threshold, upward, without
rerunning eicrecon. FAR(cut) is the number of fake candidates with
`S >= cut`, per frame. Efficiency(cut) is the fraction of real
candidates with `S >= cut`. This whole trade-off curve is computable
offline, from the stored per-candidate `S`/`S_TOF`/`cmax`/`mu` weights
alone. The external diagnostics dashboard renders this trade-off curve. It
also renders a small table of FAR and efficiency, at a few
representative cuts: the applied threshold, and the 50th, 80th, and 95th
percentile of the fake-S distribution. See Notes. Cutting harder is a
one-way, reversible knob, applied downstream of this plugin. It never
needs a new production pass.

## Cross-frame recovery

A hit near a frame edge can end up recorded in the neighbouring frame.
The slow detectors are affected most. The detector response also has a
late tail. So recovery matters in BOTH directions.
`TrkTimeAlignment_factory` deposits each frame's aligned slow
hits into `TimesliceBuffer_service`. This buffer is process-wide, keyed
by absolute frame number, and thread-safe. The unfolder fetches frames
`[N - backward_frames, N + forward_frames]`. The unfolder gates their
hits per candidate, exactly like current-frame hits.

- `backward_frames` (default 1): Frame N-1 is always already deposited
  when frame N unfolds. This recovery is free and needs no
  requirements. It catches the symmetric part of the spillover, from
  Gaussian smearing.
- `forward_frames` (default 1): The measured hit-time residuals around
  t0 have a pronounced LATE tail, similar to a Landau distribution.
  Slow particles and delayed detector response arrive after the
  collision. So spillover lands preferentially in the NEXT frame. This
  makes forward recovery at least as important as backward recovery.
  Waiting for frame N+1 requires frame N+1 to be in flight concurrently.
  This requires `nthreads >= 2` and
  `jana:max_inflight_timeslices > forward_frames`. When these conditions
  are not available, the unfolder prints a loud warning and runs with 0
  instead of aborting. Neither case adds output latency. The production
  pipeline runs eicrecon with 4 threads. So forward recovery is fully
  active there.

The buffer retains 32 frames. The buffer prunes older frames. The buffer
treats frames that never arrive, at stream edges, as absent after a
bounded wait. Correctness does not depend on thread count or processing
order.

Known seam: the TRIGGER itself still sees only one frame. So a collision
within about 50 ns of a frame edge can have its fast-hit cluster split
across two frames. Such a collision can be missed. This affects about
2.5% of events. A fix would add a forward margin to the trigger scan,
with window-start ownership to avoid double counting. This fix is not
yet implemented.

## What each PhysicsEvent contains

There are three header collections. All three use the
`edm4hep::EventHeader` type:

| Collection | Content |
|---|---|
| `EventHeader` | the child's OWN header, standard semantics only: `eventNumber = frame x 1000 + candidate`, `timeStamp` = t0 in ps, `weight = 1`. Safe for any standard tool. |
| `TimesliceHeader` | passthrough of the parent frame's header (provenance). |
| `EventBuilderInfo` | the eventbuilder payload, in its `weights` vector (below). Deliberately NOT named EventHeader so nothing mistakes these numbers for MC generator weights. |

The `EventBuilderInfo.weights` vector has this layout. It is kept in
sync with the encoding comment above `EventBuilder_factory::emitCandidate`.
That comment is the ground truth:

| index | value |
|---|---|
| 0 | t0 (ns) |
| 1 | t0sigma (ns) |
| 2 | flag (simulation only): the number of real MC collisions inside this candidate's dt window. 0 = fake, 1 = one real collision, 2+ = true pile-up. Real means flag > 0 |
| 3 | mc_t: the EARLIEST matched collision's time, in ns (simulation only; -1e9 = none). When flag >= 2 the other collisions' times are in trigger_classes (index 25+) |
| 4 | S: combined-group standardized excess (see the caveat under "Stored significance") |
| 5 | S_TOF: TOF-group standardized excess |
| 6 | p_tail = P(X >= cmax \| Poisson(mu)): the exact upper-tail background probability behind S |
| 7 | cmax: combined-group coincidence count that fired (the numerator of S) |
| 8 | cmax_TOF: TOF-group coincidence count |
| 9 | E_calo (GeV): calorimeter-coincidence energy, the summed time-aligned cluster energy with \|t_cluster − t0\| <= dt. Score only; never in the trigger decision |
| 10 | n_tracklets: vertex-pointing fast-hit pairs in the fired window. Score only |
| 11 | n_expected: in-acceptance stable charged MC particles of the earliest matched collision. Tracking-efficiency denominator |
| 12 | trk_threshold_cfg (eventbuilder:trigger:min_tracklets) |
| 13 | trk_pass: 1 if n_tracklets >= index 12 |
| 14 | cal_threshold_cfg (eventbuilder:trigger:min_cal_energy, GeV) |
| 15 | cal_pass: 1 if E_calo >= index 14 |
| 16 | trigger bitmask: bit0 = TOF group fired (cmax_TOF >= k_TOF), bit1 = combined group fired (cmax >= k), bit2: backward-ECal region tag (in-window cluster, nhits >= 10), bit3: bit2 plus >= 1 backward endcap tracklet, bit4/bit5: same for barrel ECal, bit6/bit7: same for forward ECal, bit8: B0 tag (>= 4 in-window B0 hits), bit9: ZDC tag (in-window ZDC-ECal cluster hits sum >= 50). Inclusion tags only |
| 17 | trigger_classes_mask: bit N = class_index N (0-25) had a collision within dt of this candidate. 0 = none |
| 18 | mu: winning cell's leave-one-out background (per candidate, not frame-wide) |
| 19 | k: winning cell's applied count threshold (per candidate) |
| 20 | winning cell's theta bin |
| 21 | winning cell's phi bin |
| 22 | winning grid: 1 = plain, 2 = half-cell staggered |
| 23 | S_empirical: significance of cmax against this frame's own sideband distribution (the empirical_band method), always stored regardless of the active null model |
| 24 | fano: variance/mean of the sideband window counts. 1 = Poisson-like |
| 25+ | trigger_classes (on by default via store_coincident_list): index 25 = N, the coincident-collision count (equals flag), then N (time, stream) pairs at 26+2k / 27+2k. stream = generatorStatus/1000 of the collision's base status: 0 = native primary, >= 10 = a physics class (class_index = stream-10) |

Frame-constant configuration lives in a separate **`EventBuilderFrameInfo`**
collection (same `edm4hep::EventHeader` type), one entry per frame, copied
into each child with `eventNumber = frame*1000` — join it to candidates
(`frame*1000 + candidate`) via `eventNumber / 1000`:

| index | value |
|---|---|
| 0 | n_physics_events: physics events injected into this frame (pileup-efficiency denominator; the flag > 0 candidate count is the numerator) |
| 1 | mu_TOF: measured TOF-group background (hits per cell per window) |
| 2 | k_TOF: TOF-group count threshold applied (999 marks it disabled) |
| 3 | dt, in ns: combined-group coincidence window this frame (derived per frame from `nsigma_window`) |
| 4 | dt_tof, in ns: TOF-group coincidence window this frame |
| 5 | topology:threshold |
| 6 | topology:theta_bins |
| 7 | topology:phi_bins |
| 8 | adaptive-mode flag |
| 9 | p_fake |

Separately, each PhysicsEvent also carries a **`PrefilterScores`**
collection, an `edm4hep::EventHeader`, when the GNN prefilter runs. Its
`weights[]` hold the softmax class probabilities (class 0 = BKG; the
class count depends on the loaded model, for example 27 for
`prefilter-20270803.onnx`). This is the GNN decision, stored per event.
So the background veto and the physics-class tag become reversible
OFFLINE cuts, needing no ONNX reload. See the ONNX prefilter section.
This collection is empty when no model is configured.

Two further collections, **`PrefilterLambdaK`** and **`PrefilterSegLogits`**
(same `edm4hep::EventHeader` type), carry the model's other raw outputs
when `eventbuilder:prefilter:store_aux_outputs=1` is set and the loaded
model actually declares them (looked up by name, not position).
`PrefilterLambdaK.weights[]` is the raw `lambda_k` output; `PrefilterSegLogits.weights[]`
is the raw per-hit `seg_logits` output, flattened in the same hit order
as the GNN's own hit input (see the ONNX prefilter section), truncated
to the real hit count. Neither is trained to be physically meaningful
yet — they exist for offline inspection of the model as it develops.
Both are empty unless `store_aux_outputs` is on.

Hit collections work as follows. Fast-detector hits are podio SUBSETS:
references into the frame's collections, not copies. Slow-detector hits
are owned clones, because cross-frame hits cannot be references. Raw
hits, truth links, calorimeter clusters, and MC particles are carried
through per candidate. MC particles are gated by `mc_time_window`,
default +/-50 ns around t0.

## Missing subsystems — inventory and how to add one

Stock eicrecon registers 22 detector groups, in `src/detectors/`. The
event builder has frame-level plugins for 18 of them. Of these, 16 are
active. PFRICH and LUMISPECCAL are written but stay disabled; see Notes.
Of the 18 plugins, three were added after the original wiring. These
three are cluster carry-through only. They do not participate in the
trigger:

| subsystem | what it is | eventbuilder role | notes |
|---|---|---|---|
| `BHCAL` | barrel hadronic calorimeter (SiPM tiles) | calo cluster carry-through | `src/detectors/BHCAL/BHCAL.cc`, validated (loads clean, ~1 cluster/event carried through as `HcalBarrelClusters`) |
| `FHCAL` | forward HCal (LFHCAL + EndcapP insert) | calo cluster carry-through | `src/detectors/FHCAL/FHCAL.cc` (LFHCAL Island + HcalEndcapPInsert imaging); carries `LFHCALClusters` (~3/ev) + `HcalEndcapPInsertClusters` |
| `ZDC` | zero-degree calorimeter (far-forward neutrals) | calo cluster carry-through | `src/detectors/ZDC/ZDC.cc` (EM ZDC Island + Hcal ZDC imaging); `EcalFarForwardZDCClusters` + `HcalFarForwardZDCClusters` (populate on forward-neutron events) |

The remaining 4 subsystems have **no frame-level plugin at all**. The
table below shows what adding each one would mean:

| subsystem | what it is | natural eventbuilder role | notes |
|---|---|---|---|
| `DIRC` | barrel Cherenkov (quartz bars, MCP-PMTs) | none for triggering | good intrinsic photon timing but hits are ring images, not space points; carry-through only if PID-at-frame-level is ever wanted |
| `DRICH` | forward dual-radiator RICH | none for triggering | same reasoning as DIRC; heavy reconstruction |
| `LOWQ2` | far-backward tagger tracker (Timepix) | slow collector at most | its stock chain includes ONNX ML factories that throw on beam-electron-less events — needs care |
| `RPOTS` | Roman Pots (far-forward AC-LGAD) | interesting: genuinely FAST silicon (~30 ps) | the one missing detector with real trigger-grade timing — but it sits metres downstream in the hadron beamline; adding it to the trigger would need its own time-of-flight offset handling |

Note: under the now-removed `plugins/*.cc` layout, the eventbuilder
plugin globbed every detector's carry-through file into ONE target.
Stock eicrecon instead uses one target per detector. Because of this,
any header with non-inline free-function definitions clashed when a
second detector's file included it. Adding FHCAL and ZDC surfaced this
problem in `ImagingTopoClusterConfig.h`, in the `ELayerMode` stream
operators.
The fix was to mark them `inline`. This was a latent upstream ODR bug.
Watch for the same problem when adding further imaging-cluster
detectors.

**How to add one.** Edit the detector's own stock plugin,
`src/detectors/<DET>/<DET>.cc`. Append a Timeslice-level mirror of its
existing `InitPlugin` registrations, at the end of the same function.
See `src/detectors/BTOF/BTOF.cc` for a worked example. Do not create a
separate eventbuilder-local file. `InitJANAPlugin(app)` and `using
namespace eicrecon;` already run once, at the top of the function, and
apply to the whole function. So the mirror block only adds
`app->Add(...)` calls. Follow this checklist:

1. Copy each registration from the chain above.
2. Wrap it in `->SetLevel(JEventLevel::Timeslice)`.
3. Add a `Frame` suffix to every collection and tag name it produces or
   consumes. A JANA factory binds to one event level, so the Timeslice
   variant needs distinct names from the PhysicsEvent-level ones. See
   "How cmax_TOF is computed" for why the suffix is `Frame`, not `Digi`.
4. Set timeResolution explicitly, to match the PhysicsEvent-level
   config.
5. Append the new terminal collection to the right wiring list in
   `eventbuilder.cc` (`m_simtrackerhit_collection_names` or
   `m_simcalocluster_collection_names`), without shifting existing
   indices.
6. Update the index tables in `TimeAlignment_factory.h`.
7. Whitelist the collection in your `podio:output_collections`
   configuration.
8. Validate on every beam compact, before enabling by default.

**Why not a separate file, the way it used to work.** A second copy of
each detector's registration drifts silently from the stock one over
time. Also, having eventbuilder call a second `InitPlugin`-equivalent
function, on top of the one the standard `eicrecon` executable already
loads (every default plugin, unconditionally — see `eicrecon_cli.cc`),
double-initializes it. This double-initialization caused a reliable
segfault at plugin load. One file per detector, with one registration
call per detector, is both simpler and correct.

## ONNX prefilter (optional)

A GNN classifies each FRAME as signal-like or background-only, before
any unfolding, through `ONNXRuntime_service`. A frame passes if its
background probability is below `score_threshold`. The plugin computes
6 hit features per hit, per frame:
`[x/4000, y/4000, z/5000, t/50, clip(eDep/10,0,1), det_id/6]`. In this
formula, det_id is 2 for TOF and 1 for MPGD. This prefilter is off by
default, with an empty model path.

**Veto vs score-only.** With `prefilter:veto=1`, the default, a frame
below threshold is dropped online. This drop is irreversible. With
`prefilter:veto=0`, the GNN still runs. Its class probabilities are
stored per event, in the `PrefilterScores` collection (see "What each
PhysicsEvent contains"). But NOTHING is dropped. The background veto and
the physics-class tag become REVERSIBLE OFFLINE cuts, on the stored
scores. They need no ONNX reload or kNN rebuild downstream. Score-only
mode is the recommended setting for training-dataset production. This
follows the same rationale as keeping the trigger loose and cutting on
significance offline: never throw away data online, just annotate it.
The plugin stores the scores whenever a model runs, in either mode.

**Model contract.** The GNN takes the hit features above plus a k-NN edge
list per DGCNN layer (`edge_idx_0`, `edge_idx_1`, ...; the exact layer
count is read from the loaded model itself, not hardcoded). Its first
declared output is always read as the classifier (`probs`). A model may
declare further outputs — for example `lambda_k` (frame-level) and
`seg_logits` (per-hit) in `prefilter-20270803.onnx` — but those are
ignored by default, since neither is trained to be physically meaningful
yet. Set `prefilter:store_aux_outputs=1` to store them anyway, raw, in
the `PrefilterLambdaK`/`PrefilterSegLogits` collections (see "What each
PhysicsEvent contains"), for offline inspection as the model develops.
Both are looked up by name, so their absence in an older/simpler model is
not an error.

```sh
-Peventbuilder:prefilter:onnx_model=/mnt/local/share/models/prefilter-20270803.onnx
-Peventbuilder:prefilter:score_threshold=0.5
-Peventbuilder:prefilter:max_hits=2500        # fixed input size models
-Peventbuilder:prefilter:veto=0               # score-only: store scores, drop nothing
-Peventbuilder:prefilter:store_aux_outputs=1  # also store lambda_k/seg_logits, if the model has them
```

## Outlook / TODO

This list rests on these measured facts: every count-based score
saturates at an AUC around 0.77; barrel tracklet pointing reaches AUC
0.864 with no Kalman filter; and background is over-dispersed (Fano 2-5
on raw input), because primaries produce multi-hit secondaries, so
absolute Poisson tails are nominal, not exact.

- **Tracklet pointing (`n_tracklets`, weights[10])** — pairs of gated fast
  hits, layer-aware `dphi` window (must absorb solenoid curvature, ~0.05-0.1
  rad), straight-line `z0` extrapolation to the beamline (barrel) or
  origin-consistency angle (endcap two-plane). Score-only, computed on fired
  windows; never in the trigger decision, so efficiency is untouched.
- **Calorimeter features** — cluster *counting* is measured useless
  (AUC 0.50: the slow gate floods with background). Test offline first:
  windowed energy sum / max-cluster E / E-weighted forward-backward
  asymmetry; cluster-time minus t0 residual; later tracklet-to-cluster
  matching.
- **Energy coincidence** — weight the time coincidence by hit eDep (hits
  already carry it; see `min_edep`) instead of pure counts; and/or require
  synchronous calorimeter energy. Test offline before any C++.
- **Honest significance under over-dispersion** — (i) empirical tail
  calibration: measure P(S >= s) directly on background-only frames
  (nonparametric, model-free; the number to quote); (ii) Fano-corrected
  z-score `(c-mu)/sqrt(F*mu)`, legitimate at 1x1 where mu ~ 15; (iii)
  negative-binomial null as the parametric option. Poisson tails and the
  adaptive `p_fake` target remain nominal labels, not rates.
- **Phase-resolved `mu(RF phase)` (real data)** — the simulation injects
  background uniformly in time (measured PSD: flat/white with a broadband
  2-3x burst excess, no lines), but a real machine gates everything to the
  bunch clock (~10 ns comb). With real data, estimate the background rate as
  a function of time-in-frame folded modulo the bunch period, and score each
  candidate against the rate at its own phase instead of the frame average —
  the counting-statistics analog of whitening/notching a known line.
  Common-mode for signal (collisions are bunch-locked too), so it sharpens
  calibration rather than separation, and makes out-of-time background
  (afterglow, activation) stand out.
- **Barrel-ECal clusters starved at Timeslice level** — the barrel produces
  far fewer clusters than the endcaps at Timeslice level, so E_calo
  (weights[9]) currently runs mostly on endcap information and
  barrel-only events lose their calo tag. Suspect: the merged-barrel
  cluster mirror (`EcalBarrelClusterFrame` -> TimeAlign) is missing an
  ingredient at Timeslice level; candidate fix is wiring
  `EcalBarrelScFiClusterFrame` (and Imaging) directly into the TimeAlign
  cluster list. Cluster times are fine where clusters exist.
- **Roman pots / off-momentum taggers** — RPOTS reconstruction exists
  (`src/detectors/RPOTS`, `ForwardRomanPotRecHits`) but is not wired into
  the builder; adding it (plus `ForwardOffMTrackerRecHits`) as far-forward
  tag inputs extends bit 6's coverage to lower-angle protons. NOTE the
  JOmniFactory even-variadic constraint: growing the tracker variadic to 12
  requires the calo variadic to grow to 12 as well (placeholders allowed).
  Acceptance gap caveat: protons between B0's upper edge (~20 mrad) and the
  central tracker's forward edge remain untaggable by ANY far-forward
  device — that band is a hardware gap, not a wiring one.
- **PSD monitor** — the averaged periodogram of the hit-count series (per
  region), as a standing detector-characterization panel in the audit: new
  lines = new coherent noise source; broadband steps = beam/vacuum changes.
- **Per-cell `mu_i` — IMPLEMENTED (inhomogeneous Poisson background).** The
  equal-solid-angle grid is kept, but equal solid angle does not imply equal
  expected background: the Poisson MEAN varies with detector position, each
  cell carries its own local null H_{0,i}: X_i ~ Poisson(mu_i) and (in
  adaptive mode) its own threshold k_i, with a cell -> cos-theta-band ->
  global statistics ladder (`topology:mu_mode`, `topology:mu_min_hits`).
  Estimator independence: trigger thresholds use the frame-wide estimate
  (self-contamination is small, around 2%, and irrelevant to firing); the
  stored significance is evaluated leave-one-out, with the candidate's own
  window excluded from the null it is scored against. This removes the
  dominant geometric bias at unchanged efficiency and fake rate. The
  residual p_tail non-uniformity is the signature of temporal correlations
  (burstiness), which no spatial mu map can absorb; that is the
  empirical-calibration/over-dispersed-model item above, not a defect of
  mu_i. Natural upgrade if per-frame statistics ever become limiting: a
  cross-frame running background per cell, through TimesliceBuffer_service.

## Configuration reference

Every parameter uses this canonical namespace:

```
eventbuilder:<component>:<parameter>
```

There are five components:

- `trigger`: event finding
- `topology`: the angular coincidence grid the trigger tests against
- `coincidence`: candidate-to-hit time gating
- `unfolder`: PhysicsEvent materialisation and cross-frame recovery
- `prefilter`: ONNX frame prefilter

Every parameter name is a single flat token. There is no further
`:`-nesting inside a component. So the path is always exactly three
segments. `topology` is split out of `trigger` for exactly this reason.
Its three parameters are implemented inside `EventBuilder_factory`,
alongside the rest of the trigger. But they are registered under their
own component, instead of nesting as
`eventbuilder:trigger:topology_threshold`. No bare parameters remain.
JANA may still print a spurious "parameter appears unused" warning for
perfectly-working parameters. This is a verified JANA false positive.

Trigger (`eventbuilder:trigger:*`):

| Parameter | Default | Description |
|---|---|---|
| `min_edep` | `0` | Min hit eDep (GeV) to vote in the trigger; 0 = off |
| `adaptive_threshold` | `0` | 1 = per-frame Poisson threshold (constant fake rate); 0 = fixed count |
| `p_fake` | `0.001` | adaptive mode: accidental probability per (cell x window); ~1e-5 targets ~1 fake/frame |
| `nsigma_window` | `3.0` | N-sigma confidence level for t0sigma AND for deriving the finding-stage coincidence windows `dt`/`dt_tof` per frame (`dt = nsigma_window * sqrt(2) * max hit resolution in the pool`); NOT the unfold-stage content gate, which is `eventbuilder:unfolder:nsigma_window`, a separate parameter |

Topology (`eventbuilder:topology:*`):

| Parameter | Default | Description |
|---|---|---|
| `threshold` | `2` | Min coincident hits in one angular cell to fire; also the physical floor in adaptive mode. Values `0`/`1` mean "no topology cut" and are clamped to this same floor: a single hit is not a coincidence by definition (a coincidence means two or more things happening together), and computing t0 as an inverse-variance-weighted estimate needs at least two hit times to combine. Raising it trades efficiency for a lower false-alarm rate (see "The trigger, in detail" above for the measured trade-off) |
| `theta_bins` | `12` | Topology grid: polar-angle bins of EQUAL SOLID ANGLE (uniform in cos(theta) over [-1,1], not equal 15-deg width) |
| `phi_bins` | `8` | Topology grid: azimuth bins over [0, 2pi) (45 deg cells at default); the half-cell stagger of the second grid follows automatically |
| `threshold_tof` | `2` | TOF-group floor; same semantics and clamp-to-2 as `threshold` above |

Coincidence (`eventbuilder:coincidence:*`):

| Parameter | Default | Description |
|---|---|---|
| `nsigma` | `3.0` | Gate half-width in units of each hit's own `getTimeError()` |

Unfolder (`eventbuilder:unfolder:*`):

| Parameter | Default | Description |
|---|---|---|
| `backward_frames` | `1` | Past frames merged for slow hits; 0 = off |
| `forward_frames` | `1` | Future frames merged (late-tail recovery); needs `nthreads >= 2`, degrades to 0 with a warning otherwise |
| `nsigma_window` | `3.0` | N-sigma half-width for per-detector hit acceptance (per-hit `getTimeError()`, always — no constant-override escape hatch) |
| `mc_time_window` | `50.0` | Half-width (ns) for MCParticles carried into each child; -1 = all |

Prefilter (`eventbuilder:prefilter:*`):

| Parameter | Default | Description |
|---|---|---|
| `onnx_model` | `""` | Model path; empty = pass-through |
| `score_threshold` | `0.5` | Max background probability to pass a frame (veto mode) |
| `veto` | `1` | 1 = drop below-threshold frames online; 0 = score-only (store `PrefilterScores`, drop nothing — reversible offline cut; recommended for production) |
| `max_hits` | `0` | Pad/truncate to fixed model input size; 0 = dynamic |
| `k_neighbors` | `16` | kNN graph k per DGCNN layer |
| `time_weight` | `1.0` | Time-axis scale in layer-0 spacetime kNN |
| `store_aux_outputs` | `0` | 1 = also store `PrefilterLambdaK`/`PrefilterSegLogits`, the model's raw `lambda_k`/`seg_logits` outputs, when it declares them. Neither is trained/meaningful yet |

Example:

```sh
eicrecon -Peventbuilder=true -Pplugins=eventbuilder \
    -Peventbuilder:trigger:adaptive_threshold=1 -Peventbuilder:trigger:p_fake=1e-5 \
    -Peventbuilder:topology:threshold=3 \
    -Peventbuilder:unfolder:forward_frames=1 \
    input_sim.edm4hep.root
```

Note there are NO global time-resolution parameters (a
`time_resolution_tof`/`time_resolution_mpgd` pair has circulated in older
drafts — those knobs do not exist): per-detector resolutions live only in
the digitization configs, inside each `src/detectors/<DET>/<DET>.cc`, and
every stage reads them per hit via `getTimeError()`.

## Notes

- Per-detector time resolutions live ONLY in each detector's own
  `src/detectors/<DET>/<DET>.cc` file. This file holds the digitization
  configs, for both the PhysicsEvent-level chain and its Timeslice-level
  "Frame" mirror, in the same file. Every gate and weight reads these
  resolutions per hit, through `getTimeError()`. There is no second copy
  that can drift out of sync.
- Do not let the PODIO writer probe the full standard factory set, the
  ACTS chain, on unfolded children. This aborts the process. Narrow
  `podio:output_collections` instead.

  The unfolder projects the gated `*RecHits`, `*RawHits`,
  `*RawHitLinks` (the single raw-hit truth carrier for both the tracker
  and calo families; ActsToTracks consumes them too), and TOF `Measurement2D`
  collections (`TOFBarrelClusterHits`/`TOFEndcapClusterHits`) into each
  child event, shadowing the standard factories through EulerianStore.
  The Measurement2D projection is gated around t0 by each measurement's
  own time error (`covariance.zz`). Factories whose outputs are not
  projected, the digitization intermediates, still abort if requested
  directly on children. Keep `podio:output_collections` narrowed.
- Frame-level registration exists for the central trigger trackers (TOF,
  MPGD, Si, B0), plus B0ECAL, BEMC, EEMC, FEMC, EHCAL, BHCAL, FHCAL, and
  ZDC. PFRICH and LUMISPECCAL also have frame-level code, but it stays
  disabled. PFRICH's frame-level registration aborts reconstruction if
  enabled, with the error `richgeo: ReadoutGeo is not defined for
  detector PFRICH`, on the craterlake beam compacts. This is why PFRICH
  is not in the active `m_simtrackerhit_collection_names` wiring in
  `eventbuilder.cc`, despite the code existing. DIRC, DRICH, LOWQ2, and
  RPOTS have no frame-level registration at all. To enable one, add a
  Timeslice mirror to its own `src/detectors/<DET>/<DET>.cc` file,
  following the `src/detectors/BTOF/BTOF.cc` pattern. See "Missing
  subsystems" above. None of these four detectors participate in event
  finding.
- **A JANA event-source generator competes on a hardcoded score, not on
  which one actually understands your flags.**
  `JEventSourcePODIO_generator`
  (`src/services/io/podio/JEventSourcePODIO_generator.h`) reads
  `-Peventbuilder=true` and switches the reader to Timeslice level. But
  JANA also auto-registers a generic
  `JEventSourceGeneratorT<JEventSourcePODIO>`, specialized in
  `JEventSourcePODIO.cc`. This generic generator returns a higher
  `CheckOpenable()` score for a valid PODIO file. It has no knowledge of
  `-Peventbuilder=true`. If it wins, `eicrecon` silently reads the file
  as ordinary PhysicsEvents. There is no error and no crash.
  eventbuilder's own parameters (`eventbuilder`, `eventbuilder:topology:*`,
  `eventbuilder:unfolder:*`) all just report "appears to be unused",
  since nothing ever touches them. `JEventSourcePODIO_generator::CheckOpenable()`
  is set to outscore the generic one. If this regresses again, look for
  `Creating event pool with level=PhysicsEvent` in the log, where you
  expected `level=Timeslice`.
- Diagnostics: run `make diagnostics <file|class>` in the production repo. This
  command renders purity, efficiency, applied gates, and S/S_TOF
  distributions. It also renders the offline efficiency-vs-false-alarm-rate
  curve, with an operating-point table. It renders these from any output
  file. It produces a PNG and a 2-page PDF report. It reads the
  configuration from the self-describing file.

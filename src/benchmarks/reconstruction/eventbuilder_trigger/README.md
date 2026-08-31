# `eventbuilder_trigger` — trigger performance benchmark

Measures what the software trigger actually achieved, in-process, and prints a
per-physics-class table while the job runs.

```sh
eicrecon -Peventbuilder=true \
         -Pplugins=eventbuilder,eventbuilder_trigger \
         -Peventbuilder:benchmark:report_every=250 \
         -Ppodio:output_collections="EventHeader,EventBuilderInfo,..." \
         frames.edm4hep.root
```

## Why in-process

The equivalent numbers used to come from re-reading the finished
`*.eicrecon.root` offline. That view has a structural blind spot: a frame whose
collisions the trigger missed entirely produces **no child event**, so it is
absent from the file and absent from the denominator. Efficiency measured that
way is biased high.

This plugin taps the Timeslice level, where every frame is visible — missed or
not — so the injected collisions of a zero-candidate frame still count against
the trigger.

## The two taps

| tap | level | supplies |
|---|---|---|
| `FrameBenchmark_processor` | `Timeslice` | every frame; per-class injected collisions; all STAGE 1 counters; the STAGE 2 GNN decision; the FAR numerator |
| `EventBenchmark_processor` | `PhysicsEvent` | STAGE 3 — ACTS tracks and calo clusters vs truth, per class |

Neither level can reach the other, so both fold into
`TriggerBenchmark_service`, which owns the table. Each tap registers at
`Init()` and reports at `Finish()`; the final table prints when the last one
finishes, so `Finish()` ordering does not matter.

Both use `CallbackStyle::ExpertMode` with `ProcessSequential` and
`EnableOrdering(true)`. Ordering is not needed for the totals — they are
order-free sums — but it *is* needed for the interim reports: without it,
"after 250 frames" would mean a different 250 frames at every thread count.

## Columns

Every stage prints one `eff/pur` cell.

| column | efficiency | purity |
|---|---|---|
| `time` | collisions recovered / collisions injected | real candidates / all candidates |
| `track` | real candidates passing `trk_pass` / real | `trk_pass` real / (`trk_pass` real + fake) |
| `calo` | real candidates passing `cal_pass` / real | same shape, `cal_pass` |
| `trigger` | real passing `trk_pass \|\| cal_pass` / real | same shape |
| `gnn` | frames not vetoed / frames scored | not-vetoed real / (real + false accepts) |
| `acts-trk` | majority-matched truth / in-acceptance truth | matched tracks / all tracks (ghost rate) |
| `calo-island` | matched truth groups / truth groups | matched clusters / all clusters |

`track`/`calo`/`trigger` are the **fast primitives** the trigger definition
specifies — `n_tracklets` and `E_calo` thresholded — not ACTS or CaloIsland
output. Those are the separate `acts-*` columns.

### The columns are per-stage, NOT a chain

**Each stage is denominated in the population that stage actually judged.**
`time` divides by injected collisions, because the clustering is what decides
whether a collision is recovered at all. `track`/`calo`/`trigger` divide by
real *candidates*, because the fast primitives are per-candidate quantities: a
collision the clustering never recovered has no candidate, hence no
`n_tracklets` and no `E_calo`, and cannot be said to have failed the calo term.

The consequence, and it is deliberate: **reading the row left to right does not
give you a chain.** `trigger` can read 100% next to a `time` of 93.5% — that
says the trigger accepted every candidate put in front of it, not that it
recovered every collision. Each column answers a question about its own stage,
which is what you tune a threshold against; a cumulative column would confound
the primitive's discrimination with the clustering's efficiency.

The composed figure is printed once, in the footer:

```
chain    TIME && (TRACK || CALO) = 762 / 827 = 92.1%  (end-to-end; the stage columns are per-stage)
```

That numerator counts distinct injected collisions recovered by a candidate
that also cleared the trigger. The terms are OR-ed across every candidate that
recovered a given collision rather than read off whichever one was seen first —
a collision is track-triggered if ANY candidate that found it cleared the
threshold.

The `FAKE-ACCEPT` row is how often each stage accepts a fake. Its denominator
is the fake population, so lower is better and it has no purity half.

> With the default `eventbuilder:trigger:min_cal_energy=0` the calo term is
> "off" — every candidate clears a zero threshold — so `cal_pass` is a
> tautology and `calo`, hence `trigger`, read 100% efficiency *and* 100%
> fake-accept. `min_tracklets=1` is nearly as weak on real candidates, though
> it does reject ~89% of fakes. At the defaults STAGE 1 therefore discriminates
> almost nothing real: that is the configuration reporting itself faithfully,
> not a bug. Set real thresholds to make those columns mean something, and see
> the `HcalEndcapN` note below before you pick a value.

> **`e_calo` is missing one calorimeter.** `HcalEndcapN` clusters are produced
> by EICrecon but are absent from `m_clu_in_names`/`m_clu_out_names` and from
> the time-alignment list in `eventbuilder.cc` (its only mention there is
> commented out), while its rec hits *are* wired in. Backward hadronic energy
> therefore never reaches `E_calo`, biasing it low. Harmless while the
> threshold is 0; it matters as soon as one is set.

### Per-class rows

Every stage column carries an efficiency for each class, under the same
per-stage rule as the totals — `time` over the class's injected collisions, the
rest over the class's own real candidates:

| column | per-class efficiency |
|---|---|
| `time` | collisions of this class recovered / injected |
| `track` `calo` `trigger` | real candidates carrying this class that passed `trk_pass` / `cal_pass` / either |
| `acts-trk` `calo-island` | majority-matched truth of this class / in-acceptance truth of this class |

### Per-class purity

Fakes and ghosts carry no class, so the usual real/(real+fake) purity has no
per-class form. But two contaminations ARE class-attributable, and those are
what the per-class purity halves report:

| column | per-class purity | the shortfall is |
|---|---|---|
| `time` | `found / real_cands` | duplicate candidates claiming the same collision |
| `track` `calo` `trigger` | `found_x / x_ok` | the same duplicates, restricted to candidates passing that term |
| `acts-trk` `calo-island` | `owned_by_credited / all_reconstructed` | tracks and clusters whose majority owner is NOT credited in-acceptance truth: another collision in the gate, a secondary, or a particle outside acceptance |

> Note what this purity is **not**. Every association carries a normalized
> hit fraction and one row per contributing particle, so "has an owner with
> weight > 0.5" is satisfied by essentially every track in a gated child
> event -- that ratio is 100% by construction and measures nothing. The
> denominator here is therefore ALL reconstructed objects, and the numerator
> only those owned by truth the candidate is credited with.
>
> The calo figure is low (~3%) because a child carries every gated cluster,
> most of which are low-energy junk from no credited neutral. Read it as
> cluster contamination, not as a clustering-quality metric.

The `TOTAL` row instead uses the conventional definitions -- real/(real+fake)
for the trigger stages, majority/all for ACTS -- because at that level fakes
and ghosts do count. The two are different quantities; do not read a class row
and the TOTAL row as the same measurement.

## Denominators## Denominators

- **Tracking** — stable charged MC with `pt > eventbuilder:findable:pt_min` and
  `|eta| < eventbuilder:findable:eta_abs`, mirroring the eventbuilder's own
  "findable" definition so the two agree.
- **Calorimetry** — stable neutral MC (neutrinos excluded), `|eta| < 4.0`, and
  `E > 50 MeV` for `|eta| < 1` else `30 MeV` — the Yellow Report coverage and
  minimum detectable photon energy. Targets are **angularly grouped** by
  union-find at `dR < eventbuilder:benchmark:calo_dr`: two photons merged into
  one shower cannot be resolved separately by any clustering algorithm, so
  counting them as two targets would cap efficiency below 100% even for
  perfect reconstruction.
- **Injected collisions** — frame `MCParticles` with a stable physics status,
  deduplicated by `(class, time)` within the frame's own `dt`, the same rule
  `EventBuilder_factory`'s `collision_times` loop uses.

## Standalone mode

Pointing the plugin at an already-written file works:

```sh
eicrecon -Pplugins=eventbuilder_trigger existing.eicrecon.root
```

`JEventSourcePODIO` inserts every stored collection, and JANA reuses the
matching factory's empty databundle rather than clashing, so the written
children replay as PhysicsEvents. The plugin detects the absence of a
Timeslice parent and stamps the header **`standalone (reduced denominators)`** —
in this mode the blind spot described above returns, because zero-candidate
frames were never written. STAGE 3 columns render `--` for any association
collection that was not in `podio:output_collections`.

## Why the position residual is decomposed

Do NOT use the 3D `|rec - sim|` magnitude as a tracker position resolution.
The barrel readout is a `CylindricalGridPhiZ` pinned to a FIXED radius per
segmentation -- and a module stacks several of them (`InnerPhi` at
`R - MMThinGap`, `InnerZ` at `R + MMThinGap`, `BotRad`/`TopRad` further out).
Every rec hit is therefore projected onto one of those cylinders while its sim
hit sits at its true radius, so the radial residual is QUANTIZED: measured on
both barrel MPGDs it is three discrete spikes at -1.5, 0 and +1.5 mm.

Folded into a 3D magnitude that geometry masquerades as a resolution tail. It
inflated `MPGDBarrel` to 0.32 mm and `OuterMPGDBarrel` to 0.54 mm -- a 1.7x
disagreement between two layers that are built the same way. Split out, both
read 0.26 mm and agree to better than 1%.

So `space_*` is the in-plane component (the reconstruction resolution) and
`radial_*` is kept as its own histogram, where the quantization is visible for
what it is.

## Resolutions

With `-Phistsfile=<file>` the plugin also books per-detector, per-class
residuals under `benchmark_eventbuilder/`. Each is a TH2 with the physics class
on the y axis (labelled bins), so a 19-detector x 26-class breakdown is 19
objects, and one class is a `ProjectionX` on a named bin.

| histogram | quantity | detectors |
|---|---|---|
| `time_<det>` | `t_hit - t0` [ns] | 10 trackers + 9 calo systems |
| `space_<det>` | `\|rec - sim\|` [mm] | trackers (needs `truth:simhits=1`, the default) |
| `space_<calo>` | `dR(cluster, MC)` | calo — angular, not metric |
| `energy_<calo>` | `(E_rec - E_MC) / E_MC` | calo |

### The `time` column is a SENSOR resolution

`time` fits `t_rec + |r|/c - t_sim`: the reconstructed hit time with the
event builder's propagation correction undone, against the sim hit's own time.
Both corrections are needed, and getting either wrong is what made this column
read `--` for TOF for a long time.

**Why undo `|r|/c`.** The rec hits reaching the child are TIME-ALIGNED --
`TimeAlignment_factory::correct()` subtracts `|r|/c`, the straight-line beta=1
propagation time, so hits of one collision share a common time. Comparing an
aligned hit with an unaligned sim hit measures that correction, not the sensor:
-2.95 ns on TOFBarrel, -6.35 ns on TOFEndcap, tracking flight path exactly
(tight at the fixed-z endcap, spread across the barrel's z range).

**Why not reference the collision time.** `t_hit - t_MC` on an aligned hit is
`TOF - |r|/c` -- how well the straight-line beta=1 assumption recovers t0. That
is a real quantity, but it is one-sided and a few ns wide, dominated by track
curvature in the 1.7 T field and by beta<1. There is no Gaussian core, and the
fit correctly declines rather than reporting that spread as a resolution.

Verified from the digitiser's own debug output
(`-PBTOF:TOFBarrelRawHitFrame:LogLevel=debug`), per cellID:

| check | result |
|---|---|
| `hit_time_stamp - sim_time*1000` | median **1.6 ps**, p16/p84 -20.7/+27.2 ps -- exactly the configured 25 ps smear |
| child `SimHits.time` vs the time the digitiser saw | **0.0000 ns** -- exact |
| child `RecHits.time` vs the raw timestamp | **-2.53 ns** -- the alignment correction |

So the digitisation is faithful and the sim hits are copied exactly; the whole
discrepancy was the propagation correction.

**Cost:** the residual needs a sim hit, and only hits whose MCParticle passes
the `mc_time_window` gate carry one (~12% for TOF, by design -- see
`eventbuilder.cc`). So the `entries` count for a tracker row is that subset,
not every labelled hit.

**The fit is validated against a known answer.** TOFBarrel's residual is a
pure digitisation smear of a configured 25 ps, so that row has ground truth,
and `fitCore` returns **22.5 ps** on 469 entries. It is the only row in this
table whose correct value is known independently, which makes it the one worth
re-checking after any change to the fit.

Getting there needed a fix to `fitCore` itself: it used to fit `gaus(0)+gaus(3)`
over the full axis and report whichever component came out narrower. On a clean
single Gaussian those two are degenerate -- the fit parks a narrow spike on the
top bins, lets the broad one absorb the rest, and reports the spike (TOFBarrel
came back as 9.78 ps, B0Tracker as 0.825 ns, ten times better than its silicon
siblings). It is now a single Gaussian iterated over +-2.5 sigma of the peak,
which rejects the coincidence-gate tails by not fitting them rather than by
absorbing them into a second component.

TOFEndcap still reads `--`: 179 entries against fitCore's 200 minimum.

### Time reference

Residuals are measured against **`weights[MC_T]`, the true collision time**, not
against `t0`. `t0` is the inverse-variance weighted mean of the very TOF/MPGD
hits being measured, so a `t0`-referenced residual is correlated with its own
reference and collapses toward zero when one precise hit dominates the
weighting. Fakes matched no collision and have no truth time, so they are
skipped. Note this makes the width a *system* resolution: it includes
flight-path spread and any per-layer offset the r/c alignment left behind, not
sensor jitter alone.

### TOF: two chains, two rows

TOF is digitized twice, and the table reports both:

| row | collections | what it measures |
|---|---|---|
| `TOFBarrel` / `TOFEndcap` | `*RecHits` + `*RawHitLinks` | the plain one-hit-per-deposit chain — the GNN's per-hit truth |
| `TOFBarrelShared` / `TOFEndcapShared` | `*ClusterHits` (Measurement2D) + `*SharedRawHitLinks` | the charge-sharing chain ACTS consumes, via `SiliconChargeSharing` -> `LGADHitClustering` |

They are not interchangeable: measured on the same data the plain chain gives a
3.5 mm position residual (deposit-to-rechit along the strip) and the Shared
chain 14 um (the pad resolution). `tracking.cc` makes the same distinction --
its `CentralTrackingRawHitLinks` collector takes `*SharedRawHitLinks` for TOF
and `*RawHitLinks` for everything else.

The class label is taken **per hit** on the tracker side, from the
`RawHitLinks` weight (`generatorStatus + slot*1e6`), so a hit is attributed to
the collision that produced it rather than to whatever mix its candidate was
credited with. On the calo side it comes from the matched MC particle's own
status.

> Read the time widths with care: the gated hits only span the candidate's own
> window (`dt = nsigma*sqrt(2)*max sigma`, tens of ns), so the RMS is
> gate-limited and is **not** the intrinsic detector resolution. Fit the core
> peak for that. Likewise the calo energy residual compares cluster energy
> against the full MC particle energy, with no sampling-fraction or
> containment correction applied.

Nothing is booked unless `histsfile` is set, so the default run pays nothing.

## PDF report

`-Peventbuilder:benchmark:pdf=<file>` writes a multi-page report:

1. title, headline efficiency/purity, run summary, generation time
2. software provenance (ROOT, JANA2, podio, EDM4eic, compiler) and detector config
3. per-class efficiency as a traffic-light bar chart, annotated with `(found/injected)`
4. fitted resolution per subsystem
5. the statistical table verbatim
6. residual plots, grouped by detector, four per page

The plots are booked in memory when only a PDF is asked for, so no ROOT file
is written unless `-Phistsfile` is also given.

## Parameters

| parameter | default | meaning |
|---|---|---|
| `eventbuilder:benchmark:report_every` | `250` | interim table every N frames; `0` = final only |
| `eventbuilder:benchmark:stage3` | `1` | score ACTS + calo |
| `eventbuilder:benchmark:frame_ns` | `2000` | frame length, for fakes/frame → Hz |
| `eventbuilder:benchmark:calo_dr` | `0.05` | neutral truth-group merge radius |
| `eventbuilder:benchmark:csv` | *(empty)* | write the final table as CSV |

## Conventions

Class table, `generatorStatus` decode, and all weight indices come from
[`TruthClassLabels.h`](../../../factories/eventbuilder/TruthClassLabels.h) —
the single in-repo source, so this plugin cannot drift from the factory that
writes the weights. See `src/factories/eventbuilder/README.md` for the prose.

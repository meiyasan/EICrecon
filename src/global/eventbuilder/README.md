# `eventbuilder` — streaming event building for ePIC

JANA plugin name `eventbuilder`. Layout follows the EICrecon convention —
`JOmniFactory` definitions live under `factories/`, plugin wiring under `global/`:

- `src/factories/eventbuilder/` — `EventBuilder_factory.h`, `TimeAlign_factory.h`
- `src/global/eventbuilder/` — `eventbuilder.cc` (`InitPlugin`), `EventUnfolder.h`
  (a `JEventUnfolder`, not a factory), and `digi/Digitize*Plugins.cc` (the
  per-detector digitization registration).

## Collection naming

The plugin's Timeslice-level collections carry a marker so they don't clash with
the standard PhysicsEvent-level ones (JANA keys factories by type+tag, not level):
**`<Detector><RawHit|RecHit|Cluster>Digi`** for the digitized collections (e.g.
`TOFBarrelRecHitDigi`, `TOFBarrelRawHitDigi`, `B0ECalClusterDigi`) and
**`<Detector>TimeAlign<RecHits|Clusters>`** for the time-aligned ones (e.g.
`TOFBarrelTimeAlignRecHits`).

This plugin turns a continuous **time-frame** of streaming-readout data into
discrete **physics events** in software (there is no hardware trigger in the
ePIC streaming DAQ). Enable it by reading the input at the time-frame level:

```sh
eicrecon -Peventbuilder=true -Pplugins=eventbuilder <input_sim.edm4hep.root>
#   -Peventbld=true   is accepted as an alias of -Peventbuilder=true
```

`-Peventbuilder=true` makes the PODIO source deliver each entry at
`JEventLevel::Timeslice` (the 2 µs time-frame) instead of `PhysicsEvent`. The
plugin is loaded by default but stays **idle** without this flag (its factories
sit at Timeslice level with no time-frame parent to act on), so default
reconstruction is unaffected.

## Pipeline

```
time-frame  (JANA Timeslice level)  ── -Peventbuilder=true ──┐
  └─ digitization (digi/Digitize*Plugins.cc, registered @ Timeslice)
       → *RawHitDigi / *RecHitDigi / *ClusterDigi …
  └─ TimeAlign_factory (TrkTimeAlign_factory / CalTimeAlign_factory)
       → *TimeAlignRecHits / *TimeAlignClusters
  └─ EventBuilder_factory   (finds candidate events, assigns t0 + dt0)
       → EventCandidates
  └─ EventUnfolder         (materialises one PhysicsEvent per candidate)
       → PhysicsEvent (EventHeader carries t0/dt0)
            └─ full reconstruction (ACTS, …) runs here
```

- **`EventBuilder_factory`** (`JOmniFactory` @ Timeslice): groups time-coincident
  hits from the fast trigger detectors (TOF + MPGD) with an O(N) sliding-window
  scan plus a 12×8 + 12×8 staggered θ/φ topology test, and computes the event
  time by inverse-variance weighting of the coincident hit times:
  `t0 = Σ(tᵢ/σᵢ²)/Σ(1/σᵢ²)`, `dt0 = nsigma · 1/√(Σ 1/σᵢ²)`. TOF (~30 ps)
  dominates, so **t0 is set by the fastest detector**. Analysis only.
- **`EventUnfolder`** (`JEventUnfolder`): structural only. Reads the candidate
  list and emits one `PhysicsEvent` per candidate, selecting each detector's hits
  inside its resolution-matched window `[t0 ± nsigma·σ_det]` (tight for TOF, wide
  for Si) and writing the real t0/dt0 into the child `EventHeader`.
- **`TimeAlign_factory<HitT>`**: one templated factory (aliases `TrkTimeAlign_factory`
  / `CalTimeAlign_factory`) that subtracts a per-hit time-of-flight term so all hit
  times refer to a common collision t0, then sorts by time — for tracker hits and
  calo clusters alike.
- **`digi/Digitize*Plugins.cc`**: re-register the standard ePIC digitization +
  hit-reconstruction (and calo clustering) factories at the **Timeslice** level,
  producing the `*Digi` collections the chain consumes. Only the central trigger
  trackers (TOF, MPGD, Si, B0) plus a few calorimeters are registered;
  far-forward/backward and Cherenkov detectors are intentionally excluded from
  event finding.

## Why PhysicsEvents, not fixed time-slices

A naive fixed, non-overlapping time-slice grid cannot work here: detector time
resolutions span five orders of magnitude (TOF ~30 ps … Si ~2 µs). Any single
window is either too tight (drops Si hits) or too wide (merges neighbouring
events and background). The EventBuilder instead defines an event by **hit
membership** — a precise t0 and a *per-detector* `±nsigma·σ_det` window — so each
detector contributes over its own resolution and events may overlap in time
(pileup). The input arrives as a JANA `Timeslice`; the output is `PhysicsEvent`s.

## Findings and changes (this branch)

1. **Replaced the monolithic `TimeframeSplitter`** (which fused analysis +
   unfolding and carried several bugs: destructive/no-overlap hit assignment; a
   coincidence window dominated by the slowest detector; an out-of-bounds `erase`
   using a floating-point time as an index; a `||`-vs-`&&` mistake; a missing
   staggered-bin increment; an integer timestamp counter instead of a real time;
   a stale hard-coded termination index) with `EventBuilder_factory` (analysis) +
   `EventUnfolder` (materialisation).
2. **t0/dt0 are now real and propagated.** The child `EventHeader` carries
   `timeStamp` = t0 in ps and `weights = {t0_ns, dt0_ns, phys_flag}`. Verified:
   signal-only DIS-CC → 2 frames give 2 events (t0 ≈ 0.61/0.96 ns, dt0 ≈ 20 ps);
   with background overlay → 2 frames give ~50 events (signal + background fake
   triggers spread across each frame) — the motivation for the ONNX pre-filter.
3. **Root infrastructure bug found (pre-existing): digitization registered at the
   wrong event level.** `digiBTOF/ECTOF/MPGD` and the calo digi plugins used the
   *positional* `JOmniFactoryGeneratorT` constructor, which defaults factories to
   `PhysicsEvent`; only the Si/B0 trackers used the `TypedWiring` form with
   `.level = Timeslice`. Since the trigger relies on TOF + MPGD, those *had* to be
   at Timeslice. A half-finished migration in the original branch (note the
   `// TODO: Remove me once fixed` markers in `digits/digiBEMCPlugins.cc`); the
   EventBuilder/EventUnfolder is just the first code to read these collections
   from the time-frame end-to-end, so it surfaced the latent bug. Fix: a chainable
   `JOmniFactoryGeneratorT::SetLevel(JEventLevel)` plus a `Timeslice` wrap on every
   TOF/MPGD/calo registration.
4. **JANA detail worth knowing:** `JEventUnfolder::DoUnfold` fetches *all* of an
   unfolder's inputs from the **parent** time-frame before user code runs, so any
   input collection produced at the wrong level aborts the unfold immediately.
5. **Parameter renamed:** `split_timeframes` → `eventbuilder` (alias `eventbld`).

## Removed components

- **`TimeCoincidenceFactory.h` — removed.** It was an early attempt to thin out
  noise by applying fixed time-window cuts directly on the offset-corrected,
  time-sorted hit collections. **Superseded by `EventBuilder_factory` +
  `EventUnfolder`**: a real inverse-variance t0, a *per-detector* `±nsigma·σ_det`
  acceptance window (so fast and slow detectors are each handled correctly), and a
  θ/φ topology trigger to reject random coincidences — none of which a single
  fixed time cut can do. Removed as dead code; this note records its intent.
- **`HitChecker.h` — removed.** A small debug factory that printed hit times; not
  part of the production chain.

## Known limitations / TODO

- **Calo carry-through is temporarily disabled** in `EventUnfolder` (see the
  `TODO(calo)` there), pending a clean-rebuild verification of the calo Timeslice
  migration (finding 3). The EventBuilder uses only trackers, so t0/dt0 are
  unaffected.
- **Background pre-filter (Phase B):** a shared ONNX runtime service +
  `--prefilter-onnx-model` to reject the background fake-triggers (the ~50-vs-2
  gold-old result) before reconstruction — not yet implemented.

> Note: JANA keys factories by `(type, tag)` only — *not* by event level — so a
> Timeslice and a PhysicsEvent collection of the same name would collide. That is
> why the Timeslice collections must carry the `Digi` / `TimeAlign` markers rather
> than reuse the standard names.

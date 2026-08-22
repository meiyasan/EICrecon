# Eventbuilder — structural guide

This file is the compact map of the eventbuilder. It states what each part
does, the data conventions, and the rules that must not break. The long
README (`../../global/eventbuilder/README.md`) holds derivations and
measurements. When the two disagree, fix the README.

## What the eventbuilder does

Input: one 2 µs Timeslice ("frame") of streaming-readout hits, all
detectors, no trigger. Output: zero or more PhysicsEvents ("candidates"),
each a time-windowed slice of the frame with its own reconstructed and
truth collections.

Two stages:

1. `EventBuilder_factory.h` (this directory) finds candidates. It runs at
   Timeslice level. It writes `EventCandidates` (one entry per candidate)
   and `EventBuilderFrameInfo` (one entry per frame).
2. `EventUnfolder` (`../../global/eventbuilder/eventbuilder.cc`) splits the
   frame into one child event per candidate. It copies the gated hits,
   truth, and metadata into each child.

## File map

| File | Role |
|---|---|
| `EventBuilder_factory.h` | trigger + candidate finding, weights encoding |
| `EventPrefilter_factory.h` | optional GNN frame filter (ONNX), `PrefilterScores` |
| `TimeAlignment_factory.h` | per-detector time alignment to the frame clock |
| `TimeCoincidence_factory.h` | gates slow/fast RecHits to candidate windows |
| `../../global/eventbuilder/eventbuilder.cc` | `EventUnfolder`: frame → child events |
| `../../tests/eventbuilder/EventBuilderTest.cpp` | pure-C++ unit tests (no JANA/PODIO/ROOT) |

The dataset-production pipeline (generation, hit-level noise overlay,
pile-up instance encoding, status decoding, reporting, sidecar writing)
lives in a separate production repository, outside this codebase.

## Two-stream production (how frames are made)

Physics is simulated ONCE. Noise is simulated per level and overlaid at hit
level. This is valid because Geant4 tracks never interact across
collisions; the only coupling is digitization, which runs after the merge.

```
Stage A  SBM mixes 26 physics classes (no machine bkg) -> npsim  = "vanilla"
Stage B  machine background only, per noise level      -> npsim  = "bkg-stream"
Stage C  hit-level merger overlays frame i + frame i   -> eicrecon
```

The physics content is bit-identical across noise levels by construction.

## generatorStatus — the single truth convention

| Range | Meaning |
|---|---|
| < 2000 | native codes and Geant4 secondaries (secondaries do NOT inherit a class band) |
| 2000–6999 | machine background, 1000-wide per source: 2000 SynRad, 3000 e-beam-gas, 4000 Coulomb, 5000 Touschek, 6000 p-beam-gas |
| >= 10000 | physics classes, **1000-wide per class** (LIVE): `base = 10000 + class_index*1000`, max 35999 |

Class order (index 0–25): ncdisq1, ncdisq10, ncdisq100, ncdisq1000,
ccdisq100, ccdisq1000, ddis, dvcs, ddvcs, dvmp, dempq3to10, dempq10to20,
dempq20to35, tcs, jpsi, jpsiphoto, urho, cohrho, upi0, upsilon, mesonsf,
cohphi, rho, spectroscopy, omega, photoprod.

Decode: `class_index = (gs - 10000) // 1000`, `stable = gs % 1000 == 1`.
The stable test that covers legacy trees too: `gs == 1 or (gs >= 7000 and gs % 1000 == 1)`.

**Same-class pile-up instance encoding is designed but NOT LIVE.** A wider,
10000-wide-per-class layout (`base = 10000 + class_index*10000`, instance n
stable = `base + (n-1)*1000 + 1`, decay = `+2`) exists in the production
repo's `rewrite_status.py`, meant to run between the SBM merge and npsim.
It unconditionally self-disables today because SBM writes its merge as
`.hepmc3.tree.root` (ROOT format), which the rewriter cannot parse — so it
never actually runs, and every real status is instance 1 of the 1000-wide
layout above (confirmed against `signal.sh`'s own `_CLASS_STATUS` base
values, which are 1000 apart, not 10000). Do not decode a live status as if
instance sub-bands are populated.

## Pile-up — how it is interpreted

A "collision" is one injected physics event. Collisions are deduplicated by
(class, time) within the frame's own coincidence window `dt`.

- `flag` (weights[2]) = COUNT of real collisions inside the candidate's
  window. 0 = fake. 1 = normal. 2+ = true pile-up. "Real" is `flag > 0`.
- `mc_t` (weights[3]) = the EARLIEST matched collision's time.
- `trigger_classes_mask` (weights[17]) = bitmask, bit N = class_index N
  had a collision in the window.
- `trigger_classes` (weights[25+], on by default) = count, then
  (time, stream) pairs. `stream = generatorStatus / 1000` of the collision
  base: 0 = native primary, >= 10 = physics with `class_index = stream - 10`.
- Per-hit collision slot: the unfolder bakes
  `weight = generatorStatus + slot*1e6` into truth links, where slot is the
  nearest `trigger_classes` entry to the hit's own particle time. slot must
  stay < 16 (float32 exact-integer bound; asserted at both encode sites).

## The trigger (candidate finding), in short

1. Fast hits = TOF (~30 ps) + MPGD (~10 ns) only. Silicon and calo never
   form windows.
2. `dt = nsigma_window * sqrt(2) * max(sigma)` per frame; `dt_tof` same,
   TOF pool only.
3. Sliding window over a (theta, phi) cell grid (two staggered grids). A
   cell fires when its count >= its own threshold `k_i` (fixed, or
   adaptive: smallest k with Poisson tail < `p_fake`, floored at
   `topology:threshold`).
4. TOF-only group can fire independently (`cmax_tof >= k_tof`).
5. t0 = inverse-variance weighted hit-time mean; t0sigma from the same.
6. MC matching is an inclusive count (see Pile-up above), NOT an exclusive
   nearest-match — the old 1/2/3 flag scheme is gone.

## Candidate weights (EventBuilderInfo, one row per candidate)

0 t0 | 1 t0sigma | 2 flag | 3 mc_t | 4 S | 5 S_tof | 6 p_tail | 7 cmax |
8 cmax_tof | 9 E_calo | 10 n_tracklets | 11 n_expected |
12 trk_threshold_cfg | 13 trk_pass | 14 cal_threshold_cfg | 15 cal_pass |
16 trigger bitmask | 17 trigger_classes_mask | 18 mu (winning cell) |
19 k (winning cell) | 20 win_theta | 21 win_phi | 22 win_grid |
23 S_empirical | 24 fano | 25+ trigger_classes tail (count, pairs).

`EventBuilderFrameInfo` (one row per frame, `eventNumber = frame*1000`):
0 n_physics_events | 1 mu_tof | 2 k_tof | 3 dt | 4 dt_tof |
5 topology:threshold | 6 theta_bins | 7 phi_bins | 8 adaptive_flag |
9 p_fake. Join to candidates via `eventNumber // 1000` (candidates carry
`frame*1000 + candidate`).

## Trigger bitmask (weights[16])

| Bit | Meaning |
|---|---|
| 0 | TOF-only group fired (cmax_tof >= k_tof) |
| 1 | combined TOF+MPGD group fired (cmax >= k_win) |
| 2 | backward-ECal region tag (in-window cluster, nhits >= 10) |
| 3 | bit 2 + >= 1 backward endcap tracklet |
| 4 | barrel-ECal region tag |
| 5 | bit 4 + barrel tracklet |
| 6 | forward-ECal region tag |
| 7 | bit 6 + forward endcap tracklet |
| 8 | B0 tag (>= 4 in-window B0 hits) |
| 9 | ZDC tag (cluster hits sum >= 50) |

Combined tags computed offline, never stored: cTag1 = (#set{2,4,6} >= 2)
&& bit8; cTag2 = (#set{3,5,7} >= 2) && bit8; cTag3/4 same with bit9;
cTag5 = (#set{2,4,6} == 3); cTag6 = (#set{3,5,7} >= 2).

## Truth carriers in child events

| Collection | Weight | Purpose |
|---|---|---|
| tracker `*RawHitLinks` (10 detectors) | `generatorStatus + slot*1e6` | THE per-hit class label (GNN training) |
| calo `*RawHitLinks` (11 detectors) | `generatorStatus + slot*1e6` | same, calo side; label from the highest-energy contribution in the cell |
| tracker + calo `*RawHitAssociations` | 1.0 (untouched) | internal food for ACTS truth matching only; not in SEGMENTATION output |

Both families use `podio::Link` as the single truth carrier. The calo type is
`edm4eic::MCRecoCalorimeterHitLink`
(`LinkCollection<RawCalorimeterHit, SimCalorimeterHit>`), filled by
`CalorimeterHitDigi` for all 11 detectors.

A text sidecar `<stem>.trigger_meta.txt` accompanies every eicrecon file:
`frame candidate t0_ns TRUE|FAKE class:time [class:time ...]`.

### Full truth chain (on by default)

The link's `to` side (`SimTrackerHit`) lives in the parent Timeslice frame,
which is never written into a child event — so
`RawHitLink -> SimTrackerHit -> MCParticle` does NOT resolve offline by
default. That is why the label is baked into the link weight instead: it
needs no pointer chase.

When you need the real chain (exact particle kinematics, energy deposits,
walking Geant4 secondaries back to their primary), one flag is enough:

```
-Peventbuilder:truth:simhits=1     # on by default; set 0 for GNN-only trees
```

`simhits=1` writes each candidate's `*SimHits` (10 tracker collections,
deduplicated per detector) and re-points the links at those child-local
copies.

A sim hit is written only when its MCParticle passes the same gate the
MCParticles obey, so a child event contains exactly what belongs to it and
every written sim hit resolves to a particle that is also present:
`simhit -> MCParticle` is 100% on every detector, by construction.

Do NOT add `mc_time_window=-1` for this. It is not needed and it is
expensive: it pulls the whole frame's particle list into every candidate,
and candidate windows overlap. Measured on one frame — 1121 MB of
MCParticles (96.6% of a 1161 MB file) to make 2.4% of them reachable,
against 113 MB with the gate. `-1` remains available when you genuinely
want every particle, not as a way to resolve the chain.

Ground truth does not depend on any of this. The per-hit label lives in the
RawHitLink weight and is written for every kept hit whether or not its
particle belongs to this candidate — 878323 links labelled versus 216402
sim hits written, in the same measurement.

### MCParticle gate

`eventbuilder:truth:mc_time_window` (default 50 ns) bounds which
MCParticles are cloned into a child. The effective gate is
`max(mc_time_window, dt)` — flag/trigger_classes count every collision
within +-dt of t0, so a narrower particle gate would drop the products of a
collision the candidate is credited with. dt is derived per frame from the
measured hit resolution (`3 * sqrt(2) * max_sigma`, about 42 ns for a
10 ns MPGD), so the fixed 50 ns alone would only satisfy that invariant by
coincidence for one detector mix.

## Streaming mode (one long-lived eicrecon per production stream)

`JEventSourcePODIO` (`src/services/io/podio/`) has a watch mode: after the
input file is exhausted, it polls a directory for new complete
`*.edm4hep.root` files (fEND-vs-size truncation check plus a settle delay)
and continues with them. Geometry loads once per process.

```
eicrecon ... -Ppodio:watch_directory=<sim dir> \
             -Ppodio:watch_max_files=10 \        # optional output rotation
             -Ppodio:watch_poll_ms=5000 \
             first_file.edm4hep.root
```

The first file argument stays required (it defines the schema). All new
files append to the single `podio:output_file` — the output is readable
only after the process finishes, so use `watch_max_files` to rotate when
downstream tools must follow along. `podio:watch_settle_seconds` (default
3) avoids opening files a writer has not finished.

## Rules that must not break

1. The class table and all decode arithmetic live in this file and in
   `EventBuilder_factory.h`'s comments. The external production pipeline
   keeps matching copies (its status assigner and its python decoder) —
   any change here must be coordinated with them.
2. Never re-add a duplicate of `flag` (the old weights[35]
   `n_mc_coincident` was removed for exactly that reason). Never append
   anything after the `trigger_classes` tail — it is variable-length and
   must stay last.
3. `eventNumber = frame*1000 + candidate` is the universal join key.
   The unfolder warns loudly at 1000 candidates/frame.
4. slot*1e6 packing requires slot < 16 and max generatorStatus <= 269999
   (float32 exact-integer bound, 2^24). Live max is 35999 (1000-wide
   scheme); 269999 is the dormant wide scheme's max, kept as the bound
   here so it does not need revisiting if that scheme goes live.
5. Consumers that hardcode weights indices: `EventBuilderTest.cpp` and
   `EventPrefilter_factory.h` (`build_cand_info`) in this repo, plus the
   external analysis/reporting tooling. A renumber must touch all of them
   in one coordinated pass.
6. Fakes carry no class (mask 0, empty tail): per-class purity does not
   exist; purity is per-trigger only.

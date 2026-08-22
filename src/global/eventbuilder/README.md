# `eventbuilder` — streaming event building for ePIC

JANA plugin name `eventbuilder`. Layout follows the EICrecon convention —
`JOmniFactory` definitions live under `factories/`, plugin wiring under `global/`:

- `src/factories/eventbuilder/` — `EventBuilder_factory.h`, `TimeAlignment_factory.h`, `TimeCoincidence_factory.h`
- `src/global/eventbuilder/` — `eventbuilder.cc` (`InitPlugin` + inline `EventUnfolder`
  struct, a `JEventUnfolder` not a factory), and `plugins/*.cc` (per-detector
  digitization registration, one file per detector group).

## Collection naming

The plugin's Timeslice-level collections carry a marker so they don't clash with
the standard PhysicsEvent-level ones (JANA keys factories by type+tag, not level):
**`<Detector><RawHit|RecHit|Cluster>Digi`** for the digitized collections (e.g.
`TOFBarrelRecHitDigi`, `TOFBarrelRawHitDigi`, `B0ECalClusterDigi`),
**`<Detector>TimeAlign<RecHits|Clusters>`** for the time-aligned ones (e.g.
`TOFBarrelTimeAlignRecHits` — always the current frame only),
and **`<Detector>TimeCoincRecHits`** for the candidate-t0-gated RecHit subsets
(e.g. `SiBarrelVertexTimeCoincRecHits`).

This plugin turns a continuous **time-frame** of streaming-readout data into
discrete **physics events** in software (there is no hardware trigger in the
ePIC streaming DAQ). Enable it by reading the input at the time-frame level:

```sh
eicrecon -Peventbuilder=true -Pplugins=eventbuilder <input_sim.edm4hep.root>
```

`-Peventbuilder=true` makes the PODIO source deliver each entry at
`JEventLevel::Timeslice` (the 2 µs time-frame) instead of `PhysicsEvent`. The
plugin is loaded by default but stays **idle** without this flag (its factories
sit at Timeslice level with no time-frame parent to act on), so default
reconstruction is unaffected.

## Pipeline

```
time-frame  (JANA Timeslice level)  ── -Peventbuilder=true ──┐
  └─ digitization (plugins/*.cc, registered @ Timeslice)
       → *RawHitDigi / *RecHitDigi / *ClusterDigi …
  └─ TrkTimeAlignment_factory  (tag: align)
       → *TimeAlignRecHits  [current frame only; slow (Si/B0) hits are also
                             deposited into TimesliceBuffer_service, keyed by frame number]
  └─ CalTimeAlignment_factory  → *TimeAlignClusters
  └─ EventBuilder_factory      (finds candidates, assigns t0 + t0sigma — TOF/MPGD only)
       → EventCandidates
  └─ TimeCoincidence_factory   (tag: coinc) — candidate-t0 gate on slow *TimeAlignRecHits
       → *TimeCoincRecHits   (subset of *TimeAlignRecHits for Si/B0)
  └─ EventPrefilter_factory    (optional ONNX frame-level background rejection)
       → EventCandidatesFiltered
  └─ EventUnfolder             (materialises one PhysicsEvent per candidate; merges
       → PhysicsEvent           adjacent-frame Si/B0 hits from TimesliceBuffer_service)
            └─ full reconstruction (ACTS, …) runs here
```

- **`EventBuilder_factory`** (`JOmniFactory` @ Timeslice): groups time-coincident
  hits from the fast trigger detectors (TOF + MPGD) with an O(N) sliding-window
  scan plus a 12×8 + 12×8 staggered θ/φ topology test. Computes t0 by
  inverse-variance weighting: `t0 = Σ(tᵢ/σᵢ²)/Σ(1/σᵢ²)`,
  `t0sigma = nsigma · 1/√(Σ 1/σᵢ²)`. TOF (~30 ps) dominates — **t0 is set by the
  fastest detector**. Analysis only, no reconstruction.
- **`EventUnfolder`** (`JEventUnfolder`): structural only. Reads the candidate
  list and emits one `PhysicsEvent` per candidate, selecting each detector's hits
  inside its resolution-matched window `[t0 ± nsigma·σ_det]` (tight for TOF, wide
  for Si) and writing the real t0/t0sigma into the child `EventHeader`.
- **`TrkTimeAlignment_factory`** (`TimeAlignment_factory<TrackerHit, 6>`): subtracts
  the r/c propagation term per hit and sorts by time. The output is always the
  current frame. The factory is **stateless** (JANA2 runs one instance per
  in-flight timeslice, so factory state cannot span frames); for **slow
  detectors** (indices 6-9: Si, B0) it additionally deposits the corrected hits
  into the process-wide `TimesliceBuffer_service`, keyed by absolute frame
  number, for cross-frame recovery by the EventUnfolder (see the Cross-frame
  section below).
- **`CalTimeAlignment_factory`** (`TimeAlignment_factory<Cluster>`): same r/c
  correction for calo clusters, no cross-frame deposit (all clusters are "fast").
- **`TimeCoincidence_factory`** (tag `coinc`): reads `EventCandidates` (t0 values)
  and the slow-detector `*TimeAlignRecHits` (current frame), and forwards only
  hits within `nsigma × 2000 ns` of some candidate's t0. The 2000 ns gate width
  is the MAPS integration window — a detector constant, not a free parameter.
  Output is a PODIO subset collection — no cloning. EventUnfolder reads
  `*TimeCoincRecHits` for slow detectors (indices 6-9).
  Key parameter:
  ```
  -Peventbuilder:coinc:nsigma=3.0   # N: gate = N × 2000 ns
  ```
- **`plugins/*.cc`**: re-register the standard ePIC digitization +
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

## ONNX AI integration

An optional GNN-based prefilter plugs into the pipeline via `ONNXRuntime_service`
(see `src/services/onnx/README.md` for service details and model library).

**Frame-level prefilter** (`EventPrefilter_factory`) — classifies each time-frame
as signal ("gold") or background-only before any PhysicsEvent is unfolded.
Uses a MultiClassEventGNN (4-input, probs `[1, 7]`); a frame passes if
`probs[BKG=0] < score_threshold`.

```sh
-Peventbuilder:onnx_model=/mnt/local/share/models/prefilter-gold-new-100epochs-2500hits.onnx
-Peventbuilder:score_threshold=0.5   # default: BKG probability cut
-Peventbuilder:max_hits=2500          # match model's fixed input size
-Peventbuilder:k_neighbors=16         # kNN graph k per DGCNN layer
-Peventbuilder:time_weight=1.0        # scale on time axis in layer-0 spacetime kNN
```

Hit features per frame (6 per hit, normalised):
`[x/4000, y/4000, z/5000, t/50, clip(eDep/10, 0, 1), det_id/6]`
where `det_id` is 2 for TOF (barrel + endcap) and 1 for MPGD (barrel + endcaps).

`EventBuilder_factory` always uses algorithmic inverse-variance t0:
`t0 = Σ(tᵢ/σᵢ²) / Σ(1/σᵢ²)`, `t0sigma = nsigma / √(Σ 1/σᵢ²)`.

## Configuration reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eventbuilder:backward_frames` | `1` | Past frames whose Si/B0 hits the EventUnfolder merges per candidate; `0` = off |
| `eventbuilder:forward_frames` | `0` | Future frames merged per candidate; needs `nthreads ≥ 2` and `jana:max_inflight_timeslices > forward_frames`; `0` = off |
| `eventbuilder:coinc:nsigma` | `3.0` | Gate half-width in units of σ_Si (2000 ns) for slow hit selection |
| `eventbuilder:eventBuilder:timeResolution_TOF` | `0.03` | TOF σ in ns |
| `eventbuilder:eventBuilder:timeResolution_MPGD` | `10.0` | MPGD σ in ns |
| `eventbuilder:eventBuilder:nsigma_window` | `3.0` | N-sigma half-width defining t0sigma |
| `eventbuilder:eventBuilder:coincidence_window` | `50.0` | Sliding-window width in ns for fast-hit grouping |
| `eventbuilder:onnx_model` | `""` | Frame-level GNN prefilter model path; empty = pass-through |
| `eventbuilder:score_threshold` | `0.5` | Maximum BKG class probability (`probs[0]`) to pass a frame |
| `eventbuilder:max_hits` | `0` | Fixed hit count (pad/truncate to model input size); `0` = dynamic |
| `eventbuilder:k_neighbors` | `16` | k for kNN graph construction (per DGCNN layer) |
| `eventbuilder:time_weight` | `1.0` | Scale factor on the time dimension in layer-0 spacetime kNN |

## Notes

- **Calorimeter carry-through** is enabled for B0ECAL, BEMC (EcalBarrel),
  EEMC (EcalEndcapN), and FEMC (EcalEndcapP) via `plugins/FEMC.cc` (homogeneous
  SiPM-on-tile parameters). ScFi geometry and `TrackClusterMergeSplitter` are not
  registered at Timeslice: the splitter requires ACTS `CalorimeterTrackProjections`,
  which are only available at PhysicsEvent level.
- **Cross-frame Si/B0 RecHit recovery** is on by default (backward direction;
  see the Cross-frame section below). Disable with `backward_frames=0`.
  Enable forward recovery with `forward_frames=N` (no output latency; requires
  `nthreads ≥ 2` and `jana:max_inflight_timeslices > N`).
- **ONNX prefilter** is off by default (empty `eventbuilder:onnx_model`). See
  the ONNX section above and `src/services/onnx/README.md` for usage.

---

## Cross-frame Si/B0 hit recovery

The Si vertex/tracker detectors have σ_Si ≈ 2 µs ≈ frame width, so hits from
tracks near a frame boundary can be reconstructed into a neighbouring frame.
Since each JANA `JEvent` has a single time-frame origin, such hits are invisible
to the EventUnfolder unless the pipeline explicitly merges them.

The merge is built for multithreaded JANA2: factories run in the parallel map
stage with one instance per in-flight timeslice, so no factory carries state
across frames. Instead, `TrkTimeAlignment_factory` **deposits** each frame's
corrected slow-detector hits into the `TimesliceBuffer_service`
(`src/services/eventbuilder/`), keyed by absolute frame number. The
EventUnfolder — the sequential stage of the topology — **fetches** the adjacent
frames `[N-backward_frames … N+forward_frames]` from the service and gates
their hits per candidate, exactly like current-frame hits. Correctness is
independent of thread count and of the order in which frames are processed.

**Backward recovery** (`backward_frames`, default 1): frame N-1 is always
deposited before frame N reaches the unfolder (at worst the fetch briefly
waits for an in-flight neighbour). Hits smeared *backward* (from frame N into
frame N-1) are recovered with no latency and no extra configuration.

**Forward recovery** (`forward_frames`, default 0, optional): the fetch of
frame N+1 waits until that frame's parallel stage has run. This adds **no
output latency** — it only stalls the unfolder until an already-in-flight
frame is deposited — but it requires enough concurrency for that frame to
be in flight: `nthreads ≥ 2` and `jana:max_inflight_timeslices >
forward_frames` (enforced at startup). At end of stream the final frames have
no future neighbour; their candidates are still built (after a bounded wait),
only without forward hits — no events are lost.

With σ_Si = 2 µs and a 2 µs frame width, roughly 16 % of Si hits straddle a
boundary — ~8 % smeared backward, ~8 % forward. Default (`backward_frames=1`)
recovers the backward half. Full bidirectional coverage with ±1 frame:

```
-Peventbuilder:backward_frames=1 -Peventbuilder:forward_frames=1
```

---

> Note: JANA keys factories by `(type, tag)` only — *not* by event level — so a
> Timeslice and a PhysicsEvent collection of the same name would collide. That is
> why the Timeslice collections must carry the `Digi` / `TimeAlign` markers rather
> than reuse the standard names.

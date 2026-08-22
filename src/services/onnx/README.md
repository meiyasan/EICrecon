# `ONNXRuntime_service`

A singleton `JService` that owns one `Ort::Env` per process and caches
`Ort::Session` handles keyed by model path.

## Why a shared service?

ONNX Runtime's documentation recommends **one `Ort::Env` per process**.
Creating multiple `Ort::Env` instances each spin up their own internal thread
pool and logging backend, wasting resources and occasionally causing conflicts.

The existing `ONNXInference` algorithm (in `src/algorithms/onnx/`) creates its
own `Ort::Env` per factory instance — correct for isolated use cases, but not
ideal when multiple factories load models during the same run. This service
centralises the environment and session cache so any JANA factory can share
sessions without duplicating initialisation overhead.

`Ort::Session::Run()` is **thread-safe** (const under the hood), so all JANA
worker threads may share a single session handle for the same model.

## Usage in a factory

```cpp
#include "services/onnx/ONNXRuntime_service.h"

struct MyFactory : public JOmniFactory<MyFactory> {
  Parameter<std::string> m_model_path{this, "onnx_model", "", "..."};
  std::shared_ptr<Ort::Session> m_session;
  std::string m_input_name, m_output_name;

  void Configure() {
    const std::string path = m_model_path();
    if (path.empty()) return; // optional — no model, no-op

    auto* svc  = GetApplication()->GetService<ONNXRuntime_service>();
    m_session  = svc->session(path);

    Ort::AllocatorWithDefaultOptions alloc;
    m_input_name  = m_session->GetInputNameAllocated(0, alloc).get();
    m_output_name = m_session->GetOutputNameAllocated(0, alloc).get();
  }

  void Process(int64_t, uint64_t) {
    if (!m_session) { /* no model: pass-through or default logic */ return; }
    // Build Ort::Value input tensor, call m_session->Run(...)
  }
};
```

## Session options

Sessions are created with:
- `SetIntraOpNumThreads(1)` — single-threaded within one inference call
- `SetInterOpNumThreads(1)` — no parallel graph execution

This is intentional: **JANA owns cross-event parallelism**. ORT's own threading
would compete with JANA's worker pool for CPU resources. Each JANA worker thread
calls `Session::Run()` concurrently on the same session; the session serialises
GPU memory transfers internally if needed.

## Model library

Models live under `share/models/` in the eic-shell repo, which mounts to
`/mnt/local/share/models/` inside the container.

| File | Task | Inputs | Output |
|------|------|--------|--------|
| `prefilter-gold-new-100epochs-2500hits.onnx` | Frame-level GNN background rejection | 4 inputs: `x [N,6]` float32 hit features + 3× edge index `[2,E]` int64 | `probs [1,7]` float32 softmax; `probs[0]` = BKG probability |

### Hit feature convention (GNN prefilter)

6 normalised features per fast-detector hit (TOF + MPGD, up to `max_hits`):

| Column | Field | Normalisation |
|--------|-------|---------------|
| 0 | x position | `/4000` (mm) |
| 1 | y position | `/4000` (mm) |
| 2 | z position | `/5000` (mm) |
| 3 | time | `/50` (ns) |
| 4 | energy deposit | `clip(eDep/10, 0, 1)` |
| 5 | detector id | `2/6` for TOF (barrel/endcap), `1/6` for MPGD |

## AI examples

### Frame-level GNN prefilter (background rejection before unfolding)

`EventPrefilter_factory` classifies each 2 µs time-frame using a
MultiClassEventGNN before any PhysicsEvent is unfolded. The model scores
7 classes; frames where `probs[BKG=0] < score_threshold` pass through.
Frames with a BKG score above the threshold produce no PhysicsEvents —
zero reconstruction overhead for rejected frames.

```sh
eicrecon \
  -Peventbuilder:onnx_model=/mnt/local/share/models/prefilter-gold-new-100epochs-2500hits.onnx \
  -Peventbuilder:score_threshold=0.5 \
  -Peventbuilder:max_hits=2500 \
  input_sim.edm4hep.root
```

- `onnx_model` — path to the GNN `.onnx` file; empty (default) = pass-through
- `score_threshold` — frame passes when `probs[BKG=0] < threshold`; default 0.5
- `max_hits` — truncate/pad hit count to this size (must match the model's fixed
  input; use `0` for dynamic-shape models)
- `k_neighbors` — k for kNN graph construction per DGCNN layer; default 16
- `time_weight` — scale on the time axis in layer-0 spacetime kNN; default 1.0

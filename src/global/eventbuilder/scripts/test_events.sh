#!/usr/bin/env bash
# test_events.sh — run eicrecon+eventbuilder on 1 edm4hep event per SRO class
# and write the outputs to the samples directory.
#
# Usage (called by `make test-events`):
#   test_events.sh <samples_dir> <edm4hep_root> [<eicrecon_bin>]
#
# Arguments
#   samples_dir   — output directory, e.g. .../share/samples/gold-new
#   edm4hep_root  — root of edm4hep datasets, e.g. .../vendor/sro/datasets/edm4hep
#   eicrecon_bin  — path to the eicrecon executable (default: find in PATH)
#
# Environment
#   DETECTOR_PATH — ePIC compact geometry directory (required by eicrecon)
#   DETECTOR_CONFIG — geometry config name (default: epic_craterlake_10x100)
#
# Output per class:
#   <samples_dir>/bkg.eicrecon.root
#   <samples_dir>/dis_nc.eicrecon.root
#   ...  (7 files total)
#
# Search order for input edm4hep files:
#   vanilla → gold-new → gold-old  (vanilla = pure signal, no background mixing)

set -euo pipefail

SAMPLES="${1:?Usage: test_events.sh <samples_dir> <edm4hep_root> [<eicrecon_bin>]}"
EDM_ROOT="${2:?Usage: test_events.sh <samples_dir> <edm4hep_root> [<eicrecon_bin>]}"
EIC="${3:-eicrecon}"

# Verify eicrecon is usable
if ! command -v "$EIC" &>/dev/null && [ ! -x "$EIC" ]; then
    echo "ERROR: eicrecon not found at '$EIC'" >&2
    echo "  Source the eic-shell setup or pass the path as the third argument." >&2
    exit 1
fi

# Compact geometry
CONFIG="${DETECTOR_CONFIG:-epic_craterlake_10x100}"
if [ -n "${DETECTOR_PATH:-}" ]; then
    COMPACT="$DETECTOR_PATH/${CONFIG}.xml"
else
    # try common container locations
    for d in /mnt/local/share/epic /opt/detector; do
        if [ -f "$d/${CONFIG}.xml" ]; then COMPACT="$d/${CONFIG}.xml"; break; fi
    done
fi
if [ -z "${COMPACT:-}" ] || [ ! -f "$COMPACT" ]; then
    echo "WARNING: compact geometry not found (DETECTOR_PATH not set)" >&2
    echo "  Attempting eicrecon without -Pdd4hep:xml_files — may fail." >&2
    COMPACT_FLAG=
else
    COMPACT_FLAG="-Pdd4hep:xml_files=${COMPACT}"
fi

CLASSES=(bkg dis_nc dis_cc ddis exclusive_h exclusive_d sidis)
LABELS=(vanilla gold-new gold-old)

mkdir -p "$SAMPLES"

echo "=== test-events: eventbuilder sample generation ==="
echo "    output : $SAMPLES"
echo "    edm4hep: $EDM_ROOT"
echo "    eicrecon: $(command -v "$EIC" 2>/dev/null || echo "$EIC")"
echo ""

any_run=0
for cls in "${CLASSES[@]}"; do
    out="$SAMPLES/$cls.eicrecon.root"
    if [ -f "$out" ]; then
        echo "  ↷  $cls (already exists — delete to regenerate)"
        continue
    fi

    # Search for any edm4hep sim file for this class
    infile=
    for label in "${LABELS[@]}"; do
        dir="$EDM_ROOT/$label/$cls"
        [ -d "$dir" ] || continue
        f=$(ls "$dir/"*_sim.edm4hep.root 2>/dev/null | head -1 || true)
        if [ -n "$f" ]; then infile="$f"; break; fi
    done

    if [ -z "$infile" ]; then
        echo "  ✗  $cls — no edm4hep input found under $EDM_ROOT/{${LABELS[*]}}/$cls/"
        continue
    fi

    echo "  ⚙   $cls ← $(basename "$infile")"
    "$EIC" \
        -Pjana:nevents=1 \
        -Pplugins=eventbuilder \
        -Peventbuilder=true \
        ${COMPACT_FLAG:+"$COMPACT_FLAG"} \
        -Ppodio:output_file="$out" \
        "$infile"
    echo "     ✓  → $(basename "$out")"
    any_run=1
done

echo ""
if [ "$any_run" -eq 0 ]; then
    echo "All samples already exist. Done."
else
    echo "Done. Samples written to: $SAMPLES"
fi

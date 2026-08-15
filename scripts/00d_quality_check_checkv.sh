#!/usr/bin/env bash
set -euo pipefail

# Assess viral genome quality using CheckV.
#
# Historical Yosemite workflow:
#   - CheckV end_to_end
#   - 16 threads
#   - CheckV database v1.5
#
# Software:
#   CheckV v1.0.3
#
# USER CONFIGURATION
# Set CHECKVDB to the local path of the CheckV database.
# Supply a sample identifier and VirSorter2 final-viral-combined FASTA.

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <sample_id> <final-viral-combined.fa>"
    exit 1
fi

SAMPLE="$1"
VIRAL_FASTA="$2"

THREADS=16

# Replace with the location of the CheckV database on your system.
export CHECKVDB="/path/to/checkv-db-v1.5"

OUTDIR="results/checkv/${SAMPLE}"

mkdir -p "results/checkv"

echo "[INFO] Running CheckV for ${SAMPLE}"

rm -rf "${OUTDIR}"

checkv end_to_end \
    "${VIRAL_FASTA}" \
    "${OUTDIR}" \
    -t "${THREADS}"

echo "[DONE] CheckV complete for ${SAMPLE}"
echo "[INFO] Output directory: ${OUTDIR}"

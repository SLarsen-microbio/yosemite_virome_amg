#!/usr/bin/env bash
set -euo pipefail

# Identify viral sequences in assembled metagenomic contigs using VirSorter2.
#
# Historical Yosemite workflow:
#   - input contigs were pre-filtered to >=1 kb
#   - VirSorter2 minimum viral sequence length: 5000 bp
#   - 16 threads
#   - all viral groups enabled
#
# Software:
#   VirSorter2 v2.2.4
#
# USER CONFIGURATION
# Supply a sample identifier and the FASTA file containing assembled contigs.

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <sample_id> <contigs.fasta>"
    exit 1
fi

SAMPLE="$1"
CONTIGS="$2"

THREADS=16
OUTDIR="results/virsorter2/${SAMPLE}"

mkdir -p "results/virsorter2"

echo "[INFO] Running VirSorter2 for ${SAMPLE}"

virsorter run \
    -w "${OUTDIR}" \
    -i "${CONTIGS}" \
    --min-length 5000 \
    -j "${THREADS}" \
    all

echo "[DONE] VirSorter2 complete for ${SAMPLE}"
echo "[INFO] Output directory: ${OUTDIR}"

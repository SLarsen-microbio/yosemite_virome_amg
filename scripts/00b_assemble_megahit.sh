#!/usr/bin/env bash
set -euo pipefail

# Assemble paired-end metagenomic reads with MEGAHIT.
#
# Historical Yosemite workflow:
#   - paired-end assembly
#   - 16 threads
#   - otherwise MEGAHIT default parameters
#
# Software:
#   MEGAHIT v1.1.5
#
# USER CONFIGURATION
# Replace READ1 and READ2 with the trimmed paired-end FASTQ files
# for the sample being assembled.

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <sample_id> <read1.fastq.gz> <read2.fastq.gz>"
    exit 1
fi

SAMPLE="$1"
READ1="$2"
READ2="$3"

THREADS=16
OUTDIR="results/megahit/${SAMPLE}"

mkdir -p "results/megahit"

echo "[INFO] Assembling ${SAMPLE}"

megahit \
    -1 "${READ1}" \
    -2 "${READ2}" \
    -o "${OUTDIR}" \
    -t "${THREADS}"

echo "[DONE] Assembly complete for ${SAMPLE}"
echo "[INFO] Final contigs: ${OUTDIR}/final.contigs.fa"

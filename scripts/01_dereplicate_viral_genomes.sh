#!/usr/bin/env bash
set -euo pipefail

# Dereplicate HQ/MQ viral genomes using CD-HIT-EST.
#
# Viral genomes are clustered at 95% nucleotide identity with
# alignment covering at least 85% of the shorter sequence.
#
# Original analysis:
#   Input:  909 HQ/MQ viral sequences
#   Output: 695 representative viral genome clusters
#
# Software: CD-HIT-EST v4.8.1
# USER CONFIGURATION
# Replace with the path to the FASTA containing the 909 HQ/MQ viral genomes.

INPUT="/path/to/all909_hqmc.fna"

INPUT="/path/to/all909_hqmc.fna"
OUTPUT="results/derep95_all17"

mkdir -p results

cd-hit-est \
    -i "${INPUT}" \
    -o "${OUTPUT}" \
    -c 0.95 \
    -aS 0.85 \
    -n 8 \
    -d 0 \
    -T 16 \
    -M 32000

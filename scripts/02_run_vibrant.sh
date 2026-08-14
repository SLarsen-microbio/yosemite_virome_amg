#!/usr/bin/env bash
set -euo pipefail

# Identify auxiliary metabolic genes (AMGs) in dereplicated viral genomes
# using VIBRANT.
#
# Input:
#   Viral genome representatives dereplicated at 95% nucleotide identity
#   and >=85% alignment coverage of the shorter sequence.
#
# Original analysis:
#   Input: 695 dereplicated viral genome representatives
#   VIBRANT identified 710 putative phages.
#
# AMG annotations used in downstream analyses were obtained from:
#   VIBRANT_AMG_individuals_derep95_all17.tsv
#   VIBRANT_AMG_pathways_derep95_all17.tsv
#
# Software: VIBRANT v1.2.1

INPUT="results/derep95_all17"
OUTPUT_DIR="results/vibrant_amg"

mkdir -p "${OUTPUT_DIR}"

VIBRANT_run.py \
    -i "${INPUT}" \
    -folder "${OUTPUT_DIR}" \
    -t 16 \
    -virome

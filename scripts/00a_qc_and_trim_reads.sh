#!/usr/bin/env bash
set -euo pipefail

# Quality assessment and trimming of metagenomic sequencing reads.
#
# Historical workflow used for the Yosemite lake metagenomic datasets:
#
#   1. FastQC on raw reads
#   2. fastp trimming with automatic paired-end adapter detection
#   3. FastQC on trimmed reads
#   4. MultiQC aggregation of QC reports
#
# Software:
#   FastQC v0.11.8
#   fastp v1.0.1
#   MultiQC v1.30
#
# fastp parameter retained from the original workflow:
#   --detect_adapter_for_pe
#
# USER CONFIGURATION
# Replace RAW_DIR with the location of the downloaded raw FASTQ files.

THREADS=8

RAW_DIR="/path/to/raw_fastq_files"

TRIM_DIR="results/trimmed"
QC_DIR="results/qc"
QC_RAW_DIR="${QC_DIR}/raw_fastqc"
QC_TRIM_DIR="${QC_DIR}/trimmed_fastqc"

mkdir -p \
    "${TRIM_DIR}" \
    "${QC_DIR}" \
    "${QC_RAW_DIR}" \
    "${QC_TRIM_DIR}"

# ---------------------------------------------------------------------
# 1. FastQC on raw reads
# ---------------------------------------------------------------------

echo "[STEP] Running FastQC on raw reads"

find "${RAW_DIR}" \
    -type f \
    -name "*.fastq.gz" \
    -print0 \
    | xargs -0 -n1 -P4 fastqc \
        -o "${QC_RAW_DIR}"

# ---------------------------------------------------------------------
# 2. Trim paired-end reads with fastp
# ---------------------------------------------------------------------

echo "[STEP] Trimming paired-end reads with fastp"

for R1 in "${RAW_DIR}"/*_1.fastq.gz
do
    [[ -e "${R1}" ]] || continue

    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    if [[ ! -f "${R2}" ]]; then
        echo "[WARN] Matching R2 file not found for ${R1}"
        continue
    fi

    SAMPLE="$(basename "${R1}" _1.fastq.gz)"

    OUT1="${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
    OUT2="${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"

    echo "[INFO] Trimming ${SAMPLE}"

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT1}" \
        -O "${OUT2}" \
        --detect_adapter_for_pe \
        --thread "${THREADS}" \
        -j "${QC_DIR}/${SAMPLE}.json" \
        -h "${QC_DIR}/${SAMPLE}.html"
done

# ---------------------------------------------------------------------
# 3. FastQC on trimmed reads
# ---------------------------------------------------------------------

echo "[STEP] Running FastQC on trimmed reads"

find "${TRIM_DIR}" \
    -type f \
    -name "*.trimmed.fastq.gz" \
    -print0 \
    | xargs -0 -n1 -P4 fastqc \
        -o "${QC_TRIM_DIR}"

# ---------------------------------------------------------------------
# 4. Aggregate reports with MultiQC
# ---------------------------------------------------------------------

echo "[STEP] Running MultiQC"

multiqc "${QC_DIR}" \
    -o "${QC_DIR}"

echo "[DONE] QC and trimming complete"

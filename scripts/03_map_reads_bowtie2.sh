#!/usr/bin/env bash
set -euo pipefail

# Map paired-end metagenomic reads to the dereplicated viral genome
# representatives using Bowtie2.
#
# This script maps one sample at a time. Run it once for each metagenomic
# sample listed in metadata/sample_metadata.tsv.
#
# Usage:
#   bash scripts/03_map_reads_bowtie2.sh \
#       <sample_id> \
#       <read1.fastq.gz> \
#       <read2.fastq.gz> \
#       <dereplicated_viral_genomes.fna>
#
# Example:
#   bash scripts/03_map_reads_bowtie2.sh \
#       SRR5830790 \
#       /path/to/SRR5830790_R1.fastq.gz \
#       /path/to/SRR5830790_R2.fastq.gz \
#       /path/to/dereplicated_viral_genomes.fna
#
# Inputs:
#   1. Sample identifier
#   2. Trimmed forward metagenomic reads
#   3. Trimmed reverse metagenomic reads
#   4. FASTA containing the 695 dereplicated viral genome representatives
#
# Outputs:
#   - Sorted BAM file
#   - BAM index
#   - Bowtie2 mapping log
#
# Software:
#   Bowtie2 v2.5.1
#   SAMtools

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <sample_id> <read1.fastq.gz> <read2.fastq.gz> <dereplicated_viral_genomes.fna>"
    exit 1
fi

SAMPLE="$1"
READ1="$2"
READ2="$3"
GENOMES="$4"

OUTDIR="results/read_recruitment/bam"
DB_PREFIX="results/read_recruitment/viral_db"

mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${DB_PREFIX}")"

# Build Bowtie2 index if it does not already exist.
if [[ ! -f "${DB_PREFIX}.1.bt2" && ! -f "${DB_PREFIX}.1.bt2l" ]]; then
    echo "Building Bowtie2 index"
    bowtie2-build "${GENOMES}" "${DB_PREFIX}"
fi

echo "Mapping ${SAMPLE}"

bowtie2 \
    -x "${DB_PREFIX}" \
    -1 "${READ1}" \
    -2 "${READ2}" \
    --no-unal \
    -p 8 \
    2> "${OUTDIR}/${SAMPLE}.log" \
    | samtools sort \
        -@ 8 \
        -o "${OUTDIR}/${SAMPLE}.bam"

samtools index "${OUTDIR}/${SAMPLE}.bam"

echo "Done ${SAMPLE}"

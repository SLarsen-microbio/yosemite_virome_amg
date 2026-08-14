#!/usr/bin/env bash
set -euo pipefail

# Map paired-end metagenomic reads to dereplicated viral genomes using Bowtie2.
#
# Input:
#   - Dereplicated viral genome representatives
#   - Paired-end trimmed metagenomic reads listed in metadata/sample_metadata.tsv
#
# Output:
#   - Sorted BAM file for each sample
#   - BAM index for each sample
#   - Bowtie2 mapping log for each sample
#
# Software:
#   Bowtie2 v2.5.1
#   SAMtools

GENOMES="results/derep95_all17"
DB_PREFIX="results/read_recruitment/viral_db"
METADATA="metadata/sample_metadata.tsv"
OUTDIR="results/read_recruitment/bam"

mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${DB_PREFIX}")"

# Build Bowtie2 index
bowtie2-build "${GENOMES}" "${DB_PREFIX}"

# sample_metadata.tsv must contain:
# sample    read1    read2
#
# Example:
# SRR5830790    data/SRR5830790_1.trimmed.fastq.gz    data/SRR5830790_2.trimmed.fastq.gz

while IFS=$'\t' read -r sample read1 read2; do

    # Skip header
    if [[ "${sample}" == "sample" ]]; then
        continue
    fi

    echo "Mapping ${sample}"

    bowtie2 \
        -x "${DB_PREFIX}" \
        -1 "${read1}" \
        -2 "${read2}" \
        --no-unal \
        -p 8 \
        2> "${OUTDIR}/${sample}.log" \
        | samtools sort \
            -@ 8 \
            -o "${OUTDIR}/${sample}.bam"

    samtools index "${OUTDIR}/${sample}.bam"

    echo "Done ${sample}"

done < "${METADATA}"

echo "ALL DONE"

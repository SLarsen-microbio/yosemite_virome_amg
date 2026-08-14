#!/usr/bin/env bash
set -euo pipefail

# Figure 3. Deazaguanine precursor biosynthesis gene cluster
#
# Panels A and B:
# Circular genome maps of viral contigs k141_70163 and k141_11273
# from the Upper Cathedral Lake UC1 warming tFinal sample.
#
# The two genomes were extracted from the dereplicated viral genome
# collection, annotated with Pharokka v1.9.1, and visualized using
# Pharokka's multiplotter.
#
# Panel C was assembled separately from the Pharokka annotations and
# shows the local organization of the queC, queE, queD, and folE
# deazaguanine precursor biosynthesis genes.
#
# Author: Shari Larsen


# ---------------------------------------------------------------------
# Input/output paths
# ---------------------------------------------------------------------

DEREP_FASTA="results/derep95_all17"

CANDIDATES="results/figure3/queuosine_candidates.fasta"

PHAROKKA_OUT="results/figure3/pharokka_queuosine_check"

PLOT_OUT="results/figure3/pharokka_queuosine_plots_final"

mkdir -p "results/figure3"


# ---------------------------------------------------------------------
# Extract the two viral genomes used in Figure 3
# ---------------------------------------------------------------------

# Genome identifiers:
#
# Y86_UC1_W_tf__k141_70163||full
# Y86_UC1_W_tf__k141_11273||full
#
# seqkit is used here to retrieve the two sequences from the
# dereplicated viral genome FASTA.

seqkit grep \
    -p 'Y86_UC1_W_tf__k141_70163||full' \
    -p 'Y86_UC1_W_tf__k141_11273||full' \
    "${DEREP_FASTA}" \
    > "${CANDIDATES}"


# ---------------------------------------------------------------------
# Annotate genomes with Pharokka
# ---------------------------------------------------------------------

pharokka.py \
    -i "${CANDIDATES}" \
    -o "${PHAROKKA_OUT}" \
    -t 8 \
    --force


# ---------------------------------------------------------------------
# Generate circular genome maps for Figure 3A and 3B
# ---------------------------------------------------------------------

pharokka_multiplotter.py \
    -g "${PHAROKKA_OUT}/pharokka.gbk" \
    -o "${PLOT_OUT}" \
    --dpi 300 \
    -t "Deazaguanine_pathway_cluster"


echo
echo "Figure 3A-B genome maps generated."
echo "Pharokka annotations: ${PHAROKKA_OUT}"
echo "Circular genome maps: ${PLOT_OUT}"

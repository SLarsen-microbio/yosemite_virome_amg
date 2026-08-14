#!/usr/bin/env python3

"""
Figure 3C. Deazaguanine precursor biosynthesis gene neighborhoods.

Creates a linear gene-neighborhood plot for the two viral genomes
shown in Figure 3:

    Y86_UC1_W_tf__k141_70163||full
    Y86_UC1_W_tf__k141_11273||full

Pharokka annotations are used to identify the deazaguanine-associated
genes QueC, QueE, QueD, and GTP cyclohydrolase (FolE).

Target genes are highlighted in red and surrounding genes are shown
in gray.

Author: Shari Larsen
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

ANNOTATION_PATH = Path(
    "results/figure3/pharokka_queuosine_check/"
    "pharokka_cds_final_merged_output.tsv"
)

OUTDIR = Path("results/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUTDIR / "figure3_panelC_gene_cluster.png"
OUT_SVG = OUTDIR / "figure3_panelC_gene_cluster.svg"


# ---------------------------------------------------------------------
# Genome names
# ---------------------------------------------------------------------

GENOMES = [
    "Y86_UC1_W_tf__k141_70163||full",
    "Y86_UC1_W_tf__k141_11273||full",
]


# ---------------------------------------------------------------------
# Load Pharokka CDS annotations
# ---------------------------------------------------------------------

df = pd.read_csv(
    ANNOTATION_PATH,
    sep="\t",
)


# Pharokka output columns used here:
#
# start              = CDS start coordinate
# stop               = CDS end coordinate
# strand             = + or -
# contig              = sequence identifier
# product             = predicted gene/product annotation
#
# Column names can differ slightly among Pharokka versions, so identify
# them from the known structure of the output if needed.

print("Columns in Pharokka annotation table:")
print(list(df.columns))


# ---------------------------------------------------------------------
# Identify required columns
# ---------------------------------------------------------------------

def find_column(possible_names):
    """Return first matching column name."""
    for name in possible_names:
        if name in df.columns:
            return name

    raise ValueError(
        "Could not identify required column. "
        f"Tried: {possible_names}"
    )


START_COL = find_column(
    ["start", "Start", "gene_start"]
)

END_COL = find_column(
    ["stop", "end", "End", "gene_end"]
)

STRAND_COL = find_column(
    ["strand", "Strand"]
)

CONTIG_COL = find_column(
    ["contig", "sequence", "seqid", "contig_id"]
)

PRODUCT_COL = find_column(
    [
        "product",
        "Product",
        "function",
        "annotation",
        "gene",
    ]
)


# ---------------------------------------------------------------------
# Classify target genes
# ---------------------------------------------------------------------

def target_gene(product):
    """
    Return manuscript gene label for deazaguanine pathway genes.
    """

    text = str(product).lower()

    if "quec" in text:
        return "queC"

    if "quee" in text:
        return "queE"

    if "qued" in text:
        return "queD"

    if "gtp cyclohydrolase" in text:
        return "folE"

    return None


df["target_gene"] = (
    df[PRODUCT_COL]
    .apply(target_gene)
)


# ---------------------------------------------------------------------
# Extract neighborhoods
# ---------------------------------------------------------------------

neighborhoods = {}

for genome in GENOMES:

    genome_df = df[
        df[CONTIG_COL] == genome
    ].copy()

    if genome_df.empty:
        raise ValueError(
            f"No Pharokka annotations found for {genome}"
        )

    targets = genome_df[
        genome_df["target_gene"].notna()
    ]

    print(f"\n{genome}")
    print(
        targets[
            [
                START_COL,
                END_COL,
                STRAND_COL,
                PRODUCT_COL,
                "target_gene",
            ]
        ].to_string(index=False)
    )

    if targets.empty:
        raise ValueError(
            f"No queC/queE/queD/folE genes found for {genome}"
        )

    # Display the full region spanning the target genes plus 2 kb
    # of flanking sequence on each side.

    region_start = max(
        0,
        min(
            targets[START_COL].min(),
            targets[END_COL].min(),
        )
        - 2000,
    )

    region_end = (
        max(
            targets[START_COL].max(),
            targets[END_COL].max(),
        )
        + 2000
    )

    region = genome_df[
        (
            genome_df[[START_COL, END_COL]]
            .max(axis=1)
            >= region_start
        )
        &
        (
            genome_df[[START_COL, END_COL]]
            .min(axis=1)
            <= region_end
        )
    ].copy()

    neighborhoods[genome] = (
        region,
        region_start,
        region_end,
    )


# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------

fig, axes = plt.subplots(
    2,
    1,
    figsize=(13, 6),
)

genome_labels = [
    "k141_70163",
    "k141_11273",
]


for ax, genome, genome_label in zip(
    axes,
    GENOMES,
    genome_labels,
):

    region, region_start, region_end = (
        neighborhoods[genome]
    )

    # baseline
    ax.plot(
        [region_start, region_end],
        [0, 0],
        linewidth=1,
    )

    for _, row in region.iterrows():

        start = int(row[START_COL])
        end = int(row[END_COL])
        strand = str(row[STRAND_COL])

        left = min(start, end)
        right = max(start, end)

        width = right - left

        label = row["target_gene"]

        # AMG pathway genes red; flanking genes gray
        if label is not None:
            color = "red"
        else:
            color = "lightgray"

        if strand == "+":
            arrow_start = left
            arrow_length = width
        else:
            arrow_start = right
            arrow_length = -width

        # Arrowhead size scaled conservatively so small genes remain visible
        head_length = min(
            250,
            max(60, width * 0.25),
        )

        ax.add_patch(
            FancyArrow(
                arrow_start,
                -0.13,
                arrow_length,
                0,
                width=0.26,
                head_width=0.42,
                head_length=head_length,
                length_includes_head=True,
                facecolor=color,
                edgecolor="black",
                linewidth=0.6,
            )
        )

        if label is not None:

            midpoint = (
                left + right
            ) / 2

            ax.text(
                midpoint,
                0.38,
                label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontstyle="italic",
                rotation=45,
            )

    ax.set_xlim(
        region_start,
        region_end,
    )

    ax.set_ylim(
        -0.7,
        1.0,
    )

    ax.set_yticks([])

    ax.set_ylabel(
        genome_label,
        rotation=0,
        ha="right",
        va="center",
        labelpad=50,
        fontsize=10,
    )

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.set_xlabel(
        "Genome position (bp)"
    )


fig.suptitle(
    "C. Deazaguanine precursor biosynthesis gene neighborhoods",
    fontsize=12,
)

plt.tight_layout()


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------

plt.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight",
)

plt.savefig(
    OUT_SVG,
    bbox_inches="tight",
)

plt.close()


print(f"\nSaved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")

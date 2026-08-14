#!/usr/bin/env python3

"""
Figure 2. Metabolic categories of auxiliary metabolic genes (AMGs)
identified in Yosemite alpine lake viral genomes.

AMG instances are taken from the VIBRANT AMG individual-gene output
and grouped by KEGG metabolic category using the VIBRANT AMG pathway
table.

Expected total:
    789 AMG instances across 695 dereplicated viral genomes.

Author: Shari Larsen
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

AMG_IND_PATH = Path(
    "results/vibrant_amg/"
    "VIBRANT_derep95_all17/"
    "VIBRANT_results_derep95_all17/"
    "VIBRANT_AMG_individuals_derep95_all17.tsv"
)

AMG_PATH_PATH = Path(
    "results/vibrant_amg/"
    "VIBRANT_derep95_all17/"
    "VIBRANT_results_derep95_all17/"
    "VIBRANT_AMG_pathways_derep95_all17.tsv"
)

OUTDIR = Path("results/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUTDIR / "figure2_amg_categories.png"
OUT_SVG = OUTDIR / "figure2_amg_categories.svg"
OUT_TSV = OUTDIR / "figure2_amg_category_counts.tsv"


# ---------------------------------------------------------------------
# Load VIBRANT AMG outputs
# ---------------------------------------------------------------------

amg_ind = pd.read_csv(
    AMG_IND_PATH,
    sep="\t",
)

amg_path = pd.read_csv(
    AMG_PATH_PATH,
    sep="\t",
)


# ---------------------------------------------------------------------
# Build KO-to-category lookup
# ---------------------------------------------------------------------

ko_to_category = {}

for _, row in amg_path.iterrows():

    category = row["Metabolism"]

    kos = str(
        row["Present AMG KOs"]
    ).split(",")

    for ko in kos:

        ko = ko.strip()

        if ko and ko != "nan":
            ko_to_category[ko] = category


# ---------------------------------------------------------------------
# Assign each AMG instance to a KEGG metabolic category
# ---------------------------------------------------------------------

amg_ind["category"] = (
    amg_ind["AMG KO"]
    .map(ko_to_category)
    .fillna("Other")
)


# ---------------------------------------------------------------------
# Count AMG instances
# ---------------------------------------------------------------------

category_counts = (
    amg_ind["category"]
    .value_counts()
    .rename_axis("metabolic_category")
    .reset_index(name="amg_count")
)

category_counts = (
    category_counts
    .sort_values(
        "amg_count",
        ascending=True,
    )
)


# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------

total_amgs = int(
    category_counts["amg_count"].sum()
)

print(
    f"Total AMG instances: {total_amgs}"
)

print(
    f"Number of metabolic categories: "
    f"{len(category_counts)}"
)

print("\nAMG counts by category:")
print(
    category_counts
    .sort_values(
        "amg_count",
        ascending=False,
    )
    .to_string(index=False)
)


if total_amgs != 789:
    print(
        "\nWARNING: Expected 789 AMG instances "
        f"but found {total_amgs}."
    )


# ---------------------------------------------------------------------
# Save category counts
# ---------------------------------------------------------------------

category_counts.to_csv(
    OUT_TSV,
    sep="\t",
    index=False,
)


# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(10, 7)
)

ax.barh(
    category_counts["metabolic_category"],
    category_counts["amg_count"],
)

ax.set_xlabel(
    "Number of AMG instances"
)

ax.set_ylabel(
    "KEGG metabolic category"
)

ax.set_title(
    "Metabolic categories of auxiliary metabolic genes"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()


# ---------------------------------------------------------------------
# Save figure
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


print(
    f"\nSaved: {OUT_PNG}"
)

print(
    f"Saved: {OUT_SVG}"
)

print(
    f"Saved: {OUT_TSV}"
)

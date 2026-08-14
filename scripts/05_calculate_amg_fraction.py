#!/usr/bin/env python3

"""
Calculate the fraction of viral community abundance attributable to
AMG-carrying viral genomes.

Inputs:
    1. Breadth-filtered viral genome RPKM matrix
    2. VIBRANT AMG individual-gene output

Logic:
    - Viral genomes appearing in the VIBRANT AMG output are classified
      as AMG-carrying genomes.
    - VIBRANT fragment names ending in "_fragment_N" are resolved to
      their parent viral genome names before classification.
    - For each sample:

        AMG-carrying RPKM
            = sum of RPKM from AMG-carrying viral genomes

        Total viral RPKM
            = sum of RPKM from all viral genomes

        AMG-carrying fraction
            = AMG-carrying RPKM / Total viral RPKM

    - Non-AMG-carrying RPKM is calculated as:

        Total viral RPKM - AMG-carrying RPKM

This implements the numerator/denominator logic used in the Yosemite
virome AMG analysis.
"""

from pathlib import Path
import re
import pandas as pd


RPKM_PATH = Path(
    "results/read_recruitment/rpkm_filtered/"
    "rpkm_matrix_breadth75.tsv"
)

AMG_PATH = Path(
    "results/vibrant_amg/"
    "VIBRANT_derep95_all17/"
    "VIBRANT_results_derep95_all17/"
    "VIBRANT_AMG_individuals_derep95_all17.tsv"
)

OUTPUT_PATH = Path(
    "results/amg/amg_carrying_fraction.tsv"
)


def resolve_genome_name(name):
    """
    Resolve VIBRANT fragment names to their parent viral genome.

    Example:
        genome||full_fragment_2
    becomes:
        genome||full
    """
    return re.sub(r"_fragment_\d+$", "", str(name))


def main():

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Load abundance matrix
    rpkm = pd.read_csv(
        RPKM_PATH,
        sep="\t",
        index_col=0,
    )

    # Load VIBRANT AMG calls
    amg = pd.read_csv(
        AMG_PATH,
        sep="\t",
    )

    # Resolve VIBRANT fragment names
    amg["genome_resolved"] = (
        amg["scaffold"]
        .apply(resolve_genome_name)
    )

    # A genome is AMG-carrying if at least one AMG was identified
    amg_genomes = set(
        amg["genome_resolved"].unique()
    )

    # Identify AMG-carrying genomes represented in the RPKM matrix
    amg_mask = rpkm.index.isin(amg_genomes)

    results = []

    for sample in rpkm.columns:

        total_rpkm = rpkm[sample].sum()

        amg_rpkm = rpkm.loc[
            amg_mask,
            sample,
        ].sum()

        non_amg_rpkm = (
            total_rpkm - amg_rpkm
        )

        if total_rpkm > 0:
            amg_fraction = (
                amg_rpkm / total_rpkm
            )
        else:
            amg_fraction = 0.0

        results.append(
            {
                "sample": sample,
                "total_viral_rpkm": total_rpkm,
                "amg_carrying_rpkm": amg_rpkm,
                "non_amg_rpkm": non_amg_rpkm,
                "amg_carrying_fraction": amg_fraction,
                "amg_carrying_percent": amg_fraction * 100,
            }
        )

    results_df = pd.DataFrame(results)

    results_df.to_csv(
        OUTPUT_PATH,
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    print(
        f"AMG-carrying viral genomes identified: "
        f"{len(amg_genomes)}"
    )

    print(
        f"AMG-carrying genomes present in RPKM matrix: "
        f"{amg_mask.sum()}"
    )

    print("\nAMG-carrying fraction by sample:")
    print(
        results_df[
            [
                "sample",
                "amg_carrying_rpkm",
                "total_viral_rpkm",
                "amg_carrying_fraction",
            ]
        ].to_string(index=False)
    )

    print(
        f"\nSaved: {OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

"""
Calculate viral genome RPKM after applying a sample-specific
breadth-of-coverage filter.

Workflow:
1. Read each sample BAM.
2. Calculate breadth of coverage for each viral genome:
       breadth = covered bases / genome length
3. Retain genome/sample detections with breadth >= 0.75.
4. Sum mapped reads across only the genomes passing that threshold.
5. Calculate RPKM for each passing genome:
       RPKM = mapped_reads * 1e9 / (genome_length * passing_mapped_reads)
6. Set failing genome/sample combinations to 0.

This reproduces the normalization logic recovered from the archived
Yosemite virome analysis. For example, in sample SRR5830790:
    73 genomes passed the 75% breadth threshold
    736,493 mapped reads were assigned to passing genomes

Those 736,493 reads were used as the RPKM denominator.

Software:
    Python 3
    pandas
    SAMtools
"""

from pathlib import Path
import subprocess
import pandas as pd


BAM_DIR = Path("results/read_recruitment/bam")
OUT_DIR = Path("results/read_recruitment/rpkm_filtered")
METADATA = Path("metadata/sample_metadata.tsv")

BREADTH_THRESHOLD = 0.75

OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_command(command):
    """Run a shell command and return stdout."""
    result = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def get_idxstats(bam_path):
    """
    Return genome length and mapped-read counts from samtools idxstats.
    """
    output = run_command(
        ["samtools", "idxstats", str(bam_path)]
    )

    records = []

    for line in output.strip().splitlines():
        genome, length, mapped, unmapped = line.split("\t")

        # Ignore the special '*' record
        if genome == "*":
            continue

        records.append(
            {
                "genome": genome,
                "length": int(length),
                "mapped": int(mapped),
            }
        )

    return pd.DataFrame(records).set_index("genome")


def get_breadth(bam_path, genome, genome_length):
    """
    Calculate fraction of genome positions covered by >=1 read.

    samtools depth -aa is used so zero-coverage positions are included.
    """

    output = run_command(
        [
            "samtools",
            "depth",
            "-aa",
            "-r",
            genome,
            str(bam_path),
        ]
    )

    covered = 0
    positions = 0

    for line in output.splitlines():
        fields = line.split("\t")

        if len(fields) < 3:
            continue

        depth = int(fields[2])

        positions += 1

        if depth > 0:
            covered += 1

    # Fallback in case no positions are returned
    if positions == 0:
        return 0.0

    # positions should equal genome_length when -aa is used
    return covered / genome_length


def process_sample(sample):
    """
    Calculate breadth-filtered RPKM values for one sample.
    """

    bam_path = BAM_DIR / f"{sample}.bam"

    print(f"Processing {sample}")

    stats = get_idxstats(bam_path)

    breadth_values = {}

    for genome, row in stats.iterrows():

        # Genomes with no mapped reads cannot pass breadth filtering
        if row["mapped"] == 0:
            breadth_values[genome] = 0.0
            continue

        breadth_values[genome] = get_breadth(
            bam_path,
            genome,
            row["length"],
        )

    stats["breadth"] = pd.Series(breadth_values)

    stats["passes_breadth"] = (
        stats["breadth"] >= BREADTH_THRESHOLD
    )

    passing_read_total = stats.loc[
        stats["passes_breadth"],
        "mapped",
    ].sum()

    print(
        f"  Genomes passing breadth filter: "
        f"{stats['passes_breadth'].sum()}"
    )

    print(
        f"  Mapped reads on passing genomes: "
        f"{passing_read_total}"
    )

    stats["rpkm"] = 0.0

    if passing_read_total > 0:

        passing = stats["passes_breadth"]

        stats.loc[passing, "rpkm"] = (
            stats.loc[passing, "mapped"]
            * 1e9
            / (
                stats.loc[passing, "length"]
                * passing_read_total
            )
        )

    return stats["rpkm"]


def main():

    metadata = pd.read_csv(METADATA, sep="\t")

    samples = metadata["sample_id"].tolist()

    matrix = {}

    for sample in samples:
        matrix[sample] = process_sample(sample)

    rpkm_matrix = pd.DataFrame(matrix)

    rpkm_matrix.index.name = "genome"

    output_path = OUT_DIR / "rpkm_matrix_breadth75.tsv"

    rpkm_matrix.to_csv(
        output_path,
        sep="\t",
        float_format="%.4f",
    )

    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()

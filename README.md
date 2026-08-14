# Yosemite Alpine Lake Virome and Auxiliary Metabolic Gene Analysis

Analysis code and metadata supporting the manuscript:

**Viral community structure and auxiliary metabolic genes in Yosemite alpine lakes**

## Overview

This repository contains scripts used to analyze viral communities and virus-encoded auxiliary metabolic genes (AMGs) in metagenomes from alpine lakes in Yosemite National Park, California.

The analysis includes 17 metagenomic samples representing baseline lake communities and experimental treatments. Metatranscriptomic datasets associated with the broader study were excluded from the analyses presented here.

The workflow includes viral genome identification and quality assessment, genome dereplication, AMG identification, viral read recruitment, abundance estimation, and generation of manuscript figures.

Raw sequencing reads, BAM files, and other large intermediate files are not stored in this repository.

## Repository structure

```text
yosemite_virome_amg/
├── README.md
├── metadata/
│   └── sample_metadata.tsv
└── scripts/
    ├── README.md
    ├── 01_dereplicate_viral_genomes.sh
    ├── 02_run_vibrant.sh
    ├── 03_map_reads_bowtie2.sh
    ├── 04_calculate_rpkm_breadth.py
    ├── 05_calculate_amg_fraction.py
    └── figures/
        ├── figure1_amg_abundance.py
        ├── figure2_amg_categories.py
        ├── figure3_deazaguanine_cluster.sh
        └── figure3_panelC_gene_cluster.py
```

### `metadata/`

`sample_metadata.tsv` contains metadata for the 17 metagenomes included in the analysis, including lake, year, experiment, treatment, timepoint, repository, and available public accession information.

### `scripts/`

Contains the core analysis workflow. See [`scripts/README.md`](scripts/README.md) for descriptions of individual scripts, required inputs, outputs, and analysis order.

### `scripts/figures/`

Contains scripts used to generate the AMG-related manuscript figures.

## Analysis workflow

The principal workflow was:

1. Quality trimming of metagenomic reads with Trimmomatic.
2. Metagenome assembly with MEGAHIT.
3. Viral sequence identification with VirSorter2.
4. Viral genome quality assessment with CheckV.
5. Dereplication of high- and medium-quality viral genomes using CD-HIT-EST.
6. Identification and annotation of auxiliary metabolic genes using VIBRANT.
7. Recruitment of metagenomic reads to dereplicated viral genomes using Bowtie2.
8. Viral abundance calculation using a 75% breadth-of-coverage filtered RPKM approach.
9. Quantification and characterization of AMG-carrying viral genomes.
10. Generation of manuscript figures.

Detailed commands and parameters for the downstream viral and AMG analyses are provided in the scripts.

## Viral genome dereplication

High- and medium-quality viral genomes were dereplicated using CD-HIT-EST v4.8.1 at:

- 95% nucleotide identity (`-c 0.95`)
- 85% minimum alignment coverage of the shorter sequence (`-aS 0.85`)

This resulted in **695 dereplicated viral genome representatives** used for downstream analyses.

## Viral abundance

Metagenomic reads were recruited to the dereplicated viral genomes using Bowtie2 v2.5.1.

Viral genome abundance was calculated using a **75% breadth-of-coverage filtered RPKM approach**. Genome detections with less than 75% breadth were excluded. For each sample, mapped reads associated with genomes passing the breadth threshold were used as the sample-specific denominator for RPKM calculation.

The implementation is provided in:

```text
scripts/04_calculate_rpkm_breadth.py
```

## Auxiliary metabolic genes

AMGs were identified using VIBRANT v1.2.1. The analysis identified **789 AMG instances** among the dereplicated viral genomes.

Scripts are provided for calculating the abundance of AMG-carrying viral genomes and summarizing AMG metabolic categories.

## Data availability

Raw sequencing reads are not duplicated in this GitHub repository.

Sample identities, experimental metadata, data repositories, and available accession information are provided in:

```text
metadata/sample_metadata.tsv
```

Previously published metagenomic datasets are available through the public repositories identified in the metadata table.

Five newly generated metagenomes (Y86, Y139, Y153, Y164, and Y166) are associated with NCBI BioProject **PRJNA1477654** and BioSample accessions listed in `sample_metadata.tsv`. These sequence data are currently under embargo and are scheduled for public release on December 1, 2026.

Users reproducing the analyses should download or otherwise obtain the appropriate sequencing data and provide local input-file paths when running the analysis scripts.

## Software

Major software versions used in the analysis were:

- Trimmomatic v0.39
- MEGAHIT v1.1.5
- VirSorter2 v2.2.4
- CheckV v1.0.3
- CD-HIT-EST v4.8.1
- VIBRANT v1.2.1
- Bowtie2 v2.5.1
- Pharokka v1.9.1

Additional Python and command-line dependencies are documented within the relevant scripts.

## Reproducibility

This repository is intended to document the computational analyses underlying the manuscript rather than distribute large sequencing or intermediate data files.

Users should:

1. Obtain the sequencing data identified in `metadata/sample_metadata.tsv`.
2. Supply the appropriate local input-file paths to the analysis scripts.
3. Run the analysis workflow using the software versions and parameters documented here and in `scripts/README.md`.

Intermediate files such as FASTQs, BAMs, assemblies, and large annotation outputs are intentionally not tracked in this repository.

## Citation

A permanent archived version of this repository and citation information will be provided through Zenodo upon release of the manuscript-associated version.

## License

This repository is released under the MIT License. See [`LICENSE`](LICENSE) for details.

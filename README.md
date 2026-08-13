# Yosemite Alpine Lake Virome and Auxiliary Metabolic Gene Analysis

Code supporting the manuscript Novel and diverse viral communities in high-elevation alpine lakes of Yosemite National Park reveal extreme endemism and unexpected taxonomic diversity


## Overview

This repository contains scripts used for viral identification, genome quality
assessment, auxiliary metabolic gene identification, read mapping, abundance
normalization, AMG-carrying viral fraction calculations, and generation of
Figures 1–3.

## Workflow

```text
Raw metagenomic reads
    ↓
Trimmomatic
    ↓
MEGAHIT assembly
    ↓
VirSorter2
    ↓
CheckV
    ↓
CD-HIT-EST
    ↓
VIBRANT AMG annotation
    ↓
Bowtie2 read mapping
    ↓
≥75% genome breadth filter
    ↓
RPKM normalization
    ↓
AMG-carrying viral fraction
    ↓
Figures and statistical summaries
```

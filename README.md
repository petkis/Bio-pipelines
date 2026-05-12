# Genome and Telomere Analysis Pipelines

This repository contains PBS scripts for genome alignment analysis and Oxford Nanopore telomere length analysis on an HPC cluster.

The repository is split into two separate pipelines:

1. `alignment/` – alignment simulation, alignment processing, and alignment result plotting.
2. `basecall/` – ONT basecalling, telomere length estimation, and basecaller/tool comparison.

Jobs are submitted with `qsub`. The scripts use `PBS_O_WORKDIR` as the folder where the job was submitted, and `SCRATCHDIR` for temporary job files. `PBS_O_WORKDIR` is the current working directory of the `qsub` command, so jobs must be submitted from the correct pipeline folder.

## Why this repository exists

Telomeric regions are repetitive and difficult to analyze with standard sequencing and alignment workflows. This repository provides reproducible HPC workflows for testing how different alignment references, aligners, basecallers, and telomere detection tools affect telomere-related results.

The alignment pipeline compares mapping behavior against three genome versions:

- genome containing telomeric regions
- genome without telomeric regions
- genome with selected regions masked

The basecall pipeline compares telomere length estimates produced from different ONT basecalling and telomere detection configurations.

## Repository structure

```text
alignment/
├── alignment_pipeline.pbs
├── run_analysis.pbs
├── alignment_tools.sh
├── alignment_processing.sh
├── collect_results.py
├── plot_results.R
├── dataset/
├── meryl_clean/
└── results/

basecall/
├── basecall_pipeline.pbs
├── telo_compare.pbs
├── extract.py
├── telomere_graphs.R
├── basecall_histograms.R
├── basecall_comparison.R
├── dataset/
└── results/

## Alignment pipeline

alignment/dataset/
    genome_masked.fa
    genome_no_telomeres.fa
    genome_telomeres.fa
        │
        ▼
alignment_pipeline.pbs
    - extracts a 250 kb chromosome-arm region
    - prepares Telo / NoTelo / Masked reference regions
    - simulates reads for selected read lengths
    - aligns reads with minimap2 or winnowmap
    - processes alignment results
        │
        ▼
alignment/results/{TOOL}_{ARM}_PROC_RESULTS_{CHROM}/

alignment/results/plots/

## Basecaller pipeline

basecall/dataset/
    Subset.pod5      → Dorado
    fast5_subset/    → Guppy
        │
        ▼
basecall_pipeline.pbs
    - runs Dorado or Guppy
    - creates reads.fastq
    - runs NCRF, TeloBP, or both
    - creates telomere length tables and per-run plots
        │
        ▼
basecall/results/{BASECALLER}_v{VERSION}_{MODEL_TAG}/
        │
        ▼
telo_compare.pbs
    - compares telomere length outputs
    - creates histograms, summary files, and violin plots
        │
        ▼
basecall/results/histograms/

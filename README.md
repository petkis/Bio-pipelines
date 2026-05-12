# Bio-pipelines

This repository contains HPC-based pipelines used for long-read sequencing analysis in telomeric and chromosome-level regions.

The project contains two main workflows:

1. **Alignment pipeline** – simulation and alignment of reads from selected chromosome-arm regions.
2. **Basecalling pipeline** – Oxford Nanopore basecalling and telomere length estimation.

The scripts are designed for execution on a PBS-based HPC cluster using `qsub`. Most computation is performed in `$SCRATCHDIR`, and final outputs are copied back to the corresponding `results/` folders.

---

## Purpose

Telomeric and other repetitive genomic regions are difficult to analyze because they can affect both basecalling and read alignment. This repository provides reproducible scripts for evaluating how different tools and reference genome variants influence the analysis of long-read sequencing data.

The repository was created as part of a bachelor thesis focused on long-read sequencing, telomere-related analysis, basecaller comparison, and alignment behavior in chromosome-level regions.

---

## Repository structure

```text
Bio-pipelines/
├── alignment/
│   ├── README.md
│   ├── alignment_job.pbs
│   ├── alignment_tools.sh
│   ├── alignment_processing.sh
│   ├── collect_results.py
│   ├── plot_results.R
│   ├── run_analysis.pbs
│   ├── multi_sub_job.sh
│   ├── chromosome_list.txt
│   ├── dataset/
│   ├── meryl_clean/
│   └── results/
│
└── basecall/
    ├── README.md
    ├── basecall_pipeline.pbs
    ├── visualise_results.pbs
    ├── extract.py
    ├── telomere_graphs.R
    ├── basecall_histograms.R
    ├── basecall_comparison.R
    ├── dataset/
    └── results/
```

---

## Pipelines

### Alignment pipeline

The alignment pipeline simulates reads from selected chromosome-arm regions and aligns them to different reference genome variants.

It compares three reference types:

| Reference type | Description |
|---|---|
| `Telo` | Reference genome containing telomeric regions |
| `NoTelo` | Reference genome without telomeric regions |
| `Masked` | Reference genome with selected regions masked |

The pipeline supports two aligners:

| Tool | Purpose |
|---|---|
| `minimap2` | Standard long-read aligner used as a baseline |
| `winnowmap` | Long-read aligner designed to improve mapping in repetitive regions |

Detailed instructions are available in:

```text
alignment/README.md
```

---

### Basecalling pipeline

The basecalling pipeline processes Oxford Nanopore raw signal data, generates long reads, estimates telomere lengths, and creates comparison plots.

It supports two basecallers:

| Basecaller | Input format |
|---|---|
| Dorado | POD5 |
| Guppy | FAST5 |

It supports telomere length estimation using:

| Tool | Purpose |
|---|---|
| NCRF | Detection of telomeric repeat regions in reads |
| TeloBP | Telomere length estimation from long-read sequencing data |

Detailed instructions are available in:

```text
basecall/README.md
```

---

## Dataset availability

Large input datasets are not stored directly in this Git repository.

Expected dataset locations:

```text
alignment/dataset/
basecall/dataset/
```

A permanent archive of the dataset and selected generated results should be provided through Zenodo or another long-term storage service.

```text
Zenodo DOI: TODO
```

---

## Current results

Generated results are stored in:

```text
alignment/results/
basecall/results/
```

The main result types include:

- alignment statistics,
- chromosome-level mapping summaries,
- basecaller comparison tables,
- telomere length estimates,
- histograms,
- violin plots,
- runtime and computational requirement summaries.

---

## Running the workflows

All major workflows are executed as PBS jobs.

Example alignment job:

```bash
cd alignment
qsub -v CHROM=chr21_PATERNAL,ARM=p,TOOL=winnowmap alignment_job.pbs
```

Example basecalling job:

```bash
cd basecall
qsub -v BASECALLER=dorado,VERSION=1.4.0,MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0,TELO_TOOL=both basecall_pipeline.pbs
```

---

## General notes

- Jobs should be submitted from the correct pipeline directory.
- The scripts rely on `$PBS_O_WORKDIR` and `$SCRATCHDIR`.
- Large sequencing files and large generated outputs should not be committed to Git.
- The repository is intended to document and reproduce the computational workflows used in the bachelor thesis.

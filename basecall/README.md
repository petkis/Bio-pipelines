# Basecalling Pipeline

This folder contains the Oxford Nanopore basecalling and telomere length estimation workflow.

The pipeline is designed for execution on a PBS-based HPC cluster with GPU support.

---

## Purpose

The goal of this pipeline is to compare different Oxford Nanopore basecalling configurations and evaluate their influence on telomere length estimation.

The workflow processes raw nanopore signal data, generates long reads, estimates telomere lengths, and creates summary plots for comparison between basecallers, versions, models, and telomere detection tools.

---

## Supported basecallers

| Basecaller | Input format | Description |
|---|---|---|
| Dorado | POD5 | Current Oxford Nanopore basecaller |
| Guppy | FAST5 | Older Oxford Nanopore basecaller |

---

## Supported telomere detection tools

| Tool | Description |
|---|---|
| NCRF | Detects telomeric repeat regions in reads |
| TeloBP | Estimates telomere lengths from long-read sequencing data |

---

## Pipeline overview

```text
POD5 / FAST5 input data
    ↓
Basecalling with Dorado or Guppy
    ↓
FASTQ reads
    ↓
FASTQ to FASTA conversion
    ↓
Telomere detection using NCRF and/or TeloBP
    ↓
telomere_lengths.tsv
    ↓
Histograms and violin plots
```

---

## Input files

Input files must be placed in:

```text
basecall/dataset/
```

Expected structure:

```text
basecall/dataset/
├── Subset.pod5
└── fast5_subset/
```

Dorado uses:

```text
Subset.pod5
```

Guppy uses:

```text
fast5_subset/
```

Large raw sequencing files should not be committed to Git.

---

## Main scripts

| Script | Purpose |
|---|---|
| `basecall_pipeline.pbs` | Main PBS job. Runs Dorado or Guppy and then performs telomere length estimation. |
| `extract.py` | Parses NCRF output and creates a standardized `telomere_lengths.tsv` file. |
| `telomere_graphs.R` | Creates a telomere length histogram for a single run. |
| `visualise_results.pbs` | Runs visualization scripts for collected basecalling results. |
| `basecall_histograms.R` | Creates grouped histogram plots from multiple runs. |
| `basecall_comparison.R` | Creates violin plots comparing basecallers, versions, models, and telomere detection tools. |

---

## Running basecalling

Submit the job from the `basecall/` directory.

Example Dorado run:

```bash
cd basecall
qsub -v BASECALLER=dorado,VERSION=1.4.0,MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0,TELO_TOOL=both basecall_pipeline.pbs
```

Example Guppy run:

```bash
cd basecall
qsub -v BASECALLER=guppy,VERSION=6.5.7,MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg,TELO_TOOL=both basecall_pipeline.pbs
```

Available variables:

| Variable | Example | Description |
|---|---|---|
| `BASECALLER` | `dorado`, `guppy` | Basecaller to run |
| `VERSION` | `1.4.0`, `6.5.7` | Basecaller version |
| `MODEL` | Dorado model name or Guppy config file | Basecalling model/configuration |
| `TELO_TOOL` | `ncrf`, `telobp`, `both` | Telomere detection tool |

If variables are not provided, default values defined inside `basecall_pipeline.pbs` are used.

---

## Generating comparison plots

After multiple basecalling runs are finished, run:

```bash
cd basecall
qsub visualise_results.pbs
```

This creates summary histograms and violin plots from collected `telomere_lengths.tsv` files.

---

## Output files

Each basecalling run creates a separate output directory under:

```text
basecall/results/
```

Example output directory:

```text
basecall/results/dorado_v1.4.0_dna_r10_4_1_e8_2_400bps_sup_v5_2_0/
```

Important output files include:

```text
reads.fastq
reads.fasta
ncrf/telomere_lengths.tsv
ncrf/*.ncrf
ncrf/*.pdf
telobp/telomere_lengths.tsv
telobp/*.csv
telobp/*.pdf
```

Comparison outputs are written to:

```text
basecall/results/histograms/
```

Important comparison files include:

```text
histogram_summary_stats.tsv
summary_stats.tsv
compare_basecaller_combined.pdf
compare_basecaller_combined_log.pdf
compare_versions_sup_telobp_violin.pdf
compare_versions_sup_ncrf_violin.pdf
```

---

## Output interpretation

The main output of the pipeline is a standardized table of estimated telomere lengths:

```text
telomere_lengths.tsv
```

This file is used to compare telomere length distributions between:

- basecallers,
- basecaller versions,
- basecalling models,
- telomere detection tools.

The generated histograms and violin plots are used to visualize differences in telomere length estimates between individual pipeline configurations.

---

## Dependencies

The pipeline uses:

```text
PBS/qsub
Dorado
Guppy
samtools
Python 3
R
ggplot2
dplyr
stringr
purrr
readr
patchwork
NCRF
TeloBP
CUDA-compatible GPU
```

The exact environment depends on the HPC cluster configuration.

---

## Notes

- Jobs must be submitted from the `basecall/` folder.
- Dorado expects POD5 input.
- Guppy expects FAST5 input.
- GPU basecalling requires a compatible CUDA environment.
- Temporary computation is performed in `$SCRATCHDIR`.
- Large FASTQ, BAM, POD5, FAST5, and result files should not be committed to Git.

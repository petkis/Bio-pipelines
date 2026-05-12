# Alignment Pipeline

This folder contains the alignment workflow used to test how simulated long reads from selected chromosome-arm regions map to different reference genome variants.

The pipeline is designed for execution on a PBS-based HPC cluster.

---

## Purpose

The goal of this pipeline is to evaluate how the presence, absence, or masking of telomeric regions in the reference genome affects long-read alignment.

The workflow simulates reads from selected chromosome-arm regions and aligns them to three reference genome variants:

| Reference type | Description |
|---|---|
| `Telo` | Reference genome containing telomeric regions |
| `NoTelo` | Reference genome without telomeric regions |
| `Masked` | Reference genome with selected regions masked |

Two aligners can be used:

| Aligner | Description |
|---|---|
| `minimap2` | Baseline long-read aligner |
| `winnowmap` | Long-read aligner with improved handling of repetitive regions |

---

## Pipeline overview

```text
Reference genomes
    ↓
Select chromosome and chromosome arm
    ↓
Extract 250 kb chromosome-arm region
    ↓
Prepare Telo / NoTelo / Masked region FASTA files
    ↓
Simulate reads using wgsim
    ↓
Align reads using minimap2 or winnowmap
    ↓
Convert SAM to sorted and indexed BAM
    ↓
Generate alignment statistics
    ↓
Collect results into CSV tables
    ↓
Create chromosome-level plots
```

---

## Input files

Reference genomes must be placed in:

```text
alignment/dataset/
```

Expected files:

```text
genome_masked.fa
genome_no_telomeres.fa
genome_telomeres.fa
```

For Winnowmap, repetitive k-mer files are expected in:

```text
alignment/meryl_clean/
```

Expected files:

```text
repetitive_k15_telo.txt
repetitive_k15_notelo.txt
repetitive_k15_masked.txt
```

---

## Main scripts

| Script | Purpose |
|---|---|
| `alignment_job.pbs` | Main PBS job. Runs one complete alignment workflow for a selected chromosome, arm, and aligner. |
| `alignment_tools.sh` | Simulates reads and aligns them using `minimap2` or `winnowmap`. |
| `alignment_processing.sh` | Converts SAM to BAM, sorts and indexes BAM files, and generates mapping statistics. |
| `collect_results.py` | Collects chromosome-level alignment summaries into a combined CSV file. |
| `plot_results.R` | Creates plots from the collected alignment results. |
| `run_analysis.pbs` | Runs result aggregation and plotting on the HPC cluster. |
| `multi_sub_job.sh` | Submits multiple alignment jobs using chromosomes listed in `chromosome_list.txt`. |
| `chromosome_list.txt` | List of chromosomes used for batch submission. |

---

## Running one alignment job

Submit the job from the `alignment/` directory:

```bash
cd alignment
qsub -v CHROM=chr21_PATERNAL,ARM=p,TOOL=winnowmap alignment_job.pbs
```

Available variables:

| Variable | Example | Description |
|---|---|---|
| `CHROM` | `chr21_PATERNAL` | Chromosome to analyze |
| `ARM` | `p` or `q` | Chromosome arm |
| `TOOL` | `minimap2` or `winnowmap` | Aligner used for mapping |

If variables are not provided, default values defined inside `alignment_job.pbs` are used.

---

## Running multiple alignment jobs

To submit jobs for multiple chromosomes, edit:

```text
chromosome_list.txt
```

Then run:

```bash
cd alignment
bash multi_sub_job.sh
```

Each chromosome listed in `chromosome_list.txt` is submitted as a separate PBS job.

---

## Generating summary tables and plots

After the alignment jobs finish, run:

```bash
cd alignment
qsub -v MODE=same run_analysis.pbs
```

The `MODE` variable controls how correct mapping is evaluated.

| Mode | Description |
|---|---|
| `same` | Correct mapping is evaluated against the same chromosome or haplotype. |
| `cross` | Correct mapping is evaluated against the opposite chromosome or haplotype pair. |

---

## Output files

Main outputs are written to:

```text
alignment/results/
```

Summary plots and collected result tables are written to:

```text
alignment/results/plots/
```

Important output files include:

```text
aligned.sorted.bam
aligned.sorted.bam.bai
idxstats.txt
bam_summary_stats.txt
per_chrom_all.txt
per_chrom_mapq20.txt
mapq_hist.txt
combined_results.csv
All_chromosomes_<tool>.pdf
```

---

## Output interpretation

The most important output for the thesis is the chromosome-level mapping summary.

The pipeline tracks how many simulated reads map to the expected chromosome and how this changes depending on:

- reference genome type,
- read length,
- chromosome arm,
- aligner,
- mapping quality threshold.

These results are used to compare the behavior of `minimap2` and `winnowmap`, especially in difficult or repetitive chromosome regions.

---

## Dependencies

The pipeline uses:

```text
PBS/qsub
samtools
wgsim
minimap2
winnowmap
mamba or conda
Python 3
R
tidyverse
rsync
```

The exact environment depends on the HPC cluster configuration.

---

## Notes

- Jobs must be submitted from the `alignment/` folder.
- The scripts use `$PBS_O_WORKDIR` to locate the working directory.
- Temporary computation is performed in `$SCRATCHDIR`.
- Large SAM/BAM files should not be committed to Git.
- The selected read lengths are defined inside `alignment_job.pbs`.
- The number of simulated reads is defined inside `alignment_tools.sh`.

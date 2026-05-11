# Genome and Telomere Analysis Pipelines

This repository contains PBS scripts for genome/telomere analysis on an HPC cluster.

Main workflows:

1. alignment simulation and processing
2. alignment result plotting
3. ONT basecalling and telomere length estimation
4. basecalling/telomere result comparison

Jobs are submitted with `qsub`. The scripts use `$PBS_O_WORKDIR` as the folder where the job was submitted and `$SCRATCHDIR` for temporary job files. `PBS_O_WORKDIR` is the working directory of the `qsub` command.

## Dataset note

The `dataset/` folders are intentionally empty in GitHub because the real genome, POD5, FAST5, and sequencing files are too large for normal GitHub storage. GitHub warns for files larger than 50 MiB and blocks files larger than 100 MiB.

Before running anything, copy the required data into the correct `dataset/` folder. See the README files inside the dataset folders for details.

## Folder structure

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
```

## Required input files

### Alignment

Place these files in `alignment/dataset/`:

```text
genome_masked.fa
genome_no_telomeres.fa
genome_telomeres.fa
```

Purpose:

| File | Meaning |
|---|---|
| `genome_masked.fa` | genome with selected regions masked |
| `genome_no_telomeres.fa` | genome without telomeric regions |
| `genome_telomeres.fa` | genome containing telomeric regions |

The alignment pipeline also uses:

```text
alignment/meryl_clean/repetitive_k15_telo.txt
alignment/meryl_clean/repetitive_k15_notelo.txt
alignment/meryl_clean/repetitive_k15_masked.txt
```

### Basecalling

Place these files in `basecall/dataset/`:

```text
Subset.pod5
fast5_subset/
```

Use `Subset.pod5` for Dorado and `fast5_subset/` for Guppy.

## Workflow 1: alignment pipeline

Run from the `alignment/` folder:

```bash
cd alignment
qsub alignment_pipeline.pbs
```

The script extracts a 250 kb region from the selected chromosome and runs alignment against the three genome versions: `Telo`, `NoTelo`, and `Masked`.

Default parameters:

```bash
TOOL=winnowmap
CHROM=chr21_PATERNAL
ARM=p
```

Parameters:

| Parameter | Values | Meaning |
|---|---|---|
| `TOOL` | `winnowmap`, `minimap2` | aligner |
| `CHROM` | chromosome name | chromosome/contig to process |
| `ARM` | `p`, `q` | start or end of chromosome |

Examples:

```bash
qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs
```

Read lengths are changed manually inside `alignment_pipeline.pbs`:

```bash
LENGTHS=(10000 20000 30000 40000)
```

Output:

```text
alignment/results/{TOOL}_{ARM}_PROC_RESULTS_{CHROM}/
```

## Workflow 2: plot alignment results

Run after the alignment jobs finish:

```bash
cd alignment
qsub run_analysis.pbs
```

This script collects `per_chrom_all.txt` files and creates combined plots.

Output:

```text
alignment/results/plots/
```

## Workflow 3: basecalling and telomere detection

Run from the `basecall/` folder:

```bash
cd basecall
qsub basecall_pipeline.pbs
```

This script runs Dorado or Guppy, creates FASTQ reads, and detects telomeres with NCRF, TeloBP, or both.

Default parameters:

```bash
BASECALLER=dorado
VERSION=1.0.0
MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0
TELO_TOOL=both
```

Parameters:

| Parameter | Values | Meaning |
|---|---|---|
| `BASECALLER` | `dorado`, `guppy` | basecaller |
| `VERSION` | e.g. `1.0.0` | basecaller version |
| `MODEL` | model/config name | Dorado model or Guppy config |
| `TELO_TOOL` | `ncrf`, `telobp`, `both` | telomere detection tool |

Examples:

```bash
qsub basecall_pipeline.pbs
qsub -v BASECALLER=dorado,TELO_TOOL=ncrf basecall_pipeline.pbs
qsub -v BASECALLER=dorado,TELO_TOOL=telobp basecall_pipeline.pbs
qsub -v BASECALLER=guppy,VERSION=6.5.7,MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg,TELO_TOOL=both basecall_pipeline.pbs
```

Output:

```text
basecall/results/{BASECALLER}_v{VERSION}_{MODEL_TAG}/
```

Main outputs:

```text
reads.fastq
ncrf/telomere_lengths.tsv
telobp/telomere_lengths.tsv
```

## Workflow 4: compare basecalling results

Run after one or more basecalling jobs finish:

```bash
cd basecall
qsub telo_compare.pbs
```

This script compares telomere length results and creates histograms/summary files.

Output:

```text
basecall/results/histograms/
```

## Typical run order

```bash
cd alignment

qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs

qsub run_analysis.pbs
```

```bash
cd basecall

qsub -v BASECALLER=dorado,TELO_TOOL=both basecall_pipeline.pbs
qsub telo_compare.pbs
```

## Logs

Logs are written to:

```text
job_logs/
```

Some jobs also write logs directly into the workflow folder.

## Notes

- Submit jobs from the correct workflow folder.
- Put large input data into `dataset/` manually.
- Do not commit large sequencing/genome files to GitHub.
- Alignment outputs are stored in `alignment/results/`.
- Basecalling outputs are stored in `basecall/results/`.
- Temporary files are created in `$SCRATCHDIR` and cleaned after the job.

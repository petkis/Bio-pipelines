# Genome and Telomere Analysis Pipelines

This repository contains PBS job scripts and helper scripts for running genome/telomere analysis workflows on an HPC cluster.

The repository is organized around four main workflows:

1. alignment simulation and processing
2. aggregation and plotting of alignment results
3. ONT basecalling and telomere length estimation
4. comparison and visualization of basecalling/telomere outputs

The pipelines are designed to run from the repository root using `qsub`. PBS variables such as `$PBS_O_WORKDIR`, `$PBS_JOBID`, and `$SCRATCHDIR` are used throughout the scripts. `$PBS_O_WORKDIR` is the directory from which the job was submitted, and `$PBS_JOBID` is the scheduler-assigned job ID. :contentReference[oaicite:1]{index=1}

## Important note about datasets

The `dataset/` folders in this GitHub repository are intentionally empty or contain only small placeholder/README files.

Large genome FASTA files, POD5 files, FAST5 folders, and other sequencing data are not included because of GitHub file-size limitations. GitHub warns about files larger than 50 MiB and blocks files larger than 100 MiB. :contentReference[oaicite:2]{index=2}

Before running the scripts, place the required input data into the correct `dataset/` locations. Dataset-specific README files inside the `dataset/` folders describe what should be copied there.

## General usage

Submit jobs from the repository root:

```bash
cd /path/to/repository
qsub script_name.pbs
```

Many parameters can be changed at submission time with `qsub -v`:

```bash
qsub -v VARIABLE=value script_name.pbs
```

For multiple variables:

```bash
qsub -v VARIABLE1=value1,VARIABLE2=value2 script_name.pbs
```

The scripts assume that they are submitted from the repository root. They use:

```bash
$PBS_O_WORKDIR
```

as the working directory and:

```bash
$SCRATCHDIR
```

as temporary scratch storage on the compute node.

## Required software

The PBS scripts load software through environment modules. Depending on the workflow, the following modules or tools are used:

```text
samtools
mambaforge
metabase/1
python36-modules-gcc
python/3.11.11-gcc-10.2.1-555dlyc
r-tidyverse
r-ggplot2
r/4.4.0-gcc-10.2.1-oxdi5pz
```

External tools used by the workflows include:

```text
winnowmap
minimap2
wgsim
Dorado
Guppy
NCRF
TeloBP
samtools
R
Python
```

`samtools faidx` is used to index FASTA files and extract chromosome regions. The `.fai` index allows efficient retrieval of specific FASTA/FASTQ regions. :contentReference[oaicite:3]{index=3}

## Input files

### Genome FASTA files

The alignment pipeline expects three versions of the genome:

```text
dataset/genome_masked.fa
dataset/genome_no_telomeres.fa
dataset/genome_telomeres.fa
```

These represent:

| File | Purpose |
|---|---|
| `genome_masked.fa` | Genome version with selected repetitive/telomeric regions masked. Used to test how masking affects alignment. |
| `genome_no_telomeres.fa` | Genome version with telomeric sequence removed. Used as a no-telomere comparison. |
| `genome_telomeres.fa` | Genome version containing telomeric sequence. Used as the telomere-positive reference. |

The alignment script extracts the same chromosome arm region from all three FASTA files and compares how simulated reads align to each version.

### Repetitive k-mer files

The alignment pipeline also expects:

```text
meryl_clean/repetitive_k15_telo.txt
meryl_clean/repetitive_k15_notelo.txt
meryl_clean/repetitive_k15_masked.txt
```

These files are used by the alignment tools, especially for workflows involving Winnowmap and repetitive k-mer filtering.

### ONT basecalling inputs

The basecalling pipeline expects either POD5 or FAST5 input, depending on the selected basecaller.

For Dorado:

```text
dataset/Subset.pod5
```

For Guppy:

```text
dataset/fast5_subset/
```

These files are not included in the GitHub repository. Place them manually according to the README files inside the dataset folders.

### Example submissions

(Run the commands from the alignment folder)
Run the default alignment workflow:

```bash
qsub alignment_pipeline.pbs
```

Run Winnowmap on the p arm of `chr21_PATERNAL`:

```bash
qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
```

Run Winnowmap on the q arm:

```bash
qsub -v TOOL=winnowmap,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs
```

Run Minimap2 on the p arm:

```bash
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=p alignment_pipeline.pbs
```

Run Minimap2 on the q arm:

```bash
qsub -v TOOL=minimap2,CHROM=chr21_PATERNAL,ARM=q alignment_pipeline.pbs
```
Used Lengths are changed manually inside the alignment_pipeline.pbs script

## Workflow 2: collect and plot alignment results

PBS job:

```text
genome_analysis_viz.pbs
```

PBS job name:

```text
genome_analysis_viz
```

Should be run after your desired alignment_pipeline.pbs runs.

## Workflow 3: basecalling and telomere length estimation
Should be run from basecall folder
PBS job:

```text
basecall_pipeline.pbs
```

PBS job name:

```text
basecall_pipeline
```

### Purpose

This workflow basecalls ONT reads and then estimates telomere lengths using NCRF, TeloBP, or both.

It supports:

```text
Dorado
Guppy
NCRF
TeloBP
```
### Parameters

Default values:

```bash
BASECALLER=dorado
VERSION=1.0.0
MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0
TELO_TOOL=both
```

Parameter descriptions:

| Parameter | Allowed/example values | Meaning |
|---|---|---|
| `BASECALLER` | `dorado`, `guppy` | Which ONT basecaller to use. |
| `VERSION` | e.g. `1.0.0` | Basecaller version to download/use. |
| `MODEL` | Dorado model name or Guppy config file | Basecalling model/configuration. |
| `TELO_TOOL` | `ncrf`, `telobp`, `both` | Telomere detection method. |

### Dorado input

For Dorado, the script expects:

```text
dataset/Subset.pod5
```

The script copies it to scratch as:

```text
input.pod5
```

Dorado produces BAM output, which is converted to FASTQ with:

```bash
samtools fastq basecalls.bam > reads.fastq
```

### Guppy input

For Guppy, the script expects:

```text
dataset/fast5_subset/
```

### Example submissions

Run Dorado with the default model and both telomere tools:

```bash
qsub basecall_pipeline.pbs
```

Run Dorado with NCRF only:

```bash
qsub -v BASECALLER=dorado,TELO_TOOL=ncrf basecall_pipeline.pbs
```

Run Dorado with TeloBP only:

```bash
qsub -v BASECALLER=dorado,TELO_TOOL=telobp basecall_pipeline.pbs
```

Run Guppy:

```bash
qsub -v BASECALLER=guppy,VERSION=<guppy_version>,MODEL=<guppy_config>.cfg,TELO_TOOL=both basecall_pipeline.pbs
```

Example Guppy-style command:

```bash
qsub -v BASECALLER=guppy,VERSION=6.5.7,MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg,TELO_TOOL=both basecall_pipeline.pbs
```
## Workflow 4: compare basecalling and telomere results

PBS job:

```text
telo_compare.pbs
```

PBS job name:

```text
telo_compare
```

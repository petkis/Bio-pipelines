# Alignment Pipeline

This folder contains the alignment workflow for testing how reads from telomeric chromosome-arm regions align to different genome references.

The workflow extracts a 250 kb region from the selected chromosome arm, prepares three reference variants, simulates reads of selected lengths, aligns them with either `winnowmap` or `minimap2`, and processes the resulting alignments.

## Purpose

The alignment pipeline is used to compare how telomeric or masked reference content affects read mapping. It runs the same analysis against three genome versions:

| Reference type | Meaning |
|---|---|
| `Telo` | Genome containing telomeric regions |
| `NoTelo` | Genome without telomeric regions |
| `Masked` | Genome with selected regions masked |

The pipeline can be used to evaluate differences in alignment quality, mapping behavior, and read placement near telomeric chromosome-arm regions.

## Required input files

Place these files in:

```text
alignment/dataset/
```

Required files:

```text
genome_masked.fa
genome_no_telomeres.fa
genome_telomeres.fa
```

The workflow expects all three reference FASTA files to be available before running.

## Main workflow

The main Snakefile is:

```text
alignment/Snakefile
```

The workflow performs the following steps:

1. Selects a chromosome arm and extracts a 250 kb region.
2. Builds reference-specific input files.
3. Simulates reads from the selected region.
4. Aligns reads to each reference genome.
5. Converts, sorts, and indexes alignment files.
6. Produces output files for downstream comparison.

## Configuration

Pipeline settings are defined in the Snakemake configuration file, usually:

```text
alignment/config.yaml
```

Typical configurable values include:

| Parameter | Description |
|---|---|
| `chromosome` | Chromosome used for region extraction |
| `arm` | Chromosome arm or region to analyze |
| `read_lengths` | Read lengths used for simulation |
| `aligner` | Aligner to use, such as `winnowmap` or `minimap2` |
| `threads` | Number of CPU threads used by selected rules |

Update the configuration file before running the workflow.

## Running the pipeline

From the repository root, run:

```bash
snakemake -s alignment/Snakefile --cores 8
```

To perform a dry run first:

```bash
snakemake -s alignment/Snakefile --cores 8 -n
```

To print the commands without executing them:

```bash
snakemake -s alignment/Snakefile --cores 8 -n -p
```

## Output

The pipeline produces alignment-related output files such as:

```text
alignment/results/
alignment/logs/
alignment/benchmark/
```

Depending on the configuration, expected outputs may include:

```text
*.sam
*.bam
*.sorted.bam
*.sorted.bam.bai
```

These files can be used to inspect mapping quality and compare alignment behavior between reference types.

## Dependencies

The pipeline requires Snakemake and the tools used by the selected workflow rules.

Common dependencies include:

```text
snakemake
samtools
minimap2
winnowmap
seqkit
```

Some dependencies may only be required for specific rules or aligner choices.

## Notes

- Make sure the reference FASTA files are correctly named and placed in `alignment/dataset/`.
- Check the configuration file before each run.
- Use a dry run with `-n` before launching large jobs.
- If running on a cluster, adapt the Snakemake command to the local cluster environment.

# Basecalling Pipeline

This folder contains the basecalling workflow for processing raw nanopore sequencing data into basecalled read files.

The workflow is designed to organize raw input data, run basecalling, and prepare output files for downstream analysis such as alignment, quality control, or assembly.

## Purpose

The basecalling pipeline converts raw nanopore signal data into nucleotide sequences. It provides a reproducible workflow for generating FASTQ files from raw sequencing output.

The resulting reads can be used as input for other workflows in this repository, including the alignment pipeline.

## Required input files

Place raw sequencing input data in:

```text
basecall/dataset/
```

Depending on the sequencing platform and basecaller configuration, input files may include:

```text
*.pod5
*.fast5
```

Use the input format expected by the selected basecalling tool.

## Main workflow

The main Snakefile is:

```text
basecall/Snakefile
```

The workflow generally performs the following steps:

1. Reads raw nanopore signal files from the dataset folder.
2. Runs the selected basecaller.
3. Produces basecalled reads in FASTQ format.
4. Optionally compresses or organizes output files.
5. Stores logs and benchmark information for reproducibility.

## Configuration

Pipeline settings are defined in the Snakemake configuration file, usually:

```text
basecall/config.yaml
```

Typical configurable values include:

| Parameter | Description |
|---|---|
| `input_dir` | Directory containing raw signal files |
| `output_dir` | Directory for basecalled output files |
| `model` | Basecalling model to use |
| `device` | CPU or GPU device setting |
| `threads` | Number of CPU threads |
| `basecaller` | Basecalling tool used by the workflow |

Update the configuration file before running the workflow.

## Running the pipeline

From the repository root, run:

```bash
snakemake -s basecall/Snakefile --cores 8
```

To perform a dry run first:

```bash
snakemake -s basecall/Snakefile --cores 8 -n
```

To print the commands without executing them:

```bash
snakemake -s basecall/Snakefile --cores 8 -n -p
```

## Output

The pipeline produces basecalled read files and supporting output, usually under:

```text
basecall/results/
basecall/logs/
basecall/benchmark/
```

Expected read output may include:

```text
*.fastq
*.fastq.gz
```

These files can be used for downstream workflows such as read alignment or quality assessment.

## Dependencies

The pipeline requires Snakemake and the selected basecalling software.

Common dependencies may include:

```text
snakemake
dorado
guppy
samtools
```

The exact dependencies depend on which basecaller is configured in the workflow.

## Notes

- Make sure the raw input files are placed in `basecall/dataset/`.
- Check that the selected basecalling model matches the sequencing chemistry and flow cell.
- GPU basecalling may require a compatible CUDA environment.
- Use a dry run with `-n` before launching large jobs.
- If running on a cluster, adapt the Snakemake command to the local cluster environment.

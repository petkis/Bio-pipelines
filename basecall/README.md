
---

## `basecall/README.md`

```markdown
# Basecalling and Telomere Detection Pipeline

This folder contains the ONT basecalling and telomere length estimation workflow.

The pipeline can run Dorado or Guppy, convert basecalled reads to FASTQ, detect telomeric reads with NCRF and/or TeloBP, and compare telomere length estimates between basecaller/tool configurations.

## Purpose

The basecalling pipeline is used to compare how different ONT basecallers, basecaller versions, models, and telomere detection tools affect estimated telomere lengths.

It supports:

- Dorado basecalling from POD5 input
- Guppy basecalling from FAST5 input
- telomere detection with NCRF
- telomere detection with TeloBP
- histogram and violin plot comparison of telomere length distributions

## Required input files

Place these files in `basecall/dataset/`:

```text
Subset.pod5
fast5_subset/

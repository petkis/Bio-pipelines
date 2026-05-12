
## `alignment/README.md`

```markdown
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

## Required input files

Place these files in `alignment/dataset/`:

```text
genome_masked.fa
genome_no_telomeres.fa
genome_telomeres.fa

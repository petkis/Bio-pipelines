# Dataset files

Large genome FASTA files are not stored in this GitHub repository because GitHub blocks regular Git files larger than 100 MiB.

This folder should contain the reference genome FASTA files required by `alignment_job.pbs`:

- `genome_masked.fa`
- `genome_no_telomeres.fa`
- `genome_telomeres.fa`

These files are expected to be prepared externally and copied into this folder before running the alignment pipeline.

The pipeline indexes these FASTA files using:

```bash
samtools faidx genome_masked.fa
samtools faidx genome_no_telomeres.fa
samtools faidx genome_telomeres.fa

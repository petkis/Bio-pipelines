# Winnowmap repetitive k-mer files

Large/generated Meryl files are not stored in this GitHub repository.

This folder should contain the repetitive k-mer lists required by Winnowmap in `alignment_tools.sh`:

- `repetitive_k15_telo.txt`
- `repetitive_k15_notelo.txt`
- `repetitive_k15_masked.txt`

These files are used with Winnowmap's `-W` option.

## How to generate the files

Generate one repetitive 15-mer file for each reference FASTA:

```bash
meryl count k=15 output merylDB_telo ../dataset/genome_telomeres.fa
meryl print greater-than distinct=0.9998 merylDB_telo > repetitive_k15_telo.txt

meryl count k=15 output merylDB_notelo ../dataset/genome_no_telomeres.fa
meryl print greater-than distinct=0.9998 merylDB_notelo > repetitive_k15_notelo.txt

meryl count k=15 output merylDB_masked ../dataset/genome_masked.fa
meryl print greater-than distinct=0.9998 merylDB_masked > repetitive_k15_masked.txt

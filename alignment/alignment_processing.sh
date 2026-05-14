#!/bin/bash

# Process SAM alignment outputs for one chromosome and one reference subtype.
# The script converts SAM files to sorted/indexed BAM files, creates basic BAM statistics,
# and summarizes read distribution across reference sequences and mapping qualities.

set -euo pipefail

# Read input arguments from the main pipeline script.
CHROM="$1"        # e.g. chr1_MATERNAL
CHR_LENGTH="$2"   # e.g. 109226264
SUBTYPE="$3"      # Masked | NoTelo | Telo

echo "=== proc_${SUBTYPE}.sh ==="
echo "CHROM       = $CHROM"
echo "CHR_LENGTH  = $CHR_LENGTH"
echo "SUBTYPE     = $SUBTYPE"

# Select the correct reference genome according to the subtype.
case "$SUBTYPE" in
  Masked) REF_FA="genome_masked.fa" ;;
  NoTelo) REF_FA="genome_no_telomeres.fa" ;;
  Telo)   REF_FA="genome_telomeres.fa" ;;
  *) echo "ERROR: invalid SUBTYPE: $SUBTYPE"; exit 1 ;;
esac

MM2_DIR="${CHROM}_${SUBTYPE}_align"
OUTBASE="${CHROM}_${SUBTYPE}_post_alignments"

# Check that the selected reference FASTA and its index exist.
if [[ ! -f "${REF_FA}" ]]; then
    echo "ERROR: reference not found: ${REF_FA}" >&2
    exit 1
fi
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

mkdir -p "$OUTBASE"

# Make unmatched file patterns expand to nothing instead of staying as literal text.
# Process every SAM file produced by the alignment step.
shopt -s nullglob

for SAM in "${MM2_DIR}"/*/aligned_*.sam; do
    echo "Processing $SAM"

    stem="${SUBTYPE}_$(basename "$(dirname "$SAM")")"
    outdir="${OUTBASE}/${stem}"
    mkdir -p "$outdir"

    # Convert SAM to BAM using the reference index, then sort and index the BAM file.
    samtools view -bt "${REF_FA}.fai" -o "${outdir}/aligned.bam" "$SAM"
    samtools sort -o "${outdir}/aligned.sorted.bam" "${outdir}/aligned.bam"
    rm -f "${outdir}/aligned.bam"
    samtools index "${outdir}/aligned.sorted.bam"

    # Save per-reference alignment counts and general BAM summary statistics.
    samtools idxstats "${outdir}/aligned.sorted.bam" > "${outdir}/idxstats.txt"
    samtools stats "${outdir}/aligned.sorted.bam" | grep ^SN \
        > "${outdir}/bam_summary_stats.txt"

    # Store the sorted BAM path in a shorter version
    BAM="${outdir}/aligned.sorted.bam"

    # Count reads:
    #   1) mapped to each reference sequence
    #   2) mapped to each reference sequence with mapping quality at least 20.
    samtools view -F 2304 "$BAM" \
        | awk '{c[$3]++} END{for(k in c) print k, c[k]}' \
        | sort -k2,2nr \
        > "${outdir}/per_chrom_all.txt"

    samtools view -F 2304 -q 20 "$BAM" \
        | awk '{c[$3]++} END{for(k in c) print k, c[k]}' \
        | sort -k2,2nr \
        > "${outdir}/per_chrom_mapq20.txt"
done

echo "DONE."

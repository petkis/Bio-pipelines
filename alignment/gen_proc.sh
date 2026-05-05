#!/bin/bash
set -euo pipefail

CHROM="$1"        # e.g. chr1_MATERNAL
CHR_LENGTH="$2"   # e.g. 109226264
SUBTYPE="$3"      # Masked | NoTelo | Telo

echo "=== proc_${SUBTYPE}.sh ==="
echo "CHROM       = $CHROM"
echo "CHR_LENGTH  = $CHR_LENGTH"
echo "SUBTYPE     = $SUBTYPE"

# ---------------------
# Reference selection
# ---------------------
case "$SUBTYPE" in
  Masked) REF_FA="genome_masked.fa" ;;
  NoTelo) REF_FA="genome_no_telomeres.fa" ;;
  Telo)   REF_FA="genome_telomeres.fa" ;;
  *) echo "ERROR: invalid SUBTYPE: $SUBTYPE"; exit 1 ;;
esac

MM2_DIR="${CHROM}_${SUBTYPE}_align"
OUTBASE="${CHROM}_${SUBTYPE}_post_alignments"

# reference index check
if [[ ! -f "${REF_FA}" ]]; then
    echo "ERROR: reference not found: ${REF_FA}" >&2
    exit 1
fi
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

mkdir -p "$OUTBASE"
shopt -s nullglob

for SAM in "${MM2_DIR}"/*/aligned_*.sam; do
    echo "Processing $SAM"

    stem="${SUBTYPE}_$(basename "$(dirname "$SAM")")"
    outdir="${OUTBASE}/${stem}"
    mkdir -p "$outdir"

    # BAM conversion + sorting
    samtools view -bt "${REF_FA}.fai" -o "${outdir}/aligned.bam" "$SAM"
    samtools sort -o "${outdir}/aligned.sorted.bam" "${outdir}/aligned.bam"
    rm -f "${outdir}/aligned.bam"
    samtools index "${outdir}/aligned.sorted.bam"

    # idxstats and summary
    samtools idxstats "${outdir}/aligned.sorted.bam" > "${outdir}/idxstats.txt"
    samtools stats "${outdir}/aligned.sorted.bam" | grep ^SN \
        > "${outdir}/bam_summary_stats.txt"

    # per-chrom read distribution
    BAM="${outdir}/aligned.sorted.bam"
    samtools view -F 2304 "$BAM" \
        | awk '{c[$3]++} END{for(k in c) print k, c[k]}' \
        | sort -k2,2nr \
        > "${outdir}/per_chrom_all.txt"

    samtools view -F 2304 -q 20 "$BAM" \
        | awk '{c[$3]++} END{for(k in c) print k, c[k]}' \
        | sort -k2,2nr \
        > "${outdir}/per_chrom_mapq20.txt"

    samtools view "$BAM" \
        | awk '{mq[$5]++} END{for(m in mq) print m, mq[m]}' \
        | sort -n \
        > "${outdir}/mapq_hist.txt"
done

echo "DONE."

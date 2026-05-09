#!/bin/bash
set -euo pipefail

# ======================
# INPUTS
# ======================
CHROM="$1"          # e.g. chr1_MATERNAL
READLEN="$2"        # e.g. 45000
SUBTYPE="$3"        # Masked | NoTelo | Telo
TOOL="${4:-minimap2}" # minimap2 | winnowmap
TOOL="${TOOL,,}"   # convert to lowercase
# ======================

if [[ "$TOOL" != "minimap2" && "$TOOL" != "winnowmap" ]]; then
    echo "Error: TOOL must be 'minimap2' or 'winnowmap'" >&2
    exit 1
fi

echo "=== run_${TOOL}_${SUBTYPE}.sh ==="
echo "CHROM     = $CHROM"
echo "READLEN   = $READLEN"
echo "SUBTYPE   = $SUBTYPE"

# reference selection
case "$SUBTYPE" in
  Masked) REF_FA="genome_masked.fa" ;;
  NoTelo) REF_FA="genome_no_telomeres.fa" ;;
  Telo)   REF_FA="genome_telomeres.fa" ;;
  *) echo "ERROR: invalid SUBTYPE: $SUBTYPE"; exit 1 ;;
esac

# repetitive k-mer list selection (for Winnowmap)
case "$SUBTYPE" in
  Masked) REP_K15="repetitive_k15_masked.txt" ;;
  NoTelo) REP_K15="repetitive_k15_notelo.txt" ;;
  Telo)   REP_K15="repetitive_k15_telo.txt" ;;
esac

REF_IDX="${REF_FA}.mmi"
THREADS=18
MINIMAP_PRESET="map-ont"

CHROMFA_DIR="chroms"
WGSIM_BIN="/storage/brno2/home/xpetkov/sw_wgsim/bin/wgsim"

# wgsim parameters
WGS_N=10000
WGS_LEN1="$READLEN"
WGS_LEN2=250
WGS_ERR=0
WGS_MUT=0
WGS_INDEL_FRAC=0
WGS_INDEL_EXT=0
WGS_SEED=42

# ---------------------
# CHECKS
# ---------------------
if [[ ! -x "$WGSIM_BIN" ]]; then
  echo "ERROR: wgsim not executable: $WGSIM_BIN" >&2
  exit 1
fi

if [[ ! -f "$REF_FA" ]]; then
  echo "ERROR: reference $REF_FA missing" >&2
  exit 1
fi

# INDEX REFERENCE
echo "Indexing reference..."
if [[ "$TOOL" == "winnowmap" ]]; then
  "$TOOL" -W "$REP_K15" -d "$REF_IDX" "$REF_FA"
else
  "$TOOL" -d "$REF_IDX" "$REF_FA"
fi

mkdir -p "${CHROM}_${SUBTYPE}_align"
shopt -s nullglob

#Loop
for path in "${CHROMFA_DIR}"/*.fa; do
    base=$(basename "$path" .fa)
    chrom=${base%_*}      # left part before last '_'
    subtype=${base##*_}   # right part after last '_'

    infile="${CHROMFA_DIR}/${base}.fa"
    if [[ ! -f "$infile" ]]; then
        echo "WARNING: missing $infile => skipping"
        continue
    fi

    OUT_ROOT="${CHROM}_${SUBTYPE}_align"
    outdir="${OUT_ROOT}/${chrom}_${SUBTYPE}_to_${subtype}_N_${WGS_N}_L_${READLEN}"
    mkdir -p "$outdir"

    fq="${outdir}/${chrom}_${subtype}_R1.fq"
    sam="${outdir}/aligned_${subtype}.sam"

    echo "Simulating $chrom ($subtype) -> $outdir"
    "$WGSIM_BIN" \
      -N "$WGS_N" -1 "$WGS_LEN1" -2 "$WGS_LEN2" \
      -e "$WGS_ERR" -r "$WGS_MUT" \
      -R "$WGS_INDEL_FRAC" -X "$WGS_INDEL_EXT" \
      -s "$WGS_SEED" \
      "$infile" "$fq" /dev/null

    echo "Aligning to $sam"
    if [[ "$TOOL" == "winnowmap" ]]; then
        "$TOOL" -t "$THREADS" -a -x "$MINIMAP_PRESET" \
            -W "$REP_K15" \
            --split-prefix "${chrom}_${subtype}_split" \
            "$REF_FA" "$fq" > "$sam"
    else
        "$TOOL" -t "$THREADS" -a -x "$MINIMAP_PRESET" \
            --split-prefix "${chrom}_${subtype}_split" \
            "$REF_IDX" "$fq" > "$sam"
    fi
done

echo "DONE."

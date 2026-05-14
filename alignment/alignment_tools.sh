#!/bin/bash

# Simulate reads from extracted chromosome-arm FASTA files and align them
# to one selected reference genome subtype using minimap2 or winnowmap.

set -euo pipefail

# ======================
# INPUTS (arguments passed from the main pipeline script)
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

# Select the reference genome according to the requested subtype.
case "$SUBTYPE" in
  Masked) REF_FA="genome_masked.fa" ;;
  NoTelo) REF_FA="genome_no_telomeres.fa" ;;
  Telo)   REF_FA="genome_telomeres.fa" ;;
  *) echo "ERROR: invalid SUBTYPE: $SUBTYPE"; exit 1 ;;
esac

# Select the repetitive k-mer list required by Winnowmap.
case "$SUBTYPE" in
  Masked) REP_K15="repetitive_k15_masked.txt" ;;
  NoTelo) REP_K15="repetitive_k15_notelo.txt" ;;
  Telo)   REP_K15="repetitive_k15_telo.txt" ;;
esac

# Define alignment settings and input/output locations.
REF_IDX="${REF_FA}.mmi"
THREADS=18
MINIMAP_PRESET="map-ont"

CHROMFA_DIR="chroms"
WGSIM_BIN="wgsim"

# wgsim parameters
WGS_N=10000
WGS_LEN1="$READLEN"
WGS_LEN2=250
WGS_ERR=0
WGS_MUT=0
WGS_INDEL_FRAC=0
WGS_INDEL_EXT=0
WGS_SEED=42

# Check that wgsim exists and is executable.
if [[ ! -x "$(command -v "$WGSIM_BIN")" ]]; then
  echo "ERROR: wgsim not found or not executable in PATH" >&2
  exit 1
fi

# Check that the selected reference FASTA exists.
if [[ ! -f "$REF_FA" ]]; then
  echo "ERROR: reference $REF_FA missing" >&2
  exit 1
fi

# Build an index for the selected reference.
# Winnowmap uses a repetitive k-mer file during indexing.
echo "Indexing reference..."
if [[ "$TOOL" == "winnowmap" ]]; then
  "$TOOL" -W "$REP_K15" -d "$REF_IDX" "$REF_FA"
else
  "$TOOL" -d "$REF_IDX" "$REF_FA"
fi

mkdir -p "${CHROM}_${SUBTYPE}_align"

# Make unmatched *.fa patterns expand to nothing instead of staying literal.
shopt -s nullglob

# Process each extracted chromosome-arm FASTA file from the chroms directory.
for path in "${CHROMFA_DIR}"/*.fa; do
    # Get the filename without path and extension.
    base=$(basename "$path" .fa)

    # Split the filename into chromosome name and subtype.
    # Example: chr13_PATERNAL_Telo -> chrom=chr13_PATERNAL, subtype=Telo.
    chrom=${base%_*}
    subtype=${base##*_}

    infile="${CHROMFA_DIR}/${base}.fa"

    # Skip this entry if the expected FASTA file is missing.
    if [[ ! -f "$infile" ]]; then
        echo "WARNING: missing $infile => skipping"
        continue
    fi

    OUT_ROOT="${CHROM}_${SUBTYPE}_align"
    outdir="${OUT_ROOT}/${chrom}_${SUBTYPE}_to_${subtype}_N_${WGS_N}_L_${READLEN}"
    mkdir -p "$outdir"

    fq="${outdir}/${chrom}_${subtype}_R1.fq"
    sam="${outdir}/aligned_${subtype}.sam"

    # Simulate reads from the extracted chromosome-arm FASTA file.
    # Only the first read file is kept; the second read file is discarded to /dev/null.
    echo "Simulating $chrom ($subtype) -> $outdir"
    "$WGSIM_BIN" \
      -N "$WGS_N" -1 "$WGS_LEN1" -2 "$WGS_LEN2" \
      -e "$WGS_ERR" -r "$WGS_MUT" \
      -R "$WGS_INDEL_FRAC" -X "$WGS_INDEL_EXT" \
      -s "$WGS_SEED" \
      "$infile" "$fq" /dev/null

    # Align the simulated reads to the selected reference genome.
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

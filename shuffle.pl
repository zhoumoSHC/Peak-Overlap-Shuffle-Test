#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# peak_overlap_shuffle_plot.sh
# Run empirical overlap test between two BED files,
# ------------------------------------------------------------
# Usage example:
#   ./shuffle.sh \
#       -a A.bed -b B.bed -g hg38.chrom.sizes \
#       -w 150 -n 1000 -o overlap_results.txt \
#       --seed 20251015 
# ------------------------------------------------------------

usage() {
  cat <<'USAGE'
Usage:
  peak_overlap_shuffle_plot.sh -a PEAKS_A -b PEAKS_B -g GENOME_SIZE [options]

Required:
  -a    BED file A
  -b    BED file B
  -g    Genome sizes file (e.g. hg38.chrom.sizes)

Optional:
  -w    Window size (bp, default: 150)
  -n    Iterations (default: 100)
  -o    Output file (default: overlap_null_counts.txt)
  -m    Shuffle mode: shuffle_a | shuffle_b | shuffle_both (default: shuffle_a)
  --seed INT        Random seed (default: random)
  --blacklist FILE  Exclude regions during shuffle
USAGE
}

# ---------- defaults ----------
WINDOW=150
ITERATIONS=100
OUTPUT_FILE="overlap_result.txt"
MODE="shuffle_a"
BLACKLIST=""
SEED=""
PLOT=false
TITLE="Distribution of Random Overlaps vs Observed"
OUTPREFIX=""

# ---------- parse args ----------
LONGOPTS=seed:,blacklist:,plot,title:,outprefix:
PARSED=$(getopt -o a:b:g:w:n:o:m:h --long "$LONGOPTS" -n "$0" -- "$@") || { usage; exit 1; }
eval set -- "$PARSED"

PEAKSA=""
PEAKSB=""
GENOME_SIZE=""

while true; do
  case "$1" in
    -a) PEAKSA="$2"; shift 2 ;;
    -b) PEAKSB="$2"; shift 2 ;;
    -g) GENOME_SIZE="$2"; shift 2 ;;
    -w) WINDOW="$2"; shift 2 ;;
    -n) ITERATIONS="$2"; shift 2 ;;
    -o) OUTPUT_FILE="$2"; shift 2 ;;
    -m) MODE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --plot) PLOT=true; shift ;;
    --title) TITLE="$2"; shift 2 ;;
    --outprefix) OUTPREFIX="$2"; shift 2 ;;
    -h) usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -s "${PEAKSA:-}" ]] || { echo "ERROR: -a BED file required."; exit 1; }
[[ -s "${PEAKSB:-}" ]] || { echo "ERROR: -b BED file required."; exit 1; }
[[ -s "${GENOME_SIZE:-}" ]] || { echo "ERROR: -g genome sizes file required."; exit 1; }

# ---------- random seed ----------
BASE_SEED="${SEED:-$((RANDOM + $$))}"

# ---------- temp dir ----------
TMPDIR=$(mktemp -d -t peakshuffle.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

# ---------- actual overlaps ----------
echo "[INFO] Calculating actual overlaps..."
actual_overlaps=$(bedtools window -w "$WINDOW" -a "$PEAKSA" -b "$PEAKSB" | wc -l | awk '{print $1}')
echo "[INFO] Actual overlaps = $actual_overlaps"

# ---------- randomization ----------
: > "$OUTPUT_FILE"
echo "[INFO] Running $ITERATIONS randomizations (mode=$MODE)..."
for ((i=1; i<=ITERATIONS; i++)); do
  s=$((BASE_SEED + i))
  SHUF_OPTS=(-g "$GENOME_SIZE" -seed "$s" -chrom)
  [[ -n "$BLACKLIST" ]] && SHUF_OPTS+=(-excl "$BLACKLIST")

  case "$MODE" in
    shuffle_a)
      bedtools shuffle -i "$PEAKSA" "${SHUF_OPTS[@]}" > "$TMPDIR/a.bed"
      bedtools window -w "$WINDOW" -a "$TMPDIR/a.bed" -b "$PEAKSB" | wc -l | awk '{print $1}' >> "$OUTPUT_FILE"
      ;;
    shuffle_b)
      bedtools shuffle -i "$PEAKSB" "${SHUF_OPTS[@]}" > "$TMPDIR/b.bed"
      bedtools window -w "$WINDOW" -a "$PEAKSA" -b "$TMPDIR/b.bed" | wc -l | awk '{print $1}' >> "$OUTPUT_FILE"
      ;;
    shuffle_both)
      bedtools shuffle -i "$PEAKSA" "${SHUF_OPTS[@]}" > "$TMPDIR/a.bed"
      bedtools shuffle -i "$PEAKSB" "${SHUF_OPTS[@]}" > "$TMPDIR/b.bed"
      bedtools window -w "$WINDOW" -a "$TMPDIR/a.bed" -b "$TMPDIR/b.bed" | wc -l | awk '{print $1}' >> "$OUTPUT_FILE"
      ;;
  esac
  (( i % 50 == 0 )) && echo "[INFO] Iteration $i/$ITERATIONS"
done

echo "[INFO] Randomization complete → $OUTPUT_FILE"

# ---------- compute p-value ----------
greater_equal=$(awk -v a="$actual_overlaps" '$1>=a{c++} END{print c+0}' "$OUTPUT_FILE")
p_value=$(awk -v c="$greater_equal" -v n="$ITERATIONS" 'BEGIN{printf("%.6g",(c+1)/(n+1))}')
mean=$(awk '{x+=$1} END{print x/NR}' "$OUTPUT_FILE")
sd=$(awk -v m="$mean" '{x+=($1-m)^2} END{print sqrt(x/NR)}' "$OUTPUT_FILE")
z=$(awk -v a="$actual_overlaps" -v m="$mean" -v s="$sd" 'BEGIN{if(s>0) printf("%.3f",(a-m)/s); else print "NA"}')

echo "----------- SUMMARY -----------"
echo "Actual overlaps : $actual_overlaps"
echo "Null mean (sd)  : $mean ($sd)"
echo "Z-score         : $z"
echo "P-value         : $p_value"
echo "Iterations      : $ITERATIONS"
echo "--------------------------------"


#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT_DIR="$ROOT_DIR/test/output"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/dna_result" \
	--alphabet dna \
	--motif-length 8 \
	--max-updates 1024 \
	--score-threshold 0.85 \
	--seed 1 \
	--threads 4 \
	> "$OUTPUT_DIR/dna_result.stdout.txt"

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/dna_result_explicit_one" \
	--alphabet dna \
	--motif-length 8 \
	--motif-count 1 \
	--max-updates 1024 \
	--score-threshold 0.85 \
	--seed 1 \
	--threads 4 \
	> "$OUTPUT_DIR/dna_result_explicit_one.stdout.txt"

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/dna_result_multiple" \
	--alphabet dna \
	--motif-length 8 \
	--motif-count 5 \
	--max-updates 1024 \
	--score-threshold 0.85 \
	--seed 1 \
	--threads 4 \
	> "$OUTPUT_DIR/dna_result_multiple.stdout.txt"

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/PROTEIN_TEST.fasta" \
	--output "$OUTPUT_DIR/protein_result" \
	--alphabet protein \
	--motif-length 6 \
	--max-updates 1024 \
	--score-threshold 0.80 \
	--seed 1 \
	--threads 4 \
	> "$OUTPUT_DIR/protein_result.stdout.txt"

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DUPLICATE_DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/duplicate_result" \
	--alphabet dna \
	--motif-length 8 \
	--motif-count 100 \
	--max-updates 1 \
	--score-threshold 0.80 \
	--seed 1 \
	--threads 4 \
	> "$OUTPUT_DIR/duplicate_result.stdout.txt"

grep -q "Consensus Subsequence    : ACGTGCAA" "$OUTPUT_DIR/dna_result.summary.txt"
grep -q "Consensus Subsequence    : MKLDPA" "$OUTPUT_DIR/protein_result.summary.txt"
test -s "$OUTPUT_DIR/dna_result.meme"
test -s "$OUTPUT_DIR/protein_result.meme"

for suffix in fasta offsets.tsv pwm.tsv meme; do
	diff -u "$OUTPUT_DIR/dna_result.$suffix" "$OUTPUT_DIR/dna_result_explicit_one.$suffix"
done
sed '/Elapsed Time/d' "$OUTPUT_DIR/dna_result.summary.txt" > "$OUTPUT_DIR/dna_result.summary.filtered.txt"
sed '/Elapsed Time/d' "$OUTPUT_DIR/dna_result_explicit_one.summary.txt" \
	> "$OUTPUT_DIR/dna_result_explicit_one.summary.filtered.txt"
diff -u "$OUTPUT_DIR/dna_result.summary.filtered.txt" \
	"$OUTPUT_DIR/dna_result_explicit_one.summary.filtered.txt"

grep '^\[STEP 3\] Pipeline ' "$OUTPUT_DIR/dna_result.stdout.txt" \
	> "$OUTPUT_DIR/dna_result.pipeline.txt"
grep '^\[STEP 3\] Pipeline ' "$OUTPUT_DIR/dna_result_multiple.stdout.txt" \
	> "$OUTPUT_DIR/dna_result_multiple.pipeline.txt"
diff -u "$OUTPUT_DIR/dna_result.pipeline.txt" "$OUTPUT_DIR/dna_result_multiple.pipeline.txt"

grep -q '^Requested Motif Number   : 5$' "$OUTPUT_DIR/dna_result_multiple.summary.txt"
grep -q '^Reported Motif Number    : 5$' "$OUTPUT_DIR/dna_result_multiple.summary.txt"
test "$(grep -c '^MOTIF MorbiusPlus_' "$OUTPUT_DIR/dna_result_multiple.meme")" -eq 5
test "$(grep -c '^>' "$OUTPUT_DIR/dna_result_multiple.fasta")" -eq 160
test "$(wc -l < "$OUTPUT_DIR/dna_result_multiple.offsets.tsv")" -eq 161
test "$(wc -l < "$OUTPUT_DIR/dna_result_multiple.pwm.tsv")" -eq 41

tail -n +2 "$OUTPUT_DIR/dna_result.offsets.tsv" > "$OUTPUT_DIR/dna_result.offsets.data.txt"
awk -F '\t' 'NR > 1 && $1 == 1 { print $3 "\t" $4 "\t" $5 "\t" $6 }' \
	"$OUTPUT_DIR/dna_result_multiple.offsets.tsv" > "$OUTPUT_DIR/dna_result_multiple.first.offsets.data.txt"
diff -u "$OUTPUT_DIR/dna_result.offsets.data.txt" \
	"$OUTPUT_DIR/dna_result_multiple.first.offsets.data.txt"
tail -n +2 "$OUTPUT_DIR/dna_result.pwm.tsv" > "$OUTPUT_DIR/dna_result.pwm.data.txt"
awk -F '\t' 'NR > 1 && $1 == 1 { print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' \
	"$OUTPUT_DIR/dna_result_multiple.pwm.tsv" > "$OUTPUT_DIR/dna_result_multiple.first.pwm.data.txt"
diff -u "$OUTPUT_DIR/dna_result.pwm.data.txt" "$OUTPUT_DIR/dna_result_multiple.first.pwm.data.txt"

grep -q '^Requested Motif Number   : 100$' "$OUTPUT_DIR/duplicate_result.summary.txt"
grep -q '^Reported Motif Number    : 1$' "$OUTPUT_DIR/duplicate_result.summary.txt"
test "$(grep -c '^MOTIF ' "$OUTPUT_DIR/duplicate_result.meme")" -eq 1
test "$(grep -c '^>' "$OUTPUT_DIR/duplicate_result.fasta")" -eq 4

if "$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/invalid_zero" \
	--alphabet dna \
	--motif-length 8 \
	--motif-count 0 \
	> "$OUTPUT_DIR/invalid_zero.stdout.txt" 2>&1; then
	printf "Zero motif count was accepted unexpectedly.\n"
	exit 1
fi
grep -q 'The motif count must be at least 1.' "$OUTPUT_DIR/invalid_zero.stdout.txt"

if "$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/DNA_TEST.fasta" \
	--output "$OUTPUT_DIR/invalid_overflow" \
	--alphabet dna \
	--motif-length 8 \
	--motif-count 18446744073709551616 \
	> "$OUTPUT_DIR/invalid_overflow.stdout.txt" 2>&1; then
	printf "Overflowing motif count was accepted unexpectedly.\n"
	exit 1
fi
grep -q 'Invalid value for --motif-count' "$OUTPUT_DIR/invalid_overflow.stdout.txt"

printf "All Morbius+ vanilla tests passed.\n"

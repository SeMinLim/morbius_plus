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
	--pipelines 4 \
	--max-updates 1024 \
	--score-threshold 0.85 \
	--seed 1 \
	--threads 4

"$ROOT_DIR/morbius_plus_vanila" \
	--input "$ROOT_DIR/test/PROTEIN_TEST.fasta" \
	--output "$OUTPUT_DIR/protein_result" \
	--alphabet protein \
	--motif-length 6 \
	--pipelines 4 \
	--max-updates 1024 \
	--score-threshold 0.80 \
	--seed 1 \
	--threads 4

grep -q "Consensus Subsequence    : ACGTGCAA" "$OUTPUT_DIR/dna_result.summary.txt"
grep -q "Consensus Subsequence    : MKLDPA" "$OUTPUT_DIR/protein_result.summary.txt"
test -s "$OUTPUT_DIR/dna_result.meme"
test -s "$OUTPUT_DIR/protein_result.meme"

printf "All Morbius+ vanilla tests passed.\n"

#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
rm -rf test/input test/output
mkdir -p test/input test/output

python3 - <<'PY'
from pathlib import Path

output = Path("test/input")
with (output / "UPSTREAM_TEST.fasta").open("w") as file:
    bases = "ACGT"
    for sequence_idx in range(40):
        sequence = "".join(bases[(position + sequence_idx) % 4] for position in range(1000))
        patch = ("GATTACA" * 20)[:120]
        patch_start = (sequence_idx * 17) % 800
        sequence = sequence[:patch_start] + patch + sequence[patch_start + len(patch):]
        file.write(f">background_{sequence_idx}\n{sequence}\n")

binding_sites = [
    "ATGGCTGCTGCTGCT", "ATGAAACCCGGGTTT", "GCTGACTTTCCCAAA", "CCCGGGAAATTTGCT",
    "ATGCGTGCTAAACCC", "GCTGCTGACGACGAC", "AAACCCGGGAAACCC", "TTTGCTGCTCCCGGG",
    "ATGACCGCTGACAAA", "GCTAAAGGGCCCTTT", "CCCGACTTTGCTAAA", "GGGAAACCCGCTGAC",
    "ATGGCTCCCAAAGGG", "GACGACGCTTTTCCC", "AAAGCTGCTGGGCCC", "TTTCCCAAAGCTGAC",
]
with (output / "BINDING_SITE_TEST.fasta").open("w") as file:
    for site_idx, sequence in enumerate(binding_sites):
        file.write(f">site_{site_idx}\n{sequence}\n")
PY

./synthethic_dataset_generator \
    --background test/input/UPSTREAM_TEST.fasta \
    --binding-sites test/input/BINDING_SITE_TEST.fasta \
    --output test/output \
    --seed 7 \
    --test-mode

for dataset_id in 1 2 3; do
    case "$dataset_id" in
        1) sequence_num=8 ;;
        2) sequence_num=16 ;;
        3) sequence_num=32 ;;
    esac

    dna_fasta="test/output/DATASET_DNA_${dataset_id}.fasta"
    protein_fasta="test/output/DATASET_PROTEIN_${dataset_id}.fasta"
    dna_truth="test/output/DATASET_DNA_${dataset_id}_ground_truth.tsv"
    protein_truth="test/output/DATASET_PROTEIN_${dataset_id}_ground_truth.tsv"

    test "$(grep -c '^>' "$dna_fasta")" -eq "$sequence_num"
    test "$(grep -c '^>' "$protein_fasta")" -eq "$sequence_num"
    test "$(($(wc -l < "$dna_truth") - 1))" -eq "$sequence_num"
    test "$(($(wc -l < "$protein_truth") - 1))" -eq "$sequence_num"

    awk 'BEGIN { good=1 } !/^>/ && length($0) != 1000 { good=0 } END { exit(good ? 0 : 1) }' "$dna_fasta"
    awk 'BEGIN { good=1 } !/^>/ && length($0) != 300 { good=0 } END { exit(good ? 0 : 1) }' "$protein_fasta"
done

printf "All synthetic dataset generator tests passed!\n"

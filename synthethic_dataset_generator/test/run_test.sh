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
        patch = ("GATTACATAATAGTGA" * 8)[:120]
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

python3 - <<'PY'
from pathlib import Path
import csv

CODON_TABLE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def read_fasta(filename):
    records = {}
    name = None
    sequence = []
    with filename.open() as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence)
                name = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence)
    return records


def translate(dna, source_offset, source_length):
    protein = []
    affected = []
    translated_offset = -1
    source_end = source_offset + source_length

    for dna_offset in range(0, len(dna) - 2, 3):
        amino_acid = CODON_TABLE[dna[dna_offset:dna_offset + 3]]
        if amino_acid == "*":
            continue

        protein_offset = len(protein)
        protein.append(amino_acid)
        if dna_offset < source_end and dna_offset + 3 > source_offset:
            if translated_offset < 0:
                translated_offset = protein_offset
            affected.append(amino_acid)

        if len(protein) == 300:
            break

    assert len(protein) == 300
    return "".join(protein), translated_offset, "".join(affected)


for dataset_id, sequence_num in ((1, 8), (2, 16), (3, 32)):
    output = Path("test/output")
    dna = read_fasta(output / f"DATASET_DNA_{dataset_id}.fasta")
    protein = read_fasta(output / f"DATASET_PROTEIN_{dataset_id}.fasta")

    assert len(dna) == sequence_num
    assert len(protein) == sequence_num
    assert all(len(sequence) == 1000 for sequence in dna.values())
    assert all(len(sequence) == 300 for sequence in protein.values())
    assert all(set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY") for sequence in protein.values())

    with (output / f"DATASET_DNA_{dataset_id}_ground_truth.tsv").open() as file:
        dna_truth = {row["sequence_id"]: row for row in csv.DictReader(file, delimiter="\t")}
    with (output / f"DATASET_PROTEIN_{dataset_id}_ground_truth.tsv").open() as file:
        protein_truth = {row["sequence_id"]: row for row in csv.DictReader(file, delimiter="\t")}

    assert len(dna_truth) == sequence_num
    assert len(protein_truth) == sequence_num

    for sequence_id in dna:
        row = dna_truth[sequence_id]
        source_offset = int(row["implanted_offset"])
        source_length = int(row["implanted_length"])
        expected, translated_offset, translated_sequence = translate(
            dna[sequence_id],
            source_offset,
            source_length,
        )
        assert protein[sequence_id] == expected

        protein_row = protein_truth[sequence_id]
        assert int(protein_row["source_dna_offset"]) == source_offset
        assert int(protein_row["source_dna_length"]) == source_length
        assert int(protein_row["translated_offset"]) == translated_offset
        assert int(protein_row["translated_length"]) == len(translated_sequence)
        assert protein_row["translated_sequence"] == translated_sequence

print("All synthetic dataset generator tests passed!")
PY

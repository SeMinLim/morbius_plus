# Synthetic Dataset Generator

This directory generates the six synthetic datasets used by Morbius+:

- `DATASET_DNA_1.fasta`, `DATASET_DNA_2.fasta`, and `DATASET_DNA_3.fasta`
- `DATASET_PROTEIN_1.fasta`, `DATASET_PROTEIN_2.fasta`, and `DATASET_PROTEIN_3.fasta`
- One exact ground-truth TSV file for each generated FASTA file

The production configuration follows the dataset dimensions described in the Original Morbius paper:

| Dataset | Sequence Number | Sequence Length | Binding-Site Pool |
|---|---:|---:|---:|
| DNA1 | 32,768 | 1,000 bp | 1,024 |
| DNA2 | 65,536 | 1,000 bp | 2,048 |
| DNA3 | 131,072 | 1,000 bp | 4,096 |
| Protein1 | 32,768 | 300 aa | Derived from DNA1 |
| Protein2 | 65,536 | 300 aa | Derived from DNA2 |
| Protein3 | 131,072 | 300 aa | Derived from DNA3 |

The generator splits `upstream5000.fa` into non-overlapping 1,000-base sequences, samples background sequences without replacement, selects the dataset-specific MA0007.2 binding-site pool without replacement, and implants one randomly selected binding site into every DNA sequence.

Each protein sequence is then generated directly from the corresponding synthesized DNA sequence. The generator reads codons from DNA offset 0, translates them with the standard DNA codon table, skips the three stop codons, and retains the first 300 translated amino acids. This preserves the 20-amino-acid alphabet while following the Original Morbius paper's DNA-derived protein dataset construction.

The regenerated datasets are deterministic for a given random seed. They are not intended to reproduce the original FASTA files byte-for-byte because the original random state and implantation positions were not archived.

## Input Files

- `upstream5000.fa`: upstream DNA background sequences
- `MA0007.2_sites.fasta`: MA0007.2 binding-site sequences

All binding-site sequences must have the same length.

## Build

```bash
make
```

## Generate Production Datasets

```bash
./synthethic_dataset_generator \
    --background /path/to/upstream5000.fa \
    --binding-sites /path/to/MA0007.2_sites.fasta \
    --output /path/to/output \
    --seed 1
```

## Ground-Truth Files

Each DNA ground-truth file contains:

```text
sequence_id  background_id  binding_site_id  implanted_offset  implanted_length  implanted_sequence
```

Each protein ground-truth file records the DNA source interval and the translated amino acids affected by that interval:

```text
sequence_id  binding_site_id  source_dna_offset  source_dna_length  translated_offset  translated_length  translated_sequence
```

All offsets are zero-based.

## Test

```bash
make test
```

The test generates small DNA and protein datasets and verifies that every protein sequence is the standard-codon translation of its corresponding synthesized DNA sequence.

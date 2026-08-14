# Synthetic Dataset Generator

This directory generates the six synthetic datasets used by Morbius+:

- `DATASET_DNA_1.fasta`, `DATASET_DNA_2.fasta`, and `DATASET_DNA_3.fasta`
- `DATASET_PROTEIN_1.fasta`, `DATASET_PROTEIN_2.fasta`, and `DATASET_PROTEIN_3.fasta`
- One exact ground-truth TSV file for each generated FASTA file

The production configuration follows the dataset dimensions described for Original Morbius:

| Dataset | Sequence Number | Sequence Length | Binding-Site Pool |
|---|---:|---:|---:|
| DNA1 | 32,768 | 1,000 bp | 1,024 |
| DNA2 | 65,536 | 1,000 bp | 2,048 |
| DNA3 | 131,072 | 1,000 bp | 4,096 |
| Protein1 | 32,768 | 300 aa | 1,024 translated sites |
| Protein2 | 65,536 | 300 aa | 2,048 translated sites |
| Protein3 | 131,072 | 300 aa | 4,096 translated sites |

The generator splits the background FASTA into non-overlapping 1,000-base sequences, samples background sequences without replacement, selects the dataset-specific binding-site pool without replacement, and implants one randomly selected site into every output sequence. Protein backgrounds are generated uniformly from the standard 20-amino-acid alphabet. Protein motifs are obtained by translating complete codons from the selected DNA binding sites. Because the source sites are non-coding DNA while the Morbius+ protein alphabet contains only 20 amino acids, stop codons are deterministically converted by the nearest one-base non-stop substitution: TAA/TAG to Q and TGA to W.

The regenerated datasets are deterministic for a given random seed. They are not intended to reproduce the original FASTA files byte-for-byte because the original random state and implantation positions were not archived.

## Input Files

The production generator expects the source files used by Original Morbius:

- `upstream5000.fa`: upstream DNA background sequences
- A FASTA file containing the MA0007.2 binding-site sequences

All binding-site sequences must have the same length. The Original Morbius source contained 11,206 binding-site sequences of 115 bp.

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

Each protein ground-truth file contains:

```text
sequence_id  binding_site_id  implanted_offset  implanted_length  implanted_sequence
```

Offsets are zero-based and refer to the generated FASTA sequence.

## Test

```bash
make test
```

The test executes the same generator with small dataset sizes of 8, 16, and 32 sequences.

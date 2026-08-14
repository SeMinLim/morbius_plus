# Morbius_Plus_Vanila

`Morbius_Plus_Vanila` is the pure C++ implementation of Morbius+. It supports forward-strand DNA motif discovery and protein motif discovery without FPGA hardware.

## Implemented Morbius+ Components

- Fixed-Rate Support-Guided Initialization with Markov-Adjusted Seed Ranking
  - DNA: 8-mer, third-order Markov model, Shannon-entropy filtering
  - Protein: 3-mer, first-order Markov model
  - `SampleNum = min(256, SeqNum)`
  - one selected seed, up to eight anchor candidates, exact matching, and fixed guidance rate 0.5
- Persistent Base-Pair Matrix (BPM) with pseudocount 1
- Log Probability Matrix (LPM)
- Log-domain candidate probability calculation
- Segment-based hierarchical inverse-CDF sampling
- xoshiro128+ random generation
- Overall Consensus Agreement Score
- Pipeline-local best-score management
- Independent pipeline termination and final Max Filter
- FASTA, offsets, PWM, MEME, and summary outputs

## Build

```bash
make
```

## Run

```bash
./morbius_plus_vanila \
    --input <dataset.fasta> \
    --output <result_prefix> \
    --alphabet dna \
    --motif-length 16 \
    --pipelines 16 \
    --score-threshold 0.90 \
    --seed 1
```

Protein input uses:

```bash
--alphabet protein
```

All input sequences must have the same length. DNA accepts `A`, `C`, `G`, and `T`. Protein accepts the 20 standard amino-acid symbols `ACDEFGHIKLMNPQRSTVWY`.

## Main Options

- `--pipelines <N>`: number of independent Gibbs pipelines
- `--max-updates <N>`: maximum sequence updates per pipeline
- `--score-threshold <F>`: normalized Overall Consensus Agreement Score threshold in `[0, 1]`
- `--seed <N>`: random seed
- `--threads <N>`: number of concurrent CPU threads

## Outputs

For output prefix `result/morbius_plus`, the program generates:

- `result/morbius_plus.fasta`
- `result/morbius_plus.offsets.tsv`
- `result/morbius_plus.pwm.tsv`
- `result/morbius_plus.meme`
- `result/morbius_plus.summary.txt`

## Test

```bash
make test
```

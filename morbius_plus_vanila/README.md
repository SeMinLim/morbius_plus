# Morbius_Plus_Vanila

`Morbius_Plus_Vanila` is the pure C++ implementation of Morbius+. It supports forward-strand DNA motif discovery and protein motif discovery without FPGA hardware.

## Implemented Morbius+ Components

- Fixed-Rate Support-Guided Initialization with Markov-Adjusted Seed Ranking
  - DNA: 8-mer, third-order Markov model, Shannon-entropy filtering
  - Protein: 3-mer, first-order Markov model
  - `SampleNum = min(256, SeqNum)`
  - top 16 distinct seeds, one independently selected anchor per seed from up to eight candidates, exact matching, and fixed guidance rate 0.5
- Persistent Base-Pair Matrix (BPM) with pseudocount 1
- Log Probability Matrix (LPM)
- Log-domain candidate probability calculation
- Segment-based hierarchical inverse-CDF sampling
- xoshiro128+ random generation
- Overall Consensus Agreement Score
- Pipeline-local best-score management
- Independent pipeline termination and score-ranked unique motif output
- FASTA, offsets, PWM, MEME, and summary outputs

Sequence-level seed support and Markov statistics are computed once and reused by all 16 pipelines. Eligible observed seeds retain the existing order: Markov-adjusted rank descending, support descending, then seed code ascending. Pipeline 0 receives the highest-ranked seed, pipeline 1 receives the next seed, and so on. Each seed uses the existing anchor evaluation independently; equal anchor numbers across pipelines are allowed. If fewer than 16 eligible seeds exist, the remaining pipelines use the existing random-only fallback. A selected seed with no legal sampled anchor also uses that fallback; seeds are never reused to fill pipeline slots.

Only seed assignment is diversified. The guidance rate, per-pipeline random-seed derivation, Gibbs updates, threshold and termination rules, and output ranking/deduplication are unchanged. Distinct seed words do not guarantee distinct initial offset arrays or final PWMs.

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
    --motif-count 5 \
    --score-threshold 0.80 \
    --seed 1
```

Protein input uses:

```bash
--alphabet protein
```

All input sequences must have the same length. DNA accepts `A`, `C`, `G`, and `T`. Protein accepts the 20 standard amino-acid symbols `ACDEFGHIKLMNPQRSTVWY`.

## Main Options

- `--motif-count <N>`: number of top unique motifs to output (default: `1`)
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
- `result/morbius_plus.seeds.tsv`
- `result/morbius_plus.initial_offsets.tsv`

With `--motif-count 1`, the original single-motif output formats and `MorbiusPlus` MEME ID are preserved. For larger values, candidates already produced by the 16 pipelines are ranked by score and exact duplicate PWMs are removed. The combined files identify motifs by 1-based rank, and the MEME IDs are `MorbiusPlus_1`, `MorbiusPlus_2`, and so on. If fewer unique candidates are available than requested, only the available motifs are written and the requested and reported counts are recorded in the summary.

The seed fields in the summary and console describe the first reported motif's pipeline, identified by `Seed Pipeline`. The `seeds.tsv` file records all 16 assignments, support/rank values, and anchors; `SeedValid=No` indicates random-only initialization and `NA` indicates an unavailable seed or anchor. The `initial_offsets.tsv` file records the actual offset arrays captured before the first Gibbs update, with one sequence per row and one pipeline per column. Pipeline indices, sequence indices, and offsets are zero-based. These diagnostics include all pipelines regardless of the requested motif count and allow offset-array equality to be checked without enforcing diversity or drawing additional random numbers. Capturing these arrays adds one offset copy per pipeline and diagnostic output I/O.

## Test

```bash
make test
```

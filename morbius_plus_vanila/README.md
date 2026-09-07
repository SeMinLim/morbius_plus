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
    --score-threshold 0.90 \
    --seed 1
```

Protein input uses:

```bash
--alphabet protein
```

All input sequences must have the same length. DNA accepts `A`, `C`, `G`, and `T`. Protein accepts the 20 standard amino-acid symbols `ACDEFGHIKLMNPQRSTVWY`.

## Main Options

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

## Synthetic DNA Motif Quality Evaluation

We evaluated the top-ranked motif from STREME, ProSampler, Morbius, and Morbius+ on three synthetic DNA datasets: DNA1 (32,768 sequences), DNA2 (65,536), and DNA3 (131,072). Every sequence is 1,000 bp long and contains one 115-bp implanted region. The implanted MA0007.2 core is 15 bp long. Because `implanted_offset` in each ground-truth TSV is zero-based and the core starts 50 bp into the implanted region, the true site is `[implanted_offset + 50, implanted_offset + 65)`. All 229,376 implanted regions were checked against the corresponding FASTA sequences before evaluation.

The two evaluation metrics are:

- **MA0007.2 Tomtom q-value (lower is better):** Tomtom 5.5.9 compared each discovered motif with all 205 motifs in JASPAR CORE 2014 vertebrates using Euclidean distance (`-dist ed`). A run recovers the target only if MA0007.2 is the best database hit and has `q <= 0.05`. Tomtom was run with `-thresh 1` only to retain the MA0007.2 row when it was not significant.
- **Site F1 (higher is better):** all coordinates were converted to zero-based, half-open intervals. A predicted site is a true positive when it is in the correct sequence and overlaps at least 8 bp of the 15-bp true site. Maximum one-to-one matching was used, followed by micro-averaging: `2TP / (2TP + FP + FN)`.

Each tool was run once per dataset and only its top motif was evaluated. Seed 1 was used wherever the tool exposed or accepted a seed. The main settings were:

- STREME 5.5.9: forward-only A/C/G/T alphabet, `--order 2 --hofract 0.85 --w 15 --nmotifs 1 --seed 1`.
- ProSampler 1.0.0: `-b 3 -d 1 -m 1 -f 100 -k 8 -l 6 -r 2 -p 1 -t 8 -c 1 -z 1.96`, with the default `-s 1.80` and `-w 4.50`.
- Morbius: 15-bp DNA motif, fixed seed 1, and logging-only instrumentation to record the selected offsets without changing its sampling or scoring logic.
- Morbius+: `--alphabet dna --motif-length 15 --score-threshold 0.90 --seed 1 --threads 3` with 16 independent pipelines.

STREME's default configuration exceeded the 20-GiB evaluation environment on these full datasets. The reported STREME runs therefore used the same memory-bounded configuration on all three datasets: forward-only search with 15% of sequences used for training and 85% held out. These are not default-configuration STREME results. STREME and ProSampler can report multiple sites per sequence, whereas Morbius and Morbius+ report one selected site per sequence; Site F1 includes each tool's native site-calling behavior.

Before Tomtom analysis, the STREME PWM's custom alphabet declaration was converted to standard DNA without changing its matrix or background probabilities. This makes the query compatible with the JASPAR database and allows Tomtom to evaluate both orientations of each target motif.

| Dataset | Method | MA0007.2 Tomtom q-value ↓ | Site F1 ↑ |
|---|---|---:|---:|
| DNA1 | STREME | 1.000000 | 0.000000 |
| DNA1 | ProSampler | 0.985415 | 0.001120 |
| DNA1 | Morbius | 1.000000 | 0.016113 |
| DNA1 | Morbius+ | 0.997671 | 0.020996 |
| DNA2 | STREME | 1.000000 | 0.000000 |
| DNA2 | ProSampler | 0.995220 | 0.002716 |
| DNA2 | Morbius | 1.000000 | 0.016266 |
| DNA2 | Morbius+ | 0.997671 | 0.017563 |
| DNA3 | STREME | 1.000000 | 0.000000 |
| DNA3 | ProSampler | 0.995220 | 0.004976 |
| DNA3 | Morbius | 1.000000 | 0.015602 |
| DNA3 | Morbius+ | 0.997671 | 0.018700 |

None of the 12 runs recovered MA0007.2: the target was never the best Tomtom hit and every target q-value was greater than 0.05. The evaluation pipeline itself passed both positive controls: the ground-truth-derived PWMs recovered MA0007.2 as the best hit with `q <= 5.44526e-39`, and supplying the true sites produced Site F1 = 1.0.

The dominant STREME motif was `GCTGGGATTACAGGC`. It occurred 4,150, 8,095, and 16,660 times in the complete DNA1, DNA2, and DNA3 sequences, respectively, but only 31, 28, and 98 times inside the implanted 115-bp regions. This shows that a repeated background pattern, rather than the implanted MA0007.2 motif, dominated discovery. Morbius and Morbius+ Site F1 values were also close to the uniform-start reference of approximately `15 / 986 = 0.0152`. Consequently, these measurements do not support a positive motif-recovery claim for any method on the current synthetic datasets; the background construction or preprocessing must be corrected before using this experiment as positive evidence in a paper.

## Test

```bash
make test
```

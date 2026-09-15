# Morbius_Plus_Vanila

`Morbius_Plus_Vanila` is the pure C++ implementation of Morbius+. It supports both-strand DNA motif discovery and protein motif discovery without FPGA hardware.

## Implemented Morbius+ Components

- Fixed-Rate Support-Guided Initialization with Markov-Adjusted Seed Ranking
  - DNA: 8-mer, third-order Markov model, Shannon-entropy filtering
  - Protein: 3-mer, first-order Markov model
  - `SampleNum = min(256, SeqNum)`
  - one highest-ranked seed and one anchor selected from up to eight candidates, shared across all 16 pipelines, exact matching, and fixed guidance rate 0.5
- Persistent Base-Pair Matrix (BPM) with pseudocount 1
- Log Probability Matrix (LPM)
- Log-domain candidate probability calculation
- Joint sampling of DNA `(offset, strand)` candidates
- Segment-based hierarchical inverse-CDF sampling
- xoshiro128+ random generation
- Overall Consensus Agreement Score
- Pipeline-local best-score management
- Independent pipeline termination and score-ranked unique motif output
- FASTA, offsets, PWM, MEME, and summary outputs

Sequence-level seed support and Markov statistics are computed once. The original single-seed selection chooses the highest Markov-adjusted rank, then higher support, then lower seed code to break ties. Its anchor is evaluated once with the existing scoring rule, and legal guided offsets are collected once. This complete seed model is assigned to all 16 pipelines. If the selected seed is unavailable or has no legal sampled anchor, all pipelines retain their existing random-only fallback.

Every pipeline receives the same seed word and anchor, while its original pipeline-specific RNG still chooses its initial offsets and subsequent Gibbs samples. Initial offset arrays and final PWMs may therefore differ; sharing a seed does not clone a pipeline's random state or offsets.

For DNA, each Gibbs update evaluates the forward and reverse-complement window at every legal offset and samples one site from their joint probability distribution. Each original sequence still contributes exactly one site to the BPM. The selected strand is applied when adding and removing sites and is saved with the offsets in each pipeline's best state. Final site sequences and PWM counts use these selected orientations. Protein discovery remains forward-only.

Seed-ranking and anchor-scoring formulas are unchanged; the highest-ranked seed and its anchor are shared again. The offset-initialization algorithm is unchanged and uses that shared model with each pipeline's own RNG. All initial DNA strands are forward, without drawing additional random numbers. The guidance rate, per-pipeline random-seed derivation, motif length, pipeline count, update limit, threshold and termination rules, and output ranking/deduplication are unchanged. DNA candidate counts increase from `L-W+1` to `2(L-W+1)` for sequence length `L` and motif length `W`; evaluating the extra candidates increases computation, so elapsed-time effects require measurement.

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
- `--control <FASTA>`: enable the DNA refinement described below

## Optional DNA refinement

```bash
./morbius_plus_vanila \
    --input primary.fasta --control control.fasta \
    --output result/refined --alphabet dna \
    --motif-length 10 --motif-count 16 \
    --score-threshold 0.80 --seed 1
```

This is a separate, STREME-inspired refinement of the completed Gibbs candidates.
It does not add an absent-site state to the Gibbs sampler. Enabling refinement
does not change the shared seed and anchor, initial offsets, per-pipeline RNG
streams, joint strand sampling, agreement score or termination rules. The default Gibbs agreement threshold is
still **0.80**. Without `--control`, the existing DNA and protein behavior is
preserved. Control-based refinement currently requires DNA and the same fixed
sequence length in Primary and Control; their sequence counts may differ.

The following operations are applied to every pipeline candidate before output
ranking, deduplication or the `--motif-count` limit:

1. Build one second-order Markov background from all Control sequences and reuse
   it across candidates. Tuple counts include both orientations; the prior for
   each length-`k` tuple is `1 / 4^k`. A window uses only its own preceding bases
   as context, with order 0 and 1 at its first two positions.
2. Scan every Primary and Control sequence on both strands for its best
   `log2(P(window | PWM) / P(window | background))` site. Ties prefer the forward
   strand, then the lowest original offset. Each original sequence contributes
   at most one hit to the enrichment table.
3. Search positive site-score cutoffs, admitting all Primary and Control scores
   tied at a cutoff together. A passing site has `score >= cutoff`. Among
   Primary-enriched cutoffs, choose the smallest one-sided Fisher exact p-value;
   ties retain the higher cutoff. This site-score cutoff is independent of the
   Gibbs agreement threshold. No reference motifs are used.
4. Form a proposed PWM using only the passing Primary sites, with counts in their
   selected orientations. Control sites do not enter the PWM. For `M` fitting
   sites, each DNA probability is `(count + 1) / (M + 4)`.
5. Rescan the proposed PWM and optimize its enrichment cutoff again. Accept it
   only if its enrichment p-value strictly improves the current best PWM. Rescan
   all sequences on each iteration, allowing previously excluded sequences to
   return. Stop on non-improvement, lack of support, or 20 proposals.

Refinement candidates run independently within the existing `--threads` budget
(default: 16, capped at 16 and the candidate count). Workers take the next
unprocessed candidate when ready, but results remain in pipeline-index order.
Gibbs sampling finishes before refinement workers start.

Before processing candidates, the program caches each Primary and Control
window's background log-probability on both strands and builds the Fisher
log-factorial table once. All candidates and iterations share these read-only
tables. Background scores remain `double`, Fisher arithmetic remains
`long double`, and the original accumulation and score-tie rules are unchanged.
The background cache uses two doubles per legal offset across Primary and
Control, plus sequence-start indices; it is allocated once, not per worker.
Cache construction is included in the existing refinement timer.

A proposed PWM that is exactly equal to the current PWM skips the redundant
rescan and stops with the existing `no_improvement` reason and proposal count.
There is no approximate equality threshold or candidate pruning. All 16
candidates remain eligible, and the maximum remains **20 proposals**.

The accepted PWM and the sites that actually constructed it are saved together.
Its subsequently evaluated best sites can differ from those fitting sites, so
the two cutoffs, hit counts and p-values are recorded separately. If no proposal
improves a supported initial Gibbs PWM, that original fit is retained and marked
as such. A candidate without an enriched positive-score cutoff is not reported;
if none of the 16 candidates has support, output files contain no motif and the
summary reports zero. No all-sequence fallback is used for unsupported candidates.

Final ordering continues to use the original Gibbs best score and pipeline-index
tie break. Exact duplicate **smoothed PWMs** are removed, comparing denominators
as well as counts when fitting site numbers differ. The Gibbs score is not
renormalized using the smaller fitting set. FASTA and offsets output include
only fitting sites, and PWM/MEME denominators and MEME `nsites` use their actual
number. Original sequence coordinates and selected strand are retained.

Additional files record all candidates, including ones removed by the final
output limit or deduplication:

- `.refinement.tsv`: evaluated and fitting statistics, iteration counts and stop reasons
- `.refinement_sites.tsv`: per-candidate, per-Primary presence flags and fitting coordinates
- `.refinement_candidates.meme`: all supported final candidate PWMs
- `.refinement_scoring.meme`: preceding PWMs used to select the accepted fitting sites
- `.refinement_background.tsv`: Control-derived conditional background probabilities
- `.refinement_summary.txt`: refinement settings, timing and interpretation

The retained initial Gibbs fit has no preceding refinement fitting PWM. The
metadata distinguishes that case from accepted refinement. The original
`.seeds.tsv`, `.initial_offsets.tsv` and summary's `[All Pipelines]` table continue
to describe Gibbs alone.

All MEME exports retain the same uniform letter-frequency header as the existing
evaluation output, so default Tomtom scoring uses the same convention for a motif
in `.meme` and `.refinement_candidates.meme`. Refinement itself uses the full
Control Markov model recorded in `.refinement_background.tsv`.

The method follows the sequence-level enrichment principle of
[STREME refinement](https://meme-suite.org/meme/doc/streme.html), while retaining
Morbius+'s pseudocount 1, pipeline candidates and output ranking; it is not a
reimplementation of all STREME estimation and search procedures. Refinement adds
background construction, repeated scanning and output work. With `--control`,
the reported elapsed time includes input reading, Gibbs, refinement and result
file writes up to the final summary timestamp. For a complete end-to-end benchmark,
measure the whole process externally, including its final summary/console writes;
Tomtom evaluation is a separate process and is excluded.

## Outputs

For output prefix `result/morbius_plus`, the program generates:

- `result/morbius_plus.fasta`
- `result/morbius_plus.offsets.tsv`
- `result/morbius_plus.pwm.tsv`
- `result/morbius_plus.meme`
- `result/morbius_plus.summary.txt`
- `result/morbius_plus.seeds.tsv`
- `result/morbius_plus.initial_offsets.tsv`

With `--motif-count 1`, the MEME ID is `MorbiusPlus`. For larger values, candidates already produced by the 16 pipelines are ranked by score and exact duplicate PWMs are removed. The combined files identify motifs by 1-based rank, and the MEME IDs are `MorbiusPlus_1`, `MorbiusPlus_2`, and so on. If fewer unique candidates are available than requested, only the available motifs are written and the requested and reported counts are recorded in the summary.

DNA FASTA headers include `strand=+` or `strand=-`, and DNA `offsets.tsv` includes a `Strand` column after `Offset`. Both report the zero-based leftmost coordinate of the window in the original input sequence, even for reverse-complement sites. FASTA sequences and the `Motif` column are oriented to the selected strand. DNA MEME files declare `strands: + -`; PWM and MEME matrices count each original sequence once. Protein site and matrix output formats are unchanged.

The DNA summary includes an `[All Pipelines]` table with each pipeline's best score, update count, threshold-reached flag, and termination reason (`threshold` or `max_updates`), including pipelines excluded from the final motif output. If the threshold is reached on the final allowed update, the reason is `threshold`. The summary also records the joint strand search, all-forward initial strands, original offset coordinates, and configured update limit.

The seed fields in the summary and console describe the first reported motif's pipeline, identified by `Seed Pipeline`. The `seeds.tsv` file records all 16 assignments, support/rank values, and anchors; `SeedValid=No` indicates random-only initialization and `NA` indicates an unavailable seed or anchor. The `initial_offsets.tsv` file records the actual offset arrays captured before the first Gibbs update, with one sequence per row and one pipeline per column. Pipeline indices, sequence indices, and offsets are zero-based. These diagnostics include all pipelines regardless of the requested motif count and allow offset-array equality to be checked without enforcing diversity or drawing additional random numbers. Capturing these arrays adds one offset copy per pipeline and diagnostic output I/O.

## Test

Tests require Python 3 in addition to the C++ build tools. The refinement
optimization regression compares cached, parallel execution with the original
uncached serial procedure, including exact PWM counts, site masks, positions,
strands, enrichment values, and termination diagnostics.

```bash
make test
```

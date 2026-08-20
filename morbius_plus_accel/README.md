# Morbius+ CPU-FPGA Heterogeneous Prototype

`morbius_plus_accel` implements the finalized Morbius+ design for the Alveo U50 in the same `hw/` and `sw/` organization used by blueVitis.

- `hw/`: Bluespec Gibbs pipelines, blueVitis AXI/XRT wrapper, U50 configuration, packaging scripts, and standalone simulation.
- `sw/`: Markov-adjusted Seed Initializer, Morbius Orchestrator, XRT backend, bit-accurate kernel model, result writers, and tests.

The zeroth iteration transfers sequence 0 and pipeline-specific leave-one-out BPM/LPM states.
Subsequent iterations transfer one sequence and one pipeline-specific tentative-motif offset.
Each hardware pipeline performs column-masked BPM/LPM maintenance, two-phase pipelined profiling, probability-proportional hierarchical sampling, new tentative-motif insertion, and pipeline-local score calculation.

## Quick validation without Vitis

```bash
cd hw
make sim

cd ../sw
make test
```

## blueVitis/U50 build

```bash
cd hw
make all TARGET=hw_emu
make run TARGET=hw_emu MORBIUS_ARGS="--input /path/to/input.fa --output output/result --alphabet dna --motif-length 16 --pipelines 4 --max-updates 32768 --score-threshold 1.0 --seed 1"
```

Use `TARGET=hw` for the physical U50 after enabling XRT host memory.

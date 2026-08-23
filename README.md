# Morbius+

* **A CPU-FPGA heterogeneous system for scalable probability-proportional sequence motif discovery.**
* `Morbius+` combines Markov-adjusted seed initialization with independent Gibbs pipelines, persistent BPM/LPM state, and hierarchical probability-proportional sampling.

## File structure
```text
morbius_plus/
├── morbius_plus_vanila/          pure C++ reference implementation
├── morbius_plus_accel/
│   ├── kernel_morbius_plus/      Bluespec kernel and blueVitis/U50 build flow
│   └── host_morbius_plus/        host orchestrator, XRT backend, and kernel model
└── synthethic_dataset_generator/ synthetic DNA and protein dataset generator
```

## Software version
```bash
cd morbius_plus_vanila
make
make test
```

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

Use `--alphabet protein` for protein sequences.

## Alveo U50 prototype
`morbius_plus_accel` follows the [`blueVitis`](https://github.com/SeMinLim/bluevitis) hardware/software structure and targets AMD Vitis 2025.2, XRT, and Bluespec SystemVerilog.

```bash
cd morbius_plus_accel
make test
```

```bash
cd kernel_morbius_plus
make all TARGET=hw_emu
make run TARGET=hw_emu MORBIUS_ARGS="--input <dataset.fasta> --output <result_prefix> --alphabet dna --motif-length 16 --pipelines 4 --score-threshold 0.90 --seed 1"
```

Use `TARGET=hw` for the physical Alveo U50.

## Notes
* Maintained by Se-Min Lim.
* The default accelerator platform is `xilinx_u50_gen3x16_xdma_5_202210_1`.

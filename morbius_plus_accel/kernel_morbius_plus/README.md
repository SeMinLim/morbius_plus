# Morbius+ U50 Hardware

This directory is a blueVitis-compatible RTL kernel for the Alveo U50.
The kernel uses the direct `HOST[0]` memory connection from `kernel_example_add_host` and retains sixteen independent Gibbs pipeline states on chip.

## Prototype configuration

- `NumPipeline = 16`
- `NumPE_Profiler = 16`
- `NumPE_LPM = 4`
- physical dual-mode PWL lanes: `max(NumPE_Profiler, NumPE_LPM) = 16`
- maximum sequence length: 1024 symbols
- maximum motif length: 128 symbols
- maximum alphabet size: 20 symbols
- target clock: 250 MHz

## Standalone BSV test

```bash
make sim
```

## U50 build

```bash
make all TARGET=hw_emu
```

For the physical U50:

```bash
xbutil configure --device <BDF> --host-mem --size 1G enable
make all TARGET=hw
```

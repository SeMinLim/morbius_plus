# Morbius+ Prototype Validation

Validation date: 2026-08-20

## Host-Side Model and Protocol

The following command completed successfully:

```bash
make -C morbius_plus_accel/sw test
```

Validated items:

- 16-interval Log2 PWL approximation
- 16-interval Exp2 PWL approximation
- 512-bit accelerator command/result protocol
- DNA accelerator-model execution
- Protein accelerator-model execution

Observed numerical errors against the corresponding double-precision functions:

- Maximum Log2 PWL absolute error: `0.0010731059`
- Maximum Exp2 PWL absolute error: `0.0002312184`

Small functional-test consensus results:

- DNA: `ACGTGCAA`
- Protein: `MKLDPA`

The protocol and model tests completed without AddressSanitizer or UndefinedBehaviorSanitizer failures in an additional sanitizer build.

## Bluespec Pipeline Simulation

The following command completed successfully with the official B-Lang Bluespec compiler release 2026.01:

```bash
make -C morbius_plus_accel/hw sim
```

Observed standalone Gibbs-pipeline results:

```text
Bootstrap offset=12 score=36 update=1 bestUpdate=1 terminated=0 cycle=30
Update offset=4 score=36 update=2 bestUpdate=0 terminated=1 cycle=61
All GibbsPipeline tests passed.
```

The complete blueVitis `KernelTop.bsv` hierarchy was also compiled to synthesizable Verilog with BSC. The generated kernel contained no simulation-only `$finish` statement.

## DATASET_DNA_1.fa

The bit-accurate accelerator model was executed twice on the complete production-scale DNA1 dataset.

Configuration:

- Sequence number: `32768`
- Sequence length: `1000`
- Motif length: `16`
- Pipeline number: `4`
- Maximum updates per pipeline: `32768`
- Score threshold: `1.0`
- Random seed: `1`

Observed result:

```text
Selected Seed         : GTAATCCC
Anchor Offset         : 3
Pipeline 0 Best Score : 186525
Pipeline 1 Best Score : 186063
Pipeline 2 Best Score : 186839
Pipeline 3 Best Score : 185791
Selected Pipeline     : 2
Best Score            : 186839
Normalized Best Score : 0.35633659
Consensus Subsequence : AATTTAATAAAAAAAA
```

The two runs produced byte-identical offset, FASTA, PWM, and MEME outputs. The offset-file SHA-256 digest was:

```text
c7e9fa902de9cef6a9e376a1e8ee6e04b06aa16f364341f5d7978d02b46eba62
```

Post-run invariants were checked:

- all `32768` sequences produced one motif site;
- every offset was in `[0, 984]`;
- every reported motif exactly matched the source sequence at the reported offset;
- the score recomputed from the output sites equaled `186839`;
- the recomputed consensus equaled `AATTTAATAAAAAAAA`;
- every PWM row summed to one.

## Current Validation Boundary

The available validation environment did not contain AMD Vitis, XRT, or an Alveo U50. Therefore, the following have not yet been executed here:

- XRT host-binary compilation;
- `v++` kernel linking;
- `hw_emu` execution;
- physical U50 synthesis, timing closure, or board execution.

The repository includes the blueVitis-compatible XRT host source, AXI shell, U50 configuration, packaging scripts, and build targets required to perform those steps in the target blueVitis environment.

# Morbius+ Host-Kernel Protocol

All transfers use 512-bit little-endian beats. The 16-pipeline configuration uses protocol version 2.

## Command header

| Byte | Field |
|---:|---|
| 0-3 | magic `0x4d504c53` |
| 4-5 | protocol version 2 |
| 6 | command: 0 bootstrap, 1 batch |
| 7 | alphabet size |
| 8-9 | sequence length |
| 10 | motif length |
| 11 | active pipeline count |
| 16-19 | raw score threshold |
| 20-23 | maximum updates |
| 24-27 | batch size |
| 28-31 | 512-bit beats per sequence |

Bootstrap carries sequence 0 followed by sixteen physical pipeline states. Each state contains one RNG-state beat, one BPM-column beat per motif position, and one LPM-column beat per motif position. The host has already removed sequence 0's initial tentative motif from both matrices.

A subsequent batch carries repeated pairs of one encoded sequence and one 512-bit offset beat containing sixteen 32-bit tentative-motif offsets.

## Result group

Each processed sequence produces four 512-bit result beats. Byte 6 identifies the result-beat index from 0 to 3. Each beat contains four 96-bit pipeline records, for a total of sixteen records per sequence. A record contains the new tentative-motif offset, `bestUpdate`, termination and active flags, best score, and update count.

The fixed-position summary beat is written at `batchSize * 256` bytes and reports the processed count, global completion, selected pipeline, and kernel cycle count.

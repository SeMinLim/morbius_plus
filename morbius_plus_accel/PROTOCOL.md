# Morbius+ Host-Kernel Protocol

All transfers use 512-bit little-endian beats.

## Command header

| Byte | Field |
|---:|---|
| 0-3 | magic `0x4d504c53` |
| 4-5 | protocol version 1 |
| 6 | command: 0 bootstrap, 1 batch |
| 7 | alphabet size |
| 8-9 | sequence length |
| 10 | motif length |
| 11 | active pipeline count |
| 16-19 | raw score threshold |
| 20-23 | maximum updates |
| 24-27 | batch size |
| 28-31 | 512-bit beats per sequence |

Bootstrap carries sequence 0 followed by four physical pipeline states.
Each state contains one RNG-state beat, one BPM column beat per motif position, and one LPM column beat per motif position.
The host has already removed sequence 0's initial tentative motif from both matrices.

A subsequent batch carries repeated pairs of one encoded sequence and one beat containing four 32-bit tentative-motif offsets.

## Result beat

Each result beat contains four 96-bit pipeline records.
A record contains the new tentative-motif offset, `bestUpdate`, termination and active flags, best score, and update count.
The summary beat is always written at `batchSize * 64` bytes and reports the processed count, global completion, selected pipeline, and pipeline-local best states.

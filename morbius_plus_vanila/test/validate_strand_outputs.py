"""Reconstruct DNA output, including refined membership, from forward sites."""
import csv
import io
from pathlib import Path
import re
import sys


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            records[-1][1] += line.strip()
    return records


source = fasta(sys.argv[1])
prefix = sys.argv[2]
emitted = fasta(prefix + ".fasta")
with open(prefix + ".offsets.tsv") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == len(emitted)
sites = {}
sequence_indices = {}
strands = {}
complement = str.maketrans("ACGT", "TGCA")
for row, (header, site) in zip(rows, emitted):
    rank = int(row.get("MotifRank", 1))
    idx = int(row["SequenceIdx"])
    offset = int(row["Offset"])
    strand = row["Strand"]
    assert strand == "+"
    assert row["SequenceName"] == source[idx][0]
    assert 0 <= offset <= len(source[idx][1]) - len(site)
    expected = source[idx][1][offset:offset + len(site)]
    if strand == "-":
        expected = expected.translate(complement)[::-1]
    assert row["Motif"] == site == expected
    assert f"offset={offset} strand={strand}" in header
    sites.setdefault(rank, []).append(site)
    sequence_indices.setdefault(rank, []).append(idx)
    strands.setdefault(rank, set()).add(strand)

matrices = {}
refinement = Path(prefix + ".refinement.tsv")
refined_sites = {}
if refinement.exists():
    with refinement.open() as handle:
        selected = {int(row["OutputRank"]): int(row["Pipeline"])
                    for row in csv.DictReader(handle, delimiter="\t") if row["OutputRank"] != "unreported"}
    with open(prefix + ".refinement_sites.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["SitePresent"] == "1":
                refined_sites.setdefault(int(row["Pipeline"]), []).append(int(row["SequenceIdx"]))
    assert set(selected) == set(sites)
for rank, motif_sites in sites.items():
    assert len(set(sequence_indices[rank])) == len(sequence_indices[rank])
    expected_indices = refined_sites[selected[rank]] if refinement.exists() else list(range(len(source)))
    assert sorted(sequence_indices[rank]) == expected_indices
    matrices[rank] = [
        [(sum(site[column] == base for site in motif_sites) + 1) / (len(motif_sites) + 4)
         for base in "ACGT"]
        for column in range(len(motif_sites[0]))
    ]

with open(prefix + ".pwm.tsv") as handle:
    pwm_rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(pwm_rows) == sum(map(len, matrices.values()))
for row in pwm_rows:
    expected = matrices[int(row.get("MotifRank", 1))][int(row["Position"])]
    assert all(abs(float(row[base]) - value) < 1e-6 for base, value in zip("ACGT", expected))

meme = Path(prefix + ".meme").read_text()
assert "strands: +\n" in meme
blocks = re.findall(r"letter-probability matrix:.*\n((?:[0-9. eE+-]+\n)+)", meme)
assert len(blocks) == len(matrices)
for rank, block in enumerate(blocks, 1):
    actual = [[float(value) for value in line.split()] for line in block.strip().splitlines()]
    assert len(actual) == len(matrices[rank])
    for actual_row, expected in zip(actual, matrices[rank]):
        assert len(actual_row) == 4
        assert all(abs(a - b) < 1e-6 for a, b in zip(actual_row, expected))
assert list(map(int, re.findall(r"nsites= (\d+)", meme))) == [len(sites[rank]) for rank in sorted(sites)]

summary = Path(prefix + ".summary.txt").read_text()
table = summary.split("[All Pipelines]\n", 1)[1].split("\n\n", 1)[0]
pipelines = list(csv.DictReader(io.StringIO(table), delimiter="\t"))
assert [int(row["Pipeline"]) for row in pipelines] == list(range(16))
limit = int(re.search(r"Maximum Pipeline Updates\s+: (\d+)", summary).group(1))
for row in pipelines:
    reached = row["ThresholdReached"] == "Yes"
    assert row["TerminationReason"] == ("threshold" if reached else "max_updates")
    assert 0 <= int(row["Updates"]) <= limit
    if not reached:
        assert int(row["Updates"]) == limit

if "--exact-pair" in sys.argv:
    assert all(orientations == {"+"} for orientations in strands.values())
    for rank, motif_sites in sites.items():
        assert motif_sites == [source[index][1] for index in sequence_indices[rank]]
    assert all(row["ThresholdReached"] == "No" for row in pipelines)
    assert all(int(row["BestScore"]) == 14 for row in pipelines)
print("Validated forward FASTA, offsets, PWM/MEME, and all pipeline stop reasons:", prefix)

"""Run small independent integration checks for optional post-Gibbs refinement."""
import csv
import io
import math
from pathlib import Path
import random
import re
import subprocess
import sys


ALPHABET = "ACGT"
COMPLEMENT = str.maketrans(ALPHABET, "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line.strip():
            records[-1][1] += line.strip()
    return records


def write_fasta(path, sequences):
    Path(path).write_text("".join(f">sequence_{idx}\n{sequence}\n" for idx, sequence in enumerate(sequences)))


def tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def meme(path):
    result = {}
    text = Path(path).read_text()
    for block in re.split(r"^MOTIF ", text, flags=re.M)[1:]:
        motif_id = block.split()[0]
        header = re.search(r"letter-probability matrix:.*w=\s*(\d+).*nsites=\s*(\d+).*\n", block)
        assert header, (path, motif_id)
        width, site_num = map(int, header.groups())
        matrix = [list(map(float, line.split())) for line in block[header.end():].splitlines()[:width]]
        assert len(matrix) == width and all(len(row) == 4 for row in matrix)
        assert all(abs(sum(row) - 1) < 2e-6 for row in matrix)
        result[motif_id] = (site_num, matrix)
    return result


def background(control):
    """Independent concrete tuple counts, symmetric strands and one total prior."""
    tables = []
    for length in (1, 2, 3):
        counts = {}
        from itertools import product
        for word in map("".join, product(ALPHABET, repeat=length)):
            counts[word] = 1 / 4 ** length
        for sequence in control:
            for oriented in (sequence, reverse_complement(sequence)):
                for offset in range(len(sequence) - length + 1):
                    counts[oriented[offset:offset + length]] += 1
        tables.append({word: math.log2(count / sum(counts[word[:-1] + base] for base in ALPHABET))
                       for word, count in counts.items()})
    return tables


def site_score(site, matrix, bg):
    return sum(math.log2(matrix[column][ALPHABET.index(base)]) -
               bg[min(2, column)][site[max(0, column - 2):column + 1]]
               for column, base in enumerate(site))


def best_site(sequence, matrix, bg):
    width = len(matrix)
    best = (-math.inf, 0, "+")
    for strand in ("+", "-"):
        for offset in range(len(sequence) - width + 1):
            site = sequence[offset:offset + width]
            if strand == "-":
                site = reverse_complement(site)
            score = site_score(site, matrix, bg)
            if score > best[0]:
                best = score, offset, strand
    return best


def fisher(primary_hits, control_hits, primary_num, control_num):
    total_hits = primary_hits + control_hits
    return sum(math.comb(primary_num, count) * math.comb(control_num, total_hits - count)
               for count in range(primary_hits, min(primary_num, total_hits) + 1)
               if 0 <= total_hits - count <= control_num) / math.comb(primary_num + control_num, total_hits)


def check_enrichment(row, field_prefix, primary_scores, control_scores):
    threshold = float(row[field_prefix + "ScoreThreshold"])
    assert threshold > 0
    # Printed matrices can round a site lying exactly on the retained threshold.
    primary_hits = sum(score >= threshold - 2e-7 for score in primary_scores)
    control_hits = sum(score >= threshold - 2e-7 for score in control_scores)
    assert primary_hits == int(row[field_prefix + "PrimaryHits"]), row
    assert control_hits == int(row[field_prefix + "ControlHits"]), row
    reported = float(row["FitSelectionLogPvalue" if field_prefix == "Fit" else "EvaluationLogPvalue"])
    assert abs(math.log(fisher(primary_hits, control_hits, len(primary_scores), len(control_scores))) - reported) < 2e-8
    assert primary_hits / len(primary_scores) > control_hits / len(control_scores)


def all_pipeline_table(prefix):
    text = Path(str(prefix) + ".summary.txt").read_text()
    table = text.split("[All Pipelines]\n", 1)[1].split("\n\n", 1)[0]
    rows = list(csv.DictReader(io.StringIO(table), delimiter="\t"))
    assert [int(row["Pipeline"]) for row in rows] == list(range(16))
    return rows


def check_outputs(prefix, primary, control, expected_limit):
    prefix = str(prefix)
    diagnostics = tsv(prefix + ".refinement.tsv")
    assert [int(row["Pipeline"]) for row in diagnostics] == list(range(16))
    fit_rows = tsv(prefix + ".refinement_sites.tsv")
    assert len(fit_rows) == 16 * len(primary)
    candidates = meme(prefix + ".refinement_candidates.meme")
    scoring = meme(prefix + ".refinement_scoring.meme")
    for suffix in (".meme", ".refinement_candidates.meme", ".refinement_scoring.meme"):
        text = Path(prefix + suffix).read_text()
        frequencies = re.search(r"Background letter frequencies[^\n]*\n([^\n]+)", text).group(1).split()
        assert dict(zip(frequencies[::2], map(float, frequencies[1::2]))) == dict.fromkeys(ALPHABET, 0.25)
    bg = background([sequence for _, sequence in control])
    by_pipeline = {}
    for row in fit_rows:
        by_pipeline.setdefault(int(row["Pipeline"]), []).append(row)
    selected_by_rank = {}
    accepted = 0
    for row in diagnostics:
        pipeline = int(row["Pipeline"])
        rows = by_pipeline[pipeline]
        assert [int(site["SequenceIdx"]) for site in rows] == list(range(len(primary)))
        iterations = int(row["AcceptedIterations"])
        assert 0 <= iterations <= int(row["AttemptedIterations"]) <= 20
        motif_id = f"MorbiusPlus_pipeline_{pipeline}"
        present = [site for site in rows if site["SitePresent"] == "1"]
        assert len(present) == int(row["PrimarySiteNum"])
        if row["FitSource"] == "unsupported":
            assert not present and motif_id not in candidates and motif_id not in scoring
            continue
        site_num, matrix = candidates[motif_id]
        width = len(matrix)
        assert site_num == len(present) > 0
        counts = [[0] * 4 for _ in range(width)]
        for site in rows:
            index = int(site["SequenceIdx"])
            assert site["SequenceName"] == primary[index][0]
            if site["SitePresent"] == "0":
                assert site["Offset"] == site["Strand"] == "NA"
                continue
            assert site["SitePresent"] == "1"
            offset, strand = int(site["Offset"]), site["Strand"]
            assert strand in {"+", "-"} and 0 <= offset <= len(primary[index][1]) - width
            sequence = primary[index][1][offset:offset + width]
            if strand == "-":
                sequence = reverse_complement(sequence)
            for column, base in enumerate(sequence):
                counts[column][ALPHABET.index(base)] += 1
        expected = [[(count + 1) / (site_num + 4) for count in column] for column in counts]
        assert all(abs(a - b) < 2e-6 for column, expected_column in zip(matrix, expected)
                   for a, b in zip(column, expected_column))
        primary_best = [best_site(sequence, matrix, bg) for _, sequence in primary]
        control_best = [best_site(sequence, matrix, bg) for _, sequence in control]
        check_enrichment(row, "Evaluation", [site[0] for site in primary_best], [site[0] for site in control_best])
        if iterations:
            accepted += 1
            assert row["FitSource"] == "refinement"
            assert float(row["EvaluationLogPvalue"]) < float(row["InitialLogPvalue"])
            _, fit_matrix = scoring[motif_id]
            primary_fit = [best_site(sequence, fit_matrix, bg) for _, sequence in primary]
            control_fit = [best_site(sequence, fit_matrix, bg) for _, sequence in control]
            check_enrichment(row, "Fit", [site[0] for site in primary_fit], [site[0] for site in control_fit])
            threshold = float(row["FitScoreThreshold"])
            for record, (score, offset, strand) in zip(rows, primary_fit):
                assert (record["SitePresent"] == "1") == (score >= threshold - 2e-7)
                if record["SitePresent"] == "1":
                    # Rounding can alter an exact site tie; verify chosen score independently.
                    selected = primary[int(record["SequenceIdx"])][1][int(record["Offset"]):int(record["Offset"]) + width]
                    if record["Strand"] == "-":
                        selected = reverse_complement(selected)
                    selected_score = site_score(selected, fit_matrix, bg)
                    assert abs(selected_score - score) < 2e-6
        else:
            assert row["FitSource"] == "original_gibbs" and site_num == len(primary)
            assert motif_id not in scoring
        if row["OutputRank"] != "unreported":
            selected_by_rank[int(row["OutputRank"])] = (pipeline, expected, site_num)
    assert len(selected_by_rank) <= expected_limit
    assert sorted(selected_by_rank) == list(range(1, len(selected_by_rank) + 1))
    emitted = fasta(prefix + ".fasta")
    positions = tsv(prefix + ".offsets.tsv")
    assert len(emitted) == len(positions) == sum(item[2] for item in selected_by_rank.values())
    seen = set()
    for (header, sequence), row in zip(emitted, positions):
        rank = int(row.get("MotifRank", 1))
        pipeline = selected_by_rank[rank][0]
        index, offset = int(row["SequenceIdx"]), int(row["Offset"])
        assert (rank, index) not in seen
        seen.add((rank, index))
        diagnostic = by_pipeline[pipeline][index]
        assert diagnostic["SitePresent"] == "1"
        assert row["Offset"] == diagnostic["Offset"] and row["Strand"] == diagnostic["Strand"]
        expected = primary[index][1][offset:offset + len(sequence)]
        if row["Strand"] == "-":
            expected = reverse_complement(expected)
        assert row["SequenceName"] == primary[index][0]
        assert expected == sequence == row["Motif"]
        assert f"offset={offset} strand={row['Strand']}" in header
    final_motifs = meme(prefix + ".meme")
    assert len(final_motifs) == len(selected_by_rank)
    for rank, (site_num, matrix) in enumerate(final_motifs.values(), 1):
        _, expected, count = selected_by_rank[rank]
        assert site_num == count
        assert all(abs(a - b) < 2e-6 for column, expected_column in zip(matrix, expected)
                   for a, b in zip(column, expected_column))
    pwm_rows = tsv(prefix + ".pwm.tsv")
    assert len(pwm_rows) == sum(len(item[1]) for item in selected_by_rank.values())
    for row in pwm_rows:
        expected = selected_by_rank[int(row.get("MotifRank", 1))][1][int(row["Position"])]
        assert all(abs(float(row[base]) - value) < 2e-6 for base, value in zip(ALPHABET, expected))
    return accepted, diagnostics


def main():
    binary = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(74718)
    motif = "ACGATCGA"
    primary = ["".join(rng.choices(ALPHABET, k=32)) for _ in range(36)]
    control = ["".join(rng.choices(ALPHABET, k=32)) for _ in range(36)]
    for idx in range(16):
        offset = 5 + idx % 11
        site = motif if idx % 2 == 0 else reverse_complement(motif)
        primary[idx] = primary[idx][:offset] + site + primary[idx][offset + len(site):]
    write_fasta(work / "primary.fasta", primary)
    write_fasta(work / "control.fasta", control)

    def run(name, motif_count=16, control_file=None, alphabet="dna", expect_success=True):
        prefix = work / name
        command = [binary, "--input", str(work / "primary.fasta"), "--output", str(prefix),
                   "--alphabet", alphabet, "--motif-length", "8", "--motif-count", str(motif_count),
                   "--max-updates", "720", "--seed", "1", "--threads", "1"]
        if control_file is not None:
            command += ["--control", str(control_file)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        (work / (name + ".log")).write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == expect_success, (command, result.stdout, result.stderr)
        return prefix

    baseline = run("baseline")
    refined = run("refined", control_file=work / "control.fasta")
    assert all_pipeline_table(baseline) == all_pipeline_table(refined), "Control/refinement changed Gibbs scores or termination"
    for suffix in (".seeds.tsv", ".initial_offsets.tsv"):
        assert Path(str(baseline) + suffix).read_bytes() == Path(str(refined) + suffix).read_bytes(), "refinement changed initialization"
    primary_records, control_records = fasta(work / "primary.fasta"), fasta(work / "control.fasta")
    accepted, diagnostics = check_outputs(refined, primary_records, control_records, 16)
    assert accepted > 0, "fixture must exercise an accepted refinement"
    assert any(0 < int(row["PrimarySiteNum"]) < len(primary) for row in diagnostics), "fixture must exercise zero-or-one participation"
    single = run("single", motif_count=1, control_file=work / "control.fasta")
    check_outputs(single, primary_records, control_records, 1)
    assert meme(str(single) + ".refinement_candidates.meme") == meme(str(refined) + ".refinement_candidates.meme"), "all16 candidates must refine before output limit"
    unsupported = run("unsupported", control_file=work / "primary.fasta")
    _, diagnostics = check_outputs(unsupported, primary_records, primary_records, 16)
    assert all(row["FitSource"] == "unsupported" for row in diagnostics)
    assert not meme(str(unsupported) + ".meme") and not fasta(str(unsupported) + ".fasta")
    write_fasta(work / "short_control.fasta", [sequence[:-1] for sequence in control])
    run("length_mismatch", control_file=work / "short_control.fasta", expect_success=False)
    run("protein_rejected", alphabet="protein", control_file=work / "control.fasta", expect_success=False)
    print("Validated post-Gibbs independence, all16 refinements, supported sites/PWMs, zero-support output, and CLI errors.")


if __name__ == "__main__":
    main()

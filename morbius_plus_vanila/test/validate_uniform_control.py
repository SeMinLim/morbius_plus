"""Check automatic Control against an independently materialized fixed-seed stream."""
from pathlib import Path
import random
import subprocess
import sys

from validate_refinement_outputs import all_pipeline_table, check_outputs, fasta, meme, write_fasta


def uniform_control(primary):
    mask = (1 << 64) - 1
    state, word, remaining = 1, 0, 0
    records = []
    for _, sequence in primary:
        output = []
        for _ in sequence:
            if remaining == 0:
                state = (state + 0x9E3779B97F4A7C15) & mask
                value = ((state ^ (state >> 30)) * 0xBF58476D1CE4E5B9) & mask
                value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
                word, remaining = value ^ (value >> 31), 32
            output.append("ACGT"[word & 3])
            word >>= 2
            remaining -= 1
        records.append("".join(output))
    return records


def main():
    binary = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(74718)
    primary = ["".join(rng.choices("ACGT", k=33)) for _ in range(36)]
    for idx in range(16):
        offset = 5 + idx % 11
        primary[idx] = primary[idx][:offset] + "ACGATCGA" + primary[idx][offset + 8:]
    primary_path = work / "primary.fasta"
    control_path = work / "explicit_control.fasta"
    write_fasta(primary_path, primary)
    primary_records = fasta(primary_path)
    write_fasta(control_path, uniform_control(primary_records))
    control_records = fasta(control_path)
    suffixes = {
        ".fasta", ".offsets.tsv", ".pwm.tsv", ".meme", ".seeds.tsv", ".initial_offsets.tsv",
        ".summary.txt", ".refinement.tsv", ".refinement_sites.tsv", ".refinement_scoring.meme",
        ".refinement_candidates.meme", ".refinement_background.tsv", ".refinement_summary.txt",
    }

    def run(name, threads=1, seed=1, supplied=None, expect_success=True):
        directory = work / name
        directory.mkdir()
        prefix = directory / "result"
        command = [binary, "--input", str(primary_path), "--output", str(prefix),
                   "--alphabet", "dna", "--motif-length", "8", "--motif-count", "16",
                   "--max-updates", "720", "--seed", str(seed), "--threads", str(threads)]
        if supplied is not None:
            command += ["--control", str(supplied)]
        result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=60)
        (work / (name + ".log")).write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == expect_success, (command, result.stdout, result.stderr)
        if expect_success:
            assert {path.name for path in directory.iterdir()} == {"result" + suffix for suffix in suffixes}, \
                "generation must not create a Control FASTA or any intermediate file"
        return prefix

    results = []
    for threads, seed in ((1, 1), (4, 1), (4, 7)):
        auto = run(f"auto_t{threads}_s{seed}", threads, seed)
        explicit = run(f"explicit_t{threads}_s{seed}", threads, seed, control_path)
        check_outputs(auto, primary_records, control_records, 16)
        check_outputs(explicit, primary_records, control_records, 16)
        assert all_pipeline_table(auto) == all_pipeline_table(explicit), "Control generation changed Gibbs"
        for suffix in suffixes - {".summary.txt", ".refinement_summary.txt", ".refinement_background.tsv"}:
            assert Path(str(auto) + suffix).read_bytes() == Path(str(explicit) + suffix).read_bytes(), suffix
        background = lambda prefix: [line for line in Path(str(prefix) + ".refinement_background.tsv").read_text().splitlines()
                                     if not line.startswith("#")]
        assert background(auto) == background(explicit), "automatically generated Control does not match fixed stream"
        results.append(auto)
    for suffix in suffixes - {".summary.txt", ".refinement_summary.txt"}:
        assert Path(str(results[0]) + suffix).read_bytes() == Path(str(results[1]) + suffix).read_bytes(), \
            "thread count changed deterministic results"
    assert Path(str(results[0]) + ".initial_offsets.tsv").read_bytes() != \
        Path(str(results[2]) + ".initial_offsets.tsv").read_bytes(), "fixture must exercise distinct Gibbs seeds"
    assert any(meme(str(prefix) + ".meme") for prefix in results), "fixture must yield a supported auto-Control motif"
    same_as_primary = run("primary_as_control", supplied=primary_path)
    _, diagnostics = check_outputs(same_as_primary, primary_records, primary_records, 16)
    assert all(row["FitSource"] == "unsupported" for row in diagnostics), "external Control was replaced by random data"
    assert not meme(str(same_as_primary) + ".meme")
    run("missing_control", supplied=work / "does_not_exist.fasta", expect_success=False)
    assert "File not found:" in (work / "missing_control.log").read_text()
    print("Validated automatic/explicit Control equivalence, supplied-Control bypass, thread/seed isolation, and memory-only generation.")


if __name__ == "__main__":
    main()

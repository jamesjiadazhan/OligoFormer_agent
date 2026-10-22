"""Verification of cli_wrapper.compare_sirna_mismatch_kon against `python scripts/main.py -m 1`.

Upstream binding: scripts/main.py -> scripts/mismatch.py:mismatch (Args.mismatch == 1), which reads
the base and the mismatched RNA from stdin and prints calculate_kon(R1, R2). mismatch.py calls
model.eval() and the MismatchModule has no dropout, so kon is deterministic: reference values come
from executor mismatch_compare_pair (3 identical runs + a replay) and are compared exactly.
"""
from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

from conftest import (
    NOTEBOOKS,
    assert_tool_contract,
    call_tool,
    call_tool_expect_error,
)

TOOL = "compare_sirna_mismatch_kon"
PAIR_DIR = NOTEBOOKS / "mismatch_compare_pair"


def upstream_kon(case: str) -> tuple[str, str, float]:
    """Read a saved direct-upstream case: its stdin pair and the kon printed on stdout."""
    stdin_path = PAIR_DIR / "data" / f"{case}.stdin.txt"
    if case == "reference":
        stdout_path = PAIR_DIR / "outputs" / "run1" / "stdout.txt"
    else:
        stdout_path = PAIR_DIR / "outputs" / "changed_input" / case / "stdout.txt"
    base, mismatched = [l.strip() for l in stdin_path.read_text().splitlines() if l.strip()][:2]
    lines = [l.strip() for l in stdout_path.read_text().splitlines() if l.strip()]
    return base, mismatched, float(lines[-1])


ALL_CASES = ["reference", "a_identical", "b_2mm", "c_5mm", "d_other_base_1mm"]


def test_upstream_reference_is_reproducible():
    """Guard the expectation source itself: the three upstream runs + replay agree exactly."""
    runs = [PAIR_DIR / "outputs" / f"run{i}" / "stdout.txt" for i in (1, 2, 3)]
    runs += [PAIR_DIR / "outputs" / "replay" / f"run{i}" / "stdout.txt" for i in (1, 2, 3)]
    values = {float(p.read_text().strip().splitlines()[-1].strip()) for p in runs}
    assert values == {0.9940418601036072}


@pytest.mark.parametrize("case", ALL_CASES)
async def test_matches_direct_upstream_kon(mcp_server, case):
    base, mismatched, expected = upstream_kon(case)
    data = await call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched})
    assert_tool_contract(data, min_artifacts=0)
    assert data["kon"] == pytest.approx(expected, abs=1e-7), (
        f"{case}: kon {data['kon']!r} != direct upstream {expected!r}"
    )
    assert data["base_rna"] == base
    assert data["mismatched_rna"] == mismatched
    expected_positions = [i + 1 for i, (a, b) in enumerate(zip(base, mismatched)) if a != b]
    assert data["mismatch_positions"] == expected_positions
    assert data["num_mismatches"] == len(expected_positions)
    assert data["artifacts"] == []


async def test_reference_case_is_bit_identical_across_repeated_calls(mcp_server):
    base, mismatched, expected = upstream_kon("reference")
    first = await call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched})
    second = await call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched})
    assert first["kon"] == second["kon"] == expected


async def test_identical_pair_is_not_hardcoded_to_one(mcp_server):
    """The identical-pair case must return the learned score, not a synthesised 1.0."""
    base, mismatched, expected = upstream_kon("a_identical")
    assert base == mismatched
    data = await call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched})
    assert data["kon"] == pytest.approx(expected, abs=1e-7)
    assert data["kon"] != 1.0
    assert data["num_mismatches"] == 0


async def test_arguments_are_not_swapped(mcp_server):
    """kon(a, b) and kon(b, a) are different upstream, so argument order must be preserved."""
    base, mismatched, expected = upstream_kon("d_other_base_1mm")
    forward = await call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched})
    reverse = await call_tool(mcp_server, TOOL, {"base_rna": mismatched, "mismatched_rna": base})
    assert forward["kon"] == pytest.approx(expected, abs=1e-7)
    assert reverse["kon"] != pytest.approx(forward["kon"], abs=1e-6)


@pytest.mark.parametrize(
    "base,mismatched,fragment",
    [
        # upstream: "The length of mismatched RNA is not 19 nt!" (executor case e_len_19_vs_18)
        ("UCGCAAAUUAAAUUGAACC", "UCGCAAAUUAAAUUGAAC", "mismatched_rna"),
        ("UCGCAAAUUAAAUUGAAC", "UCGCAAAUUAAAUUGAACC", "base_rna"),
        # upstream: bare KeyError in one_hot_encode (executor cases f1_lowercase / f2_T_containing)
        ("ucgcaaauuaaauugaacc", "ucgcaaauuaaauugaacg", "base_rna"),
        ("UCGCAAAUUAAAUUGAACC", "TCGCAAATTAAATTGAACG", "mismatched_rna"),
        ("UCGCAAAUUAAAUUGAAC-", "UCGCAAAUUAAAUUGAACG", "base_rna"),
        ("", "UCGCAAAUUAAAUUGAACG", "base_rna"),
    ],
)
async def test_invalid_sequences_raise_mcp_errors(mcp_server, base, mismatched, fragment):
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched}
    )
    assert fragment in text
    assert "19 nt" in text


async def test_missing_required_argument_is_an_error(mcp_server):
    text = await call_tool_expect_error(mcp_server, TOOL, {"base_rna": "UCGCAAAUUAAAUUGAACC"})
    assert "mismatched_rna" in text


async def test_wrong_runtime_directory_is_reported(mcp_server, monkeypatch, tmp_path):
    empty = tmp_path / "not-oligoformer"
    empty.mkdir()
    monkeypatch.setenv("OLIGOFORMER_RUNTIME", str(empty))
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": "UCGCAAAUUAAAUUGAACC", "mismatched_rna": "UCGCAAAUUAAAUUGAACG"},
    )
    assert "not a complete OligoFormer runtime" in text
    assert "scripts/main.py" in text


async def test_missing_runtime_directory_is_reported(mcp_server, monkeypatch, tmp_path):
    monkeypatch.setenv("OLIGOFORMER_RUNTIME", str(tmp_path / "absent"))
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": "UCGCAAAUUAAAUUGAACC", "mismatched_rna": "UCGCAAAUUAAAUUGAACG"},
    )
    assert "runtime directory not found" in text


async def test_foreign_program_in_the_runtime_is_refused(mcp_server, monkeypatch, tmp_path, wrapper_module):
    """A runtime override whose scripts/main.py is a different program must be rejected."""
    real = Path(wrapper_module.PROJECT_ROOT) / "runtime" / "OligoFormer"
    fake = tmp_path / "fake-runtime"
    fake.mkdir()
    for relative in (
        "scripts/infer.py",
        "scripts/mismatch.py",
        "scripts/RNA-FM.sh",
        "model/best_model.pth",
        "model/mismatch_model.pth",
        "toxicity/cell_viability.txt",
        "RNA-FM/torch_home/hub/checkpoints/RNA-FM_pretrained.pth",
    ):
        target = fake / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("placeholder\n")
    (fake / "scripts" / "main.py").write_text("print('some other program')\n")
    assert (real / "scripts" / "main.py").is_file()
    monkeypatch.setenv("OLIGOFORMER_RUNTIME", str(fake))
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": "UCGCAAAUUAAAUUGAACC", "mismatched_rna": "UCGCAAAUUAAAUUGAACG"},
    )
    assert "does not look like the OligoFormer CLI entry point" in text


async def test_timeout_becomes_an_mcp_error(mcp_server, wrapper_module, monkeypatch):
    """A hung child process must surface as a tool error, not a truncated success."""
    monkeypatch.setattr(wrapper_module, "_MISMATCH_PAIR_TIMEOUT", 1)
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": "UCGCAAAUUAAAUUGAACC", "mismatched_rna": "UCGCAAAUUAAAUUGAACG"},
    )
    assert "timed out after 1s" in text


async def test_nonzero_exit_is_not_swallowed(mcp_server, wrapper_module, monkeypatch):
    """If the CLI exits nonzero the wrapper must raise instead of parsing partial stdout."""
    original = wrapper_module._run_cli

    def failing(runtime, args, **kwargs):
        return original(runtime, ["--infer", "99", "--unknown-flag"], **kwargs)

    monkeypatch.setattr(wrapper_module, "_run_cli", failing)
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": "UCGCAAAUUAAAUUGAACC", "mismatched_rna": "UCGCAAAUUAAAUUGAACG"},
    )
    assert "failed with exit status" in text


async def test_concurrent_calls_are_serialized_and_both_succeed(mcp_server):
    """The runtime lock must let two simultaneous calls through with correct, identical results."""
    base, mismatched, expected = upstream_kon("reference")
    _, other_mm, other_expected = upstream_kon("c_5mm")
    first, second = await asyncio.gather(
        call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": mismatched}),
        call_tool(mcp_server, TOOL, {"base_rna": base, "mismatched_rna": other_mm}),
    )
    assert first["kon"] == pytest.approx(expected, abs=1e-7)
    assert second["kon"] == pytest.approx(other_expected, abs=1e-7)

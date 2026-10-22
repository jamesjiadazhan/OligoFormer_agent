"""Verification of cli_wrapper.scan_sirna_mismatch_variants against `main.py -m 2|-m 3`.

Upstream binding: scripts/main.py -> scripts/mismatch.py:mismatch (Args.mismatch in (2, 3)), which
reads one base RNA from stdin, enumerates generate_mismatch_sequences(R1, k) for k = 1..2 (-m 2) or
1..3 (-m 3) and writes Args.output_dir + '/result/mismatch.txt'. The module is deterministic
(model.eval(), no dropout), so every kon is compared exactly against the executor's saved tables.
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
    read_table,
    sha256,
)

TOOL = "scan_sirna_mismatch_variants"
SCAN = NOTEBOOKS / "mismatch_scan_variants"
REFERENCE_BASE = (SCAN / "data" / "base_input.txt").read_text().strip()
CHANGED_BASE = (SCAN / "data" / "changed_input.txt").read_text().strip()


def upstream_rows(case: str) -> dict[str, tuple[str, str]]:
    header, rows = read_table(SCAN / "outputs" / case / "result" / "mismatch.txt")
    assert header == ["seq", "kon", "num_mismatch"], header
    table = {r["seq"]: (r["kon"], r["num_mismatch"]) for r in rows}
    assert len(table) == len(rows), f"{case}: duplicate seq keys in the upstream table"
    return table


def assert_table_matches_upstream(actual: Path, case: str) -> dict:
    expected = upstream_rows(case)
    header, rows = read_table(actual)
    assert header == ["seq", "kon", "num_mismatch"]
    got = {r["seq"]: (r["kon"], r["num_mismatch"]) for r in rows}
    assert len(got) == len(rows)
    assert set(got) == set(expected), (
        f"{case}: variant set differs (only-wrapper {len(set(got) - set(expected))}, "
        f"only-upstream {len(set(expected) - set(got))})"
    )
    worst = 0.0
    for seq, (kon, num) in expected.items():
        got_kon, got_num = got[seq]
        assert got_num == num, f"{case}: {seq} num_mismatch {got_num} != {num}"
        delta = abs(float(got_kon) - float(kon))
        worst = max(worst, delta)
        assert delta <= 1e-7, f"{case}: {seq} kon {got_kon} != upstream {kon}"
    return {"rows": len(rows), "max_abs_kon_delta": worst}


def test_upstream_reference_is_reproducible():
    """Guard the expectation source: the three upstream runs and the replay are identical."""
    hashes = {
        sha256(SCAN / "outputs" / c / "result" / "mismatch.txt")
        for c in ("run1", "run2", "run3", "replay")
    }
    assert len(hashes) == 1


async def test_reference_depth2_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(out_root)}
    )
    assert_tool_contract(data, min_artifacts=1, within=out_root)
    assert data["base_rna"] == REFERENCE_BASE
    assert data["max_mismatches"] == 2
    assert data["num_rows"] == 1596
    assert data["columns"] == ["seq", "kon", "num_mismatch"]
    assert data["rows_per_num_mismatch"] == {"1": 57, "2": 1539}
    table = Path(data["artifacts"][0]["path"])
    assert table.name == "mismatch.txt"
    stats = assert_table_matches_upstream(table, "run1")
    assert stats["rows"] == 1596


async def test_changed_base_matches_direct_upstream(mcp_server, out_root):
    """A different base sequence must reproduce the executor's changed-input run exactly."""
    data = await call_tool(
        mcp_server, TOOL, {"base_rna": CHANGED_BASE, "output_dir": str(out_root)}
    )
    assert data["base_rna"] == CHANGED_BASE
    assert data["num_rows"] == 1596
    assert_table_matches_upstream(Path(data["artifacts"][0]["path"]), "changed_input")
    # and it must differ from the reference base's table (no hardcoded example)
    reference = upstream_rows("run1")
    assert set(upstream_rows("changed_input")) != set(reference)


@pytest.mark.slow
async def test_depth3_parameter_is_honoured(mcp_server, out_root):
    """max_mismatches=3 must reach the upstream -m 3 branch (27759 variants, ~95 s)."""
    data = await call_tool(
        mcp_server,
        TOOL,
        {"base_rna": REFERENCE_BASE, "max_mismatches": 3, "output_dir": str(out_root)},
    )
    assert data["max_mismatches"] == 3
    assert data["num_rows"] == 27759
    assert data["rows_per_num_mismatch"] == {"1": 57, "2": 1539, "3": 26163}
    assert_table_matches_upstream(Path(data["artifacts"][0]["path"]), "depth3")


async def test_default_depth_is_two(mcp_server, out_root):
    """The documented default must stay 2; -m 3 must not be selected implicitly."""
    schema = None
    from fastmcp import Client

    async with Client(mcp_server) as client:
        for tool in await client.list_tools():
            if tool.name == TOOL:
                schema = tool.inputSchema
    assert schema is not None
    assert schema["properties"]["max_mismatches"]["default"] == 2
    assert schema["required"] == ["base_rna"]
    assert schema["properties"]["max_mismatches"]["enum"] == [2, 3]
    for name, prop in schema["properties"].items():
        assert prop.get("description"), f"parameter {name} has no Annotated description"


@pytest.mark.parametrize("depth", [1, 4, 0, -1])
async def test_unsupported_depth_is_rejected(mcp_server, out_root, depth):
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {"base_rna": REFERENCE_BASE, "max_mismatches": depth, "output_dir": str(out_root)},
    )
    assert "max_mismatches" in text


@pytest.mark.parametrize(
    "base",
    [
        "UCGCAAAUUAAAUUGAAC",  # 18 nt -> upstream Exception (executor case err_18nt)
        "TCGCAAATTAAATTGAACC",  # DNA T -> upstream KeyError (executor case err_T)
        "ucgcaaauuaaauugaacc",
        "UCGCAAAUUAAAUUGAACCU",
        "",
    ],
)
async def test_invalid_base_raises_mcp_error(mcp_server, out_root, base):
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"base_rna": base, "output_dir": str(out_root)}
    )
    assert "base_rna" in text and "19 nt" in text
    # nothing may be written for a rejected input
    assert not list(out_root.rglob("mismatch.txt"))


async def test_repeated_calls_use_distinct_output_directories(mcp_server, out_root):
    first = await call_tool(
        mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(out_root)}
    )
    second = await call_tool(
        mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(out_root)}
    )
    assert first["output_dir"] != second["output_dir"]
    assert Path(first["artifacts"][0]["path"]) != Path(second["artifacts"][0]["path"])
    assert Path(first["artifacts"][0]["path"]).is_file()
    assert Path(second["artifacts"][0]["path"]).is_file()
    assert sha256(Path(first["artifacts"][0]["path"])) == sha256(
        Path(second["artifacts"][0]["path"])
    ), "the deterministic mismatch table must be identical across repeated calls"


async def test_output_directory_with_spaces(mcp_server, out_root):
    spaced = out_root / "dir with spaces" / "second level"
    spaced.mkdir(parents=True, exist_ok=True)
    data = await call_tool(
        mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(spaced)}
    )
    table = Path(data["artifacts"][0]["path"])
    assert "dir with spaces" in str(table)
    assert_table_matches_upstream(table, "run1")


async def test_concurrent_scans_are_serialized_with_distinct_outputs(mcp_server, out_root):
    """Two simultaneous scans must both succeed in separate directories (runtime lock)."""
    first, second = await asyncio.gather(
        call_tool(mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(out_root)}),
        call_tool(mcp_server, TOOL, {"base_rna": CHANGED_BASE, "output_dir": str(out_root)}),
    )
    assert first["output_dir"] != second["output_dir"]
    assert_table_matches_upstream(Path(first["artifacts"][0]["path"]), "run1")
    assert_table_matches_upstream(Path(second["artifacts"][0]["path"]), "changed_input")


async def test_timeout_becomes_an_mcp_error(mcp_server, wrapper_module, monkeypatch, out_root):
    monkeypatch.setattr(wrapper_module, "_MISMATCH_SCAN_TIMEOUT", 1)
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"base_rna": REFERENCE_BASE, "output_dir": str(out_root)}
    )
    assert "timed out after 1s" in text

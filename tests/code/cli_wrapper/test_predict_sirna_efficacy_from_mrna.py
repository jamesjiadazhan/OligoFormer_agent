"""Verification of cli_wrapper.predict_sirna_efficacy_from_mrna against `main.py --infer 1 -i1 ...`.

Upstream binding: scripts/main.py -> scripts/infer.py:infer with Args.infer == 1 and no
Args.infer_siRNA_fasta, i.e. the branch that builds antiRNA(mRNA[i:i+19]) for every window.
Reference tables come from direct upstream runs:
  * default flags                 notebooks/infer_mrna_traverse/outputs/run1 (+run2/run3/replay)
  * -tox                          notebooks/infer_toxicity_filter/outputs/run1
  * -tox --toxicity_threshold 80  notebooks/infer_toxicity_filter/outputs/changed_params/thr80
  * -tox -nf                      notebooks/infer_toxicity_filter/outputs/changed_params/nf
  * -nf                           tests/results/cli_wrapper/upstream_nf_only/run1 (this verifier)
  * changed input (HuMRNA3)       notebooks/infer_mrna_traverse/outputs/changed_input
"""
from __future__ import annotations

from pathlib import Path

import pytest

from conftest import (
    DATA,
    NOTEBOOKS,
    RESULTS,
    assert_tool_contract,
    call_tool,
    call_tool_expect_error,
    compare_result_tables,
    read_table,
    sha256,
)

TOOL = "predict_sirna_efficacy_from_mrna"
TRAVERSE = NOTEBOOKS / "infer_mrna_traverse" / "outputs"
TOXICITY = NOTEBOOKS / "infer_toxicity_filter" / "outputs"
NF_ONLY = RESULTS / "upstream_nf_only" / "run1"

EXAMPLE = DATA / "example.fa"
CHANGED = DATA / "changed_input_Hu_mRNA_RNA3.fa"


def result_table(data: dict, record: str, suffix: str = "") -> Path:
    path = Path(data["output_dir"]) / f"{record}{suffix}.txt"
    assert path.is_file(), f"{path} missing; artifacts were {data['artifacts']}"
    return path


def test_upstream_reference_rows_are_stable_across_runs():
    """Guard the expectation source: the deterministic columns agree across the 3 upstream runs."""
    for record in ("RNA0", "RNA1"):
        base = TRAVERSE / "run1" / f"{record}.txt"
        for other in ("run2", "run3", "replay"):
            compare_result_tables(
                TRAVERSE / other / f"{record}.txt", base, label=f"{other}/{record}"
            )


async def test_reference_case_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)})
    assert_tool_contract(data, min_artifacts=6, within=out_root)
    assert data["records"] == ["RNA0", "RNA1"]
    assert data["windows_per_record"] == {"RNA0": 39, "RNA1": 39}
    assert data["stochastic_efficacy"] is True
    for record in ("RNA0", "RNA1"):
        summary = data["per_record_summary"][record]
        assert summary["rows"] == 39
        assert summary["columns"] == ["pos", "sense", "siRNA", "efficacy", "func_filter", "filter"]
        stats = compare_result_tables(
            result_table(data, record), TRAVERSE / "run1" / f"{record}.txt", label=record
        )
        assert stats["rows"] == 39
    # func_filter distribution and filter pass counts are deterministic
    assert data["per_record_summary"]["RNA0"]["func_filter_counts"] == {"0": 16, "1": 22, "4": 1}
    assert data["per_record_summary"]["RNA1"]["func_filter_counts"] == {"0": 24, "1": 15}
    assert data["per_record_summary"]["RNA0"]["rows_passing_filter"] == 16
    assert data["per_record_summary"]["RNA1"]["rows_passing_filter"] == 24


async def test_ranked_and_filtered_tables_follow_upstream_contract(mcp_server, out_root):
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)})
    for record in ("RNA0", "RNA1"):
        header, rows = read_table(result_table(data, record))
        _, ranked = read_table(result_table(data, record, "_ranked"))
        _, filtered = read_table(result_table(data, record, "_ranked_filtered"))
        assert len(ranked) == len(rows)
        # ranked is the same row multiset sorted by efficacy descending
        assert sorted(r["pos"] for r in ranked) == sorted(r["pos"] for r in rows)
        efficacies = [float(r["efficacy"]) for r in ranked]
        assert efficacies == sorted(efficacies, reverse=True)
        assert [r for r in rows if r["filter"] == "0"].__len__() == len(filtered)
        assert all(r["filter"] == "0" for r in filtered)
        assert data["per_record_summary"][record]["ranked_filtered_rows"] == len(filtered)


async def test_window_fastas_match_direct_upstream_exactly(mcp_server, out_root):
    """The per-window siRNA/mRNA FASTAs handed to RNA-FM are deterministic; compare byte-wise."""
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)})
    for record in ("RNA0", "RNA1"):
        for name in ("siRNA.fa", "mRNA.fa"):
            produced = Path(data["output_dir"]) / "windows" / record / name
            expected = TRAVERSE / "run1" / "data_infer" / record / name
            assert produced.is_file(), produced
            assert sha256(produced) == sha256(expected), f"{record}/{name} differs from upstream"


async def test_toxicity_filter_flag_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {"mrna_fasta_path": str(EXAMPLE), "include_toxicity": True, "output_dir": str(out_root)},
    )
    for record in ("RNA0", "RNA1"):
        summary = data["per_record_summary"][record]
        assert summary["columns"] == [
            "pos", "sense", "siRNA", "efficacy", "func_filter",
            "Seed", "cell_viability", "toxicity_filter", "filter",
        ]
        compare_result_tables(
            result_table(data, record), TOXICITY / "run1" / f"{record}.txt", label=f"tox/{record}"
        )
    assert data["per_record_summary"]["RNA0"]["rows_passing_filter"] == 10
    assert data["per_record_summary"]["RNA1"]["rows_passing_filter"] == 18


async def test_toxicity_threshold_80_matches_direct_upstream(mcp_server, out_root):
    """A changed threshold must reach --toxicity_threshold and change toxicity_filter upstream."""
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "include_toxicity": True,
            "toxicity_threshold": 80.0,
            "output_dir": str(out_root),
        },
    )
    for record in ("RNA0", "RNA1"):
        compare_result_tables(
            result_table(data, record),
            TOXICITY / "changed_params" / "thr80" / f"{record}.txt",
            label=f"thr80/{record}",
        )
    # the parameter really changed the outcome relative to the default threshold
    _, thr50 = read_table(TOXICITY / "run1" / "RNA0.txt")
    _, thr80 = read_table(result_table(data, "RNA0"))
    assert [r["toxicity_filter"] for r in thr50] != [r["toxicity_filter"] for r in thr80]


async def test_functionality_filter_off_with_toxicity_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "apply_functionality_filter": False,
            "include_toxicity": True,
            "output_dir": str(out_root),
        },
    )
    for record in ("RNA0", "RNA1"):
        summary = data["per_record_summary"][record]
        assert "func_filter" not in summary["columns"]
        compare_result_tables(
            result_table(data, record),
            TOXICITY / "changed_params" / "nf" / f"{record}.txt",
            label=f"nf_tox/{record}",
        )


async def test_functionality_filter_off_alone_matches_direct_upstream(mcp_server, out_root):
    """-nf without -tox: compared against this verifier's own direct upstream run."""
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "apply_functionality_filter": False,
            "output_dir": str(out_root),
        },
    )
    for record in ("RNA0", "RNA1"):
        summary = data["per_record_summary"][record]
        assert summary["columns"] == ["pos", "sense", "siRNA", "efficacy", "filter"]
        assert summary["rows_passing_filter"] == 39
        compare_result_tables(
            result_table(data, record), NF_ONLY / f"{record}.txt", label=f"nf/{record}"
        )
    assert "filter off" in data["message"]


async def test_changed_input_record_matches_direct_upstream(mcp_server, out_root):
    """A different transcript (Hu_mRNA.fa RNA3) must reproduce the executor's changed-input run."""
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(CHANGED), "output_dir": str(out_root)})
    assert data["records"] == ["HuMRNA3"]
    assert data["windows_per_record"] == {"HuMRNA3": 39}
    compare_result_tables(
        result_table(data, "HuMRNA3"),
        TRAVERSE / "changed_input" / "HuMRNA3.txt",
        label="HuMRNA3",
    )
    # not the reference example's result
    _, changed_rows = read_table(result_table(data, "HuMRNA3"))
    _, reference_rows = read_table(TRAVERSE / "run1" / "RNA0.txt")
    assert [r["siRNA"] for r in changed_rows] != [r["siRNA"] for r in reference_rows]


async def test_single_record_fasta_matches_the_two_record_run(mcp_server, out_root):
    """Record handling must be per-record: RNA0 alone equals RNA0 inside example.fa."""
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(out_root)})
    assert data["records"] == ["RNA0"]
    compare_result_tables(
        result_table(data, "RNA0"), TRAVERSE / "run1" / "RNA0.txt", label="only_RNA0"
    )


async def test_user_input_file_is_not_modified(mcp_server, out_root):
    before = sha256(EXAMPLE)
    stat_before = EXAMPLE.stat().st_mtime_ns
    await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)})
    assert sha256(EXAMPLE) == before
    assert EXAMPLE.stat().st_mtime_ns == stat_before


async def test_repeated_calls_produce_distinct_artifact_directories(mcp_server, out_root):
    first = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(out_root)})
    second = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(out_root)})
    assert first["output_dir"] != second["output_dir"]
    assert set(a["path"] for a in first["artifacts"]).isdisjoint(a["path"] for a in second["artifacts"])
    for data in (first, second):
        compare_result_tables(
            result_table(data, "RNA0"), TRAVERSE / "run1" / "RNA0.txt", label="repeat"
        )


async def test_output_directory_with_spaces(mcp_server, out_root):
    spaced = out_root / "run outputs with spaces"
    spaced.mkdir(parents=True, exist_ok=True)
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(spaced)})
    assert "run outputs with spaces" in data["output_dir"]
    compare_result_tables(
        result_table(data, "RNA0"), TRAVERSE / "run1" / "RNA0.txt", label="spaced"
    )


async def test_input_path_with_spaces(mcp_server, out_root, tmp_path):
    spaced_input = tmp_path / "fasta inputs" / "only RNA0.fa"
    spaced_input.parent.mkdir(parents=True)
    spaced_input.write_bytes((DATA / "only_RNA0.fa").read_bytes())
    data = await call_tool(mcp_server, TOOL, {"mrna_fasta_path": str(spaced_input), "output_dir": str(out_root)})
    compare_result_tables(
        result_table(data, "RNA0"), TRAVERSE / "run1" / "RNA0.txt", label="spaced_input"
    )


@pytest.mark.parametrize(
    "fixture,fragment",
    [
        ("missing.fa", "mRNA FASTA not found"),
        ("bad_empty.fa", "No FASTA records found"),
        ("bad_header_space.fa", "Unsupported FASTA record name"),
        ("bad_short_record.fa", "requires at least 19 nt"),
        ("bad_alphabet.fa", "unsupported characters"),
        ("bad_duplicate_header.fa", "Duplicate FASTA record name"),
        ("bad_header_parent.fa", "single path component"),
    ],
)
async def test_invalid_inputs_raise_mcp_errors(mcp_server, out_root, fixture, fragment):
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(DATA / fixture), "output_dir": str(out_root)}
    )
    assert fragment in text


async def test_parent_directory_record_name_cannot_touch_the_runtime(mcp_server, wrapper_module):
    """A record named '..' would make upstream's ./data/infer/<name>/ scratch path - and the
    wrapper's cleanup of it - resolve to the runtime's whole data/ directory, so it is rejected
    before anything runs."""
    runtime_data = Path(wrapper_module.PROJECT_ROOT) / "runtime" / "OligoFormer" / "data"
    before = sorted(p.name for p in runtime_data.iterdir())
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "bad_header_parent.fa")}
    )
    assert "single path component" in text
    assert sorted(p.name for p in runtime_data.iterdir()) == before
    assert (runtime_data / "example.fa").is_file()


async def test_schema_documents_every_parameter(mcp_server):
    from fastmcp import Client

    async with Client(mcp_server) as client:
        tools = {t.name: t for t in await client.list_tools()}
    schema = tools[TOOL].inputSchema
    assert schema["required"] == ["mrna_fasta_path"]
    for name, prop in schema["properties"].items():
        assert prop.get("description"), f"parameter {name} has no Annotated description"
    assert schema["properties"]["apply_functionality_filter"]["default"] is True
    assert schema["properties"]["include_toxicity"]["default"] is False
    assert schema["properties"]["toxicity_threshold"]["default"] == 50.0
    assert tools[TOOL].description.strip().count("\n") == 1


async def test_timeout_becomes_an_mcp_error(mcp_server, wrapper_module, monkeypatch, out_root):
    monkeypatch.setattr(wrapper_module, "_INFER_TIMEOUT_BASE", 1)
    monkeypatch.setattr(wrapper_module, "_INFER_TIMEOUT_PER_WINDOW", 0)
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(out_root)}
    )
    assert "timed out after 1s" in text


async def test_missing_result_table_is_detected(mcp_server, wrapper_module, monkeypatch, out_root):
    """A silently incomplete run (upstream ignores its sub-process exit codes) must error."""
    original = wrapper_module._run_cli

    def truncated(runtime, args, **kwargs):
        proc = original(runtime, args, **kwargs)
        for table in Path(out_root).rglob("RNA0.txt"):
            table.unlink()
        return proc

    monkeypatch.setattr(wrapper_module, "_run_cli", truncated)
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(DATA / "only_RNA0.fa"), "output_dir": str(out_root)}
    )
    assert "was not produced for record" in text

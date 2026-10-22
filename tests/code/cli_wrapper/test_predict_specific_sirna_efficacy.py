"""Verification of cli_wrapper.predict_specific_sirna_efficacy against `-i1 ... -i2 ...`.

Upstream binding: scripts/main.py -> scripts/infer.py:infer with Args.infer == 1 and
Args.infer_siRNA_fasta set, i.e. the branch that keeps only those user siRNAs whose
antiRNA(siRNA) occurs in the (upper-cased, T->U) mRNA, using the first match's start offset.

Reference tables come from direct upstream runs:
  * example.fa + example_siRNA.fa   notebooks/infer_specific_sirnas/outputs/run1 (+run2/3/replay)
  * example.fa + changed_siRNA.fa   notebooks/infer_specific_sirnas/outputs/changed_input
  * the record pre-filtering adaptation is checked against this verifier's own upstream run of the
    matching record alone: tests/results/cli_wrapper/upstream_only_RNA0/run{1,2}
"""
from __future__ import annotations

import re
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

TOOL = "predict_specific_sirna_efficacy"
SPECIFIC = NOTEBOOKS / "infer_specific_sirnas" / "outputs"
ONLY_RNA0 = RESULTS / "upstream_only_RNA0"

EXAMPLE = DATA / "example.fa"
EXAMPLE_SIRNA = DATA / "example_siRNA.fa"
CHANGED_SIRNA = DATA / "changed_siRNA.fa"
ZERO_MATCH_SIRNA = DATA / "zero_match_RNA1_siRNA.fa"


def result_table(data: dict, record: str, suffix: str = "") -> Path:
    path = Path(data["output_dir"]) / f"{record}{suffix}.txt"
    assert path.is_file(), f"{path} missing; artifacts were {data['artifacts']}"
    return path


def upstream_sirna_order(case: str, record: str) -> list[str]:
    _, rows = read_table(SPECIFIC / case / f"{record}.txt")
    return [r["siRNA"] for r in rows]


# ------------------------------------------------------------------ source-level review of the
# pre-scan: it must use infer.py's own matching rule and only drop records, never siRNAs.


def test_prescan_matching_rule_equals_infer_py(wrapper_module):
    """_anti_rna must agree with infer.py:antiRNA on every supported input."""
    import importlib.util

    upstream_path = (
        Path(wrapper_module.PROJECT_ROOT) / "repo" / "OligoFormer" / "scripts" / "infer.py"
    )
    source = upstream_path.read_text()
    start = source.index("def antiRNA(RNA):")
    end = source.index("def Calculate_DGH")
    namespace: dict = {}
    exec(compile(source[start:end], str(upstream_path), "exec"), namespace)
    upstream_anti = namespace["antiRNA"]

    import itertools
    import random

    rng = random.Random(0)
    cases = ["".join(p) for p in itertools.product("ACGU", repeat=4)]
    cases += ["".join(rng.choice("ACGU") for _ in range(19)) for _ in range(500)]
    for seq in cases:
        assert wrapper_module._anti_rna(seq) == upstream_anti(seq), seq


def test_prescan_only_selects_records(wrapper_module):
    """The pre-scan must pass every user siRNA to upstream and only filter whole records."""
    source = (Path(wrapper_module.__file__)).read_text()
    body = source[source.index("def predict_specific_sirna_efficacy") : source.index("def screen_sirna_offtarget")]
    # the siRNA FASTA handed to upstream is always the user's file, never a rewritten one
    assert '"-i2", str(sirna_path)' in body
    assert "sirna_path" in body and "filtered_fasta.write_text" in body
    # the rewritten FASTA only ever contains mRNA records
    assert 'f">{name}\\n{seq}\\n" for name, seq in kept' in body
    # and the pre-scan uses re.search on the reverse complement, like infer.py
    assert "re.search(_anti_rna(s), norm)" in body


def test_prescan_agrees_with_upstream_row_order(wrapper_module):
    """The pre-scan's matched lists must equal the siRNA column of the direct upstream tables."""
    _, records = wrapper_module._load_mrna_fasta(str(EXAMPLE))
    for sirna_file, case in ((EXAMPLE_SIRNA, "run1"), (CHANGED_SIRNA, "changed_input")):
        _, sirnas = wrapper_module._load_sirna_fasta(str(sirna_file))
        for name, seq in records:
            norm = seq.upper().replace("T", "U")
            hits = [s for s in sirnas if re.search(wrapper_module._anti_rna(s), norm)]
            assert hits == upstream_sirna_order(case, name), (name, case)


# ------------------------------------------------------------------------- scientific comparisons


def test_upstream_reference_rows_are_stable_across_runs():
    for record in ("RNA0", "RNA1"):
        base = SPECIFIC / "run1" / f"{record}.txt"
        for other in ("run2", "run3", "replay"):
            compare_result_tables(SPECIFIC / other / f"{record}.txt", base, label=f"{other}/{record}")


async def test_reference_case_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(EXAMPLE_SIRNA),
            "output_dir": str(out_root),
        },
    )
    assert_tool_contract(data, min_artifacts=6, within=out_root)
    assert data["records_scored"] == ["RNA0", "RNA1"]
    assert data["unmatched_records"] == []
    assert data["unmatched_sirnas"] == []
    assert data["matched_sirnas_per_record"] == {
        "RNA0": ["UCGCAAAUUAAAUUGAACC"],
        "RNA1": ["UCUCUUUCGCACAUUAAAU"],
    }
    for record in ("RNA0", "RNA1"):
        assert data["per_record_summary"][record]["rows"] == 1
        compare_result_tables(
            result_table(data, record), SPECIFIC / "run1" / f"{record}.txt", label=record
        )
    # no rewritten FASTA artifact when every record matched
    assert not any("actually passed upstream" in a["description"] for a in data["artifacts"])


async def test_changed_sirna_set_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(CHANGED_SIRNA),
            "output_dir": str(out_root),
        },
    )
    assert data["records_scored"] == ["RNA0", "RNA1"]
    assert data["unmatched_records"] == []
    assert data["unmatched_sirnas"] == ["GCUAGCUAGCUAGCUAGCU"]
    assert data["matched_sirnas_per_record"] == {
        "RNA0": ["GGAGGUCCACUGGCGCCCA", "GGCGCCCACGUCCGUAAGG"],
        "RNA1": ["GGCGCCCACGUCCGUAAGG", "CACAUUAAAUUGAACCUAU"],
    }
    for record in ("RNA0", "RNA1"):
        assert data["per_record_summary"][record]["rows"] == 2
        compare_result_tables(
            result_table(data, record),
            SPECIFIC / "changed_input" / f"{record}.txt",
            label=f"changed/{record}",
        )
        # every row of the changed set is filtered out upstream -> empty ranked_filtered table
        assert data["per_record_summary"][record]["ranked_filtered_rows"] == 0


async def test_record_prefiltering_reproduces_the_single_record_upstream_run(mcp_server, out_root):
    """The adaptation for upstream's zero-match crash must give the same table as running the
    matching record alone directly upstream (tests/results/cli_wrapper/upstream_only_RNA0)."""
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(ZERO_MATCH_SIRNA),
            "output_dir": str(out_root),
        },
    )
    assert data["records_scored"] == ["RNA0"]
    assert data["unmatched_records"] == ["RNA1"]
    assert data["unmatched_sirnas"] == []
    compare_result_tables(
        result_table(data, "RNA0"), ONLY_RNA0 / "run1" / "RNA0.txt", label="prefiltered"
    )
    compare_result_tables(
        result_table(data, "RNA0"), ONLY_RNA0 / "run2" / "RNA0.txt", label="prefiltered_run2"
    )
    # and it also equals the rows upstream managed to write before it crashed on RNA1
    compare_result_tables(
        result_table(data, "RNA0"),
        SPECIFIC / "error_cases" / "zero_match_RNA1" / "RNA0.txt",
        label="upstream_partial",
    )
    assert not (Path(data["output_dir"]) / "RNA1.txt").exists()
    # the FASTA actually handed upstream is reported and contains only the retained record
    rewritten = [a["path"] for a in data["artifacts"] if "actually passed upstream" in a["description"]]
    assert len(rewritten) == 1
    assert Path(rewritten[0]).read_text().count(">") == 1
    assert ">RNA0" in Path(rewritten[0]).read_text()


async def test_upstream_crashes_on_the_same_input_without_prefiltering():
    """Evidence that the adaptation is necessary, not a scientific change (executor exit 1)."""
    commands = (NOTEBOOKS / "infer_specific_sirnas" / "logs" / "commands.tsv").read_text()
    row = [l for l in commands.splitlines() if l.startswith("error_zero_match_RNA1")][0]
    assert row.split("\t")[1] == "1"
    stderr = (NOTEBOOKS / "infer_specific_sirnas" / "logs" / "error_zero_match_RNA1.stderr").read_text()
    assert "FileNotFoundError" in stderr and "./data/infer/RNA1/siRNA.fa" in stderr


async def test_no_matching_sirna_at_all_is_an_mcp_error(mcp_server, out_root, tmp_path):
    nomatch = tmp_path / "nomatch_siRNA.fa"
    nomatch.write_text(">nomatch\nGCUAGCUAGCUAGCUAGCU\n")
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(nomatch),
            "output_dir": str(out_root),
        },
    )
    assert "match any record" in text


async def test_user_input_files_are_not_modified(mcp_server, out_root):
    before = (sha256(EXAMPLE), sha256(EXAMPLE_SIRNA))
    await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(EXAMPLE_SIRNA),
            "output_dir": str(out_root),
        },
    )
    assert (sha256(EXAMPLE), sha256(EXAMPLE_SIRNA)) == before


async def test_repeated_calls_produce_distinct_artifact_directories(mcp_server, out_root):
    args = {
        "mrna_fasta_path": str(EXAMPLE),
        "sirna_fasta_path": str(EXAMPLE_SIRNA),
        "output_dir": str(out_root),
    }
    first = await call_tool(mcp_server, TOOL, args)
    second = await call_tool(mcp_server, TOOL, args)
    assert first["output_dir"] != second["output_dir"]
    assert set(a["path"] for a in first["artifacts"]).isdisjoint(a["path"] for a in second["artifacts"])


async def test_paths_with_spaces(mcp_server, out_root, tmp_path):
    spaced_dir = tmp_path / "inputs with spaces"
    spaced_dir.mkdir()
    mrna = spaced_dir / "example copy.fa"
    sirna = spaced_dir / "siRNA copy.fa"
    mrna.write_bytes(EXAMPLE.read_bytes())
    sirna.write_bytes(EXAMPLE_SIRNA.read_bytes())
    out = out_root / "results with spaces"
    out.mkdir(parents=True, exist_ok=True)
    data = await call_tool(
        mcp_server,
        TOOL,
        {"mrna_fasta_path": str(mrna), "sirna_fasta_path": str(sirna), "output_dir": str(out)},
    )
    assert "results with spaces" in data["output_dir"]
    for record in ("RNA0", "RNA1"):
        compare_result_tables(
            result_table(data, record), SPECIFIC / "run1" / f"{record}.txt", label=f"spaced/{record}"
        )


def test_sirna_alphabet_restriction_is_source_backed():
    """Direct upstream runs by this verifier: T-containing / lowercase siRNAs crash with a bare
    KeyError in infer.py:Calculate_end_diff, so rejecting them is a clarified upstream failure."""
    logs = RESULTS.parent.parent / "logs" / "cli_wrapper"
    with_t = (logs / "upstream_bad_sirna_with_T.stderr").read_text()
    lowercase = (logs / "upstream_bad_sirna_lowercase.stderr").read_text()
    assert "KeyError: 'TC'" in with_t
    assert "KeyError: 'uc'" in lowercase
    for case in ("bad_sirna_with_T", "bad_sirna_lowercase"):
        assert not list((RESULTS / "upstream_invalid_sirna" / case).glob("*.txt"))


@pytest.mark.parametrize(
    "sirna_fixture,fragment",
    [
        # upstream: Exception("The length of some siRNA is not 19 nt!") (executor error_len20)
        ("error_20nt_siRNA.fa", "Invalid siRNA on line"),
        ("bad_sirna_with_T.fa", "Invalid siRNA on line"),
        ("bad_sirna_lowercase.fa", "Invalid siRNA on line"),
        ("bad_sirna_blank_line.fa", "Invalid siRNA on line"),
        ("missing_sirna.fa", "siRNA FASTA not found"),
        ("bad_empty.fa", "No siRNA sequences found"),
    ],
)
async def test_invalid_sirna_input_raises_mcp_error(mcp_server, out_root, sirna_fixture, fragment):
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "sirna_fasta_path": str(DATA / sirna_fixture),
            "output_dir": str(out_root),
        },
    )
    assert fragment in text


@pytest.mark.parametrize(
    "mrna_fixture,fragment",
    [
        ("missing.fa", "mRNA FASTA not found"),
        ("bad_header_space.fa", "Unsupported FASTA record name"),
        ("bad_short_record.fa", "requires at least 19 nt"),
    ],
)
async def test_invalid_mrna_input_raises_mcp_error(mcp_server, out_root, mrna_fixture, fragment):
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(DATA / mrna_fixture),
            "sirna_fasta_path": str(EXAMPLE_SIRNA),
            "output_dir": str(out_root),
        },
    )
    assert fragment in text


async def test_schema_documents_every_parameter(mcp_server):
    from fastmcp import Client

    async with Client(mcp_server) as client:
        tools = {t.name: t for t in await client.list_tools()}
    schema = tools[TOOL].inputSchema
    assert schema["required"] == ["mrna_fasta_path", "sirna_fasta_path"]
    for name, prop in schema["properties"].items():
        assert prop.get("description"), f"parameter {name} has no Annotated description"
    assert tools[TOOL].description.strip().count("\n") == 1


async def test_missing_required_argument_is_an_error(mcp_server, out_root):
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)}
    )
    assert "sirna_fasta_path" in text

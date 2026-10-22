"""Verification of cli_wrapper.screen_sirna_offtarget against `main.py --infer 1 -i1 ... -off`.

Upstream binding: scripts/main.py -> scripts/infer.py:infer with Args.off_target and
Args.top_n == -1, i.e. the branch that runs `bash scripts/pita.sh` and `bash scripts/targetscan.sh`
over data/infer/<record>/siRNA.fa and merges the per-siRNA minimum PITA ddG / maximum TargetScan
context score back onto the result table.

Reference tables come from direct upstream runs (executor infer_offtarget_toxicity):
  * defaults                                       outputs/run1 (+run2/run3/replay, identical)
  * --pita_threshold -5 --targetscan_threshold 0.2 outputs/changed_params

Each successful off-target run costs several minutes, so exactly two are executed here; they are
marked `slow` but are not skipped.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from conftest import (
    DATA,
    NOTEBOOKS,
    assert_tool_contract,
    call_tool,
    call_tool_expect_error,
    compare_result_tables,
    read_table,
    sha256,
)

TOOL = "screen_sirna_offtarget"
OFF = NOTEBOOKS / "infer_offtarget_toxicity" / "outputs"
EXAMPLE = DATA / "example.fa"
RECORDS = ("RNA0", "RNA1")

EXPECTED_COLUMNS = [
    "pos", "sense", "siRNA", "efficacy", "func_filter",
    "pita_score", "targetscan_score", "off_target_filter",
    "Seed", "cell_viability", "toxicity_filter", "filter",
]
# row counts of the deterministic PITA / TargetScan tables (executor outputs/comparison.json)
EXPECTED_PITA_ROWS = {"RNA0": 1771, "RNA1": 1732}


def result_table(data: dict, record: str, suffix: str = "") -> Path:
    path = Path(data["output_dir"]) / f"{record}{suffix}.txt"
    assert path.is_file(), f"{path} missing; artifacts were {data['artifacts']}"
    return path


def test_upstream_reference_is_reproducible():
    """Guard the expectation source: PITA/TargetScan tables are byte-identical across runs."""
    for record in RECORDS:
        for name in ("pita.tab", "targetscan.tab"):
            hashes = {
                sha256(OFF / case / "infer" / record / name)
                for case in ("run1", "run2", "run3", "replay")
            }
            assert len(hashes) == 1, f"{record}/{name} is not reproducible upstream"
        base = OFF / "run1" / f"{record}.txt"
        for case in ("run2", "run3", "replay"):
            compare_result_tables(OFF / case / f"{record}.txt", base, label=f"{case}/{record}")


def test_top_n_and_all_human_are_not_exposed(wrapper_module):
    """Both upstream options are documented exclusions; they must not reappear as parameters."""
    import inspect

    target = wrapper_module.screen_sirna_offtarget
    signature = inspect.signature(getattr(target, "fn", target))
    assert "top_n" not in signature.parameters
    assert "all_human" not in signature.parameters
    source = Path(wrapper_module.__file__).read_text()
    assert "-top" not in source and "--top_n" not in source
    assert "all_human" not in source and '"-a"' not in source


@pytest.mark.slow
async def test_reference_case_matches_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "include_toxicity": True,
            "output_dir": str(out_root),
        },
    )
    # 3 result tables + pita.tab + targetscan.tab per record
    assert_tool_contract(data, min_artifacts=10, within=out_root)
    assert data["records"] == list(RECORDS)
    assert Path(data["utr_reference"]).name == "100_human_UTR.txt"
    assert Path(data["orf_reference"]).name == "100_human_ORF.txt"
    for record in RECORDS:
        summary = data["per_record_summary"][record]
        assert summary["rows"] == 39
        assert summary["columns"] == EXPECTED_COLUMNS
        assert summary["pita_rows"] == EXPECTED_PITA_ROWS[record]
        assert summary["pita_score_missing"] == 0
        assert summary["targetscan_score_missing"] == 0
        compare_result_tables(
            result_table(data, record), OFF / "run1" / f"{record}.txt", label=record
        )
        # the raw off-target tables must be byte-identical to the direct upstream ones
        for name in ("pita.tab", "targetscan.tab"):
            produced = Path(data["output_dir"]) / "offtarget" / record / name
            assert sha256(produced) == sha256(OFF / "run1" / "infer" / record / name), (
                f"{record}/{name} differs from the direct upstream table"
            )
    assert any(a["description"].endswith("raw PITA ddG table") for a in data["artifacts"])
    assert any(
        a["description"].endswith("raw TargetScan context score table") for a in data["artifacts"]
    )


@pytest.mark.slow
async def test_changed_thresholds_match_direct_upstream(mcp_server, out_root):
    data = await call_tool(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "include_toxicity": True,
            "pita_threshold": -5.0,
            "targetscan_threshold": 0.2,
            "output_dir": str(out_root),
        },
    )
    assert "-5.0" in data["message"] and "0.2" in data["message"]
    for record in RECORDS:
        compare_result_tables(
            result_table(data, record),
            OFF / "changed_params" / f"{record}.txt",
            label=f"changed/{record}",
        )
        # the thresholds genuinely changed the flag relative to the defaults
        _, default_rows = read_table(OFF / "run1" / f"{record}.txt")
        _, changed_rows = read_table(result_table(data, record))
        assert [r["off_target_filter"] for r in default_rows] != [
            r["off_target_filter"] for r in changed_rows
        ]
        # and the underlying scores did not change: only the thresholds did
        assert [r["pita_score"] for r in default_rows] == [r["pita_score"] for r in changed_rows]
        assert [r["targetscan_score"] for r in default_rows] == [
            r["targetscan_score"] for r in changed_rows
        ]
        # off_target_filter is exactly upstream's rule on those scores
        for row in changed_rows:
            expected = 1 if (float(row["pita_score"]) < -5.0 or float(row["targetscan_score"]) > 0.2) else 0
            assert int(row["off_target_filter"]) == expected


# ------------------------------------------------------------------ cheap failure paths (no run)


@pytest.mark.parametrize("which", ["utr_reference_path", "orf_reference_path"])
async def test_empty_reference_file_is_rejected(mcp_server, out_root, which):
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            which: str(DATA / "broken_utr_reference.txt"),
            "output_dir": str(out_root),
        },
    )
    assert "missing or empty" in text


@pytest.mark.parametrize("which", ["utr_reference_path", "orf_reference_path"])
async def test_missing_reference_file_is_rejected(mcp_server, out_root, which, tmp_path):
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            which: str(tmp_path / "absent.txt"),
            "output_dir": str(out_root),
        },
    )
    assert "missing or empty" in text


async def test_reference_path_with_spaces_is_rejected_with_an_explanation(mcp_server, out_root, tmp_path):
    """infer.py interpolates --utr/--orf into os.system, so a space must be refused, not run."""
    spaced = tmp_path / "ref dir" / "100_human_UTR.txt"
    spaced.parent.mkdir()
    spaced.write_text(">x\nACGU\n")
    text = await call_tool_expect_error(
        mcp_server,
        TOOL,
        {
            "mrna_fasta_path": str(EXAMPLE),
            "utr_reference_path": str(spaced),
            "output_dir": str(out_root),
        },
    )
    assert "passed unquoted to a shell" in text


@pytest.mark.parametrize(
    "fixture,fragment",
    [
        ("missing.fa", "mRNA FASTA not found"),
        ("bad_header_space.fa", "Unsupported FASTA record name"),
        ("bad_short_record.fa", "requires at least 19 nt"),
        ("bad_alphabet.fa", "unsupported characters"),
    ],
)
async def test_invalid_mrna_input_raises_mcp_error(mcp_server, out_root, fixture, fragment):
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(DATA / fixture), "output_dir": str(out_root)}
    )
    assert fragment in text


async def test_missing_offtarget_toolchain_is_reported(mcp_server, out_root, monkeypatch, tmp_path):
    monkeypatch.setenv("OLIGOFORMER_OFFTARGET_ENV", str(tmp_path / "no-conda-env"))
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)}
    )
    assert "Off-target toolchain not found" in text
    assert "RNAplfold" in text


async def test_wrong_runtime_directory_is_reported(mcp_server, out_root, monkeypatch, tmp_path):
    empty = tmp_path / "empty-runtime"
    empty.mkdir()
    monkeypatch.setenv("OLIGOFORMER_RUNTIME", str(empty))
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "output_dir": str(out_root)}
    )
    assert "not a complete OligoFormer runtime" in text


def _inject_offtarget_run(monkeypatch, wrapper_module, *, pita: str, targetscan: str, blank_scores: bool):
    """Replace the (multi-minute) CLI call with fabricated upstream-shaped output, so that the
    wrapper's post-run validation of the silently-failing off-target scripts can be exercised."""
    import subprocess

    def fake_run(runtime, args, **kwargs):
        out = Path(args[args.index("--output_dir") + 1])
        for record in RECORDS:
            for suffix in ("", "_ranked", "_ranked_filtered"):
                text = (OFF / "run1" / f"{record}{suffix}.txt").read_text()
                if blank_scores:
                    lines = text.splitlines()
                    header = lines[0].split("\t")
                    index = header.index("pita_score")
                    rewritten = [lines[0]]
                    for line in lines[1:]:
                        cells = line.split("\t")
                        cells[index] = ""
                        rewritten.append("\t".join(cells))
                    text = "\n".join(rewritten) + "\n"
                (out / f"{record}{suffix}.txt").write_text(text)
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_collect(runtime, names, dest, files):
        copied = []
        for name in names:
            directory = dest / name
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "pita.tab").write_text(pita)
            (directory / "targetscan.tab").write_text(targetscan)
            copied += [directory / "pita.tab", directory / "targetscan.tab"]
        return copied

    monkeypatch.setattr(wrapper_module, "_run_cli", fake_run)
    monkeypatch.setattr(wrapper_module, "_collect_infer_scratch", fake_collect)
    monkeypatch.setattr(wrapper_module, "_cleanup_infer_scratch", lambda runtime, names: None)


async def test_header_only_pita_table_is_detected(mcp_server, out_root, monkeypatch, wrapper_module):
    """pita.sh runs through os.system; an empty PITA table must not pass as a successful screen."""
    _inject_offtarget_run(
        monkeypatch,
        wrapper_module,
        pita="RefSeq\tmicroRNA\tSites\tScore\n",
        targetscan="ENST1\tRNA0\t0.1\n",
        blank_scores=False,
    )
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "include_toxicity": True, "output_dir": str(out_root)}
    )
    assert "is not usable" in text


async def test_empty_targetscan_table_is_detected(mcp_server, out_root, monkeypatch, wrapper_module):
    _inject_offtarget_run(
        monkeypatch,
        wrapper_module,
        pita="RefSeq\tmicroRNA\tSites\tScore\nENST1\tRNA0\t3\t-12.0\n",
        targetscan="",
        blank_scores=False,
    )
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "include_toxicity": True, "output_dir": str(out_root)}
    )
    assert "is empty" in text and "TargetScan" in text


async def test_all_nan_scores_are_detected(mcp_server, out_root, monkeypatch, wrapper_module):
    """All-NaN scores let every siRNA pass the off-target filter; that must be an error."""
    _inject_offtarget_run(
        monkeypatch,
        wrapper_module,
        pita="RefSeq\tmicroRNA\tSites\tScore\nENST1\tRNA0\t3\t-12.0\n",
        targetscan="ENST1\tRNA0\t0.1\n",
        blank_scores=True,
    )
    text = await call_tool_expect_error(
        mcp_server, TOOL, {"mrna_fasta_path": str(EXAMPLE), "include_toxicity": True, "output_dir": str(out_root)}
    )
    assert "Every pita_score value is missing" in text


async def test_schema_documents_every_parameter(mcp_server):
    from fastmcp import Client

    async with Client(mcp_server) as client:
        tools = {t.name: t for t in await client.list_tools()}
    schema = tools[TOOL].inputSchema
    assert schema["required"] == ["mrna_fasta_path"]
    for name, prop in schema["properties"].items():
        assert prop.get("description"), f"parameter {name} has no Annotated description"
    assert schema["properties"]["pita_threshold"]["default"] == -10.0
    assert schema["properties"]["targetscan_threshold"]["default"] == 1.0
    assert tools[TOOL].description.strip().count("\n") == 1

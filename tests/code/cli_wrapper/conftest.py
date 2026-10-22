"""Shared fixtures/comparators for the independent verification of src/tools/cli_wrapper.py.

Expectations come from direct upstream executions only:
  * the stage-2 executors' saved outputs under notebooks/<id>/outputs/
  * this verifier's own direct upstream run under tests/results/cli_wrapper/upstream_only_RNA0/
Nothing here imports the wrapper to derive an expectation.
"""
from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DATA = PROJECT_ROOT / "tests" / "data" / "cli_wrapper"
NOTEBOOKS = PROJECT_ROOT / "notebooks"
RESULTS = PROJECT_ROOT / "tests" / "results" / "cli_wrapper"
VERIFY_OUT = PROJECT_ROOT / "tmp" / "outputs" / "verify"

# Columns that upstream computes deterministically (no dropout path): compared exactly.
DETERMINISTIC_COLUMNS = (
    "pos",
    "sense",
    "siRNA",
    "func_filter",
    "Seed",
    "cell_viability",
    "toxicity_filter",
    "pita_score",
    "targetscan_score",
    "off_target_filter",
    "filter",
)

# scripts/infer.py never calls model.eval() and never seeds torch, so nn.Dropout(0.02) stays
# active during inference and the efficacy column is stochastic. Executor evidence (3 runs per
# assignment, plus a replay) bounded the run-to-run difference at ~0.18 absolute with Spearman
# >= 0.93, so a per-row tolerance of 0.25 and a per-record Spearman floor of 0.9 are used.
EFFICACY_ABS_TOL = 0.25
EFFICACY_SPEARMAN_MIN = 0.9
SPEARMAN_MIN_ROWS = 5


def read_table(path: Path) -> tuple[list[str], list[dict]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return [], []
    header = rows[0]
    return header, [dict(zip(header, r)) for r in rows[1:] if r]


def _as_float(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _cell_equal(actual: str, expected: str) -> bool:
    if actual == expected:
        return True
    fa, fe = _as_float(actual), _as_float(expected)
    if fa is None or fe is None:
        return False
    if math.isnan(fa) and math.isnan(fe):
        return True
    return abs(fa - fe) <= 1e-9


def _spearman(xs: list[float], ys: list[float]) -> float:
    from scipy.stats import spearmanr

    return float(spearmanr(xs, ys).statistic)


def compare_result_tables(
    actual_path: Path,
    expected_path: Path,
    *,
    efficacy_abs_tol: float = EFFICACY_ABS_TOL,
    spearman_min: float = EFFICACY_SPEARMAN_MIN,
    label: str = "",
) -> dict:
    """Compare a wrapper result table against a direct upstream reference table.

    Header, row count, row identity and every deterministic column must match exactly;
    the stochastic efficacy column is compared with a tolerance and a rank correlation.
    """
    tag = label or f"{actual_path.name}"
    ha, ra = read_table(actual_path)
    he, re_rows = read_table(expected_path)
    assert he, f"{tag}: reference table {expected_path} is empty"
    assert ha == he, f"{tag}: column names differ: {ha} != {he}"
    assert len(ra) == len(re_rows), f"{tag}: row count {len(ra)} != reference {len(re_rows)}"

    key = lambda row: (row["pos"], row["siRNA"])
    assert [key(r) for r in ra] == [key(r) for r in re_rows], (
        f"{tag}: (pos, siRNA) rows differ from the upstream reference, in this order"
    )

    compared = [c for c in DETERMINISTIC_COLUMNS if c in ha]
    for i, (row_a, row_e) in enumerate(zip(ra, re_rows), start=1):
        for column in compared:
            assert _cell_equal(row_a[column], row_e[column]), (
                f"{tag}: row {i} column {column!r}: {row_a[column]!r} != reference {row_e[column]!r}"
            )

    stats = {"rows": len(ra), "deterministic_columns": compared}
    if "efficacy" in ha:
        ea = [float(r["efficacy"]) for r in ra]
        ee = [float(r["efficacy"]) for r in re_rows]
        assert all(math.isfinite(v) for v in ea), f"{tag}: non-finite efficacy in {actual_path}"
        deltas = [abs(x - y) for x, y in zip(ea, ee)]
        stats["max_abs_efficacy_delta"] = max(deltas)
        assert max(deltas) <= efficacy_abs_tol, (
            f"{tag}: max |efficacy delta| {max(deltas):.4f} exceeds the justified stochastic "
            f"tolerance {efficacy_abs_tol}"
        )
        if len(ea) >= SPEARMAN_MIN_ROWS:
            rho = _spearman(ea, ee)
            stats["efficacy_spearman"] = rho
            assert rho >= spearman_min, (
                f"{tag}: efficacy Spearman {rho:.4f} below the justified floor {spearman_min}"
            )
    return stats


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope="session")
def wrapper_module():
    from src.tools import cli_wrapper

    return cli_wrapper


@pytest.fixture(scope="session")
def mcp_server(wrapper_module):
    return wrapper_module.cli_wrapper_mcp


@pytest.fixture
def out_root(request, tmp_path_factory) -> Path:
    """A per-test output root under the verifier-owned tmp/outputs/verify tree."""
    safe = request.node.name.replace("/", "_").replace("[", "_").replace("]", "")
    root = VERIFY_OUT / safe
    root.mkdir(parents=True, exist_ok=True)
    return root


async def call_tool(server, name: str, arguments: dict):
    """Exercise a tool the way the MCP client does and return result.data."""
    from fastmcp import Client

    async with Client(server) as client:
        result = await client.call_tool(name, arguments)
    assert not result.is_error, f"{name} returned an MCP error: {result.content}"
    return result.data


async def call_tool_expect_error(server, name: str, arguments: dict) -> str:
    """Assert the tool call surfaces as an MCP tool error and return its text."""
    from fastmcp import Client
    from fastmcp.exceptions import ToolError

    async with Client(server) as client:
        try:
            result = await client.call_tool(name, arguments)
        except ToolError as exc:
            return str(exc)
        except Exception as exc:  # pragma: no cover - any transport-level error is also a failure
            return f"{type(exc).__name__}: {exc}"
    assert result.is_error, (
        f"{name} unexpectedly succeeded with {arguments!r}; an MCP tool error was required "
        f"(got data={getattr(result, 'data', None)!r})"
    )
    return str(result.content)


def artifact_paths(data: dict) -> list[Path]:
    return [Path(a["path"]) for a in data["artifacts"]]


def assert_tool_contract(data: dict, *, min_artifacts: int, within: Path | None = None) -> None:
    """Every tool must return the documented message/reference/artifacts contract."""
    assert isinstance(data, dict)
    for field in ("message", "reference", "artifacts"):
        assert field in data, f"missing return field {field!r}"
    assert data["message"].strip()
    assert data["reference"].startswith("https://github.com/lulab/OligoFormer/blob/")
    assert isinstance(data["artifacts"], list)
    assert len(data["artifacts"]) >= min_artifacts, (
        f"expected >= {min_artifacts} artifacts, got {len(data['artifacts'])}"
    )
    seen = set()
    for entry in data["artifacts"]:
        assert set(entry) >= {"description", "path"}, entry
        path = Path(entry["path"])
        assert path.is_absolute(), f"artifact path is not absolute: {path}"
        assert path.exists(), f"artifact path does not exist: {path}"
        assert str(path) not in seen, f"duplicate artifact path {path}"
        seen.add(str(path))
        if within is not None:
            assert within.resolve() in path.resolve().parents or path.resolve() == within.resolve(), (
                f"artifact {path} escaped the requested output root {within}"
            )

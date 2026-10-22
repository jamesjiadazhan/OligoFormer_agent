"""Tools extracted from the OligoFormer command-line interface (repo/OligoFormer/scripts/main.py)."""
from __future__ import annotations

import csv
import fcntl
import os
import re
import shutil
import subprocess
import sys
import threading
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Literal

from fastmcp import FastMCP
from fastmcp.exceptions import ToolError

COMMIT = "e2f53ad63387bbe166bf123949151e2bc9bf6ec3"
_BLOB = f"https://github.com/lulab/OligoFormer/blob/{COMMIT}"
REFERENCE = f"{_BLOB}/scripts/main.py"
REFERENCE_INFER = f"{_BLOB}/scripts/infer.py"
REFERENCE_MISMATCH = f"{_BLOB}/scripts/mismatch.py"

cli_wrapper_mcp = FastMCP(name="cli_wrapper")

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# scripts/infer.py runs the model without model.eval() and with nn.Dropout(0.02) layers active,
# and never seeds torch, so the efficacy column is not reproducible run to run (stage-2 executions
# observed up to 0.18 absolute / 31% relative differences). Upstream behaviour, not repaired here.
STOCHASTIC_NOTE = (
    "Upstream OligoFormer inference leaves dropout active (scripts/infer.py never calls "
    "model.eval() and does not seed torch), so the efficacy column is stochastic: repeated runs "
    "on the same input differ by up to ~0.18 absolute (~30% relative) and the ranked order "
    "changes. All other columns (pos, sense, siRNA, func_filter, pita/targetscan scores, Seed, "
    "cell_viability, toxicity_filter, filter) are deterministic."
)

_MRNA_HEADER_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_SIRNA_RE = re.compile(r"^[ACGU]{19}$")
_MRNA_ALPHABET = set("ACGTU")
_PITA_HEADER = ["RefSeq", "microRNA", "Sites", "Score"]

_PROCESS_LOCK = threading.Lock()
_LOG_TAIL_CHARS = 4000

# Timeouts (seconds). Reference runs: ~45 s for 78 windows without off-target, ~7 min with
# off-target, ~95 s for the 3-mismatch scan. The per-window terms keep longer transcripts inside
# the limit while still bounding a hung child process.
_INFER_TIMEOUT_BASE = 900
_INFER_TIMEOUT_PER_WINDOW = 30
_OFFTARGET_TIMEOUT_BASE = 3600
_OFFTARGET_TIMEOUT_PER_WINDOW = 60
_MISMATCH_PAIR_TIMEOUT = 600
_MISMATCH_SCAN_TIMEOUT = 3600


# --------------------------------------------------------------------------------------- runtime


def _runtime_root() -> Path:
    """Resolve and cheaply identify the OligoFormer runtime copy the child process runs in."""
    raw = os.environ.get("OLIGOFORMER_RUNTIME") or str(PROJECT_ROOT / "runtime" / "OligoFormer")
    root = Path(raw).expanduser()
    if not root.is_dir():
        raise ToolError(
            f"OligoFormer runtime directory not found: {root}. Build one with "
            "deps/setup_oligoformer_runtime.sh <target> or set OLIGOFORMER_RUNTIME."
        )
    root = root.resolve()
    main_py = root / "scripts" / "main.py"
    required = [
        main_py,
        root / "scripts" / "infer.py",
        root / "scripts" / "mismatch.py",
        root / "scripts" / "RNA-FM.sh",
        root / "model" / "best_model.pth",
        root / "model" / "mismatch_model.pth",
        root / "toxicity" / "cell_viability.txt",
        root / "RNA-FM" / "torch_home" / "hub" / "checkpoints" / "RNA-FM_pretrained.pth",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise ToolError(
            f"{root} is not a complete OligoFormer runtime; missing: {', '.join(missing)}. "
            "Rebuild it with deps/setup_oligoformer_runtime.sh."
        )
    head = main_py.read_text(errors="replace")[:8000]
    if "OligoFormer" not in head or "--infer" not in head:
        raise ToolError(
            f"{main_py} does not look like the OligoFormer CLI entry point (no OligoFormer "
            "argparse description / --infer option); refusing to run an unknown program."
        )
    return root


def _offtarget_bin() -> Path:
    """Resolve the conda environment that provides perl/RNAplfold/GNU grep for the off-target run."""
    raw = os.environ.get("OLIGOFORMER_OFFTARGET_ENV") or str(PROJECT_ROOT / "deps" / "offtarget-env")
    env_root = Path(raw).expanduser()
    bin_dir = env_root / "bin"
    if not (bin_dir / "perl").exists() or not (bin_dir / "RNAplfold").exists():
        raise ToolError(
            f"Off-target toolchain not found in {bin_dir} (needs perl with BioPerl/Statistics::Lite "
            "and RNAplfold). Create it with deps/create_offtarget_env.sh or set "
            "OLIGOFORMER_OFFTARGET_ENV."
        )
    return bin_dir.resolve()


def _child_env(runtime: Path, with_offtarget: bool) -> dict:
    """Build the child environment: project python first on PATH, TORCH_HOME set, PERL5LIB dropped."""
    env = dict(os.environ)
    # Not resolve(): the venv's bin/python is a symlink to the base interpreter, whose directory
    # does not contain the installed torch/RNA-FM dependencies that scripts/RNA-FM.sh needs.
    parts = [os.path.dirname(os.path.abspath(sys.executable))]
    if with_offtarget:
        parts.append(str(_offtarget_bin()))
    existing = env.get("PATH", "")
    env["PATH"] = os.pathsep.join(parts + ([existing] if existing else []))
    # scripts/RNA-FM.sh loads the checkpoint through torch.hub's cache; no download happens.
    env["TORCH_HOME"] = str(runtime / "RNA-FM" / "torch_home")
    # An inherited PERL5LIB can shadow the conda perl's BioPerl/Statistics::Lite.
    env.pop("PERL5LIB", None)
    return env


@contextmanager
def _runtime_lock(runtime: Path):
    """Serialize runs: upstream writes fixed paths (./data/infer, ./result, off-target/tmp)."""
    lock_path = runtime / ".cli_wrapper.lock"
    with _PROCESS_LOCK:
        handle = open(lock_path, "a+")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


# ------------------------------------------------------------------------------------ utilities


def _fresh_output_dir(output_dir: str | None, prefix: str) -> Path:
    """Create a collision-resistant directory for this invocation and return it."""
    base = Path(output_dir).expanduser().resolve() if output_dir else PROJECT_ROOT / "tmp" / "outputs" / "cli_wrapper"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir = Path(base) / f"{prefix}_{stamp}_{uuid.uuid4().hex[:8]}"
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        raise ToolError(f"Could not create output directory {run_dir}: {exc}") from exc
    return run_dir


def _tail(text: str) -> str:
    text = (text or "").strip()
    return text if len(text) <= _LOG_TAIL_CHARS else "...(truncated)...\n" + text[-_LOG_TAIL_CHARS:]


def _run_cli(
    runtime: Path,
    args: list[str],
    *,
    timeout: int,
    with_offtarget: bool = False,
    stdin_text: str | None = None,
    log_dir: Path | None = None,
    label: str = "main",
) -> subprocess.CompletedProcess:
    """Run `python scripts/main.py ...` from the runtime root with shell=False and check the status."""
    cmd = [sys.executable, "scripts/main.py", *args]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(runtime),
            env=_child_env(runtime, with_offtarget),
            input=stdin_text if stdin_text is not None else "",
            capture_output=True,
            text=True,
            timeout=timeout,
            shell=False,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolError(
            f"OligoFormer CLI timed out after {timeout}s: {' '.join(cmd)}\n"
            f"stdout tail:\n{_tail(exc.stdout.decode() if isinstance(exc.stdout, bytes) else exc.stdout or '')}\n"
            f"stderr tail:\n{_tail(exc.stderr.decode() if isinstance(exc.stderr, bytes) else exc.stderr or '')}"
        ) from exc
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        (log_dir / f"{label}.stdout.log").write_text(proc.stdout or "")
        (log_dir / f"{label}.stderr.log").write_text(proc.stderr or "")
        (log_dir / f"{label}.cmd").write_text(" ".join(cmd) + f"\n# cwd={runtime}\n")
    if proc.returncode != 0:
        raise ToolError(
            f"OligoFormer CLI failed with exit status {proc.returncode}: {' '.join(cmd)}\n"
            f"stdout tail:\n{_tail(proc.stdout)}\nstderr tail:\n{_tail(proc.stderr)}"
        )
    return proc


def _anti_rna(seq: str) -> str:
    """Reverse complement used only to locate a user siRNA's site, mirroring infer.py:antiRNA."""
    table = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return "".join(table[c] for c in reversed(seq))


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:]
            chunks = []
        elif line:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _load_mrna_fasta(mrna_fasta_path: str) -> tuple[Path, list[tuple[str, str]]]:
    """Validate the mRNA FASTA against the upstream contract and shell-safety requirements."""
    path = Path(mrna_fasta_path).expanduser()
    if not path.is_file():
        raise ToolError(f"mRNA FASTA not found: {path}")
    path = path.resolve()
    records = _read_fasta(path)
    if not records:
        raise ToolError(f"No FASTA records found in {path} (expected '>name' headers with sequence lines).")
    seen: set[str] = set()
    for name, seq in records:
        if not _MRNA_HEADER_RE.match(name) or name in {".", ".."}:
            raise ToolError(
                f"Unsupported FASTA record name {name!r} in {path}. Upstream interpolates record "
                "names unquoted into shell commands and mangles spaces, so names must match "
                "^[A-Za-z0-9_.-]+$ (no spaces or other characters) and must be a single path "
                "component ('.' and '..' are rejected: upstream writes ./data/infer/<name>/, so "
                "they would escape the per-record scratch directory)."
            )
        if name in seen:
            raise ToolError(f"Duplicate FASTA record name {name!r} in {path}; names must be unique.")
        seen.add(name)
        norm = seq.upper()
        if len(norm) < 19:
            raise ToolError(
                f"Record {name!r} in {path} is {len(norm)} nt; OligoFormer requires at least 19 nt."
            )
        bad = sorted(set(norm) - _MRNA_ALPHABET)
        if bad:
            raise ToolError(
                f"Record {name!r} in {path} contains unsupported characters {bad}. Only A/C/G/T/U "
                "are supported (upstream silently drops other characters, which would shift the "
                "window sequences)."
            )
    return path, records


def _load_sirna_fasta(sirna_fasta_path: str) -> tuple[Path, list[str]]:
    """Validate the siRNA FASTA: upstream requires exactly 19 nt per line and no case/T conversion."""
    path = Path(sirna_fasta_path).expanduser()
    if not path.is_file():
        raise ToolError(f"siRNA FASTA not found: {path}")
    path = path.resolve()
    sirnas: list[str] = []
    lines = path.read_text().splitlines()
    for lineno, raw in enumerate(lines, start=1):
        if raw.startswith(">"):
            continue
        # Blank lines are not skipped: upstream reads every non-header line and raises the 19-nt
        # exception on an empty one, so they are reported here with a clearer message.
        if not _SIRNA_RE.match(raw.strip()):
            raise ToolError(
                f"Invalid siRNA on line {lineno} of {path}: {raw.strip()!r}. Each sequence line must "
                "be exactly 19 nt of A/C/G/U. Upstream raises on any other length and does not "
                "upper-case or convert T to U for siRNA input, so lowercase, DNA 'T' and blank "
                "lines are rejected here instead of being silently altered."
            )
        sirnas.append(raw.strip())
    if not sirnas:
        raise ToolError(f"No siRNA sequences found in {path}.")
    return path, sirnas


def _validate_19nt_rna(value: str, label: str) -> str:
    if not _SIRNA_RE.match(value or ""):
        raise ToolError(
            f"{label} must be exactly 19 nt of A/C/G/U (got {value!r}). The upstream mismatch model "
            "one-hot encoder raises a bare KeyError for lowercase letters, DNA 'T' or any other "
            "character, so the input is validated rather than converted."
        )
    return value


def _read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], [r for r in rows[1:] if r]


def _table_summary(result_path: Path, ranked_filtered_path: Path) -> dict:
    header, rows = _read_table(result_path)
    summary: dict = {"rows": len(rows), "columns": header}
    index = {name: i for i, name in enumerate(header)}
    if "func_filter" in index:
        counts: dict[str, int] = {}
        for row in rows:
            key = row[index["func_filter"]]
            counts[key] = counts.get(key, 0) + 1
        summary["func_filter_counts"] = dict(sorted(counts.items()))
    if "filter" in index:
        summary["rows_passing_filter"] = sum(1 for row in rows if row[index["filter"]] == "0")
    for column in ("pita_score", "targetscan_score", "cell_viability"):
        if column in index:
            values = [row[index[column]] for row in rows]
            summary[f"{column}_missing"] = sum(1 for v in values if v.strip() == "" or v.strip().lower() == "nan")
    if ranked_filtered_path.exists():
        _, filtered_rows = _read_table(ranked_filtered_path)
        summary["ranked_filtered_rows"] = len(filtered_rows)
    return summary


def _record_artifacts(run_dir: Path, names: list[str]) -> tuple[list[dict], dict]:
    """Check that the upstream result tables were written and collect them as artifacts."""
    artifacts: list[dict] = []
    summaries: dict = {}
    for name in names:
        result_path = run_dir / f"{name}.txt"
        ranked_path = run_dir / f"{name}_ranked.txt"
        filtered_path = run_dir / f"{name}_ranked_filtered.txt"
        if not result_path.is_file():
            raise ToolError(
                f"Expected result table {result_path} was not produced for record {name!r}. The "
                "upstream CLI ignores the exit status of its RNA-FM and off-target sub-processes, "
                "so a zero exit does not prove the run completed."
            )
        header, rows = _read_table(result_path)
        if not rows:
            raise ToolError(f"Result table {result_path} for record {name!r} has no data rows.")
        artifacts.append({"description": f"{name}: predictions in positional order", "path": str(result_path)})
        if ranked_path.is_file():
            artifacts.append({"description": f"{name}: predictions ranked by efficacy", "path": str(ranked_path)})
        if filtered_path.is_file():
            artifacts.append(
                {"description": f"{name}: ranked predictions passing all filters", "path": str(filtered_path)}
            )
        summaries[name] = _table_summary(result_path, filtered_path)
    return artifacts, summaries


def _collect_infer_scratch(runtime: Path, names: list[str], dest: Path, files: tuple[str, ...]) -> list[Path]:
    """Copy the per-record intermediates out of the fixed ./data/infer/<record> scratch directory."""
    copied: list[Path] = []
    for name in names:
        src_dir = runtime / "data" / "infer" / name
        for filename in files:
            src = src_dir / filename
            if src.is_file():
                out_dir = dest / name
                out_dir.mkdir(parents=True, exist_ok=True)
                target = out_dir / filename
                shutil.copy2(src, target)
                copied.append(target)
    return copied


def _cleanup_infer_scratch(runtime: Path, names: list[str]) -> None:
    for name in names:
        scratch = runtime / "data" / "infer" / name
        if scratch.is_dir():
            shutil.rmtree(scratch, ignore_errors=True)


def _infer_timeout(records: list[tuple[str, str]], off_target: bool) -> int:
    windows = sum(max(0, len(seq) - 19 + 1) for _, seq in records)
    if off_target:
        return int(_OFFTARGET_TIMEOUT_BASE + _OFFTARGET_TIMEOUT_PER_WINDOW * windows)
    return int(_INFER_TIMEOUT_BASE + _INFER_TIMEOUT_PER_WINDOW * windows)


def _filter_flags(apply_functionality_filter: bool, include_toxicity: bool, toxicity_threshold: float) -> list[str]:
    args: list[str] = []
    if not apply_functionality_filter:
        args.append("-nf")
    if include_toxicity:
        args += ["-tox", "--toxicity_threshold", repr(float(toxicity_threshold))]
    return args


# ---------------------------------------------------------------------------------------- tools


@cli_wrapper_mcp.tool()
def predict_sirna_efficacy_from_mrna(
    mrna_fasta_path: Annotated[str, "FASTA of target mRNA records (each >= 19 nt, header matching ^[A-Za-z0-9_.-]+$)"],
    apply_functionality_filter: Annotated[bool, "Apply the known functionality criteria filter (off adds -nf)"] = True,
    include_toxicity: Annotated[bool, "Add the seed-based cell-viability toxicity filter (-tox)"] = False,
    toxicity_threshold: Annotated[float, "Cell viability percent below which a seed is flagged toxic"] = 50.0,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created per call"] = None,
) -> dict:
    """Predict siRNA efficacy for every 19-nt antisense window of each mRNA record.
    Input is an mRNA FASTA; output is three tab-separated tables per record (all, ranked, ranked and filtered).
    """
    runtime = _runtime_root()
    fasta_path, records = _load_mrna_fasta(mrna_fasta_path)
    names = [name for name, _ in records]
    run_dir = _fresh_output_dir(output_dir, "predict_mrna")
    args = ["--infer", "1", "-i1", str(fasta_path), "--output_dir", f"{run_dir}{os.sep}"]
    args += _filter_flags(apply_functionality_filter, include_toxicity, toxicity_threshold)
    with _runtime_lock(runtime):
        try:
            _run_cli(runtime, args, timeout=_infer_timeout(records, False), log_dir=run_dir / "logs")
            _collect_infer_scratch(runtime, names, run_dir / "windows", ("siRNA.fa", "mRNA.fa"))
        finally:
            _cleanup_infer_scratch(runtime, names)
    artifacts, summaries = _record_artifacts(run_dir, names)
    windows = {name: max(0, len(seq) - 19 + 1) for name, seq in records}
    return {
        "message": (
            f"Scored {sum(windows.values())} 19-nt windows across {len(names)} record(s) with "
            f"OligoFormer (functionality filter {'on' if apply_functionality_filter else 'off'}, "
            f"toxicity filter {'on' if include_toxicity else 'off'}). {STOCHASTIC_NOTE}"
        ),
        "reference": REFERENCE_INFER,
        "artifacts": artifacts,
        "records": names,
        "windows_per_record": windows,
        "per_record_summary": summaries,
        "output_dir": str(run_dir),
        "stochastic_efficacy": True,
    }


@cli_wrapper_mcp.tool()
def predict_specific_sirna_efficacy(
    mrna_fasta_path: Annotated[str, "FASTA of target mRNA records (each >= 19 nt, header matching ^[A-Za-z0-9_.-]+$)"],
    sirna_fasta_path: Annotated[str, "FASTA of candidate siRNAs, each exactly 19 nt of A/C/G/U"],
    apply_functionality_filter: Annotated[bool, "Apply the known functionality criteria filter (off adds -nf)"] = True,
    include_toxicity: Annotated[bool, "Add the seed-based cell-viability toxicity filter (-tox)"] = False,
    toxicity_threshold: Annotated[float, "Cell viability percent below which a seed is flagged toxic"] = 50.0,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created per call"] = None,
) -> dict:
    """Predict efficacy only for user-supplied 19-nt siRNAs whose site occurs in the target mRNA.
    Input is an mRNA FASTA plus an siRNA FASTA; output is the same three tables per matched record.
    """
    runtime = _runtime_root()
    fasta_path, records = _load_mrna_fasta(mrna_fasta_path)
    sirna_path, sirnas = _load_sirna_fasta(sirna_fasta_path)

    # Upstream crashes (FileNotFoundError in loader.py) on a record with zero matching siRNAs and
    # aborts every later record, so pre-scan matches with the same rule as infer.py
    # (re.search(antiRNA(siRNA), mRNA) after .upper().replace('T','U')) and pass only records that
    # have at least one match. No scientific change: matching and scoring stay upstream's.
    matched: dict[str, list[str]] = {}
    for name, seq in records:
        norm = seq.upper().replace("T", "U")
        hits = [s for s in sirnas if re.search(_anti_rna(s), norm) is not None]
        matched[name] = hits
    kept = [(name, seq) for name, seq in records if matched[name]]
    unmatched_records = [name for name, _ in records if not matched[name]]
    unmatched_sirnas = [s for s in sirnas if not any(s in hits for hits in matched.values())]
    if not kept:
        raise ToolError(
            f"None of the {len(sirnas)} siRNA(s) in {sirna_path} match any record in {fasta_path}. "
            "A match requires the reverse complement of the siRNA to occur in the mRNA sequence."
        )

    run_dir = _fresh_output_dir(output_dir, "predict_sirnas")
    if unmatched_records:
        # Fresh FASTA with the retained records only; sequences are copied verbatim so upstream
        # applies its own normalisation.
        filtered_fasta = run_dir / "matched_records.fa"
        filtered_fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in kept))
        infer_fasta = filtered_fasta
    else:
        infer_fasta = fasta_path
    names = [name for name, _ in kept]
    args = [
        "--infer", "1",
        "-i1", str(infer_fasta),
        "-i2", str(sirna_path),
        "--output_dir", f"{run_dir}{os.sep}",
    ]
    args += _filter_flags(apply_functionality_filter, include_toxicity, toxicity_threshold)
    with _runtime_lock(runtime):
        try:
            _run_cli(runtime, args, timeout=_infer_timeout(kept, False), log_dir=run_dir / "logs")
            _collect_infer_scratch(runtime, names, run_dir / "windows", ("siRNA.fa", "mRNA.fa"))
        finally:
            _cleanup_infer_scratch(runtime, names)
    artifacts, summaries = _record_artifacts(run_dir, names)
    if infer_fasta != fasta_path:
        artifacts.append(
            {"description": "mRNA FASTA actually passed upstream (records with >= 1 siRNA match)",
             "path": str(infer_fasta)}
        )
    return {
        "message": (
            f"Scored {len(sirnas)} candidate siRNA(s) against {len(names)} matched record(s); "
            f"{len(unmatched_records)} record(s) and {len(unmatched_sirnas)} siRNA(s) had no match "
            "and were excluded (upstream drops unmatched siRNAs silently and crashes on a record "
            "with no match). Note: the 'pos' column is a row counter over matched siRNAs, not a "
            f"transcript position, and only the first site per record is scored. {STOCHASTIC_NOTE}"
        ),
        "reference": REFERENCE_INFER,
        "artifacts": artifacts,
        "records_scored": names,
        "matched_sirnas_per_record": {name: matched[name] for name in names},
        "unmatched_records": unmatched_records,
        "unmatched_sirnas": unmatched_sirnas,
        "per_record_summary": summaries,
        "output_dir": str(run_dir),
        "stochastic_efficacy": True,
    }


@cli_wrapper_mcp.tool()
def screen_sirna_offtarget(
    mrna_fasta_path: Annotated[str, "FASTA of target mRNA records (each >= 19 nt, header matching ^[A-Za-z0-9_.-]+$)"],
    apply_functionality_filter: Annotated[bool, "Apply the known functionality criteria filter (off adds -nf)"] = True,
    include_toxicity: Annotated[bool, "Add the seed-based cell-viability toxicity filter (-tox)"] = False,
    toxicity_threshold: Annotated[float, "Cell viability percent below which a seed is flagged toxic"] = 50.0,
    pita_threshold: Annotated[float, "PITA ddG below which a siRNA is flagged as off-target"] = -10.0,
    targetscan_threshold: Annotated[float, "TargetScan context score above which a siRNA is flagged"] = 1.0,
    utr_reference_path: Annotated[str | None, "Reference UTR table; defaults to the bundled 100_human_UTR.txt"] = None,
    orf_reference_path: Annotated[str | None, "Reference ORF table; defaults to the bundled 100_human_ORF.txt"] = None,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created per call"] = None,
) -> dict:
    """Screen every 19-nt siRNA window of an mRNA for off-target binding with PITA and TargetScan.
    Input is an mRNA FASTA plus reference UTR/ORF tables; output is the efficacy tables with pita_score, targetscan_score and off_target_filter, plus the raw PITA/TargetScan tables.
    """
    runtime = _runtime_root()
    _offtarget_bin()
    fasta_path, records = _load_mrna_fasta(mrna_fasta_path)
    names = [name for name, _ in records]

    def _reference(path: str | None, default: Path, label: str) -> Path:
        resolved = Path(path).expanduser().resolve() if path else default
        # infer.py interpolates these paths unquoted into os.system('bash scripts/pita.sh ...'),
        # so a space or shell metacharacter would break (or redirect) the off-target run.
        if not re.fullmatch(r"[A-Za-z0-9_./+=:-]+", str(resolved)):
            raise ToolError(
                f"{label} reference path {resolved} contains characters that the upstream "
                "off-target scripts cannot handle: the path is passed unquoted to a shell. Use a "
                "path made of letters, digits and _ . / + = : - only (no spaces)."
            )
        if not resolved.is_file() or resolved.stat().st_size == 0:
            raise ToolError(
                f"{label} reference file {resolved} is missing or empty. The upstream off-target "
                "scripts ignore failures, which silently produces all-NaN scores and lets every "
                "siRNA pass the off-target filter."
            )
        return resolved

    utr = _reference(utr_reference_path, runtime / "off-target" / "ref" / "100_human_UTR.txt", "UTR")
    orf = _reference(orf_reference_path, runtime / "off-target" / "ref" / "100_human_ORF.txt", "ORF")

    run_dir = _fresh_output_dir(output_dir, "offtarget")
    args = [
        "--infer", "1",
        "-i1", str(fasta_path),
        "-off",
        "--utr", str(utr),
        "--orf", str(orf),
        "--pita_threshold", repr(float(pita_threshold)),
        "--targetscan_threshold", repr(float(targetscan_threshold)),
        "--output_dir", f"{run_dir}{os.sep}",
    ]
    args += _filter_flags(apply_functionality_filter, include_toxicity, toxicity_threshold)
    with _runtime_lock(runtime):
        try:
            _run_cli(
                runtime,
                args,
                timeout=_infer_timeout(records, True),
                with_offtarget=True,
                log_dir=run_dir / "logs",
            )
            _collect_infer_scratch(
                runtime, names, run_dir / "offtarget", ("pita.tab", "targetscan.tab", "siRNA.fa")
            )
        finally:
            _cleanup_infer_scratch(runtime, names)

    artifacts, summaries = _record_artifacts(run_dir, names)
    for name in names:
        pita_path = run_dir / "offtarget" / name / "pita.tab"
        ts_path = run_dir / "offtarget" / name / "targetscan.tab"
        if not pita_path.is_file() or not ts_path.is_file():
            raise ToolError(
                f"Off-target tables are missing for record {name!r} ({pita_path}, {ts_path}). "
                "pita.sh/targetscan.sh run through os.system and their failures are ignored "
                "upstream, so the CLI's exit status cannot be trusted here."
            )
        pita_header, pita_rows = _read_table(pita_path)
        if pita_header[: len(_PITA_HEADER)] != _PITA_HEADER or not pita_rows:
            raise ToolError(
                f"PITA table {pita_path} for record {name!r} is not usable (header {pita_header}, "
                f"{len(pita_rows)} data rows). Expected header {_PITA_HEADER} and at least one row; "
                "an empty table means the PITA run failed silently."
            )
        if ts_path.stat().st_size == 0:
            raise ToolError(
                f"TargetScan table {ts_path} for record {name!r} is empty; the TargetScan run "
                "failed silently."
            )
        summary = summaries[name]
        for column in ("pita_score", "targetscan_score"):
            if summary.get(f"{column}_missing") == summary["rows"]:
                raise ToolError(
                    f"Every {column} value is missing for record {name!r}. The off-target scripts "
                    "failed without reporting it, so all siRNAs would falsely pass the off-target "
                    "filter."
                )
        summary["pita_rows"] = len(pita_rows)
        artifacts.append({"description": f"{name}: raw PITA ddG table", "path": str(pita_path)})
        artifacts.append({"description": f"{name}: raw TargetScan context score table", "path": str(ts_path)})

    return {
        "message": (
            f"Off-target screening finished for {len(names)} record(s) using PITA (threshold "
            f"{pita_threshold}) and TargetScan (threshold {targetscan_threshold}) against "
            f"{utr.name} / {orf.name}; PITA and TargetScan tables were verified to be populated. "
            f"{STOCHASTIC_NOTE}"
        ),
        "reference": REFERENCE_INFER,
        "artifacts": artifacts,
        "records": names,
        "per_record_summary": summaries,
        "utr_reference": str(utr),
        "orf_reference": str(orf),
        "output_dir": str(run_dir),
        "stochastic_efficacy": True,
    }


@cli_wrapper_mcp.tool()
def compare_sirna_mismatch_kon(
    base_rna: Annotated[str, "Base 19-nt RNA sequence (A/C/G/U only)"],
    mismatched_rna: Annotated[str, "Mismatched 19-nt RNA sequence to compare (A/C/G/U only)"],
) -> dict:
    """Score the relative on-rate (kon) of a mismatched siRNA against its base sequence.
    Input is two 19-nt RNA sequences; output is the single kon value printed by the mismatch module.
    """
    runtime = _runtime_root()
    base = _validate_19nt_rna(base_rna, "base_rna")
    mismatched = _validate_19nt_rna(mismatched_rna, "mismatched_rna")
    with _runtime_lock(runtime):
        proc = _run_cli(
            runtime,
            ["-m", "1"],
            timeout=_MISMATCH_PAIR_TIMEOUT,
            stdin_text=f"{base}\n{mismatched}\n",
        )
    lines = [line.strip() for line in (proc.stdout or "").splitlines() if line.strip()]
    if not lines:
        raise ToolError("The mismatch module produced no output on stdout.")
    try:
        kon = float(lines[-1])
    except ValueError as exc:
        raise ToolError(
            f"Could not parse a kon value from the mismatch module's last stdout line: {lines[-1]!r}"
        ) from exc
    mismatch_positions = [i + 1 for i, (a, b) in enumerate(zip(base, mismatched)) if a != b]
    return {
        "message": (
            f"kon = {kon} for {len(mismatch_positions)} mismatch(es) at position(s) "
            f"{mismatch_positions or 'none'}. The mismatch model calls model.eval() and has no "
            "dropout, so this value is deterministic. It is a learned score and is not guaranteed "
            "to be highest for an identical pair."
        ),
        "reference": REFERENCE_MISMATCH,
        "artifacts": [],
        "kon": kon,
        "base_rna": base,
        "mismatched_rna": mismatched,
        "num_mismatches": len(mismatch_positions),
        "mismatch_positions": mismatch_positions,
    }


@cli_wrapper_mcp.tool()
def scan_sirna_mismatch_variants(
    base_rna: Annotated[str, "Base 19-nt RNA sequence (A/C/G/U only)"],
    max_mismatches: Annotated[Literal[2, 3], "Enumerate variants up to 2 mismatches (1596 rows) or 3 (27759 rows)"] = 2,
    output_dir: Annotated[str | None, "Base output directory; a fresh subdirectory is created per call"] = None,
) -> dict:
    """Score every mismatch variant of a 19-nt siRNA up to the requested mismatch depth.
    Input is one 19-nt RNA sequence; output is a tab-separated table of variant sequence, kon and mismatch count.
    """
    runtime = _runtime_root()
    base = _validate_19nt_rna(base_rna, "base_rna")
    if int(max_mismatches) not in (2, 3):
        raise ToolError("max_mismatches must be 2 or 3 (the upstream -m 2 / -m 3 modes).")
    run_dir = _fresh_output_dir(output_dir, "mismatch_scan")
    # mismatch.py writes Args.output_dir + '/result/mismatch.txt', so pass the directory without a
    # trailing separator and pre-create its result/ subdirectory.
    (run_dir / "result").mkdir(parents=True, exist_ok=True)
    with _runtime_lock(runtime):
        _run_cli(
            runtime,
            ["-m", str(int(max_mismatches)), "--output_dir", str(run_dir)],
            timeout=_MISMATCH_SCAN_TIMEOUT,
            stdin_text=f"{base}\n",
            log_dir=run_dir / "logs",
        )
    table = run_dir / "result" / "mismatch.txt"
    if not table.is_file():
        raise ToolError(f"The mismatch module did not write its result table {table}.")
    header, rows = _read_table(table)
    if not rows:
        raise ToolError(f"The mismatch result table {table} has no data rows.")
    index = {name: i for i, name in enumerate(header)}
    counts: dict[str, int] = {}
    for row in rows:
        key = row[index["num_mismatch"]]
        counts[key] = counts.get(key, 0) + 1
    return {
        "message": (
            f"Scored {len(rows)} mismatch variants of {base} up to {int(max_mismatches)} "
            "mismatches. The mismatch model is deterministic (model.eval(), no dropout)."
        ),
        "reference": REFERENCE_MISMATCH,
        "artifacts": [
            {"description": "variant sequence, kon and mismatch count", "path": str(table)}
        ],
        "base_rna": base,
        "max_mismatches": int(max_mismatches),
        "num_rows": len(rows),
        "columns": header,
        "rows_per_num_mismatch": dict(sorted(counts.items())),
        "output_dir": str(run_dir),
    }

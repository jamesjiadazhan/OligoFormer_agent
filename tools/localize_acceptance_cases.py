"""Rewrite the published acceptance cases with absolute paths for this checkout.

`reports/mcp-acceptance-cases.json` is published with repository-relative paths so the evidence
carries no machine-specific paths. The acceptance helper (`verify_mcp_server.py`) requires
`artifact_root` to be absolute, so resolve the paths against this checkout before running it:

    python tools/localize_acceptance_cases.py            # -> tmp/mcp-acceptance-cases.local.json
    python <skill>/scripts/verify_mcp_server.py \
        --server src/OligoFormer_mcp.py \
        --expected reports/expected-mcp-tools.json \
        --cases tmp/mcp-acceptance-cases.local.json \
        --require-all-tools --call-timeout 1800 \
        --report tmp/acceptance.json

Only path-shaped values are rewritten; every assertion (expected_subset, required_inputs,
min_artifacts, error_contains, unique_artifacts) is copied through unchanged.
"""
from __future__ import annotations

import argparse
import json
import pathlib

PATH_KEYS = {"mrna_fasta_path", "sirna_fasta_path", "output_dir", "utr_reference_path",
             "orf_reference_path", "artifact_root"}


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=str(root / "reports" / "mcp-acceptance-cases.json"))
    ap.add_argument("--out", default=str(root / "tmp" / "mcp-acceptance-cases.local.json"))
    args = ap.parse_args()

    cases = json.loads(pathlib.Path(args.cases).read_text())
    for case in cases:
        if isinstance(case.get("artifact_root"), str):
            case["artifact_root"] = str(root / case["artifact_root"])
        for key, value in case.get("arguments", {}).items():
            if key in PATH_KEYS and isinstance(value, str) and not value.startswith("/"):
                case["arguments"][key] = str(root / value)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cases, indent=2) + "\n")
    print(f"wrote {out} ({len(cases)} cases resolved against {root})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

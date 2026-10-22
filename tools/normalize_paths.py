"""Replace machine-absolute paths in the published evidence with repository-relative paths.

The reports, logs and reference outputs were produced in a build workspace, so they recorded
absolute paths from that machine. This rewrites them for publication:

    <workspace>/mcp-work/<x>  ->  <x>          (the build project root is this repository root)
    <workspace>/mcp-work      ->  .
    <workspace>/<x>           ->  <WORKSPACE>/<x>
    ~/.claude/skills/paper2agent -> <SKILL_ROOT>
    /private/tmp/.../<session>/scratchpad     -> <SCRATCH>
    any other /Users/<user>   ->  <HOME>

Run from the repository root:  python tools/normalize_paths.py [--check]

--check reports what would change and exits 1 if anything would, without writing.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

PROJECT = "/Users/james/Local/paper2agent/paper2agent_test/OligoFormer/mcp-work"
WORKSPACE = "/Users/james/Local/paper2agent/paper2agent_test/OligoFormer"
SKILL_ROOT = "/Users/james/.claude/skills/paper2agent"
CONDA = "/Users/james/miniconda"
HOME = "/Users/james"
SCRATCH_RE = re.compile(r"/private/tmp/claude-\d+/[^/\s\"']+/[^/\s\"']+/scratchpad")
SCRATCH_FALLBACK_RE = re.compile(r"/private/tmp/claude-\d+/[^/\s\"']+/[^/\s\"']+")

# Files whose bytes are scientific reference outputs or fixtures: never rewritten.
SKIP_SUFFIXES = {".jpg", ".png", ".pth", ".zip", ".npy", ".gif", ".pdf", ".xlsx"}
SKIP_DIRS = {".git", "runtime", ".venv", "__pycache__"}


def normalize(text: str) -> str:
    text = text.replace(PROJECT + "/", "").replace(PROJECT, ".")
    text = text.replace(SKILL_ROOT, "<SKILL_ROOT>")
    text = text.replace(CONDA, "<CONDA>")
    text = SCRATCH_RE.sub("<SCRATCH>", text)
    text = SCRATCH_FALLBACK_RE.sub("<SCRATCH>", text)
    text = text.replace(WORKSPACE + "/", "<WORKSPACE>/").replace(WORKSPACE, "<WORKSPACE>")
    text = text.replace(HOME + "/", "<HOME>/").replace(HOME, "<HOME>")
    return text


def candidates(root: pathlib.Path):
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.is_symlink():
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix.lower() in SKIP_SUFFIXES:
            continue
        yield p


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="report changes without writing")
    args = ap.parse_args()

    root = pathlib.Path(__file__).resolve().parents[1]
    changed = []
    for p in candidates(root):
        if p.resolve() == pathlib.Path(__file__).resolve():
            continue  # this script documents the literal prefixes it replaces
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        new = normalize(text)
        if new != text:
            changed.append(p.relative_to(root))
            if not args.check:
                p.write_text(new, encoding="utf-8")

    verb = "would change" if args.check else "rewrote"
    print(f"{verb} {len(changed)} files")
    for c in changed[:20]:
        print("  ", c)
    if len(changed) > 20:
        print(f"   ... and {len(changed) - 20} more")
    return 1 if (args.check and changed) else 0


if __name__ == "__main__":
    sys.exit(main())

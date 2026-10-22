# Published evidence has been path-normalized

The reports, logs and reference outputs in this repository were produced in a build workspace, so
they originally recorded absolute paths from the machine that produced them
(`/Users/<user>/.../mcp-work/...`, scratch directories, and the skill install location). For
publication those were rewritten by [`tools/normalize_paths.py`](../tools/normalize_paths.py):

| Original | Published as |
| --- | --- |
| `<workspace>/mcp-work/<path>` | `<path>` — the build project root is this repository root |
| `<workspace>/mcp-work` | `.` |
| `<workspace>/<path>` (outside the project) | `<WORKSPACE>/<path>` |
| `~/.claude/skills/paper2agent` | `<SKILL_ROOT>` |
| temporary session scratch directories | `<SCRATCH>` |
| any other path under the build user's home | `<HOME>/...`, `<CONDA>/...` |

Re-running `python tools/normalize_paths.py --check` from the repository root reports no remaining
changes.

## What this means for the recorded hashes

**This is a byte-level rewrite of evidence files, so the sha256 values recorded inside
`agent-runs.json` (run `report_sha256`, `produced_files`, `tested_files` and `completion.files`) no
longer match the published copies of the text files that were rewritten.** Those hashes were
computed over the original bytes during the build, and they are left as they were rather than
refreshed, because refreshing them would make a rewritten file look like the file that was actually
verified. For the same reason, running the Paper2Agent workflow checker against this repository will
report hash mismatches; that is expected here and is not a validation failure.

Unaffected, and still hash-exact:

- `dist/OligoFormer-mcp.zip` — the archive that passed delivery validation, sha256
  `dde0644144f5c38f23c826516c1031ff5e2656a7ca2d4185c11bf060f6f1655f`.
- `src/tools/cli_wrapper.py` — the verified production module, sha256
  `833e57bbd04e6e4ad02972a10c9b35999ce7454d44a51c906c2be17d54bfb3f0`, as recorded in
  `verification-cli_wrapper.json`.
- Every scientific reference output: result tables, `pita.tab`, `targetscan.tab`, `mismatch.txt`,
  FASTA inputs, figures and CSVs. Binary and tabular reference data were never rewritten (the
  normalizer skips `.jpg`, `.png`, `.pth`, `.zip`, `.npy`, `.gif`, `.pdf`, `.xlsx`), and the
  comparisons those files support are unchanged.

The authoritative, unmodified evidence — including the original absolute paths and bytes matching
every recorded hash — remains in the build workspace outside this repository.

## Running the acceptance cases

`mcp-acceptance-cases.json` now holds repository-relative paths. The acceptance helper requires an
absolute `artifact_root`, so resolve them against your checkout first:

```bash
python tools/localize_acceptance_cases.py     # writes tmp/mcp-acceptance-cases.local.json
```

That rewrites only path-shaped values; every assertion is copied through unchanged.

# Coordinator delivery repair (stage 6, after delivery-validation attempt 1)

The independent delivery verifier (agent a71e7150dc44332a3) validated archive sha256 `8e1229...c462` and
reported three defects in the packaged installer/documentation. The wrapper and server were unaffected
(`src/tools/cli_wrapper.py` stayed at the verified revision `833e57bb...`), so no scientific re-verification
of the tools was required; only the installation path changed.

Changes made by the coordinator:

1. `deps/setup_oligoformer_runtime.sh` (finding F4, high). The host PITA binaries are now prepared *before*
   the `git checkout` / `git clean` of `off-target/pita`, and a sha256 mismatch against the published-binary
   constants now rebuilds from the pinned ViennaRNA 1.6 source instead of calling `die`. Previously a
   recipient without the published arm64 binaries got a locally built RNAddG4 whose hash necessarily differs
   (ad-hoc signature / LC_UUID); the second run of the script then aborted *after* resetting the tree,
   leaving Linux x86-64 binaries and an unsubstituted `EXE_BASE_DIR`, i.e. a broken `screen_sirna_offtarget`.
   Coordinator test before repackaging: run 1 then run 2 with an empty `PITA_BIN_DIR` both exit 0, run 2 logs
   the rebuild, and the staged `RNAduplex`/`RNAddG4` are Mach-O arm64 with `EXE_BASE_DIR` substituted.
2. `deps/create_offtarget_env.sh` (finding F3, medium). A shared `smoke_ok()` now runs the perl modules,
   `RNAplfold --version` and the GNU sed/grep checks for both the "already complete" short-circuit and the
   post-`conda create` check, with one retry before failing. The old short-circuit only tested `-x` on
   RNAplfold, so a prefix whose RNAplfold aborted was reported as complete.
3. `USAGE.md` (finding F2, medium): the documented `uv venv` branch now passes `--seed` (otherwise the next
   documented `python -m pip install` fails with "No module named pip"), with a `uv pip install` alternative,
   plus a paragraph documenting the first-run PITA build and that re-running the setup script is safe.

Artifact history: `deps/setup_oligoformer_runtime.sh` and `deps/create_offtarget_env.sh` were produced by
`environment-attempt-2`. The bytes that attempt produced are retained unchanged at
`deps/archive/setup_oligoformer_runtime.attempt-2.sh` (sha256 `5f8f66036f25...`) and
`deps/archive/create_offtarget_env.attempt-2.sh` (sha256 `f9e9adcb911d...`); that run record now points at the
archives, so its historical hashes still verify. The repaired files are owned by this coordinator record.

The repair invalidated the delivery result for archive `8e1229...c462`. The package was rebuilt as
`dist/OligoFormer-mcp.zip` sha256 `dde0644144f5c38f23c826516c1031ff5e2656a7ca2d4185c11bf060f6f1655f` and
re-validated by the same independent delivery verifier (attempt 2, success=true, 18/18 acceptance cases and
16/16 scientific comparisons from a fresh extraction, with F2/F3/F4 explicitly retested). Attempt 1 and its
findings are preserved in `reports/delivery-validation.json` under `previous_attempts`.

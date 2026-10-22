# Stage 2 executor brief (OligoFormer, CLI route, python backend)

Bound values
- PROJECT_ROOT = .
- SKILL_ROOT = <SKILL_ROOT>/paper2mcp
- PROJECT_PYTHON = $PROJECT_ROOT/OligoFormer-env/bin/python (Python 3.11.12, torch 2.2.1 CPU)
- Pinned source: $PROJECT_ROOT/repo/OligoFormer (lulab/OligoFormer @ e2f53ad63387bbe166bf123949151e2bc9bf6ec3). READ-ONLY: never run commands that write inside it.
- Scanner inventory: $PROJECT_ROOT/reports/tutorial-scanner.json (your assignment entry + tool_review + wrapping_hazards).
- Environment report: $PROJECT_ROOT/reports/environment-manager_results.md (read the "Attempt 2: runtime copy" section).

Read first: $SKILL_ROOT/references/agents/cli/tutorial-executor-cli.md, $SKILL_ROOT/references/agents/tutorial-executor.md, $SKILL_ROOT/references/routes/cli.md ("Real execution and reference outputs"), $SKILL_ROOT/references/runtime.md (stage 2).

Your private runtime copy
- Several executors run concurrently and OligoFormer writes to fixed paths inside its working tree (./data/infer/<record>/, ./result/, off-target/pita tmp files), so you MUST build and use your own runtime copy:
  `OLIGOFORMER_SRC=$PROJECT_ROOT/repo/OligoFormer RNAFM_SRC=$PROJECT_ROOT/deps/downloads/rnafm_git/RNA-FM RNAFM_WEIGHTS=$PROJECT_ROOT/deps/downloads/rnafm_weights_hf/RNA-FM_pretrained.pth bash $PROJECT_ROOT/deps/setup_oligoformer_runtime.sh $PROJECT_ROOT/notebooks/<ID>/runtime/OligoFormer`
  (check the script header for exact variable semantics; if RNAFM_SRC isn't a valid git clone, omit it and let the script use its default). Then `source $PROJECT_ROOT/notebooks/<ID>/runtime/OligoFormer.env.sh` and run commands with cwd = that runtime copy. Record the runtime copy's `git rev-parse HEAD` and `git status --short | wc -l`.
- Do not install packages. If something is missing, stop and report the exact need (the environment manager owns dependencies).

Known upstream behaviour (record, do not fix)
- scripts/infer.py never calls model.eval() and never seeds torch, while scripts/model.py has nn.Dropout(0.02) layers, so efficacy predictions are stochastic (~1% run-to-run variation on CPU). Therefore: run the documented command **3 times** (fresh output dirs each, same inputs) and record per-row spread (max abs diff, Spearman of efficacy across runs, whether the ranked order / filter columns are identical). Deterministic columns (pos, sense, siRNA, func_filter, toxicity/off-target scores) must match exactly across runs. This spread defines the justified tolerance the verifier will use.
- --output_dir is string-concatenated, so it must end with '/' for --infer; the mismatch module writes <output_dir>/result/mismatch.txt and the directory must exist.

Evidence layout (you own only these)
- $PROJECT_ROOT/notebooks/<ID>/ : runtime/ (your copy), `run.sh` driver (bash, `set -euo pipefail`, records each command, exit status, wall time), `data/` (exact input files copied from the runtime copy with sha256), `outputs/run1..3/` (all native outputs: result tables, data/infer/<record>/ intermediate files you consider relevant such as pita.tab/targetscan.tab/top_n_siRNA.fa — NOT the large .npy embeddings unless small), `logs/` (stdout/stderr per run).
- $PROJECT_ROOT/reports/executed_notebook_<ID>.json : keyed by "<ID>"; include source path + commit, status (success/failed/blocked), execution_mode "native_command", execution_path (your run.sh), commands with exit statuses, runtimes, executable identity (python path/version, perl/RNAplfold versions if used), cwd, env vars set (no secrets), input/output paths with sha256, run-to-run comparison results, replay check (rerun from saved inputs in notebooks/<ID>/data and compare), observed hazards/bugs, limitations, attempts.
- Do not write reports/executed_notebooks.json or other agents' paths. At most 5 attempts. When done you may delete notebooks/<ID>/runtime/OligoFormer/RNA-FM/redevelop/pretrained weights hard link? NO — leave the runtime copy in place (weights are hard links, no extra disk).

Final reply: ≤12 lines: status, commands + exit codes, runtime, key output shapes/row counts, run-to-run spread, report path.

# OligoFormer MCP server

MCP tools for siRNA design with **OligoFormer** (Bai *et al.*, *Bioinformatics* 40(10):btae577, 2024; code: [lulab/OligoFormer](https://github.com/lulab/OligoFormer) at commit `e2f53ad63387bbe166bf123949151e2bc9bf6ec3`). Each tool runs the original OligoFormer command line (`python scripts/main.py ...`) in a local runtime copy and returns its native result tables. The wrappers don't reimplement the model, the filters or the off-target scoring.

## Tools

| Tool | When to use it | Upstream command |
| --- | --- | --- |
| `predict_sirna_efficacy_from_mrna` | Rank every possible 19-nt siRNA against one or more target mRNAs | `--infer 1 -i1 <mRNA.fa>` (+ `-nf`, `-tox`, `--toxicity_threshold`) |
| `predict_specific_sirna_efficacy` | Score only the siRNA candidates you already have | `--infer 1 -i1 <mRNA.fa> -i2 <siRNA.fa>` |
| `screen_sirna_offtarget` | Add PITA + TargetScan off-target screening (slow) to the full window scan | `--infer 1 -i1 <mRNA.fa> -off` (+ thresholds, `--utr/--orf`, `-tox`) |
| `compare_sirna_mismatch_kon` | Relative on-rate (kon) of one mismatched siRNA against its base sequence | `-m 1` (two sequences on stdin) |
| `scan_sirna_mismatch_variants` | kon for every 1–2 (or 1–3) mismatch variant of a siRNA | `-m 2` / `-m 3` (base sequence on stdin) |

All tools return a JSON object with `message`, `reference` (upstream source URL), `artifacts` (`[{description, path}]`, absolute paths in a fresh per-call directory) and compact summary fields.

### Inputs and parameters

- **mRNA FASTA** (`mrna_fasta_path`): one or more records, each at least 19 nt. DNA (T) is converted to U by upstream. Each header must be unique and match `^[A-Za-z0-9_.-]+$`, with no spaces: upstream uses the name in unquoted shell commands and output file names.
- **siRNA FASTA** (`sirna_fasta_path`, specific-siRNA tool only): every sequence exactly 19 nt of `A/C/G/U` (guide/antisense strand, 5'→3'). The tool scores only siRNAs whose reverse complement occurs in a record; it uses the same matching rule as `scripts/infer.py`. Upstream crashes when a record has no matching siRNA. The wrapper passes upstream only the records that have matches and reports `unmatched_records` and `unmatched_sirnas`. It returns an error if nothing matches.
- **Filters** (efficacy and off-target tools):
  - `apply_functionality_filter` (default `true`): GC content 30–65%, no 5-base runs, no 6 G/C in a row, no palindromes. `false` passes `-nf`.
  - `include_toxicity` (default `false`): seed (nt 2–7) cell-viability lookup.
  - `toxicity_threshold` (default `50.0`): viability % below which a seed is flagged. It only applies when `include_toxicity` is true, as upstream.
- **Off-target** (`screen_sirna_offtarget`):
  - `pita_threshold` (default `-10`): a window is flagged when its PITA ddG is below this value.
  - `targetscan_threshold` (default `1.0`): a window is flagged when its TargetScan context score is above this value.
  - `utr_reference_path` / `orf_reference_path`: default to the bundled 100-gene references (`off-target/ref/100_human_{UTR,ORF}.txt`) in the runtime copy. Paths must contain no spaces (upstream passes them unquoted to a shell).
- **Mismatch tools**: `base_rna` and `mismatched_rna` are exactly 19 nt of `A/C/G/U`. `max_mismatches` is 2 (1596 variants, ~6 s) or 3 (27759 variants, ~95 s).
- **`output_dir`** (optional): base directory; each call creates its own timestamped subdirectory. Defaults to `tmp/outputs/cli_wrapper` inside the installed package.

### Output tables

The efficacy tools write three tab-separated tables per mRNA record: `<record>.txt` (all windows), `<record>_ranked.txt` (sorted by efficacy) and `<record>_ranked_filtered.txt` (only rows where `filter == 0`). Columns: `pos` (row counter), `sense` (passenger strand), `siRNA` (guide strand), `efficacy`, and the filter columns in use (`func_filter`, `Seed`/`cell_viability`/`toxicity_filter`, `pita_score`/`targetscan_score`/`off_target_filter`, and their sum `filter`). `screen_sirna_offtarget` additionally returns the raw `pita.tab` and `targetscan.tab` per record.

## Reproducibility: efficacy scores are stochastic

Upstream `scripts/infer.py` never calls `model.eval()`, so the model's dropout layers stay active during prediction, and it never seeds torch. **Repeated runs on identical input give different efficacy values** — measured differences up to ~0.18 absolute (~30% relative), with Spearman correlation 0.93–0.99 between runs, and a different `_ranked` order each time. This is upstream behaviour and is deliberately not "fixed" here, since changing it would alter the published method's output.

Everything else is reproducible: `pos`, `sense`, `siRNA`, `func_filter`, `Seed`, `cell_viability`, `toxicity_filter`, `pita_score`, `targetscan_score`, `off_target_filter` and `filter` match exactly, and the PITA/TargetScan tables are byte-identical across runs. Treat the efficacy column as a ranking signal, not an exact number, and prefer the set of windows passing the filter over their precise order.

The mismatch tools (`compare_sirna_mismatch_kon`, `scan_sirna_mismatch_variants`) use a model in eval mode and **are** deterministic.

## Installation

Requires macOS or Linux, Python 3.11, `git`, and — for `screen_sirna_offtarget` only — `conda` (Miniconda/Miniforge) plus a C compiler (`make`, `cc`; on macOS the Xcode Command Line Tools), used to rebuild PITA's ViennaRNA 1.6 helpers for the host architecture. Tested on macOS 15 (Apple Silicon, arm64) with Python 3.11.12 and torch 2.2.1 (CPU). Installation needs roughly 3 GB of downloads (the 1.2 GB RNA-FM checkpoint, the conda prefix, the checkouts) and about 4 GB of disk.

```bash
unzip OligoFormer-mcp.zip && cd OligoFormer-mcp
chmod +x deps/*.sh                      # ZIP extraction can drop the executable bit

python3.11 -m venv .venv                # or: uv venv .venv --python 3.11 --seed
./.venv/bin/python -m pip install -r src/requirements.txt
# With uv you can skip pip entirely:
#   uv venv .venv --python 3.11 && uv pip install --python .venv/bin/python -r src/requirements.txt

# The scientific runtime: OligoFormer @ e2f53ad + RNA-FM + PITA/TargetScan.
# Downloads the 1.2 GB RNA-FM checkpoint (sha256-verified) on first run.
./deps/create_offtarget_env.sh "$PWD/deps/offtarget-env"          # needed only for screen_sirna_offtarget
PROJECT_ENV="$PWD/.venv" ./deps/setup_oligoformer_runtime.sh "$PWD/runtime/OligoFormer"
```

Install the runtime at a path **without spaces**: the upstream off-target scripts pass paths unquoted to a shell, so `screen_sirna_offtarget` fails under a path containing a space. The other four tools are unaffected.

On a host without the published arm64 PITA binaries, the script builds `RNAduplex`/`RNAddG4` from the bundled ViennaRNA 1.6 source on the first run and reuses them afterwards; a locally built binary differs from the published hashes, so later runs rebuild rather than fail. Re-running the script is safe: the host binaries are prepared before the checkout is reset, so a failed build leaves the previous working runtime in place.

`setup_oligoformer_runtime.sh` is idempotent and accepts `OLIGOFORMER_SRC`, `RNAFM_SRC`, `RNAFM_WEIGHTS` (an existing checkpoint to hard-link instead of downloading), `PROJECT_ENV`, `OFFTARGET_ENV`, `PITA_BIN_DIR`, `PITA_FORCE_BUILD=1` and `SKIP_OFFTARGET=1`. On macOS it rebuilds PITA's `RNAduplex`/`RNAddG4` for the host architecture (the upstream binaries are Linux x86-64) and runs `make install`; on Linux it keeps the shipped binaries. Moving the runtime afterwards requires re-running the script, because `make install` writes absolute paths into the PITA scripts.

## Running the server

```bash
./.venv/bin/python src/OligoFormer_mcp.py
```

The server finds its runtime at `<package>/runtime/OligoFormer` and the off-target tools at `<package>/deps/offtarget-env`. Override with `OLIGOFORMER_RUNTIME` and `OLIGOFORMER_OFFTARGET_ENV`. No credentials or network access are needed at call time.

Client registration (not performed by this package; syntax from the tested `fastmcp` 4.0.3):

```bash
./.venv/bin/fastmcp install claude-code src/OligoFormer_mcp.py:mcp   --name OligoFormer --with-requirements src/requirements.txt
```

`fastmcp install` also supports `claude-desktop`, `cursor`, `gemini-cli`, `goose` and `mcp-json`.

## Runtimes and timeouts

Measured on the test host (CPU, 20 cores) for `data/example.fa`, 2 records × 39 windows:

| Tool | Typical runtime |
| --- | --- |
| `predict_sirna_efficacy_from_mrna`, `predict_specific_sirna_efficacy` | 15–50 s |
| `screen_sirna_offtarget` | 3–11 min (perl runs under Rosetta 2 on Apple Silicon) |
| `compare_sirna_mismatch_kon` | ~2 s |
| `scan_sirna_mismatch_variants` | ~6 s (depth 2), ~95 s (depth 3) |

Runtime scales with the number of 19-nt windows, which is `len(mRNA) - 18` per record. Give MCP clients a call timeout of at least 1800 s for off-target screening. Calls are serialized by a lock, because upstream writes to fixed paths inside the runtime copy; concurrent calls queue rather than interfere.

## Validated scope

- 108 scientific tests (one file per tool) compare the tools against direct upstream runs of the same commands: deterministic columns and row sets exactly, PITA/TargetScan tables by sha256, mismatch values within 1e-7, efficacy within 0.25 absolute and Spearman ≥ 0.9 per record. 105 passed with 3 slow tests run separately.
- 18 stdio acceptance cases (every tool, plus input errors and repeated-call artifact isolation) passed in the development environment and again in a freshly installed environment with an independently rebuilt runtime.
- Only the five tools above are covered. This does not establish that the rest of the repository, or the paper's published results, reproduce.

## Not included, and why

- **Training** (`scripts/main.py` without `--infer`): needs a CUDA GPU (~60 min on a V100) and was not run.
- **Evaluation mode** (`-t`, `-s -t`): needs precomputed RNA-FM embedding files (`data/RNAFM/<dataset>_{siRNA,mRNA}`) that the repository doesn't ship and only training regenerates. It also reports metrics to a log rather than an artifact.
- **`-top N`** (off-target for the top N siRNAs only): upstream writes the top-N FASTA with headers of the top-ranked siRNAs but sequences taken in positional order, so PITA/TargetScan scores are attached to the wrong rows. Exposing it would return mislabelled results.
- **`-a/--all_human`** (genome-wide references): upstream sets the ORF reference to the UTR file, and the full references ship only as `.zip`. The run exits 0 with all-NaN scores, so every siRNA silently passes the off-target filter.
- **`--infer 2`** (interactive mRNA entry): duplicates `predict_sirna_efficacy_from_mrna` and writes to a fixed record name.

The wrapper verifies after each off-target run that `pita.tab` and `targetscan.tab` are present and non-empty and that the score columns are not entirely missing, because upstream ignores the exit status of its perl subprocesses.

## Licence

OligoFormer is free for academic and non-commercial use; commercial use requires written authorization from Tsinghua University (`ott@tsinghua.edu.cn`). See `THIRD_PARTY_NOTICES/`. RNA-FM, PITA and TargetScan carry their own licences.

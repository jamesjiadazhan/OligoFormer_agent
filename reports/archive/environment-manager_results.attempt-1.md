# Environment manager results — OligoFormer (Paper2MCP stage 1, CLI route / Python backend)

Host: macOS 24.6.0 (Darwin), arm64 Apple Silicon, 20 CPUs, 64 GB RAM, **no CUDA**.
Route: `cli` / `cli_backend: python` (`.pipeline/language.json`).

## 1. Interpreter and environment

| Item | Value |
| --- | --- |
| PROJECT_ENV | `OligoFormer-env` |
| PROJECT_PYTHON | `OligoFormer-env/bin/python` |
| Python | 3.11.12 (uv-managed `cpython-3.11.12-macos-aarch64-none`) |
| Creator | `uv 0.7.3` (`<HOME>/.local/bin/uv`) |
| Jupyter kernel | `oligoformer-mcp` ("OligoFormer MCP (py3.11)") at `OligoFormer-env/share/jupyter/kernels/oligoformer-mcp` |

### Deviation from README's Python 3.8 (recorded)

README badge/instructions say `conda create -n oligoformer python=3.8*`. **No deviation from the
scientific pins was needed**: every pin in `requirements.txt` installed unchanged on Python 3.11
alongside `fastmcp==4.0.3`. There is therefore a single environment; no separate Python-3.8
scientific runtime and no second interpreter for the `sh scripts/RNA-FM.sh` subprocess.

Resolved versions (exact repo pins preserved):

```
torch==2.2.1  torchvision==0.17.1  numpy==1.24.4  pandas==2.0.3  scipy==1.10.1
scikit-learn==1.3.1  matplotlib==3.7.5  bio==1.6.2  prefetch-generator==1.0.3
ptflops==0.7.3  pytorch-ignite==0.5.0.post2  tqdm==4.66.1  yacs==0.1.8
fastmcp==4.0.3  pytest==9.1.1  pytest-asyncio==1.4.0
papermill / nbclient / ipykernel / jupytext (notebook evidence)
```

Note: `bio==1.6.2` pulls `biopython==1.88` (unpinned transitive); `Bio.__version__` is 1.88.
Full snapshot: `reports/environment-requirements.txt` (149 packages). Version: `reports/python-version.txt`.

### Source identity

`repo/OligoFormer` is a git clone of https://github.com/lulab/OligoFormer at
`e2f53ad63387bbe166bf123949151e2bc9bf6ec3`; **working tree clean** after my work (only untracked
`RNA-FM/`, which I own). `scripts/` is not an installable package — it is run as
`python scripts/main.py ...` from the repo root, so nothing from the repo is pip-installed.
Model weights shipped in the repo:
- `model/best_model.pth` sha256 `8c2fb4bc847e20c63dfe69381c5dc5c341be71e65c011f4e9a4cd6fdb0e34eca`
- `model/mismatch_model.pth` sha256 `fb32c9a31d436113ff5271edd138cc86fe4065b7b1082e7c1cadb3563d627043`

### Device availability (measured)

```
torch.cuda.is_available() -> False
torch.backends.mps.is_built() -> True ; torch.backends.mps.is_available() -> True
torch.get_num_threads() -> 20
```
OligoFormer selects `cuda:0` if available else `cpu`, so it runs on **CPU** here (MPS is never
selected by the repo code). RNA-FM required a patch to reach its CPU path (see §2).

## 2. RNA-FM dependency

`scripts/infer.py` runs `sh scripts/RNA-FM.sh ../../data/infer/<name>`, which `cd ./RNA-FM/redevelop`
and calls `python launch/predict.py --config=pretrained/extract_embedding.yml ...` with whatever
`python` is first on PATH.

### Provenance

- **README source 1 (Tsinghua Cloud) is DEAD.** `curl -L 'https://cloud.tsinghua.edu.cn/f/46d71884ee8848b3a958/?dl=1'`
  returns HTTP 200 with a 7220-byte Seafile HTML page reading "File does not exist" (no tarball).
  Saved as `deps/downloads/rnafm/source1_response.html`, sha256
  `9f58c7619ca41f510e82e70c3c928966e087a501740b51e2ed348b7d4f4de0fd`. **No RNA-FM.tar.gz sha256 exists.**
- **README source 2 used instead:** `git clone https://github.com/ml4bio/RNA-FM.git`, commit
  `348951516e0963d22bbb33b3c9fc18c89081d38e` (2025-05-27, "Update pretrained.py"), cloned into
  `deps/downloads/rnafm_git/RNA-FM` and rsynced (minus `.git`) to `repo/OligoFormer/RNA-FM/` (untracked, 1.1 GB).
- **Weights:** the README gdrive folder was not used; the upstream code's own checkpoint URL was.
  `fm/pretrained.py` at the cloned commit points at
  `https://proj.cse.cuhk.edu.hk/rnafm/api/download?filename=RNA-FM_pretrained.pth`; the immediately
  preceding commit `224fa9e` used the HuggingFace mirror, which I downloaded:
  `https://huggingface.co/cuhkaih/rnafm/resolve/main/RNA-FM_pretrained.pth`
  → 1,194,424,423 bytes, **sha256 `5b5d7d87b37c291ef42c140ef9edf7aea29f255fa2a4fd435f776c52e93d5e99`**,
  placed at `repo/OligoFormer/RNA-FM/redevelop/pretrained/RNA-FM_pretrained.pth`
  (the path `redevelop/pretrained/extract_embedding.yml` + `fm.pretrained.rna_fm_t12()` resolve to).
  Downloaded into its own empty directory and inspected before use: it is a torch zip archive whose
  `archive/data.pkl` references only `argparse.Namespace`, `collections.OrderedDict`,
  `torch.FloatStorage`, `torch.DoubleStorage`, `torch._utils._rebuild_tensor_v2` (benign;
  `weights_only=True` refuses it only because of `argparse.Namespace`).

### Making it run on CPU

`redevelop/launch/predict.py` defaulted `--device gpu`, and `scripts/RNA-FM.sh` passes no `--device`,
so on this host it died with `AssertionError: Torch not compiled with CUDA enabled`.
One-line patch to the **untracked dependency only** (no OligoFormer tracked file touched):

```diff
-        "--device", default="gpu", choices=["cpu", "gpu"]
+        "--device", default=("gpu" if torch.cuda.is_available() else "cpu"), choices=["cpu", "gpu"]
```
File: `repo/OligoFormer/RNA-FM/redevelop/launch/predict.py`; diff saved at `deps/rnafm_cpu_default.patch`.
No separate interpreter is needed: `predict.py` does `sys.path.append('..')`, so the sibling `fm/`
package is imported from the checkout and PROJECT_PYTHON runs it as-is (`rna-fm` is **not** pip-installed).

### PATH / env requirements for the subprocesses

- `PATH` must have `OligoFormer-env/bin` first so `scripts/RNA-FM.sh`'s bare `python` is PROJECT_PYTHON.
- `TORCH_HOME=<PROJECT_ROOT>/deps/torch_home` is set in my runs (the hub checkpoint dir is symlinked
  to the local weights file). It is **not required** — the local `pretrained/RNA-FM_pretrained.pth`
  is loaded directly — but it prevents any accidental network download into `~/.cache/torch`.
- All commands must be run with cwd = `repo/OligoFormer` (relative paths `./data`, `./model`,
  `./RNA-FM`, `./off-target`, `./toxicity`, `./result`).

## 3. Commands run, with exit status

| # | Command (cwd) | Exit |
| --- | --- | --- |
| 1 | `uv venv --python 3.11 OligoFormer-env` | 0 |
| 2 | `uv pip install --python $PROJECT_PYTHON -r repo/OligoFormer/requirements.txt fastmcp==4.0.3 pytest pytest-asyncio` | 0 |
| 3 | `uv pip install --python $PROJECT_PYTHON papermill nbclient ipykernel jupytext` | 0 |
| 4 | `uv pip check --python $PROJECT_PYTHON` → "Checked 149 packages … All installed packages are compatible" | 0 |
| 5 | `$PROJECT_PYTHON -m ipykernel install --prefix OligoFormer-env --name oligoformer-mcp` | 0 |
| 6 | import check (torch, numpy, pandas, scipy, sklearn, matplotlib, torchvision, ignite, ptflops, yacs, Bio, fastmcp, tqdm, prefetch_generator) | 0 |
| 7 | `curl … cloud.tsinghua.edu.cn/f/46d71884ee8848b3a958/?dl=1` (RNA-FM tarball) | 0 but **HTML "File does not exist"** |
| 8 | `git clone https://github.com/ml4bio/RNA-FM.git` | 0 |
| 9 | `curl … huggingface.co/cuhkaih/rnafm/resolve/main/RNA-FM_pretrained.pth` (1.19 GB, 19 s) | 0 |
| 10 | RNA-FM `predict.py` probe, default device (pre-patch) | **1** (`Torch not compiled with CUDA enabled`) |
| 11 | RNA-FM `predict.py` probe `--device cpu` → `saved embedding shape: (29, 640)` | 0 |
| 12 | **`python scripts/main.py --infer 1 -i1 data/example.fa -i2 data/example_siRNA.fa --output_dir <scratch>/`** | **0** |
| 13 | `curl … cloud.tsinghua.edu.cn/f/cab2afdf951140a48fec/?dl=1` (PerlLib.zip) | 0 but **HTML "File does not exist"** |
| 14 | `conda create -p deps/offtarget-env -c conda-forge -c bioconda perl perl-bioperl perl-statistics-lite viennarna` | 0 |
| 15 | `conda install -p deps/offtarget-env -c conda-forge sed grep` (GNU sed 4.10, GNU grep 3.12) | 0 |
| 16 | `perl -e 'use Statistics::Lite; use Bio::TreeIO;'` with conda perl 5.32.1 / BioPerl 1.7.8 | 0 |
| 17 | ViennaRNA-1.6 (bundled PITA source, copied to `deps/pita-vienna`) `./configure … && make` | 2 (only `RNAforester` C++ fails; `RNAduplex` built) |
| 18 | `cc -O2 -std=gnu89 … -o RNAddG4 RNAddG4.c ../lib/libRNA.a -lm` | 0 |
| 19 | `cd off-target/pita && make install` (GNU sed first on PATH) | 0 |
| 20 | **`python scripts/main.py --infer 1 -i1 data/example.fa -i2 data/example_siRNA.fa -off -tox --output_dir …`** | **0** (first try had BSD-`grep -P` errors; clean after #15) |
| 21 | `python scripts/main.py --infer 1 … -tox --output_dir …` (toxicity only) | 0 |
| 22 | `printf '…\n…\n' \| python scripts/main.py -m 1` → `0.9940418601036072` | 0 |
| 23 | `cd off-target/pita && make clean` + `git checkout -- off-target/pita` → tree clean | 0 |

### Inference smoke test (#12) detail

```
cd repo/OligoFormer
PATH=<PROJECT_ENV>/bin:$PATH TORCH_HOME=<PROJECT_ROOT>/deps/torch_home \
  python scripts/main.py --infer 1 -i1 data/example.fa -i2 data/example_siRNA.fa \
  --output_dir <PROJECT_ROOT>/tmp/outputs/env_check/infer1/
```
Exit 0. **Runtime ~15 s wall** (script self-report: "10 seconds"), CPU only, 2 mRNAs × 1 matching
siRNA each. Outputs: `RNA0.txt`, `RNA0_ranked.txt`, `RNA0_ranked_filtered.txt` and the RNA1
equivalents in the given `--output_dir` (note: `--output_dir` is concatenated, so it **must end in `/`**).

**Side-effect writes outside `--output_dir` (recorded):**
- `repo/OligoFormer/data/infer/<seq_name>/` — `siRNA.fa`, `mRNA.fa`, and RNA-FM embedding dirs
  `siRNA/`, `mRNA/` (plus `pita.tab`, `targetscan.tab` with `-off`). `infer()` does
  `rm -rf ./data/infer/<name>` first, so a re-run wipes that subdirectory.
- `repo/OligoFormer/result/` — created unconditionally by `scripts/main.py` even when `--output_dir`
  points elsewhere (empty unless `-m 2/3`).
- With `-off`: `off-target/pita/` receives `<name>_pita_results.tab` and **overwrites the tracked
  files `off-target/pita/RNA0_pita_results.tab`, `tmp_seqfile1`, `tmp_seqfile2`**;
  `off-target/tmp/` is created and removed by `scripts/targetscan.sh`.
Both `data/infer/` and `result/` were deleted after my smoke tests, and the tracked PITA files
restored with `git checkout`; `git status` shows only untracked `RNA-FM/`.

**Reproducibility caveat for the executors:** `efficacy` is not bit-stable across runs on this CPU
host — RNA0 gave 1.029039, 1.023211, 1.016862, 1.016512 across four runs (~1% spread) with the same
seed. Reference comparisons must use a tolerance, not exact equality.

## 4. Off-target and toxicity readiness

**Off-target (`-off`): WORKING, with a documented non-default PATH and a staging step.**

Requirements found and how they were met:

| Requirement | Status |
| --- | --- |
| perl `Statistics::Lite`, `Bio::TreeIO` (TargetScan scripts) | conda prefix `deps/offtarget-env` (perl 5.32.1, BioPerl 1.7.8, perl-statistics-lite). System `/usr/bin/perl` 5.34.1 has **neither**. |
| ViennaRNA `RNAplfold` (TargetScan context scores) | `deps/offtarget-env/bin/RNAplfold`, ViennaRNA **2.7.2** (conda-forge). |
| GNU `grep -P` | `targetscan_70_context_scores.pl:1970` shells out to `grep -A7 -P`; macOS BSD grep fails ("invalid option -- P") and silently loses the site-accessibility feature — `targetscan_score` changed 0.409 → 0.414 once GNU grep 3.12 was on PATH. **GNU grep is mandatory for correct scores.** |
| GNU `sed -i` | `off-target/pita/Makefile` uses GNU-style `sed -i 's@…@…@g' file`; BSD sed misparses it. GNU sed 4.10 in the conda prefix. (The repo already ships 24 tracked `off-target/pita/lib/sed*` junk files from an earlier BSD-sed `make install` upstream.) |
| PITA `RNAduplex` + `RNAddG4` | **Shipped binaries are Linux x86-64 ELF** (`Bin/ViennaRNA/ViennaRNA-1.6/Progs/`) and cannot run here. I rebuilt both from the bundled ViennaRNA-1.6 C source out-of-tree in `deps/pita-vienna` (after deleting stale Linux `*.o`/`*.a`) and staged the arm64 Mach-O binaries in `deps/pita-bin-arm64/`: `RNAduplex` sha256 `b01a67068a5168f87df18a967c72398d0ea4d6ff6847f055ee6b743aea2c0df8`, `RNAddG4` sha256 `ada0a5af757ae3a086adda5dc5968d1b7703bcb352c706f3f61e94e64f28c607`. `RNAforester` (C++) does not build with modern clang and is not used. |
| `cd off-target/pita && make install` | Succeeds with GNU sed. |
| TargetScan reference data | present in-tree (`off-target/targetscan/`, `off-target/ref/100_human_{UTR,ORF}.txt`). |

**Documented procedure for an off-target run (must be repeated; I reverted it):**
1. `cp deps/pita-bin-arm64/{RNAduplex,RNAddG4} repo/OligoFormer/off-target/pita/Bin/ViennaRNA/ViennaRNA-1.6/Progs/` (+`chmod +x`)
2. `cd repo/OligoFormer/off-target/pita && PATH=<PROJECT_ROOT>/deps/offtarget-env/bin:$PATH make install`
3. run with `PATH=<PROJECT_ENV>/bin:<PROJECT_ROOT>/deps/offtarget-env/bin:$PATH` from `repo/OligoFormer`
4. afterwards `make clean` and `git checkout -- off-target` to restore pinned source.

**Exactly which pinned files `make install` + binary staging change** (measured with `git status`/`git diff`):
- content: 25 files — `off-target/pita/pita_prediction.pl` and `off-target/pita/lib/{add_column,body,cap,cut,extract_sequence,fasta2stab,find_potential_mirna_targets,format_number,join,libfile,libtable,load_args,merge_columns,modify_column,pita_run,RNAddG_compute,sequence_helpers,stab2length,tab2feature_gxt,transpose,uniquify,average_rows,libstats,ltqnorm}.pl` — the literal `EXE_BASE_DIR` is replaced by the absolute path of `off-target/pita` (95 insertions / 95 deletions).
- mode: ~676 further files under `off-target/pita/{lib,Bin}` from `chmod -R 755` (mode-only).
- my staging additionally replaces 2 binaries under `Bin/ViennaRNA/ViennaRNA-1.6/Progs/`.
- `make clean` reverses the content substitution; `git checkout -- off-target/pita` restored everything (verified clean). Diff stat tail kept at `tmp/outputs/env_check/pita_install_diffstat.txt`.

Verified end-to-end (#20, ~60 s wall): `RNA0.txt` carried real `pita_score = -5.47`,
`targetscan_score = 0.414`, `off_target_filter = 0` with an error-free log.

**Toxicity (`-tox`): CONFIRMED trivial.** `infer.py` only reads `./toxicity/cell_viability.txt`
(present in-tree) and left-joins on the siRNA seed (positions 2-7); no perl, no ViennaRNA, no
network. Run #21 took ~12 s and produced `Seed=CGCAAA`, `cell_viability=75.0`, `toxicity_filter=0`.
`-tox` needs **only** PROJECT_ENV on PATH (RNA-FM still runs, as always for `--infer`).

**`-a / --all_human` is broken independently of the environment:** `infer.py` sets
`Args.orf = './off-target/ref/human_UTR.txt'` (UTR, a copy-paste bug) and neither
`human_UTR.txt` nor `human_ORF.txt` exists — only the `.zip` archives are shipped, so they must be
unzipped first. Recorded as a repo bug for the tool-wrapping agents; not fixed by me.

## 5. Mismatch module (`-m`)

- `scripts/mismatch.py` loads `model/mismatch_model.pth` **at import time** (module level,
  `torch.load` without `map_location`) — i.e. on *every* `scripts/main.py` invocation, including
  `--infer`. It loads fine on CPU here.
- **It is interactive: `-m 1` calls `input()` twice** ("please input base RNA:", "please input
  mismatched RNA:"), `-m 2/3` calls `input()` once. There is no flag to pass the sequences, so an MCP
  wrapper must feed them on **stdin** (verified: `printf 'R1\nR2\n' | python scripts/main.py -m 1`
  → `0.9940418601036072`, exit 0).
- `-m 2/3` writes to `Args.output_dir + '/result/mismatch.txt'`, i.e. `result//result/mismatch.txt`
  by default — the `<output_dir>/result/` directory must exist beforehand or the write fails.
  `-m 3` enumerates 1/2/3-mismatch variants of a 19-mer and will be slow.

## 6. Blockers / unresolved

1. **Both Tsinghua Cloud downloads in the README are dead** (RNA-FM.tar.gz and PerlLib.zip) — both
   return a "File does not exist" HTML page. Worked around via ml4bio/RNA-FM + HuggingFace weights
   and a conda perl prefix. Record in the server docs; no tarball hash can ever be reported.
2. **RNA-FM's `predict.py` cannot run CPU-only unpatched**; the one-line patch above lives in the
   untracked `RNA-FM/` dependency and must ship with the environment (or `--device cpu` must be
   injected, which would require editing the tracked `scripts/RNA-FM.sh`).
3. **PITA ships Linux x86-64 binaries**; the arm64 rebuild in `deps/pita-bin-arm64/` must be staged
   into the repo for every off-target run, and `make install` dirties 25 pinned perl files +
   ~676 file modes. Any source-identity check must run after `make clean` + `git checkout`.
4. **No CUDA**; training (`~60 min on a V100` per README) is out of scope on this host. Inference,
   testing, off-target, toxicity and mismatch all have working CPU paths.
5. **`-a/--all_human` is broken** (wrong ORF path + unextracted `human_{UTR,ORF}.txt.zip`).
6. Efficacy values vary ~1% run to run on CPU (see §3) — reference outputs need tolerances.

## 7. Paths I own

- `OligoFormer-env/` (PROJECT_ENV, kernel `oligoformer-mcp`)
- `deps/offtarget-env/` (conda: perl+BioPerl+Statistics::Lite+ViennaRNA 2.7.2+GNU sed/grep)
- `deps/pita-vienna/` (out-of-tree ViennaRNA-1.6 build), `deps/pita-bin-arm64/` (RNAduplex, RNAddG4)
- `deps/torch_home/` (hub checkpoint symlink), `deps/rnafm_cpu_default.patch`
- `deps/downloads/{rnafm,rnafm_git,rnafm_weights_hf,perllib}/` (untrusted downloads, each isolated)
- `repo/OligoFormer/RNA-FM/` (untracked dependency + weights)
- `tmp/outputs/env_check/` (probe/smoke logs and outputs)
- `reports/environment-manager_results.md`, `reports/environment-requirements.txt`, `reports/python-version.txt`

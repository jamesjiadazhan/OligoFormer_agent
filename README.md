# OligoFormer_agent

An MCP server and a paper skill generated from **OligoFormer** — Bai Y, Zhong H, Wang T, Lu ZJ, "OligoFormer: an accurate and robust prediction method for siRNA design", *Bioinformatics* 40(10):btae577 (2024) — and its code repository [lulab/OligoFormer](https://github.com/lulab/OligoFormer) at commit `e2f53ad63387bbe166bf123949151e2bc9bf6ec3`.

Both parts were built by Paper2Agent: five MCP tools that call OligoFormer's own command line, and a reading package of the paper plus its supplementary material.

## What's here

| Path | What it is |
| --- | --- |
| `src/` | The MCP server (`OligoFormer_mcp.py`), the tool module (`tools/cli_wrapper.py`) and pinned runtime requirements |
| `deps/` | Scripts that build the scientific runtime: OligoFormer + RNA-FM + the PITA/TargetScan toolchain |
| `USAGE.md` | **Full installation and usage guide** — read this first |
| `skill/oligoformer-paper/` | The paper skill: paper and supplement as searchable Markdown, 42 figure/equation crops, 16 tables as CSV |
| `dist/OligoFormer-mcp.zip` | The exact archive that passed delivery validation (sha256 `dde0644144f5c38f23c826516c1031ff5e2656a7ca2d4185c11bf060f6f1655f`) |
| `tests/`, `notebooks/`, `reports/` | Validation evidence: the test suite, upstream reference outputs, and every verification report |

## The five tools

| Tool | What it does |
| --- | --- |
| `predict_sirna_efficacy_from_mrna` | Ranks every 19-nt siRNA window of a target mRNA |
| `predict_specific_sirna_efficacy` | Scores siRNA candidates you already have |
| `screen_sirna_offtarget` | Adds PITA + TargetScan off-target screening (slow: ~7 min) |
| `compare_sirna_mismatch_kon` | Relative on-rate of a mismatched siRNA vs its base sequence |
| `scan_sirna_mismatch_variants` | Scores every 1–2 (or 1–3) mismatch variant of an siRNA |

Each one runs OligoFormer's `scripts/main.py` and returns its native tables. No scientific logic was reimplemented.

## Quick start

```bash
git clone https://github.com/jamesjiadazhan/OligoFormer_agent.git
cd OligoFormer_agent
chmod +x deps/*.sh

uv venv .venv --python 3.11 --seed      # or: python3.11 -m venv .venv
./.venv/bin/python -m pip install -r src/requirements.txt

# Builds runtime/OligoFormer: clones OligoFormer + RNA-FM at pinned commits and
# downloads the 1.2 GB RNA-FM checkpoint (sha256-verified).
./deps/create_offtarget_env.sh "$PWD/deps/offtarget-env"   # only needed for screen_sirna_offtarget
PROJECT_ENV="$PWD/.venv" ./deps/setup_oligoformer_runtime.sh "$PWD/runtime/OligoFormer"

./.venv/bin/python src/OligoFormer_mcp.py                   # start the server
```

Install at a path **without spaces**: OligoFormer's off-target scripts pass paths to an unquoted shell, so `screen_sirna_offtarget` fails under a path containing a space. The other four tools are unaffected.

See [`USAGE.md`](USAGE.md) for parameters, output columns, timeouts, client registration and the full prerequisite list.

## Two things to know before you trust a number

**Efficacy scores are not reproducible.** OligoFormer's `scripts/infer.py` never calls `model.eval()`, so dropout stays active during prediction, and it never seeds torch. Repeated runs on identical input differ by up to ~0.18 absolute (~30% relative) and the ranked order changes. This is upstream behaviour and was deliberately left alone — "fixing" it would change the published method. Every other column is exactly reproducible, and the PITA/TargetScan tables are byte-identical between runs. Treat efficacy as a ranking signal, not an exact value.

**Two upstream options are deliberately not exposed.** `-top N` writes its FASTA with the top-ranked headers but positionally-ordered sequences, so off-target scores land on the wrong rows. `-a/--all_human` exits successfully while producing all-NaN scores, silently passing every siRNA. Training and evaluation modes are also excluded: one needs a CUDA GPU, the other needs embedding files the repository doesn't ship. `USAGE.md` has the details.

## Reproducing the validation

```bash
./.venv/bin/python -m pytest tests/code/cli_wrapper -m "not slow"   # 105 fast tests
./.venv/bin/python -m pytest tests/code/cli_wrapper                 # add the 3 slow ones
```

The tests resolve paths relative to this repository root and compare the tools against the saved upstream reference outputs in `notebooks/`. They need the runtime built first. Reference logs under `tests/logs/` and `reports/` contain absolute paths from the machine that produced them; that is recorded evidence and is left unmodified.

What was validated: 108 tests against direct upstream runs of the same commands; 18 stdio acceptance cases passing in three separate environments (development, a clean install, and the extracted ZIP); PITA/TargetScan tables byte-identical to upstream; all deterministic columns exactly equal. This covers the five tools above — not the whole repository, and not the paper's published results.

Tested only on macOS 15 (Apple Silicon, arm64) with Python 3.11.12 and torch 2.2.1 on CPU. Linux is untested and takes a different PITA code path. The off-target toolchain is an osx-64 conda prefix, so it runs under Rosetta 2 on Apple Silicon.

## Licensing and attribution

This repository contains **no OligoFormer source code or model weights** — `deps/setup_oligoformer_runtime.sh` fetches those from the upstream repository at the pinned commit when you install.

- **OligoFormer** is free for academic research, non-profit and personal use; commercial use needs written authorization from Tsinghua University (`ott@tsinghua.edu.cn`). See `THIRD_PARTY_NOTICES/OligoFormer-LICENSE`. That restriction governs anything you do with these tools.
- **The paper content** in `skill/oligoformer-paper/` is reproduced from the open-access article under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). © The Author(s) 2024, published by Oxford University Press. Cite the original article.
- **RNA-FM**, **PITA** and **TargetScan** carry their own licenses; RNA-FM's is in `THIRD_PARTY_NOTICES/`.

```bibtex
@article{bai2024oligoformer,
  title={OligoFormer: an accurate and robust prediction method for siRNA design},
  author={Bai, Yilan and Zhong, Haochen and Wang, Taiwei and Lu, Zhi John},
  journal={Bioinformatics}, volume={40}, number={10}, pages={btae577}, year={2024},
  doi={10.1093/bioinformatics/btae577}
}
```

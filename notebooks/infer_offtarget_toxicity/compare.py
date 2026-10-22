"""Post-run checks for infer_offtarget_toxicity (reads saved outputs only; no algorithm re-implementation).
Usage: python compare.py <notebook_dir>  -> prints JSON summary and writes outputs/comparison.json"""
import json, sys, re, hashlib
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

N = Path(sys.argv[1])
OUT = N / "outputs"; LOGS = N / "logs"
RECS = ["RNA0", "RNA1"]
RUNS = ["run1", "run2", "run3", "replay"]
DET = ["pos", "sense", "siRNA", "func_filter", "pita_score", "targetscan_score", "off_target_filter",
       "Seed", "cell_viability", "toxicity_filter", "filter"]
PITA_HDR = "RefSeq\tmicroRNA\tSites\tScore"  # header of PITA <prefix>_pita_results_targets.tab
ERR_PAT = re.compile(r"(error|fatal|no such file|cannot open|usage:|invalid option|can't locate|died|segmentation)", re.I)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def log_scan(label):
    hits = []
    for kind in ("stdout", "stderr"):
        p = LOGS / f"{label}.{kind}.log"
        for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
            if ERR_PAT.search(line):
                hits.append(f"{kind}:{i}: {line.strip()[:200]}")
    return hits


def tab_check(label, rec):
    d = OUT / label / "infer" / rec
    p, t = d / "pita.tab", d / "targetscan.tab"
    res = {}
    if p.exists():
        lines = p.read_text().splitlines()
        res["pita_tab"] = {"rows": len(lines) - 1, "header_ok": bool(lines) and lines[0] == PITA_HDR,
                           "header": lines[0] if lines else None, "sha256": sha(p)}
    if t.exists():
        lines = [l for l in t.read_text().splitlines() if l.strip()]
        res["targetscan_tab"] = {"rows": len(lines), "cols_ok": all(len(l.split("\t")) == 3 for l in lines) and len(lines) > 0,
                                 "sha256": sha(t)}
    return res


summary = {"log_scan": {}, "tabs": {}, "tables": {}, "comparisons": {}}
for label in RUNS + ["changed_params", "readme_all_human"]:
    summary["log_scan"][label] = log_scan(label)
    summary["tabs"][label] = {r: tab_check(label, r) for r in RECS}
    summary["tables"][label] = {}
    for r in RECS:
        df = pd.read_csv(OUT / label / f"{r}.txt", sep="\t")
        summary["tables"][label][r] = {
            "shape": list(df.shape), "columns": list(df.columns),
            "n_filter0": int((df["filter"] == 0).sum()),
            "n_off_target_filter1": int((df["off_target_filter"] == 1).sum()),
            "n_pita_nan": int(df["pita_score"].isna().sum()), "n_ts_nan": int(df["targetscan_score"].isna().sum()),
            "first_row": df.iloc[0].to_dict(),
        }

base = "run1"
for other in ["run2", "run3", "replay"]:
    comp = {}
    for r in RECS:
        a = pd.read_csv(OUT / base / f"{r}.txt", sep="\t"); b = pd.read_csv(OUT / other / f"{r}.txt", sep="\t")
        det_equal = {c: bool(a[c].equals(b[c])) for c in DET}
        ra = pd.read_csv(OUT / base / f"{r}_ranked.txt", sep="\t"); rb = pd.read_csv(OUT / other / f"{r}_ranked.txt", sep="\t")
        fa = pd.read_csv(OUT / base / f"{r}_ranked_filtered.txt", sep="\t"); fb = pd.read_csv(OUT / other / f"{r}_ranked_filtered.txt", sep="\t")
        tabs_equal = {f: sha(OUT / base / "infer" / r / f) == sha(OUT / other / "infer" / r / f)
                      for f in ("pita.tab", "targetscan.tab", "siRNA.fa")}
        comp[r] = {
            "deterministic_columns_equal": det_equal,
            "efficacy_max_abs_diff": float((a["efficacy"] - b["efficacy"]).abs().max()),
            "efficacy_max_rel_diff": float(((a["efficacy"] - b["efficacy"]).abs() / a["efficacy"].abs()).max()),
            "efficacy_spearman": float(spearmanr(a["efficacy"], b["efficacy"]).correlation),
            "ranked_order_identical": bool(list(ra["pos"]) == list(rb["pos"])),
            "ranked_filtered_set_identical": bool(set(fa["pos"]) == set(fb["pos"])),
            "ranked_filtered_order_identical": bool(list(fa["pos"]) == list(fb["pos"])),
            "intermediate_files_sha256_equal": tabs_equal,
        }
    summary["comparisons"][f"{base}_vs_{other}"] = comp

# changed params vs run1: deterministic pita/targetscan scores equal, filters recomputed
ch = {}
for r in RECS:
    a = pd.read_csv(OUT / "run1" / f"{r}.txt", sep="\t"); c = pd.read_csv(OUT / "changed_params" / f"{r}.txt", sep="\t")
    ch[r] = {
        "scores_equal": bool(a["pita_score"].equals(c["pita_score"]) and a["targetscan_score"].equals(c["targetscan_score"])),
        "off_target_filter1_run1": int((a["off_target_filter"] == 1).sum()),
        "off_target_filter1_changed": int((c["off_target_filter"] == 1).sum()),
        "rows_off_target_filter_changed": int((a["off_target_filter"] != c["off_target_filter"]).sum()),
        "filter_consistent_with_thresholds": bool(all(
            int(o) == int((p < -5) or (t > 0.2)) for o, p, t in zip(c["off_target_filter"], c["pita_score"], c["targetscan_score"]))),
        "filter_sum_consistent": bool((c["filter"] == c["func_filter"] + c["off_target_filter"] + c["toxicity_filter"]).all()),
    }
summary["comparisons"]["run1_vs_changed_params"] = ch
(OUT / "comparison.json").write_text(json.dumps(summary, indent=1, default=str))
print(json.dumps(summary, indent=1, default=str))

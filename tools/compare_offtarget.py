"""Coordinator check: compare a wrapper off-target output dir with the upstream reference run1."""
import hashlib, sys, pathlib
import pandas as pd
out, ref = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
DET = ["pos", "sense", "siRNA", "func_filter", "pita_score", "targetscan_score", "off_target_filter",
       "Seed", "cell_viability", "toxicity_filter", "filter"]
ok = True
for rec in ("RNA0", "RNA1"):
    for tab in ("pita.tab", "targetscan.tab"):
        same = sha(out / "offtarget" / rec / tab) == sha(ref / "infer" / rec / tab)
        ok &= same; print(rec, tab, "identical" if same else "DIFFERENT")
    w = pd.read_csv(next(out.rglob(f"{rec}.txt")), sep="\t"); r = pd.read_csv(ref / f"{rec}.txt", sep="\t")
    cols = [c for c in DET if c in w.columns and c in r.columns]
    eq = w[cols].equals(r[cols]); ok &= eq
    d = (w["efficacy"] - r["efficacy"]).abs().max(); rho = w["efficacy"].corr(r["efficacy"], method="spearman")
    print(rec, "rows", len(w), "deterministic cols", cols, "equal" if eq else "DIFFER", f"efficacy max|d|={d:.3f} spearman={rho:.3f}")
    ok &= d <= 0.25 and rho >= 0.9
print("PASS" if ok else "FAIL")

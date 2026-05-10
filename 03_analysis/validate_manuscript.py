#!/usr/bin/env python3
"""
manuscript validation script (robust)
covers:
  V1. unassigned ARG bias (mobility + drug-class string compare)
  V2. Enterobacteriaceae effluent - ARG-bearing vs mobile distinction (sample-level + pooled)
  V3. Gini pooled vs per-sample (mobile + ARG-bearing)
  V4. geNomad threshold sensitivity (find plasmid/virus score tables; summarize score distribution)
  V5. contig_length for assigned vs unassigned (if contig_length exists OR optional external contig length table)

run:
  module load anaconda3/2023.09-0
  source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
  python3 validate_manuscript.py > validation_report.txt 2>&1
"""

import os, re, glob
import pandas as pd
import numpy as np
from scipy import stats

RGI_CSV  = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv"
META_TSV = "/QRISdata/Q6636/sra_ww_mobilization/results/sample_map_complete.tsv"

# If you have a global contig-length table, set it here (optional)
# expected columns: Contig (or seq_name) + length
CONTIG_LEN_TSV = ""  # e.g. "/QRISdata/Q6636/sra_ww_mobilization/results/contig_length.tsv"

# geNomad search base (may vary across your project)
GENOMAD_BASE = "/QRISdata/Q6636"

SEP = "\n" + "="*70 + "\n"

ENTERO = {"Escherichia","Klebsiella","Enterobacter","Citrobacter",
          "Salmonella","Serratia","Proteus","Morganella","Providencia"}

KW8 = ["glycopeptide","tetracycline","macrolide","aminoglycoside",
       "beta-lactam","fluoroquinolone","sulfonamide","phenicol"]


# ----------------------------- helpers
def pick_col(df, cands):
    for c in cands:
        if c in df.columns:
            return c
    return None

def gini(x):
    x = np.array(x, dtype=float)
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    x.sort()
    n = len(x)
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))

def pooled_gini(sub, col="genus"):
    vc = sub[col].value_counts()
    return gini(vc.values), len(sub), len(vc)

def per_sample_gini_df(sub, col="genus"):
    rows = []
    for sid, g in sub.groupby("sample_id"):
        vc = g[col].value_counts()
        rows.append({
            "sample_id": sid,
            "gini": gini(vc.values),
            "n_hits": len(g),
            "n_genera": len(vc)
        })
    return pd.DataFrame(rows)

def binom_upper_95_zero_success(n):
    # one-sided upper bound if k=0 among n trials
    if n <= 0:
        return np.nan
    return 1 - (0.05 ** (1 / n))

def fmt_pct(x, nd=2):
    if pd.isna(x):
        return "NA"
    return f"{x*100:.{nd}f}%"

def safe_str_series(s):
    return s.fillna("").astype(str)

def find_genomad_score_tables(base):
    """
    Try multiple common geNomad outputs; return candidate tsv files.
    """
    patterns = [
        os.path.join(base, "**", "*plasmid*summary*.tsv"),
        os.path.join(base, "**", "*virus*summary*.tsv"),
        os.path.join(base, "**", "*contigs_summary*.tsv"),
        os.path.join(base, "**", "*scores*.tsv"),
        os.path.join(base, "**", "*genomad*", "*.tsv"),
    ]
    files = set()
    for pat in patterns:
        for f in glob.glob(pat, recursive=True):
            bn = os.path.basename(f).lower()
            # keep likely summary/score tables but avoid huge unrelated
            if any(k in bn for k in ["plasmid","virus","contig","summary","score","scores"]) and bn.endswith(".tsv"):
                files.add(f)
    return sorted(files)


# ----------------------------- load
print("Loading data...")
df = pd.read_csv(RGI_CSV, low_memory=False)
meta = pd.read_csv(META_TSV, sep="\t")

# ensure sample_id + category
if "sample_id" not in df.columns:
    raise SystemExit("ERROR: rgi table missing sample_id column.")
if "category" not in df.columns:
    df = df.merge(meta[["sample_id","category"]], on="sample_id", how="left")

# mobility
if ("on_plasmid" not in df.columns) or ("on_virus" not in df.columns):
    raise SystemExit("ERROR: rgi table missing on_plasmid / on_virus columns.")
df["is_mobile"] = (df["on_plasmid"].fillna(False) | df["on_virus"].fillna(False)).astype(int)

# de-dup keys (robust)
keys = [c for c in ["sample_id","Contig","Best_Hit_ARO"] if c in df.columns]
before = len(df)
if keys:
    df = df.drop_duplicates(subset=keys)
print(f"De-dup: {before} -> {len(df)} rows (keys={keys if keys else 'NONE'})")

# genus assigned subset
if "genus" not in df.columns:
    raise SystemExit("ERROR: rgi table missing genus column.")
df_assigned   = df[df["genus"].notna()].copy()
df_unassigned = df[df["genus"].isna()].copy()

# detect drug class column robustly
drug_col = pick_col(df, ["Drug Class","drug_class","Drug Class (Summary)","Drug Class (ARO)"])
if drug_col is None:
    print("WARNING: drug class column not found; V1 drug-class comparison will be skipped.")


# ================================================================== V1
print(SEP + "V1: UNASSIGNED ARG BIAS CHECK (assigned vs unassigned)")

a_mob = df_assigned["is_mobile"].mean()
u_mob = df_unassigned["is_mobile"].mean()
print(f"Assigned   mobility: {a_mob*100:.2f}%  (n={len(df_assigned)})")
print(f"Unassigned mobility: {u_mob*100:.2f}%  (n={len(df_unassigned)})")
print(f"Delta (unassigned-assigned): {(u_mob - a_mob)*100:+.2f}%")

# chi-square test (assigned vs mobile)
ct = pd.crosstab(df["genus"].notna(), df["is_mobile"])
chi2, p, dof, exp = stats.chi2_contingency(ct)
print(f"Chi2 test (host-assigned vs mobile): chi2={chi2:.3f}, p={p:.4g}")
print("CONCLUSION:", "NO significant bias (mobility independent of assignment)" if p > 0.05 else "WARNING: significant bias detected")

# drug class keyword compare (note: strings can contain multiple classes)
if drug_col:
    print("\nDrug class keyword comparison (strings; fractions can sum >100%):")
    print(f"{'keyword':<16} {'assigned%':>10} {'unassigned%':>12}  note")
    for kw in KW8:
        a_pct = safe_str_series(df_assigned[drug_col]).str.lower().str.contains(kw).mean() * 100
        u_pct = safe_str_series(df_unassigned[drug_col]).str.lower().str.contains(kw).mean() * 100
        note = "<-- differ >5pp" if abs(a_pct - u_pct) > 5 else ""
        print(f"{kw:<16} {a_pct:>9.2f}% {u_pct:>11.2f}%  {note}")
else:
    print("\n[SKIP] No drug class column detected.")


# ================================================================== V5 (contig length)
print(SEP + "V5: CONTIG LENGTH (assigned vs unassigned)")

len_col = None
for c in ["contig_length","Contig_length","contig_len","Contig_len","sequence_length","Sequence length","length"]:
    if c in df.columns:
        len_col = c
        break

if len_col:
    for label, sub in [("Assigned", df_assigned), ("Unassigned", df_unassigned)]:
        x = pd.to_numeric(sub[len_col], errors="coerce").dropna()
        if len(x) == 0:
            print(f"{label}: contig length column present but no numeric values.")
        else:
            print(f"{label}: n={len(x)} median={np.median(x):.0f} IQR=({np.quantile(x,0.25):.0f}-{np.quantile(x,0.75):.0f}) mean={np.mean(x):.0f}")
else:
    print("No contig_length column in rgi table.")
    if CONTIG_LEN_TSV and os.path.exists(CONTIG_LEN_TSV):
        print(f"External contig-length table detected: {CONTIG_LEN_TSV}")
        contig_col = pick_col(df, ["Contig","contig","contig_id","seq_name"])
        if contig_col is None:
            print("WARNING: cannot join external length table because contig id column not found.")
        else:
            clen = pd.read_csv(CONTIG_LEN_TSV, sep="\t")
            clen_id = pick_col(clen, ["Contig","contig","contig_id","seq_name"])
            clen_len = pick_col(clen, ["length","contig_length","len"])
            if clen_id and clen_len:
                tmp = df[[contig_col,"genus"]].merge(clen[[clen_id, clen_len]], left_on=contig_col, right_on=clen_id, how="left")
                tmp["assigned"] = tmp["genus"].notna()
                for label, sub in [("Assigned", tmp[tmp["assigned"]]), ("Unassigned", tmp[~tmp["assigned"]])]:
                    x = pd.to_numeric(sub[clen_len], errors="coerce").dropna()
                    print(f"{label}: n={len(x)} median={np.median(x):.0f} IQR=({np.quantile(x,0.25):.0f}-{np.quantile(x,0.75):.0f}) mean={np.mean(x):.0f}")
            else:
                print("WARNING: external length table missing contig id or length columns.")
    else:
        print("No external contig length source provided; V5 limited.")


# ================================================================== V2
print(SEP + "V2: ENTEROBACTERIACEAE EFFLUENT - ARG-BEARING vs MOBILE")

# Note: use host-assigned subset for genus-based Entero check
eff = df_assigned[df_assigned["category"] == "ww_effluent_municipal"].copy()
inf = df_assigned[df_assigned["category"] == "ww_influent_municipal"].copy()

def entero_block(label, sub_all_assigned, raw_df):
    # raw sample counts (including non-assigned)
    raw_cat = raw_df[raw_df["category"] == ("ww_effluent_municipal" if "EFFLUENT" in label else "ww_influent_municipal")]
    raw_samples = raw_cat["sample_id"].nunique()

    n_samp = sub_all_assigned["sample_id"].nunique()
    n_rows = len(sub_all_assigned)
    mob_rows = int(sub_all_assigned["is_mobile"].sum())

    ent_argb = sub_all_assigned[sub_all_assigned["genus"].isin(ENTERO)]
    ent_mob  = sub_all_assigned[sub_all_assigned["genus"].isin(ENTERO) & (sub_all_assigned["is_mobile"] == 1)]

    samp_with_ent_argb = sub_all_assigned.groupby("sample_id").apply(lambda x: x["genus"].isin(ENTERO).any()).sum()
    samp_with_ent_mob  = sub_all_assigned.groupby("sample_id").apply(lambda x: (x["genus"].isin(ENTERO) & (x["is_mobile"]==1)).any()).sum()

    print(f"\n--- {label} ---")
    print(f"Raw samples in category (all data):         {raw_samples}")
    print(f"Samples with host-assigned ARG rows:        {n_samp}")
    print(f"ARG-host rows (host-assigned):              {n_rows}")
    print(f"Total mobile hits (within host-assigned):   {mob_rows}")
    print(f"Entero ARG-bearing rows:                    {len(ent_argb)} ({len(ent_argb)/n_rows*100:.2f}% of ARG-host rows)" if n_rows else "Entero ARG-bearing rows: NA")
    print(f"Entero MOBILE rows:                         {len(ent_mob)} ({len(ent_mob)/mob_rows*100:.2f}% of mobile hits)" if mob_rows else "Entero MOBILE rows: NA")
    print(f"Samples with any Entero ARG-bearing:        {int(samp_with_ent_argb)}/{n_samp}")
    print(f"Samples with any Entero mobile:             {int(samp_with_ent_mob)}/{n_samp}")

    if mob_rows > 0 and len(ent_mob) == 0:
        p_upper = binom_upper_95_zero_success(mob_rows)
        print(f"[Binomial 95% upper bound] Entero mobile proportion < {p_upper*100:.2f}% (0/{mob_rows})")

    # mobility rate within Entero vs non-Entero (ARG-bearing denominator)
    if len(ent_argb) > 0:
        print(f"Entero mobility rate (within Entero ARG rows): {len(ent_mob)/len(ent_argb)*100:.2f}%")
    nonent = sub_all_assigned[~sub_all_assigned["genus"].isin(ENTERO)]
    if len(nonent) > 0:
        nonent_mob = nonent[nonent["is_mobile"]==1]
        print(f"Non-Entero mobility rate:                     {len(nonent_mob)/len(nonent)*100:.2f}%")

entero_block("INFLUENT (municipal)", inf, df)
entero_block("EFFLUENT (municipal)", eff, df)

print("\nINTERPRETATION GUIDE:")
print("  - If effluent has Entero ARG-bearing >0 but Entero mobile=0: strong claim is 'mobilome-level absence'.")
print("  - If effluent Entero ARG-bearing is near 0: then '0 mobile' is largely driven by low presence; soften wording.")


# ================================================================== V3
print(SEP + "V3: GINI - POOLED vs PER-SAMPLE (pooled is what paper currently uses)")

blocks = [
    ("INFLUENT mobile", inf[inf["is_mobile"]==1]),
    ("EFFLUENT mobile", eff[eff["is_mobile"]==1]),
    ("INFLUENT ARG-bearing", inf),
    ("EFFLUENT ARG-bearing", eff),
]

for label, sub in blocks:
    pg, n_hits, n_genera = pooled_gini(sub, col="genus")
    ps_df = per_sample_gini_df(sub, col="genus")
    ps_df = ps_df[ps_df["n_hits"] > 0]
    if len(ps_df) == 0:
        print(f"\n[{label}] no hits")
        continue

    mean_unw = ps_df["gini"].mean()
    mean_w   = np.average(ps_df["gini"], weights=np.maximum(ps_df["n_hits"].values, 1))
    med      = ps_df["gini"].median()
    q25, q75 = ps_df["gini"].quantile([0.25, 0.75])
    pct_zero = (ps_df["gini"] == 0).mean() * 100

    print(f"\n[{label}]")
    print(f"Pooled Gini:                 {pg:.3f} (n_hits={n_hits}, genera={n_genera})")
    print(f"Per-sample mean (unweighted): {mean_unw:.3f}")
    print(f"Per-sample mean (weighted):   {mean_w:.3f}")
    print(f"Per-sample median:            {med:.3f}  IQR=[{q25:.3f},{q75:.3f}]")
    print(f"Samples with Gini=0:          {pct_zero:.1f}%")

print("\nNOTE FOR METHODS:")
print("  Pooled Gini = concentration of the aggregated mobilome across the dataset.")
print("  Per-sample Gini = typical within-sample concentration; report in Supplement if needed.")


# ================================================================== V4
print(SEP + "V4: geNomad THRESHOLD SENSITIVITY (score distribution only unless join exists)")

cand = find_genomad_score_tables(GENOMAD_BASE)
print(f"Candidate geNomad-related TSV files found: {len(cand)}")
if len(cand) == 0:
    print("WARNING: No geNomad TSV found under base. Check GENOMAD_BASE.")
else:
    print("Examples (first 8):")
    for f in cand[:8]:
        print("  ", f)

    # Try to load plasmid/virus score tables if columns exist
    loaded = []
    for f in cand:
        bn = os.path.basename(f).lower()
        if ("plasmid" in bn or "virus" in bn or "score" in bn) and os.path.getsize(f) < 200_000_000:
            try:
                t = pd.read_csv(f, sep="\t")
            except Exception:
                continue
            # find score columns
            pl = pick_col(t, ["plasmid_score","plasmidScore","score_plasmid"])
            vi = pick_col(t, ["virus_score","virusScore","score_virus"])
            sid = pick_col(t, ["sample_id","sample","run","Run","SRA"])
            seq = pick_col(t, ["seq_name","contig","contig_id","sequence","name"])
            if pl or vi:
                keep = []
                if seq: keep.append(seq)
                if sid: keep.append(sid)
                if pl: keep.append(pl)
                if vi: keep.append(vi)
                tt = t[keep].copy()
                tt["__file__"] = f
                loaded.append(tt)
            if len(loaded) >= 30:  # cap for speed
                break

    if not loaded:
        print("Could not parse any score tables with plasmid_score/virus_score columns.")
    else:
        gdf = pd.concat(loaded, ignore_index=True)
        print(f"\nLoaded score rows (sampled): {len(gdf)} from {len(loaded)} files")

        if "plasmid_score" in gdf.columns:
            s = pd.to_numeric(gdf["plasmid_score"], errors="coerce").dropna()
            if len(s)>0:
                print(f"[plasmid_score] mean={s.mean():.3f}, median={s.median():.3f}, q90={s.quantile(0.90):.3f}")
                for thr in [0.5,0.6,0.7,0.8]:
                    n = int((s>=thr).sum())
                    print(f"  plasmid_score >= {thr:.1f}: {n}/{len(s)} ({n/len(s)*100:.1f}%)")
        if "virus_score" in gdf.columns:
            s = pd.to_numeric(gdf["virus_score"], errors="coerce").dropna()
            if len(s)>0:
                print(f"[virus_score] mean={s.mean():.3f}, median={s.median():.3f}, q90={s.quantile(0.90):.3f}")
                for thr in [0.5,0.6,0.7,0.8]:
                    n = int((s>=thr).sum())
                    print(f"  virus_score >= {thr:.1f}: {n}/{len(s)} ({n/len(s)*100:.1f}%)")

        print("\nNOTE:")
        print("  This V4 currently summarizes score distributions only.")
        print("  To compute how mobility rate changes at thresholds (0.5/0.6/0.7),")
        print("  we need a consistent contig-id join between RGI contigs and geNomad contig scores for all samples.")


# ================================================================== SUMMARY
print(SEP + "SUMMARY FOR MANUSCRIPT (template statements; fill with numbers above)")

print("""
V1 UNASSIGNED BIAS:
  -> If chi-square p > 0.05 AND mobility% similar:
     'Unassigned ARG hits showed comparable mobility rates to host-assigned hits,
      suggesting limited systematic bias in host-assigned mobilome conclusions.'
  -> If p <= 0.05 OR big delta:
     'Unassigned ARG hits differed in mobility, indicating potential host-assignment bias;
      we therefore interpret host-resolved mobilome patterns cautiously and report sensitivity analyses.'

V2 ENTEROBACTERIACEAE IN EFFLUENT:
  -> If effluent Entero ARG-bearing >0 but Entero mobile=0:
     'Enterobacteriaceae remained detectable among ARG-bearing hosts, but their mobile ARGs were absent
      (0/k mobile hits; 95% upper bound < X%), indicating mobilome-level loss rather than complete host removal.'
  -> If effluent Entero ARG-bearing ~0:
     'Enterobacteriaceae were rare among effluent ARG-bearing hosts in this dataset, consistent with the lack of
      Enterobacteriaceae-associated mobile ARGs.'

V3 GINI METHODOLOGY:
  -> Add to Methods: 'Gini coefficients were computed on pooled, sample-aggregated genus counts, capturing
     global mobilome concentration across the dataset. Per-sample Gini distributions are provided in Supplement.'

V4 GENOMAD THRESHOLD:
  -> Use as Supplement if you can complete the contig-level join:
     report mobility% under score>=0.5/0.6/0.7 and show rank robustness for MBI.

[Validation complete]
""")

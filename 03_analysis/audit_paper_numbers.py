#!/usr/bin/env python3
"""
audit_paper_numbers.py
================================================================
Bunya 上的论文数据核对脚本 - 一次性核对所有关键数字

依赖文件:
  /QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv
  /QRISdata/Q6636/sra_ww_mobilization/results/sample_map_complete.tsv
  /QRISdata/Q6636/sra_ww_mobilization/results/bracken_genus_abundance.csv

运行方式:
  module load anaconda3/2023.09-0
  source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
  python3 audit_paper_numbers.py 2>&1 | tee ~/audit_log.txt

输出:
  ~/paper_audit_results.tsv   - 详细对照表 (Excel 可打开)
  ~/paper_audit_summary.txt   - 不匹配项汇总
  ~/audit_log.txt             - 完整运行日志
================================================================
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
from collections import Counter

# ======================== CONFIG ========================
RESULTS_DIR   = "/QRISdata/Q6636/sra_ww_mobilization/results"
RGI_HOST_FILE = os.path.join(RESULTS_DIR, "arg_analysis", "rgi_with_mag_taxonomy.csv")
META_FILE     = os.path.join(RESULTS_DIR, "sample_map_complete.tsv")
BRACKEN_FILE  = os.path.join(RESULTS_DIR, "bracken_genus_abundance.csv")
OUTPUT_TSV    = os.path.expanduser("~/paper_audit_results.tsv")
OUTPUT_TXT    = os.path.expanduser("~/paper_audit_summary.txt")

MUNI_INF = "ww_influent_municipal"
MUNI_EFF = "ww_effluent_municipal"

# Enterobacteriaceae 属代理（如 family 列不存在时使用）
ENTERO_GENERA = {"Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                 "Salmonella", "Serratia", "Proteus", "Morganella",
                 "Providencia", "Cronobacter", "Pantoea", "Yersinia"}

# ======================== HELPERS ========================
def gini(x):
    """Gini 系数；输入为计数或频率向量"""
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) == 0: return 0.0
    x = np.sort(x); n = len(x); idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))

def bray_curtis(p, q):
    """Bray-Curtis 距离；两向量需 index 对齐"""
    p, q = np.asarray(p, float), np.asarray(q, float)
    s = p.sum() + q.sum()
    return np.sum(np.abs(p - q)) / s if s > 0 else 0.0

def upper_bound_zero(n, alpha=0.05):
    """0 successes 时二项比例的 95% 单侧上界"""
    return 1 - alpha ** (1.0 / n) if n > 0 else 1.0

# ======================== RECORDING ========================
RECORDS = []

def record(num, section, label, paper, computed, status="?", notes=""):
    sym_map = {"OK": "[OK]    ", "DIFF": "[DIFF]  ", "?": "[CHECK] ", "!": "[WARN]  "}
    sym = sym_map.get(status, status)
    RECORDS.append({
        "num": num, "section": section, "label": label,
        "paper_claim": str(paper), "computed": str(computed),
        "status": status, "notes": notes
    })
    print(f"  [{num}] {label}")
    print(f"        paper:    {paper}")
    print(f"        computed: {computed}    {sym}")
    if notes:
        print(f"        note:     {notes}")
    print()

def section(title):
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72 + "\n")

# ======================== LOAD DATA ========================
print("=" * 72)
print("  Loading data files...")
print("=" * 72)

try:
    df_raw = pd.read_csv(RGI_HOST_FILE, low_memory=False)
    print(f"  RGI-host: {df_raw.shape[0]:,} rows x {df_raw.shape[1]} cols")
except Exception as e:
    print(f"FATAL: {e}"); sys.exit(1)

try:
    meta = pd.read_csv(META_FILE, sep="\t")
    print(f"  Metadata: {meta.shape[0]:,} samples")
except Exception as e:
    print(f"FATAL: {e}"); sys.exit(1)

try:
    bracken = pd.read_csv(BRACKEN_FILE, index_col=0)
    print(f"  Bracken:  {bracken.shape[0]:,} genera x {bracken.shape[1]:,} samples")
except Exception as e:
    print(f"WARN: Cannot read Bracken: {e}")
    bracken = None

# 检测 category 列
cat_col = None
for c in ["category", "category_final"]:
    if c in df_raw.columns:
        cat_col = c
        break
if cat_col is None and "category" in meta.columns:
    df_raw = df_raw.merge(meta[["sample_id", "category"]], on="sample_id", how="left")
    cat_col = "category"
print(f"  Category column: '{cat_col}'")

# 必需列检查
print(f"\n  All columns: {sorted(df_raw.columns.tolist())}")
required = ["sample_id", "on_plasmid", "on_virus", "genus", "Best_Hit_ARO", "Drug Class"]
missing = [c for c in required if c not in df_raw.columns]
if missing:
    print(f"  WARN: Missing columns: {missing}")

# 移动性标志
df_raw["is_mobile"] = (df_raw["on_plasmid"].fillna(False).astype(bool) |
                       df_raw["on_virus"].fillna(False).astype(bool)).astype(int)

# 去重
dedup_keys = [c for c in ["sample_id", "Contig", "Best_Hit_ARO"] if c in df_raw.columns]
df = df_raw.drop_duplicates(subset=dedup_keys) if dedup_keys else df_raw.copy()
removed = len(df_raw) - len(df)
print(f"  After dedup ({dedup_keys}): {df.shape[0]:,} rows ({removed:,} removed)")

# 子集
inf = df[df[cat_col] == MUNI_INF]
eff = df[df[cat_col] == MUNI_EFF]
muni = df[df[cat_col].isin([MUNI_INF, MUNI_EFF])]
host_inf = inf[inf["genus"].notna()]
host_eff = eff[eff["genus"].notna()]
host_muni = muni[muni["genus"].notna()]
print(f"\n  Subsets:")
print(f"    Influent:  {len(inf):>6,} rows ({len(host_inf):>6,} host-assigned)")
print(f"    Effluent:  {len(eff):>6,} rows ({len(host_eff):>6,} host-assigned)")
print(f"    Municipal: {len(muni):>6,} rows ({len(host_muni):>6,} host-assigned)")

# ============================================================
# §2.1 / §3.1  Sample counts
# ============================================================
section("Section 2.1 / 3.1  Sample counts")

cat_counts = meta["category"].value_counts() if "category" in meta.columns else pd.Series()

record("S1.1", "2.1", "Total samples", 185, len(meta),
       "OK" if len(meta) == 185 else "DIFF")

for paper_n, cat_name in [(108, MUNI_INF), (61, MUNI_EFF), (9, "ww_sludge"),
                            (3, "ww_influent_hospital"), (2, "ww_effluent_hospital")]:
    actual = int(cat_counts.get(cat_name, 0))
    record(f"S1.{cat_name[3:8]}", "2.1", f"Samples in {cat_name}", paper_n, actual,
           "OK" if actual == paper_n else "DIFF")

au_n = int(cat_counts.get("ww_AU_municipal", 0))
record("S1.AU", "2.1", "ww_AU_municipal (NOT in paper)",
       "(unmentioned)", au_n,
       "DIFF" if au_n > 0 else "OK",
       "Paper omits this category; if non-zero, must explain")

muni_total = int(cat_counts.get(MUNI_INF, 0)) + int(cat_counts.get(MUNI_EFF, 0))
record("S1.muni", "2.1/3.1", "Municipal total (inf+eff)", 169, muni_total,
       "OK" if muni_total == 169 else "DIFF")

# 大洲数
cont_col = next((c for c in ["continent", "Continent"] if c in meta.columns), None)
country_col = next((c for c in ["country_std", "country", "Country"] if c in meta.columns), None)
if cont_col:
    n_cont = meta[cont_col].dropna().nunique()
    record("S1.cont", "2.1/Fig1A", "Number of continents",
           "5 (paper §2.1) vs 4 (Fig 1A)", n_cont,
           "OK" if n_cont == 4 else "DIFF",
           "Paper §2.1 says 5 'continents'; Fig 1A says 4 'continents' (project: 4)")
if country_col:
    n_country = meta[country_col].dropna().nunique()
    record("S1.cnt", "2.1/3.1", "Number of countries", 29, n_country,
           "OK" if n_country == 29 else "DIFF")

# ============================================================
# §2.5 / §3.1  ARG hit counts
# ============================================================
section("Section 2.5 / 3.1  ARG hit counts")

record("S2.1", "2.5.1", "Total ARG hits (185 samples, raw)", 38555, len(df_raw),
       "OK" if len(df_raw) == 38555 else "?")
record("S2.1b", "2.5.1", "Total ARG hits (after dedup)", "(implied 38,555)", len(df), "?",
       f"Removed {removed:,} duplicates")

if "Drug Class" in df.columns:
    n_dc_raw = df["Drug Class"].dropna().nunique()
    record("S2.2", "2.5.1 / 3.1", "Distinct drug class strings",
           "28 (§2.5.1) vs 66 (§3.1)", n_dc_raw,
           "DIFF" if n_dc_raw not in (28, 66) else "?",
           "Paper §2.5.1 and §3.1 contradict each other")

record("S2.3", "3.1", "ARG hits in 169 municipal", 13332, len(muni),
       "OK" if len(muni) == 13332 else "DIFF",
       f"inf={len(inf):,}, eff={len(eff):,}")

# 每样本统计
sample_counts = muni.groupby("sample_id").size()
all_muni = meta[meta["category"].isin([MUNI_INF, MUNI_EFF])]["sample_id"].tolist() \
           if "category" in meta.columns else []
zero_samples = [s for s in all_muni if s not in sample_counts.index]
if zero_samples:
    sample_counts = pd.concat([
        sample_counts,
        pd.Series([0] * len(zero_samples), index=zero_samples)
    ])
if len(sample_counts) > 0:
    record("S2.4", "3.1", "Mean ARG hits per municipal sample", 82,
           round(sample_counts.mean(), 1), "?")
    record("S2.5", "3.1", "Median ARG hits", 72, int(sample_counts.median()), "?")
    record("S2.6", "3.1", "Range (min-max)", "1-383",
           f"{int(sample_counts.min())}-{int(sample_counts.max())}", "?")

# ============================================================
# §2.6 / §3.1  Mobility
# ============================================================
section("Section 2.6 / 3.1  Mobility")

n_muni_pl = int(muni["on_plasmid"].fillna(False).astype(bool).sum())
n_muni_ph = int(muni["on_virus"].fillna(False).astype(bool).sum())
n_muni_mob = int(muni["is_mobile"].sum())
mr_muni = 100.0 * n_muni_mob / len(muni) if len(muni) > 0 else 0

record("S3.1", "2.6.1", "Plasmid hits in 169 municipal", 1283, n_muni_pl,
       "OK" if n_muni_pl == 1283 else "DIFF")
record("S3.2", "2.6.1", "Phage hits in 169 municipal", 126, n_muni_ph,
       "OK" if n_muni_ph == 126 else "DIFF")
record("S3.3", "2.6.1", "Total mobile (1409=1283+126)", 1409, n_muni_mob,
       "OK" if n_muni_mob == 1409 else "DIFF")
record("S3.4", "2.6.1", "Municipal mobility rate %", 10.6, round(mr_muni, 2),
       "OK" if abs(mr_muni - 10.6) < 0.5 else "DIFF")

n_all_pl = int(df["on_plasmid"].fillna(False).astype(bool).sum())
n_all_ph = int(df["on_virus"].fillna(False).astype(bool).sum())
n_all_mob = int(df["is_mobile"].sum())
mr_all = 100.0 * n_all_mob / len(df)
record("S3.5", "2.6.1", "Plasmid hits across 185", 3926, n_all_pl,
       "OK" if n_all_pl == 3926 else "DIFF")
record("S3.6", "2.6.1", "Phage hits across 185", 392, n_all_ph,
       "OK" if n_all_ph == 392 else "DIFF")
record("S3.7", "2.6.1", "Global mobility rate %", 11.2, round(mr_all, 2),
       "OK" if abs(mr_all - 11.2) < 0.5 else "DIFF")

# ============================================================
# §2.8  Host assignment
# ============================================================
section("Section 2.8  Host assignment")

record("S4.1", "2.8", "Host-assigned ARG in 169 municipal", 12446, len(host_muni),
       "OK" if len(host_muni) == 12446 else "DIFF")
record("S4.1b", "2.8", "Host assignment %", 93.4,
       round(100.0 * len(host_muni) / len(muni), 1), "?")

# ============================================================
# §3.2  ARG-bearing host distribution
# ============================================================
section("Section 3.2  ARG-bearing host distribution")

record("S5.1", "3.2", "Influent host-assigned (paper 9,387)", 9387, len(host_inf),
       "OK" if len(host_inf) == 9387 else "DIFF")
record("S5.2", "3.2", "Effluent host-assigned (paper 2,983)", 2983, len(host_eff),
       "OK" if len(host_eff) == 2983 else "DIFF")
sum_inf_eff = len(host_inf) + len(host_eff)
record("S5.3", "3.2 vs 2.8", "Sum check (12,446 vs 12,370)",
       "9387+2983=12,370 vs 12,446 in §2.8",
       f"sum={sum_inf_eff}, total in §2.8={len(host_muni)}, diff={len(host_muni)-sum_inf_eff}",
       "OK" if sum_inf_eff == len(host_muni) else "DIFF")

inf_genera = host_inf["genus"].value_counts()
eff_genera = host_eff["genus"].value_counts()
record("S5.4", "3.2", "Influent: distinct genera (paper 1,060)", 1060, len(inf_genera),
       "OK" if len(inf_genera) == 1060 else "?")
record("S5.5", "3.2", "Effluent: distinct genera (paper 807)", 807, len(eff_genera),
       "OK" if len(eff_genera) == 807 else "?")

print("Top 5 influent ARG-bearing genera:")
for g, n in inf_genera.head(5).items():
    print(f"    {g:25s}: {n} ({100*n/len(host_inf):.2f}%)")
print()
print("Top 5 effluent ARG-bearing genera:")
for g, n in eff_genera.head(5).items():
    print(f"    {g:25s}: {n} ({100*n/len(host_eff):.2f}%)")
print()

esch_l2 = 100.0 * inf_genera.get("Escherichia", 0) / len(host_inf)
record("S5.6", "3.2/3.3", "Escherichia in ARG-bearing influent (paper 6.0%)", 6.0,
       f"{esch_l2:.2f}%", "OK" if abs(esch_l2 - 6.0) < 0.3 else "?")

g_inf = gini(inf_genera.values)
g_eff = gini(eff_genera.values)
record("S5.7", "3.2/Fig5B", "Influent ARG-host Gini (paper 0.681)", 0.681, round(g_inf, 3),
       "OK" if abs(g_inf - 0.681) < 0.01 else "?")
record("S5.8", "3.2/Fig5B", "Effluent ARG-host Gini (paper 0.555)", 0.555, round(g_eff, 3),
       "OK" if abs(g_eff - 0.555) < 0.01 else "?")

# Bray-Curtis ARG-host
common = sorted(set(inf_genera.index) | set(eff_genera.index))
inf_v = np.array([inf_genera.get(g, 0) for g in common])
eff_v = np.array([eff_genera.get(g, 0) for g in common])
bc_arg = bray_curtis(inf_v, eff_v)
record("S5.9", "3.2/3.5", "Bray-Curtis ARG-host (paper 0.331)", 0.331, round(bc_arg, 3),
       "OK" if abs(bc_arg - 0.331) < 0.02 else "?")

# Enterobacteriaceae 检测
fam_col = next((c for c in ["family", "Family"] if c in host_inf.columns), None)
if fam_col:
    n_entero = (host_inf[fam_col] == "Enterobacteriaceae").sum()
    method = "family-level"
else:
    n_entero = host_inf["genus"].isin(ENTERO_GENERA).sum()
    method = f"genus proxy (n={len(ENTERO_GENERA)} genera)"
pct_entero = 100.0 * n_entero / len(host_inf)
record("S5.10", "3.2/3.3", "Enterobacteriaceae % in inf ARG-bearing (paper 8.49%)",
       8.49, f"{pct_entero:.2f}%",
       "OK" if abs(pct_entero - 8.49) < 0.5 else "?",
       f"Detection: {method}")

# ============================================================
# §3.3  Mobilome concentration
# ============================================================
section("Section 3.3  Mobilome concentration")

mob_inf = host_inf[host_inf["is_mobile"] == 1]
mob_eff = host_eff[host_eff["is_mobile"] == 1]
record("S6.1", "3.3", "Mobile (host-assigned) in influent", "(implied)", len(mob_inf), "?")
record("S6.2", "3.3", "Mobile (host-assigned) in effluent", "(paper says 0/301)", len(mob_eff),
       "?", "Paper denominator is 301; check exact")

total_mob_municipal = len(mob_inf) + len(mob_eff)
record("S6.3", "3.3-3.6 KEY", "TOTAL host-assigned mobile in 169 municipal (paper 1,514)",
       1514, total_mob_municipal,
       "OK" if total_mob_municipal == 1514 else "DIFF",
       "If <1514, paper's '1,514' must include non-municipal samples")

host_all = df[df["genus"].notna()]
n_all_host_mob = int(host_all["is_mobile"].sum())
record("S6.4", "3.3-3.6 KEY", "All-185 host-assigned mobile",
       "(was '1,514' from 169 muni in paper)", n_all_host_mob,
       "?", "If 1,514 was actually all-185 host-assigned, this should match")

# Mobilome Gini
mob_inf_genera = mob_inf["genus"].value_counts()
mob_eff_genera = mob_eff["genus"].value_counts()
g_mob_inf = gini(mob_inf_genera.values)
g_mob_eff = gini(mob_eff_genera.values)
record("S6.5", "3.3/Fig3A", "Influent mobilome Gini (paper 0.538)", 0.538,
       round(g_mob_inf, 3),
       "OK" if abs(g_mob_inf - 0.538) < 0.01 else "?")
record("S6.6", "3.3/Fig3A", "Effluent mobilome Gini (paper 0.345)", 0.345,
       round(g_mob_eff, 3),
       "OK" if abs(g_mob_eff - 0.345) < 0.01 else "?")

# Bray-Curtis mobilome
common_m = sorted(set(mob_inf_genera.index) | set(mob_eff_genera.index))
mi_v = np.array([mob_inf_genera.get(g, 0) for g in common_m])
me_v = np.array([mob_eff_genera.get(g, 0) for g in common_m])
bc_mob = bray_curtis(mi_v, me_v)
record("S6.7", "3.5/Fig5A", "Bray-Curtis mobilome (paper 0.681)", 0.681,
       round(bc_mob, 3),
       "OK" if abs(bc_mob - 0.681) < 0.02 else "?")

# 排列零模型
print("Permutation null for influent mobilome Gini (1000 iter)...")
rng = np.random.RandomState(42)
gini_null = []
n_mob_inf = len(mob_inf)
inf_pool = host_inf["genus"].dropna().values
for _ in range(1000):
    if len(inf_pool) < n_mob_inf:
        idx = rng.choice(len(inf_pool), size=n_mob_inf, replace=True)
    else:
        idx = rng.choice(len(inf_pool), size=n_mob_inf, replace=False)
    samp_counts = pd.Series(inf_pool[idx]).value_counts()
    gini_null.append(gini(samp_counts.values))
gini_null = np.array(gini_null)
null_mean = gini_null.mean()
null_std = gini_null.std()
z = (g_mob_inf - null_mean) / null_std if null_std > 0 else 0
record("S6.8", "3.3", "Random expected Gini (paper 0.270)", 0.270, round(null_mean, 3),
       "OK" if abs(null_mean - 0.270) < 0.05 else "?",
       f"std={null_std:.4f}; permutation: random sample without replacement from ARG-bearing pool")
record("S6.9", "3.3", "Z-score (paper 30.7)", 30.7, round(z, 1),
       "OK" if abs(z - 30.7) < 5 else "?",
       "Different permutation schemes give slightly different Z; this is approximate")

# ============================================================
# §3.3  Three-layer enrichment
# ============================================================
section("Section 3.3  Three-layer enrichment (Bracken)")

if bracken is not None:
    inf_samples = meta[meta["category"] == MUNI_INF]["sample_id"].tolist()
    bracken_samples = [s for s in inf_samples if s in bracken.columns]
    print(f"Bracken inf samples available: {len(bracken_samples)}/{len(inf_samples)}")
    if bracken_samples:
        bracken_inf = bracken[bracken_samples]
        bracken_inf_rel = bracken_inf.div(bracken_inf.sum(axis=0), axis=1)
        bracken_inf_mean = bracken_inf_rel.mean(axis=1)
        
        print("\nTop 10 community genera (Bracken mean across inf):")
        for g, v in bracken_inf_mean.sort_values(ascending=False).head(10).items():
            print(f"    {str(g)[:25]:25s}: {100*v:.2f}%")
        print()
        
        # Escherichia 三层
        esch_l1 = bracken_inf_mean.get("Escherichia", np.nan)
        if pd.isna(esch_l1):
            for v in bracken_inf_mean.index:
                if "Escherichia" in str(v):
                    esch_l1 = bracken_inf_mean[v]
                    print(f"    [INFO] Escherichia matched as '{v}'")
                    break
        esch_l1_pct = 100 * esch_l1 if not pd.isna(esch_l1) else np.nan
        
        record("S7.1", "3.3", "Escherichia L1 community % (paper 1.15%)", 1.15,
               f"{esch_l1_pct:.2f}%" if not pd.isna(esch_l1_pct) else "NA",
               "OK" if not pd.isna(esch_l1_pct) and abs(esch_l1_pct - 1.15) < 0.2 else "?")
        
        esch_l3 = 100.0 * mob_inf_genera.get("Escherichia", 0) / len(mob_inf) if len(mob_inf) > 0 else 0
        record("S7.2", "3.3", "Escherichia L3 mobile % (paper 15.57%)", 15.57,
               f"{esch_l3:.2f}%", "OK" if abs(esch_l3 - 15.57) < 0.5 else "?")
        
        if not pd.isna(esch_l1_pct) and esch_l1_pct > 0:
            fold = esch_l3 / esch_l1_pct
            record("S7.3", "3.3", "Escherichia L1->L3 fold (paper 13.57)", 13.57,
                   round(fold, 2),
                   "OK" if abs(fold - 13.57) < 0.5 else "DIFF",
                   f"15.57/1.15 = 13.54 (direct); paper writes 13.57 (rounded compound)")
        
        # Enterobacteriaceae 三层
        entero_in_bracken = [g for g in ENTERO_GENERA if g in bracken_inf_mean.index]
        entero_l1 = sum(bracken_inf_mean.get(g, 0) for g in entero_in_bracken)
        entero_l1_pct = 100 * entero_l1
        record("S7.4", "3.3", "Enterobacteriaceae L1 community % (paper 4.01%)", 4.01,
               f"{entero_l1_pct:.2f}%",
               "OK" if abs(entero_l1_pct - 4.01) < 1 else "?",
               f"Genus-proxy detected: {entero_in_bracken[:5]}")
        
        if len(mob_inf) > 0:
            entero_l3 = 100.0 * mob_inf["genus"].isin(ENTERO_GENERA).sum() / len(mob_inf)
            record("S7.5", "3.3", "Enterobacteriaceae L3 mobile % (paper 20.4%)", 20.4,
                   f"{entero_l3:.2f}%",
                   "OK" if abs(entero_l3 - 20.4) < 1 else "?")
        
        # Acinetobacter (decreased)
        acin_l1 = 100 * bracken_inf_mean.get("Acinetobacter", 0)
        acin_l3 = 100.0 * mob_inf_genera.get("Acinetobacter", 0) / len(mob_inf) if len(mob_inf) > 0 else 0
        record("S7.6", "3.3", "Acinetobacter L1 (paper 6.19%)", 6.19, f"{acin_l1:.2f}%", "?")
        record("S7.7", "3.3", "Acinetobacter L3 (paper 1.08%)", 1.08, f"{acin_l3:.2f}%", "?")
else:
    print("Bracken not loaded; skipping three-layer enrichment")

# Effluent Enterobacteriaceae - critical zero claim
n_entero_mob_eff = mob_eff["genus"].isin(ENTERO_GENERA).sum() if not fam_col else \
                   (mob_eff[fam_col] == "Enterobacteriaceae").sum()
samples_with_entero_eff = host_eff[host_eff["genus"].isin(ENTERO_GENERA)]["sample_id"].nunique() if not fam_col else \
                           host_eff[host_eff[fam_col] == "Enterobacteriaceae"]["sample_id"].nunique()

record("S7.8", "3.3/3.5 KEY", "Enterobacteriaceae mobile in effluent (paper 0/301)",
       "0/301", f"{n_entero_mob_eff}/{len(mob_eff)}",
       "OK" if n_entero_mob_eff == 0 else "DIFF")
record("S7.9", "3.3/3.5", "Effluent samples with Enterobacteriaceae ARG (paper 29/61)",
       "29/61", f"{samples_with_entero_eff}/{host_eff['sample_id'].nunique()}",
       "OK" if samples_with_entero_eff == 29 else "?")

if len(mob_eff) > 0:
    ub = upper_bound_zero(len(mob_eff))
    record("S7.10", "3.3/3.5", "95% upper bound (paper <0.99%)", "<0.99%",
           f"{ub*100:.2f}%", "OK" if abs(ub - 0.0099) < 0.005 else "?")

# ============================================================
# §3.6  Mechanism mobility rates
# ============================================================
section("Section 3.6  Resistance mechanism mobility")

if "Resistance Mechanism" in host_muni.columns:
    print("Raw mechanism distribution (top 10, host-assigned in 169 municipal):")
    mech_raw = host_muni["Resistance Mechanism"].value_counts()
    for m, n in mech_raw.head(10).items():
        sub = host_muni[host_muni["Resistance Mechanism"] == m]
        mob = 100.0 * sub["is_mobile"].sum() / len(sub)
        print(f"    {str(m)[:50]:50s}: n={n:>6,}, mob={mob:>5.2f}%")
    print()
    
    def map_mech(s):
        if pd.isna(s): return "unknown"
        s = str(s).lower()
        if "replacement" in s: return "replacement"
        if "alteration" in s and "efflux" in s: return "alteration_efflux_mixed"
        if "alteration" in s: return "alteration"
        if "inactivation" in s: return "inactivation"
        if "protection" in s: return "protection"
        if "efflux" in s: return "efflux"
        return "other"
    
    host_muni2 = host_muni.copy()
    host_muni2["mech_g"] = host_muni2["Resistance Mechanism"].apply(map_mech)
    
    paper_mob = {"replacement": 30.2, "inactivation": 23.2, "protection": 19.7,
                 "efflux": 15.4, "alteration": 3.2}
    paper_mbi = {"replacement": 1.43, "alteration": -1.81}
    
    print("Aggregated to paper's 5 categories (host-assigned, in 169 municipal):")
    for mg in ["replacement", "inactivation", "protection", "efflux", "alteration",
               "alteration_efflux_mixed"]:
        sub = host_muni2[host_muni2["mech_g"] == mg]
        if len(sub) == 0: continue
        n = len(sub); mb = sub["is_mobile"].sum(); pct = 100.0*mb/n
        print(f"    {mg:30s}: n={n:>6,}, mobile={mb:>5,}, rate={pct:.2f}%")
        if mg in paper_mob:
            paper = paper_mob[mg]
            record(f"S8.{mg}", "3.6", f"Mobility — {mg} (paper {paper}%)", paper,
                   f"{pct:.2f}% (n={n})",
                   "OK" if abs(pct - paper) < 1.5 else "DIFF")
            if mg in paper_mbi:
                m_val = np.log2(pct / 11.2) if pct > 0 else np.nan
                record(f"S8.MBI.{mg}", "3.6", f"MBI — {mg} (paper {paper_mbi[mg]})",
                       paper_mbi[mg], round(m_val, 2),
                       "OK" if abs(m_val - paper_mbi[mg]) < 0.3 else "?")

# vanR/vanS check
if "Best_Hit_ARO" in host_muni.columns:
    vrs = host_muni[host_muni["Best_Hit_ARO"].astype(str).str.match(r"^van[RS]\b", case=False, na=False)]
    vrs_pl = vrs[vrs["on_plasmid"].fillna(False).astype(bool)]
    record("S8.van", "3.6", "vanR/vanS on plasmid (paper '2 of 14,625')",
           "2 of 14,625", f"{len(vrs_pl)} of {len(host_muni)}", "?",
           f"vanR/vanS total in host_muni: {len(vrs)}")

# ============================================================
# §3.1 / §3.6  Drug class composition
# ============================================================
section("Section 3.1 / 3.6  Drug class composition")

if "Drug Class" in muni.columns:
    keywords = ["sulfonamide", "phenicol", "macrolide", "aminoglycoside",
                "tetracycline", "fluoroquinolone", "glycopeptide", "beta-lactam"]
    
    def classify_dc(s):
        if pd.isna(s): return "other"
        s = str(s).lower()
        for k in keywords:
            if k in s: return k
            if k == "beta-lactam" and any(x in s for x in
                ["cephalosporin", "carbapenem", "penam", "monobactam", "penicillin"]):
                return "beta-lactam"
        return "other"
    
    muni_dc = muni.copy()
    muni_dc["dc"] = muni_dc["Drug Class"].apply(classify_dc)
    dist = muni_dc["dc"].value_counts(normalize=True) * 100
    
    paper_dc = {"glycopeptide": 45.8, "tetracycline": 14.5, "beta-lactam": 5.0,
                "aminoglycoside": 4.7, "macrolide": 4.6, "phenicol": 2.1,
                "sulfonamide": 0.8, "other": 20.9}
    
    print("Drug class composition (priority keyword matching, in 169 municipal):")
    for k in keywords + ["other"]:
        obs = dist.get(k, 0)
        paper = paper_dc.get(k, "(not in paper)")
        n_class = (muni_dc["dc"] == k).sum()
        print(f"    {k:18s}: paper={str(paper):>10s}, observed={obs:>6.2f}% (n={n_class:>6,})")
        if isinstance(paper, (int, float)):
            record(f"S9.{k}", "3.1", f"Drug class — {k}", paper,
                   f"{obs:.2f}%",
                   "OK" if abs(obs - paper) < 1 else "DIFF")
        else:
            record(f"S9.{k}", "3.1 MISSING", f"Drug class — {k} (NOT IN PAPER §3.1)",
                   "(unmentioned)", f"{obs:.2f}%",
                   "DIFF" if obs > 5 else "?",
                   "If substantial, paper §3.1 should add this class")

# Drug class fold range
fold_dc = 2 ** (1.88 - (-2.47))
record("S10.fold_dc", "3.6/Fig6C", "Drug class MBI range (paper 17-fold)", 17,
       round(fold_dc, 2),
       "DIFF" if abs(fold_dc - 17) > 1 else "OK",
       "From sulfonamide MBI=+1.88 to glycopeptide MBI=-2.47, actual is 20.4-fold")

# Mechanism fold range
fold_mech = 2 ** (1.43 - (-1.81))
record("S10.fold_mech", "3.6", "Mechanism MBI range (paper 9.5-fold)", 9.5,
       round(fold_mech, 2),
       "OK" if abs(fold_mech - 9.5) < 0.2 else "DIFF",
       "From replacement MBI=+1.43 to alteration MBI=-1.81")

# Escherichia enrichment fold (paper 13.57)
fold_esch = 15.57 / 1.15
record("S10.fold_esch", "3.3/Conclusions",
       "Escherichia enrichment fold (paper 13.57)", 13.57,
       round(fold_esch, 2),
       "DIFF" if abs(fold_esch - 13.57) > 0.05 else "OK",
       "15.57/1.15 = 13.539; paper rounds via 5.24x x 2.59x = 13.57")

# Enterobacteriaceae enrichment fold
fold_entero = 20.39 / 4.01
record("S10.fold_entero", "3.3", 
       "Enterobacteriaceae enrichment fold (paper 5.09)", 5.09,
       round(fold_entero, 3),
       "OK" if abs(fold_entero - 5.09) < 0.05 else "DIFF",
       "20.39/4.01 = 5.085")

# ============================================================
# 输出
# ============================================================
print("\n" + "=" * 72)
print("  Writing outputs")
print("=" * 72)

records_df = pd.DataFrame(RECORDS)
records_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print(f"\n  Detail: {OUTPUT_TSV}  ({len(records_df)} checks)")

with open(OUTPUT_TXT, "w") as fh:
    fh.write("PAPER AUDIT SUMMARY\n")
    fh.write("=" * 70 + "\n\n")
    by_status = records_df["status"].value_counts()
    fh.write(f"Total checks: {len(records_df)}\n")
    for s, label in [("OK", "matches"), ("DIFF", "mismatches"),
                      ("?", "unverified"), ("!", "warnings")]:
        fh.write(f"  {s:5s} ({label}): {by_status.get(s, 0)}\n")
    fh.write("\n")
    
    mismatches = records_df[records_df["status"] == "DIFF"]
    fh.write(f"\n--- MISMATCHES ({len(mismatches)}) — must fix ---\n\n")
    for _, r in mismatches.iterrows():
        fh.write(f"[{r['num']}] {r['section']} - {r['label']}\n")
        fh.write(f"   PAPER:    {r['paper_claim']}\n")
        fh.write(f"   COMPUTED: {r['computed']}\n")
        if r["notes"]:
            fh.write(f"   NOTE:     {r['notes']}\n")
        fh.write("\n")
    
    unverified = records_df[records_df["status"] == "?"]
    fh.write(f"\n--- UNVERIFIED ({len(unverified)}) — needs visual check ---\n\n")
    for _, r in unverified.iterrows():
        fh.write(f"[{r['num']}] {r['label']}: paper={r['paper_claim']}, computed={r['computed']}\n")

print(f"  Summary: {OUTPUT_TXT}")
print("\n" + "=" * 72)
print(f"  DONE")
print(f"  Mismatches: {by_status.get('DIFF', 0)}")
print(f"  Total checks: {len(records_df)}")
print("=" * 72)

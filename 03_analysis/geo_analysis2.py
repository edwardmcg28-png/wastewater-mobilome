import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kruskal, spearmanr
import warnings
warnings.filterwarnings('ignore')

# ── 载入所有数据源 ──────────────────────────────────────────
metrics = pd.read_csv("/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/per_sample_viz_metrics.csv")
bracken = pd.read_csv("/QRISdata/Q6636/sra_ww_mobilization/results/bracken_genus_abundance.csv", index_col=0)

# ARGs-OAP abundance (copies per 16S)
try:
    argsoap = pd.read_csv("/QRISdata/Q6636/sra_ww_mobilization/results/args_oap/ARG_abundance.txt", sep='\t', index_col=0)
    print(f"ARGs-OAP shape: {argsoap.shape}")
    argsoap_total = argsoap.sum(axis=0).reset_index()
    argsoap_total.columns = ['sample_id','ARG_copies_per16S']
except Exception as e:
    print(f"ARGs-OAP load error: {e}")
    argsoap_total = None

# Alpha diversity (Shannon) from bracken
def shannon(v):
    v = v[v>0]
    p = v / v.sum()
    return -np.sum(p * np.log(p))

shannon_div = bracken.apply(shannon, axis=0).reset_index()
shannon_div.columns = ['sample_id','shannon_diversity']

# Total bacterial reads proxy: sum of bracken abundances
bracken_total = bracken.sum(axis=0).reset_index()
bracken_total.columns = ['sample_id','total_bracken_reads']

# Richness (number of genera detected >0)
richness = (bracken > 0).sum(axis=0).reset_index()
richness.columns = ['sample_id','genus_richness']

# Merge all
df = metrics[metrics["income_group"].notna()].copy()
df = df.merge(shannon_div, on='sample_id', how='left')
df = df.merge(bracken_total, on='sample_id', how='left')
df = df.merge(richness, on='sample_id', how='left')
if argsoap_total is not None:
    df = df.merge(argsoap_total, on='sample_id', how='left')

# Add AU sample
au = {"sample_id":"24_In_E250051725_L01_561","income_group":"High",
      "income_order":1,"antibiotic_use_DDD":21.2,"koppen":"Cfa",
      "koppen_group":"Temperate","mobility_rate":10/179,
      "entero_community_pct":0.0244}
au = {k:v for k,v in au.items() if k in df.columns}
df = pd.concat([df, pd.DataFrame([au])], ignore_index=True)

print(f"\n总样本数: {len(df)}")
print(f"可用指标: {[c for c in df.columns if c not in ['sample_id','income_group','koppen','koppen_group','country','city']]}")

METRICS = {}
for col in ['mobility_rate','n_ARG_total','ARG_copies_per16S',
            'shannon_diversity','genus_richness','total_bracken_reads',
            'entero_community_pct','gini_mobilome','gini_arg_host']:
    if col in df.columns:
        METRICS[col] = col

print(f"\n纳入分析的指标: {list(METRICS.keys())}\n")

def run_spearman(ctx_col, ctx_name):
    print(f"\n  {ctx_name}")
    for mc in METRICS:
        sub = df[[ctx_col,mc]].dropna()
        if len(sub)<10: continue
        r,p = spearmanr(sub[ctx_col],sub[mc])
        sig = "**" if p<0.01 else ("*" if p<0.05 else "ns")
        print(f"    {mc:<28} rho={r:+.3f}  p={p:.4f}  {sig}  n={len(sub)}")

def run_kruskal(grp_col, grp_name):
    if grp_col not in df.columns: return
    print(f"\n  {grp_name}")
    for mc in METRICS:
        sub = df[[grp_col,mc]].dropna()
        groups = [g[mc].values for _,g in sub.groupby(grp_col) if len(g)>=3]
        if len(groups)<2: continue
        h,p = kruskal(*groups)
        sig = "**" if p<0.01 else ("*" if p<0.05 else "ns")
        # group medians
        meds = {k: round(np.median(g[mc].dropna()),3) for k,g in sub.groupby(grp_col) if len(g[mc].dropna())>=1}
        print(f"    {mc:<28} H={h:.2f}  p={p:.4f}  {sig}  medians={meds}")

print("="*65)
print("1. SPEARMAN CORRELATIONS")
print("="*65)
run_spearman("antibiotic_use_DDD", "Antibiotic use DDD")
run_spearman("income_order",       "Income order (1=High → 4=Low)")

print("\n"+"="*65)
print("2. KRUSKAL-WALLIS BY GROUP")
print("="*65)
run_kruskal("income_group",  "Income group")
run_kruskal("koppen_group",  "Köppen group")

print("\n"+"="*65)
print("3. PAIRWISE: HIGH INCOME vs LOW+LOWER-MIDDLE")
print("="*65)
df["income_bin"] = df["income_group"].map(
    lambda x: "High" if x=="High" else ("Low/LMI" if x in ("Low","Lower-middle") else "UMI") if pd.notna(x) else np.nan)
for mc in METRICS:
    hi = df[df["income_bin"]=="High"][mc].dropna()
    lo = df[df["income_bin"]=="Low/LMI"][mc].dropna()
    if len(hi)<5 or len(lo)<3: continue
    _,p = stats.mannwhitneyu(hi,lo,alternative='two-sided')
    sig = "**" if p<0.01 else ("*" if p<0.05 else "ns")
    print(f"  {mc:<28} High(n={len(hi)}) med={np.median(hi):.3f}  Low/LMI(n={len(lo)}) med={np.median(lo):.3f}  p={p:.4f} {sig}")


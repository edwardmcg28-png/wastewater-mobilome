import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kruskal, spearmanr
import ast
import warnings
warnings.filterwarnings('ignore')

metrics = pd.read_csv("/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/per_sample_viz_metrics.csv")
top3    = pd.read_csv("/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/per_sample_top3.csv")
df = metrics[metrics["income_group"].notna()].merge(top3, on="sample_id")

au_row = {"sample_id":"24_In_E250051725_L01_561","income_group":"High",
          "income_order":1,"antibiotic_use_DDD":21.2,"koppen":"Cfa",
          "koppen_group":"Temperate","mobility_rate":10/179,
          "entero_mobile_pct":0.5,"entero_community_pct":0.0244,
          "n_ARG_total":179,"n_mobile":10,
          "layer3_mobile":"{'Klebsiella':0.5,'Janthinobacterium':0.1,'Others':0.4}",
          "layer1_dc":"{'Glycopeptide':0.5754,'Rifamycin':0.1006,'Disinfectant':0.0726,'Others':0.2514}",
          "layer2_host":"{'Klebsiella':0.0542,'Others':0.9458}"}
au_row = {k:v for k,v in au_row.items() if k in df.columns}
df = pd.concat([df, pd.DataFrame([au_row])], ignore_index=True)

ENTERO = {'Escherichia','Klebsiella','Citrobacter','Enterobacter','Serratia','Rahnella','Kluyvera'}

def entero_L3(s):
    try: return sum(v for k,v in ast.literal_eval(str(s)).items() if k in ENTERO)
    except: return np.nan

def glyco_L1(s):
    try: return ast.literal_eval(str(s)).get("Glycopeptide",0)
    except: return np.nan

df["entero_L3"] = df["layer3_mobile"].apply(entero_L3)
df["glyco_L1"]  = df["layer1_dc"].apply(glyco_L1)

METRICS = {
    "Mobility rate":       "mobility_rate",
    "Entero L3 fraction":  "entero_L3",
    "Glycopeptide L1":     "glyco_L1",
    "Entero community":    "entero_community_pct",
}

print(f"n={len(df)}\n")

print("="*60)
print("1. SPEARMAN: continuous context variables")
print("="*60)
for ctx_name, ctx_col in [("Antibiotic use DDD","antibiotic_use_DDD"),("Income order 1=High","income_order")]:
    print(f"\n  {ctx_name}")
    for mn,mc in METRICS.items():
        sub = df[[ctx_col,mc]].dropna()
        if len(sub)<10: continue
        r,p = spearmanr(sub[ctx_col],sub[mc])
        print(f"    {mn:<25} rho={r:+.3f}  p={p:.4f} {'**' if p<0.01 else '*' if p<0.05 else ''}")

print("\n"+"="*60)
print("2. KRUSKAL-WALLIS by group")
print("="*60)
for grp_name,grp_col in [("Income group","income_group"),("Köppen group","koppen_group"),("Köppen code","koppen")]:
    if grp_col not in df.columns: continue
    print(f"\n  {grp_name}")
    for mn,mc in METRICS.items():
        sub = df[[grp_col,mc]].dropna()
        groups = [g[mc].values for _,g in sub.groupby(grp_col) if len(g)>=3]
        if len(groups)<2: continue
        h,p = kruskal(*groups)
        print(f"    {mn:<25} H={h:.2f}  p={p:.4f} {'**' if p<0.01 else '*' if p<0.05 else ''}")

print("\n"+"="*60)
print("3. DESCRIPTIVE by income group")
print("="*60)
ilab={1:"High",2:"Upper-mid",3:"Lower-mid",4:"Low"}
for mn,mc in METRICS.items():
    print(f"\n  {mn}")
    for io in sorted(df["income_order"].dropna().unique()):
        g=df[df["income_order"]==io][mc].dropna()
        if len(g)==0: continue
        print(f"    {ilab.get(int(io),io):<12} n={len(g):3d}  med={np.median(g):.3f}  IQR=[{np.percentile(g,25):.3f},{np.percentile(g,75):.3f}]")

print("\n"+"="*60)
print("4. DESCRIPTIVE by Köppen group")
print("="*60)
for mn,mc in METRICS.items():
    print(f"\n  {mn}")
    if "koppen_group" not in df.columns: continue
    for kg,gdf in df.groupby("koppen_group"):
        g=gdf[mc].dropna()
        if len(g)<3: continue
        print(f"    {kg:<14} n={len(g):3d}  med={np.median(g):.3f}  IQR=[{np.percentile(g,25):.3f},{np.percentile(g,75):.3f}]")

print("\n"+"="*60)
print("5. HIGH vs LOW antibiotic use (split at median)")
print("="*60)
if "antibiotic_use_DDD" in df.columns:
    med=df["antibiotic_use_DDD"].median()
    print(f"  Median DDD={med:.1f}")
    df["ddd_hi"]=df["antibiotic_use_DDD"]>=med
    for mn,mc in METRICS.items():
        hi=df[df["ddd_hi"]==True][mc].dropna()
        lo=df[df["ddd_hi"]==False][mc].dropna()
        if len(hi)<5 or len(lo)<5: continue
        _,p=stats.mannwhitneyu(hi,lo,alternative='two-sided')
        print(f"  {mn:<25} Hi med={np.median(hi):.3f}  Lo med={np.median(lo):.3f}  p={p:.4f} {'**' if p<0.01 else '*' if p<0.05 else ''}")

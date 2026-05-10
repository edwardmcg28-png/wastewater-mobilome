#!/usr/bin/env python3
"""
audit_followup.py - 跟进核对脚本
================================================================
专门解决第一轮 audit 中标记为 [CHECK] 的关键项目：

1. Bray-Curtis 用相对丰度 (relative abundance) 重新计算
   - 论文 0.331 (ARG-host) vs 0.681 (mobilome) 是核心结论之一
   
2. Permutation null 用样本内 shuffle 重新做
   - 论文 random Gini = 0.270, Z = 30.7
   - 第一轮用了简单的 pool 抽样得到 0.477
   
3. 检查 29/61 vs 35/61 Enterobacteriaceae samples 的差异
   - 是否论文用了"≥1 mobile" 或其他更严格的过滤
   
4. 1,315 / 1,414 / 1,514 三组数字的最终核实
   - 检查所有可能的 subset 组合

5. vanR/vanS 0 vs 2 的差异
   - 是否论文用了 ARO 名称变体

运行:
    python3 audit_followup.py 2>&1 | tee ~/audit_followup_log.txt
================================================================
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np

RESULTS_DIR   = "/QRISdata/Q6636/sra_ww_mobilization/results"
RGI_HOST_FILE = os.path.join(RESULTS_DIR, "arg_analysis", "rgi_with_mag_taxonomy.csv")
META_FILE     = os.path.join(RESULTS_DIR, "sample_map_complete.tsv")

MUNI_INF = "ww_influent_municipal"
MUNI_EFF = "ww_effluent_municipal"

ENTERO_GENERA = {"Escherichia", "Klebsiella", "Enterobacter", "Citrobacter",
                 "Salmonella", "Serratia", "Proteus", "Morganella",
                 "Providencia", "Cronobacter", "Pantoea", "Yersinia"}

def gini(x):
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) == 0: return 0.0
    x = np.sort(x); n = len(x); idx = np.arange(1, n+1)
    return (2*np.sum(idx*x) - (n+1)*np.sum(x)) / (n*np.sum(x))

def bc(p, q):
    p, q = np.asarray(p, float), np.asarray(q, float)
    s = p.sum() + q.sum()
    return np.sum(np.abs(p-q))/s if s > 0 else 0.0

def section(t):
    print("\n" + "="*72 + "\n  " + t + "\n" + "="*72 + "\n")

# ===== LOAD =====
print("Loading...")
df = pd.read_csv(RGI_HOST_FILE, low_memory=False)
meta = pd.read_csv(META_FILE, sep="\t")

cat_col = "category" if "category" in df.columns else "category_final"
if cat_col not in df.columns:
    df = df.merge(meta[["sample_id", "category"]], on="sample_id", how="left")
    cat_col = "category"

df["is_mobile"] = (df["on_plasmid"].fillna(False).astype(bool) |
                   df["on_virus"].fillna(False).astype(bool)).astype(int)
df = df.drop_duplicates(subset=["sample_id", "Contig", "Best_Hit_ARO"])

inf  = df[df[cat_col] == MUNI_INF]
eff  = df[df[cat_col] == MUNI_EFF]
muni = df[df[cat_col].isin([MUNI_INF, MUNI_EFF])]
host_inf  = inf[inf["genus"].notna()]
host_eff  = eff[eff["genus"].notna()]
host_muni = muni[muni["genus"].notna()]
mob_inf = host_inf[host_inf["is_mobile"]==1]
mob_eff = host_eff[host_eff["is_mobile"]==1]
fam_col = "family" if "family" in df.columns else None

print(f"  Influent:  ARG-bearing={len(host_inf):>5}, mobile={len(mob_inf):>5}")
print(f"  Effluent:  ARG-bearing={len(host_eff):>5}, mobile={len(mob_eff):>5}")

# ============================================================
# 1. Bray-Curtis 用 relative abundance 计算
# ============================================================
section("1. Bray-Curtis with relative abundance (paper version)")

def bc_relabund(group_inf, group_eff, key="genus"):
    """计算 inf vs eff 的 BC 用 relative abundance"""
    inf_counts = group_inf[key].value_counts()
    eff_counts = group_eff[key].value_counts()
    
    # 转为相对丰度
    inf_rel = inf_counts / inf_counts.sum()
    eff_rel = eff_counts / eff_counts.sum()
    
    common = sorted(set(inf_rel.index) | set(eff_rel.index))
    inf_v = np.array([inf_rel.get(g, 0) for g in common])
    eff_v = np.array([eff_rel.get(g, 0) for g in common])
    
    return bc(inf_v, eff_v)

# (a) ARG-host BC: relative abundance
bc_arg_rel = bc_relabund(host_inf, host_eff)
print(f"  ARG-host BC (relative abundance) = {bc_arg_rel:.3f}")
print(f"     [paper says 0.331, with raw counts gives 0.577]")

# (b) Mobilome BC: relative abundance
bc_mob_rel = bc_relabund(mob_inf, mob_eff)
print(f"  Mobilome BC (relative abundance) = {bc_mob_rel:.3f}")
print(f"     [paper says 0.681, with raw counts gives 0.731]")

# (c) 比值
ratio = bc_mob_rel / bc_arg_rel if bc_arg_rel > 0 else np.nan
print(f"  Mobilome / ARG-host ratio = {ratio:.2f}x")
print(f"     [paper says 2.1x: 0.681/0.331 = 2.06]")

# (d) 也试一下 sample-level mean BC
print()
print("  --- Sample-level pairwise BC (alternative method) ---")

def per_sample_genus_freq(host_df):
    """每个样本的 genus 频率向量"""
    out = {}
    for sid, g in host_df.groupby("sample_id"):
        c = g["genus"].value_counts()
        out[sid] = c / c.sum() if c.sum() > 0 else c
    return out

inf_freqs = per_sample_genus_freq(host_inf)
eff_freqs = per_sample_genus_freq(host_eff)
inf_freqs_mob = per_sample_genus_freq(mob_inf)
eff_freqs_mob = per_sample_genus_freq(mob_eff)

def mean_bc_between(d1, d2):
    """所有 d1 vs d2 样本对的平均 BC"""
    all_g = set()
    for f in list(d1.values()) + list(d2.values()):
        all_g |= set(f.index)
    all_g = sorted(all_g)
    bcs = []
    for s1, f1 in d1.items():
        v1 = np.array([f1.get(g, 0) for g in all_g])
        for s2, f2 in d2.items():
            v2 = np.array([f2.get(g, 0) for g in all_g])
            bcs.append(bc(v1, v2))
    return np.mean(bcs) if bcs else np.nan, np.std(bcs) if bcs else np.nan

bc_arg_mean, bc_arg_std = mean_bc_between(inf_freqs, eff_freqs)
print(f"  ARG-host mean pairwise BC = {bc_arg_mean:.3f} (sd={bc_arg_std:.3f})")

bc_mob_mean, bc_mob_std = mean_bc_between(inf_freqs_mob, eff_freqs_mob)
print(f"  Mobilome mean pairwise BC = {bc_mob_mean:.3f} (sd={bc_mob_std:.3f})")

print(f"  ratio = {bc_mob_mean/bc_arg_mean:.2f}x")

# ============================================================
# 2. Permutation null for Gini - within-sample shuffle
# ============================================================
section("2. Permutation null Gini - within-sample shuffle")

def gini_within_sample_shuffle(host_df_inf, mob_df_inf, n_iter=1000, seed=42):
    """
    论文方法：每个样本内 shuffle genus labels, 保持每个样本的 ARG hit 数和 mobile 数
    然后计算 pooled mobilome 的 Gini
    """
    rng = np.random.RandomState(seed)
    null_ginis = []
    
    # 每个样本：所有 ARG-bearing 的 genus list
    sample_genera = {}
    sample_n_mobile = {}
    for sid, g in host_df_inf.groupby("sample_id"):
        sample_genera[sid] = g["genus"].dropna().values
        sample_n_mobile[sid] = int(g["is_mobile"].sum())
    
    for it in range(n_iter):
        # 模拟 mobilome: 每个样本随机抽 n_mobile 个 genus
        sim_genera = []
        for sid, all_gs in sample_genera.items():
            n_mob = sample_n_mobile[sid]
            if n_mob > 0 and len(all_gs) > 0:
                # 随机抽样（不替换）
                if n_mob >= len(all_gs):
                    sim = all_gs
                else:
                    sim = rng.choice(all_gs, size=n_mob, replace=False)
                sim_genera.extend(sim)
        # 计算 pooled Gini
        if sim_genera:
            counts = pd.Series(sim_genera).value_counts()
            null_ginis.append(gini(counts.values))
    
    return np.array(null_ginis)

print("  Running 1000 within-sample permutations...")
null_ginis = gini_within_sample_shuffle(host_inf, mob_inf, n_iter=1000)
print(f"  Null Gini mean = {null_ginis.mean():.3f}")
print(f"  Null Gini std  = {null_ginis.std():.4f}")

mob_inf_genera = mob_inf["genus"].value_counts()
g_obs = gini(mob_inf_genera.values)
z = (g_obs - null_ginis.mean()) / null_ginis.std() if null_ginis.std() > 0 else 0
print(f"  Observed Gini  = {g_obs:.3f}")
print(f"  Z-score        = {z:.1f}")
print()
print(f"  [paper: null=0.270, Z=30.7]")
print(f"  [first audit (pool sampling): null=0.477, Z=6.7]")
print(f"  [this method (within-sample): null={null_ginis.mean():.3f}, Z={z:.1f}]")

# ============================================================
# 3. 29/61 vs 35/61 Enterobacteriaceae samples
# ============================================================
section("3. Enterobacteriaceae detection in effluent samples")

# Method A: family-level (current default)
if fam_col:
    samples_a = host_eff[host_eff[fam_col] == "Enterobacteriaceae"]["sample_id"].nunique()
    print(f"  Method A (family): {samples_a}/61 samples")
else:
    samples_a = host_eff[host_eff["genus"].isin(ENTERO_GENERA)]["sample_id"].nunique()
    print(f"  Method A (genus proxy): {samples_a}/61 samples")

# Method B: 仅 Escherichia + Klebsiella + Enterobacter + Citrobacter (paper highlights these)
narrow_genera = {"Escherichia", "Klebsiella", "Enterobacter", "Citrobacter"}
samples_b = host_eff[host_eff["genus"].isin(narrow_genera)]["sample_id"].nunique()
print(f"  Method B (narrow: 4 main genera): {samples_b}/61 samples")

# Method C: 要求 ≥X hits (more stringent)
for min_hits in [2, 3, 5]:
    if fam_col:
        sample_counts = host_eff[host_eff[fam_col] == "Enterobacteriaceae"]["sample_id"].value_counts()
    else:
        sample_counts = host_eff[host_eff["genus"].isin(ENTERO_GENERA)]["sample_id"].value_counts()
    n = (sample_counts >= min_hits).sum()
    print(f"  Method C (>= {min_hits} Enterobacteriaceae hits): {n}/61 samples")

# Method D: 仅 Escherichia 单属
samples_d = host_eff[host_eff["genus"] == "Escherichia"]["sample_id"].nunique()
print(f"  Method D (Escherichia only): {samples_d}/61 samples")

print()
print(f"  [paper: 29/61]")
print(f"  Look for which method gives 29/61")

# ============================================================
# 4. 1,315 / 1,414 / 1,514 完整核实
# ============================================================
section("4. Mobile ARG counts - all possible subsets")

print("Subset breakdown:")
print(f"  All 185 samples (post-dedup):")
print(f"    Total ARG hits: {len(df):>6,}")
print(f"    Mobile: {int(df['is_mobile'].sum()):>6,}")
print(f"    Host-assigned: {len(df[df['genus'].notna()]):>6,}")
print(f"    Host-assigned mobile: {int(df[df['genus'].notna()]['is_mobile'].sum()):>6,}")
print()

print(f"  169 municipal (inf+eff):")
print(f"    Total ARG hits: {len(muni):>6,}")
print(f"    Mobile: {int(muni['is_mobile'].sum()):>6,}")
print(f"    Host-assigned: {len(host_muni):>6,}")
print(f"    Host-assigned mobile: {int(host_muni['is_mobile'].sum()):>6,}")
print()

# Per category
for cat in df[cat_col].dropna().unique():
    sub = df[df[cat_col] == cat]
    sub_h = sub[sub["genus"].notna()]
    print(f"  {cat}:")
    print(f"    n_samples = {sub['sample_id'].nunique()}")
    print(f"    ARG hits: {len(sub):>6,}")
    print(f"    Mobile: {int(sub['is_mobile'].sum()):>6,}")
    print(f"    Host-assigned: {len(sub_h):>6,}")
    print(f"    Host-assigned mobile: {int(sub_h['is_mobile'].sum()):>6,}")

# Combinations
print()
print("Possible combinations to match '1,514':")
all_cats = df[cat_col].dropna().unique()
print("  Looking for cat combos giving total host-assigned mobile = 1,514...")
from itertools import combinations
for r in range(2, len(all_cats)+1):
    for combo in combinations(all_cats, r):
        sub = df[df[cat_col].isin(combo)]
        sub_h = sub[sub["genus"].notna()]
        n_mob = int(sub_h["is_mobile"].sum())
        if n_mob == 1514:
            print(f"    EXACT MATCH 1,514: {combo}")
        elif abs(n_mob - 1514) < 50:
            print(f"    Close ({n_mob}): {combo}")

# ============================================================
# 5. vanR/vanS check
# ============================================================
section("5. vanR/vanS plasmid hits")

# Various match patterns
patterns = {
    "Best_Hit_ARO == 'vanR'": df[df["Best_Hit_ARO"].astype(str).str.fullmatch("vanR", case=False, na=False)],
    "Best_Hit_ARO == 'vanS'": df[df["Best_Hit_ARO"].astype(str).str.fullmatch("vanS", case=False, na=False)],
    "Best_Hit_ARO matches '^van[RS]'": df[df["Best_Hit_ARO"].astype(str).str.match(r"^van[RS]\b", case=False, na=False)],
    "Best_Hit_ARO matches '^Van[RS]'": df[df["Best_Hit_ARO"].astype(str).str.match(r"^Van[RS]", case=False, na=False)],
    "Best_Hit_ARO contains 'vanR'": df[df["Best_Hit_ARO"].astype(str).str.contains(r"vanR", case=False, na=False)],
    "Best_Hit_ARO contains 'vanS'": df[df["Best_Hit_ARO"].astype(str).str.contains(r"vanS", case=False, na=False)],
}

for name, sub in patterns.items():
    n_total = len(sub)
    n_pl = int(sub["on_plasmid"].fillna(False).astype(bool).sum())
    print(f"  {name}: n={n_total}, on plasmid={n_pl}")

# What if it's a more complex name? Check if there's anything containing 'van' and S/R
van_all = df[df["Best_Hit_ARO"].astype(str).str.contains("van", case=False, na=False)]
print(f"\n  All 'van*' ARO hits in dataset: {len(van_all)}")
van_aro_uniq = van_all["Best_Hit_ARO"].value_counts().head(20)
print(f"  Top 20 distinct 'van*' ARO names:")
for aro, n in van_aro_uniq.items():
    print(f"    {aro}: {n}")

# ============================================================
# 6. Drug class composition - alternative classification
# ============================================================
section("6. Drug class composition - check priority ordering")

if "Drug Class" in muni.columns:
    # Show raw distribution of unique drug class strings
    print("Top 25 raw 'Drug Class' values in 169 municipal:")
    raw_dc = muni["Drug Class"].fillna("(NA)").value_counts().head(25)
    for v, n in raw_dc.items():
        print(f"  {str(v)[:80]:80s}: {n} ({100*n/len(muni):.2f}%)")
    print()
    
    # Try alternative priority orderings
    # Paper's priority (Section 2.5.1): sulfonamide, phenicol, macrolide, aminoglycoside, 
    #                                    tetracycline, fluoroquinolone, glycopeptide, beta-lactam
    priority_paper = ["sulfonamide", "phenicol", "macrolide", "aminoglycoside",
                      "tetracycline", "fluoroquinolone", "glycopeptide", "beta-lactam"]
    
    # But some classes contain multiple keywords - the priority order matters
    # E.g. "macrolide; tetracycline" -> by paper priority -> macrolide
    # E.g. "tetracycline" alone -> tetracycline
    # E.g. "fluoroquinolone; tetracycline" -> by paper priority -> tetracycline
    
    def classify_by_priority(s, priority):
        if pd.isna(s): return "other"
        s_low = str(s).lower()
        for k in priority:
            if k in s_low:
                return k
            if k == "beta-lactam" and any(x in s_low for x in 
                ["cephalosporin", "carbapenem", "penam", "monobactam", "penicillin"]):
                return "beta-lactam"
        return "other"
    
    muni2 = muni.copy()
    muni2["dc_paper"] = muni2["Drug Class"].apply(lambda s: classify_by_priority(s, priority_paper))
    
    print("\nWith paper's priority order (sulfonamide -> phenicol -> macrolide -> ...):")
    for k in priority_paper + ["other"]:
        n = (muni2["dc_paper"] == k).sum()
        print(f"  {k:18s}: n={n:>6,} ({100*n/len(muni2):.2f}%)")
    
    # Alternative: tetracycline first (if paper uses different order)
    priority_alt = ["sulfonamide", "tetracycline", "phenicol", "macrolide", 
                    "aminoglycoside", "fluoroquinolone", "glycopeptide", "beta-lactam"]
    muni2["dc_alt"] = muni2["Drug Class"].apply(lambda s: classify_by_priority(s, priority_alt))
    
    print("\nAlt priority (tetracycline early - because paper says 14.5%, our default gave 11.27%):")
    for k in priority_alt + ["other"]:
        n = (muni2["dc_alt"] == k).sum()
        print(f"  {k:18s}: n={n:>6,} ({100*n/len(muni2):.2f}%)")

print("\n" + "="*72)
print("  DONE")
print("="*72)

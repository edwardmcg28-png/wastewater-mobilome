"""
validate_all_findings.py

Systematic validation of all core manuscript findings.
Covers:
  1. Effluent Enterobacteriaceae 0/301 result
  2. Gini coefficient and permutation
  3. Three-layer enrichment numbers
  4. Bray-Curtis mobilome vs ARG-host
  5. Mechanism mobility hierarchy
  6. MBI drug class range
  7. Sample counts and data integrity
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr

RGI_PATH  = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy_fixed.csv"
META_PATH = "/QRISdata/Q6636/sra_ww_mobilization/results/sample_map_complete.tsv"
BRACKEN_PATH = "/QRISdata/Q6636/sra_ww_mobilization/results/bracken_genus_abundance.csv"

ENTERO_CORE = {"Escherichia","Klebsiella","Citrobacter","Enterobacter","Serratia"}
ENTERO_EXT  = ENTERO_CORE | {"Rahnella","Kluyvera","Salmonella","Proteus","Hafnia"}

def gini(x):
    x = np.array(x, dtype=float); x = x[x>0]
    if len(x)==0: return 0.0
    x = np.sort(x); n = len(x); idx = np.arange(1,n+1)
    return float((2*np.sum(idx*x)-(n+1)*np.sum(x))/(n*np.sum(x)))

def bray_curtis(p, q):
    p=np.array(p,dtype=float); q=np.array(q,dtype=float)
    d=p.sum()+q.sum()
    return float(np.sum(np.abs(p-q))/d) if d>0 else 0.0

def simplify_mech(x):
    x=str(x).lower()
    if "replacement" in x or "substitution" in x: return "Target replacement"
    if "alteration"  in x or "modification"  in x: return "Target alteration"
    if "inactivation" in x: return "Inactivation"
    if "efflux" in x:       return "Efflux"
    if "protection" in x:   return "Target protection"
    return None

# ── Load ──────────────────────────────────────────────────
print("Loading data...")
rgi     = pd.read_csv(RGI_PATH, low_memory=False)
meta    = pd.read_csv(META_PATH, sep="\t")
bracken = pd.read_csv(BRACKEN_PATH, index_col=0)

# Ensure category column
if "category" not in rgi.columns:
    rgi = rgi.merge(meta[["sample_id","category","country_std"]],
                    on="sample_id", how="left")
else:
    rgi = rgi.merge(meta[["sample_id","country_std"]],
                    on="sample_id", how="left")

rgi["is_mobile"] = (rgi["on_plasmid"].fillna(False) |
                    rgi["on_virus"].fillna(False))

dedup = [c for c in ["sample_id","Contig","Best_Hit_ARO"] if c in rgi.columns]
if dedup: rgi = rgi.drop_duplicates(subset=dedup)

inf_ids = set(meta[meta["category"]=="ww_influent_municipal"]["sample_id"])
eff_ids = set(meta[meta["category"]=="ww_effluent_municipal"]["sample_id"])
muni_ids = inf_ids | eff_ids

rgi_muni = rgi[rgi["sample_id"].isin(muni_ids)]
rgi_inf  = rgi[rgi["sample_id"].isin(inf_ids)]
rgi_eff  = rgi[rgi["sample_id"].isin(eff_ids)]

PASS = "✅ PASS"
FAIL = "❌ FAIL"

results = []

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 1: Dataset overview")
print("="*65)
n_muni = len(muni_ids); n_inf = len(inf_ids); n_eff = len(eff_ids)
n_hits = len(rgi_muni)
n_mobile = rgi_muni["is_mobile"].sum()
mob_rate = n_mobile/n_hits*100
print(f"  Municipal samples: {n_muni} (inf={n_inf}, eff={n_eff})")
print(f"  ARG hits:          {n_hits:,}")
print(f"  Mobile hits:       {n_mobile:,} ({mob_rate:.1f}%)")
r1 = PASS if n_muni==169 and abs(mob_rate-10.6)<0.5 else FAIL
print(f"  Expected: n=169, mob~10.6% → {r1}")
results.append(("Dataset overview", r1))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 2: Effluent Enterobacteriaceae 0/301")
print("="*65)
mob_eff = rgi_eff[rgi_eff["is_mobile"] & rgi_eff["genus"].notna()]
ent_mob_eff = rgi_eff[rgi_eff["genus"].isin(ENTERO_CORE) & rgi_eff["is_mobile"]]
ent_any_eff = rgi_eff[rgi_eff["genus"].isin(ENTERO_CORE) & rgi_eff["genus"].notna()]
n_eff_samps_with_entero = rgi_eff[rgi_eff["genus"].isin(ENTERO_CORE)]["sample_id"].nunique()

print(f"  Effluent mobile hits total:          {len(mob_eff)}")
print(f"  Entero ARG-bearing hits (effluent):  {len(ent_any_eff)} "
      f"(in {n_eff_samps_with_entero}/61 samples)")
print(f"  Entero MOBILE hits (effluent):       {len(ent_mob_eff)}")

# 95% upper bound
n_mob_eff = int(rgi_eff["is_mobile"].sum())
p_upper   = 1 - (0.05 ** (1/n_mob_eff))
print(f"  95% upper bound:                     {p_upper*100:.3f}%")

# Extended genus list
ent_ext_mob = rgi_eff[rgi_eff["genus"].isin(ENTERO_EXT) & rgi_eff["is_mobile"]]
print(f"  Extended Entero list mobile hits:    {len(ent_ext_mob)}")

# Per-country breakdown
print(f"\n  Per-country breakdown:")
for country, grp in rgi_eff[rgi_eff["genus"].notna()].groupby("country_std"):
    em = grp[grp["genus"].isin(ENTERO_CORE) & grp["is_mobile"]]
    tm = grp["is_mobile"].sum()
    ns = grp["sample_id"].nunique()
    print(f"    {country:<15} n={ns:2d}  Entero_mob={len(em)}/{tm}")

r2 = PASS if len(ent_mob_eff)==0 and n_eff_samps_with_entero>0 else FAIL
print(f"\n  0 mobile hits, {n_eff_samps_with_entero} samples with Entero ARG-bearing → {r2}")
results.append(("Effluent Entero 0/301", r2))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 3: Gini coefficient influent mobilome")
print("="*65)
mob_inf = rgi_inf[rgi_inf["is_mobile"] & rgi_inf["genus"].notna()]
obs_gini = gini(mob_inf["genus"].value_counts().values)

# Null: multinomial uniform
n_hits_mob = len(mob_inf)
n_genera   = mob_inf["genus"].nunique()
rng = np.random.default_rng(42)
null = np.array([gini(rng.multinomial(n_hits_mob, [1/n_genera]*n_genera))
                 for _ in range(1000)])
z = (obs_gini - null.mean()) / null.std()
print(f"  Observed Gini:  {obs_gini:.3f}  (paper: 0.538)")
print(f"  Null mean±SD:   {null.mean():.3f}±{null.std():.4f}  (paper: 0.270±0.009)")
print(f"  Z-score:        {z:.1f}  (paper: 30.7)")
r3 = PASS if abs(obs_gini-0.538)<0.01 and z>20 else FAIL
print(f"  → {r3}")
results.append(("Gini influent mobilome", r3))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 4: Three-layer enrichment (Escherichia)")
print("="*65)
br_inf  = bracken[[c for c in bracken.columns if c in inf_ids]]
l1_esc  = float(br_inf.loc["Escherichia"].mean()*100) if "Escherichia" in bracken.index else 0
rgi_inf_host = rgi_inf[rgi_inf["genus"].notna()]
l2_esc  = rgi_inf_host["genus"].value_counts().get("Escherichia",0)/len(rgi_inf_host)*100
l3_esc  = mob_inf["genus"].value_counts().get("Escherichia",0)/len(mob_inf)*100
total_enrich = l3_esc/l1_esc if l1_esc>0 else 0
print(f"  L1 community:    {l1_esc:.2f}%  (paper: 1.15%)")
print(f"  L2 ARG-bearing:  {l2_esc:.2f}%  (paper: 6.01%)")
print(f"  L3 mobile:       {l3_esc:.2f}%  (paper: 15.57%)")
print(f"  Total enrichment:{total_enrich:.2f}×  (paper: 13.57×)")
r4 = PASS if abs(total_enrich-13.57)<1.0 else FAIL
print(f"  → {r4}")
results.append(("Three-layer Escherichia", r4))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 5: Bray-Curtis mobilome > ARG-host")
print("="*65)
all_g   = sorted(set(rgi_inf_host["genus"]) | set(rgi_eff[rgi_eff["genus"].notna()]["genus"]))
rgi_eff_host = rgi_eff[rgi_eff["genus"].notna()]
mob_eff2 = rgi_eff[rgi_eff["is_mobile"] & rgi_eff["genus"].notna()]

def prof(df, g): 
    vc=df["genus"].value_counts(); t=max(len(df),1)
    return np.array([vc.get(x,0)/t for x in g])

bc_argb = bray_curtis(prof(rgi_inf_host,all_g), prof(rgi_eff_host,all_g))
all_gm  = sorted(set(mob_inf["genus"]) | set(mob_eff2["genus"]))
bc_mob  = bray_curtis(prof(mob_inf,all_gm), prof(mob_eff2,all_gm))
print(f"  BC ARG-bearing: {bc_argb:.3f}  (paper: 0.331)")
print(f"  BC Mobilome:    {bc_mob:.3f}  (paper: 0.681)")
print(f"  Ratio:          {bc_mob/bc_argb:.2f}×  (paper: 2.1×)")
r5 = PASS if abs(bc_argb-0.331)<0.05 and abs(bc_mob-0.681)<0.05 else FAIL
print(f"  → {r5}")
results.append(("Bray-Curtis ratio", r5))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 6: Mechanism mobility hierarchy")
print("="*65)
rgi_muni2 = rgi_muni.copy()
rgi_muni2["mech"] = rgi_muni2["Resistance Mechanism"].apply(simplify_mech)
mech_df = rgi_muni2[rgi_muni2["mech"].notna()]
mech_rates = mech_df.groupby("mech")["is_mobile"].mean().sort_values(ascending=False)
print("  Mechanism mobility rates:")
for m,v in mech_rates.items():
    print(f"    {m:<25} {v*100:.1f}%")
# Check hierarchy: replacement > inactivation > protection > efflux > alteration
mechs_ordered = ["Target replacement","Inactivation","Target protection",
                 "Efflux","Target alteration"]
rates = [float(mech_rates.get(m,0)) for m in mechs_ordered]
hierarchy_ok = all(rates[i]>rates[i+1] for i in range(len(rates)-1))
r6 = PASS if hierarchy_ok else FAIL
print(f"  Hierarchy correct: {hierarchy_ok} → {r6}")
results.append(("Mechanism hierarchy", r6))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINDING 7: Enterobacteriaceae influent mobilome 20.4%")
print("="*65)
entero_mob_inf = mob_inf[mob_inf["genus"].isin(ENTERO_EXT)]
entero_pct = len(entero_mob_inf)/len(mob_inf)*100
print(f"  Entero mobile hits: {len(entero_mob_inf)}/{len(mob_inf)} = {entero_pct:.1f}%")
print(f"  Paper value: 20.4%")
r7 = PASS if abs(entero_pct-20.4)<1.0 else FAIL
print(f"  → {r7}")
results.append(("Entero influent mobilome %", r7))

# ══════════════════════════════════════════════════════════
print("\n" + "="*65)
print("SUMMARY")
print("="*65)
for finding, status in results:
    print(f"  {status}  {finding}")
n_pass = sum(1 for _,s in results if "PASS" in s)
print(f"\n  {n_pass}/{len(results)} findings validated")


import os, glob
import pandas as pd
import numpy as np

RGI_CSV = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv"
GENOMAD_BASE = "/QRISdata/Q6636/data"
CATEGORIES = ["ww_influent_municipal","ww_effluent_municipal",
              "ww_influent_hospital","ww_effluent_hospital",
              "ww_sludge","ww_AU_municipal"]

# load RGI
rgi = pd.read_csv(RGI_CSV, usecols=["sample_id","Contig","Best_Hit_ARO","on_plasmid","on_virus"],
                  low_memory=False)
rgi = rgi.drop_duplicates(subset=["sample_id","Contig","Best_Hit_ARO"])
print(f"RGI rows: {len(rgi)}")

# load aggregated_classification scores
rows = []
for cat in CATEGORIES:
    pat = os.path.join(GENOMAD_BASE, cat, "genomad", "*",
                       "*_contigs_aggregated_classification",
                       "*_contigs_aggregated_classification.tsv")
    for f in glob.glob(pat):
        sid = f.split(os.sep)[f.split(os.sep).index("genomad")+1]
        try:
            t = pd.read_csv(f, sep="\t",
                            usecols=["seq_name","plasmid_score","virus_score"])
            t["sample_id"] = sid
            rows.append(t)
        except: pass

scores = pd.concat(rows, ignore_index=True).rename(columns={"seq_name":"Contig"})
print(f"Score rows: {len(scores)}")

# join RGI with scores
merged = rgi.merge(scores, on=["sample_id","Contig"], how="left")
matched = merged["plasmid_score"].notna().sum()
print(f"Matched: {matched}/{len(merged)} ({matched/len(merged)*100:.1f}%)")

# baseline: geNomad summary-based (current paper)
orig = (merged["on_plasmid"].fillna(False) | merged["on_virus"].fillna(False)).mean()*100
print(f"\nBaseline (geNomad summary output): {orig:.2f}%")

# recompute at different thresholds using aggregated scores
print("\nThreshold sensitivity (aggregated_classification.tsv):")
print(f"{'Threshold':<12} {'Mobility%':>10} {'Delta vs baseline':>18}")
for thr in [0.5, 0.6, 0.7, 0.8, 0.9]:
    on_pl = merged["plasmid_score"].fillna(0) >= thr
    on_vi = merged["virus_score"].fillna(0) >= thr
    mob = (on_pl | on_vi).mean() * 100
    delta = mob - orig
    print(f"  >= {thr:.1f}     {mob:>9.2f}%  {delta:>+17.2f}pp")

# check MBI rank stability at key thresholds
print("\nMechanism MBI rank at different thresholds:")
mech_col = "Resistance Mechanism"
if mech_col not in merged.columns:
    # reload with mechanism
    rgi2 = pd.read_csv(RGI_CSV,
                       usecols=["sample_id","Contig","Best_Hit_ARO",
                                "Resistance Mechanism"],
                       low_memory=False)
    rgi2 = rgi2.drop_duplicates(subset=["sample_id","Contig","Best_Hit_ARO"])
    merged = merged.merge(rgi2, on=["sample_id","Contig","Best_Hit_ARO"], how="left")

results = {}
for thr in [0.5, 0.7, 0.9]:
    on_pl = merged["plasmid_score"].fillna(0) >= thr
    on_vi = merged["virus_score"].fillna(0) >= thr
    merged[f"mob_{thr}"] = ((on_pl | on_vi)).astype(int)
    r = merged.groupby(mech_col)[f"mob_{thr}"].mean().sort_values(ascending=False)
    results[thr] = r
    print(f"\n  threshold={thr}:")
    for mech, val in r.items():
        print(f"    {mech[:50]:<50} {val*100:.1f}%")

# check rank correlation between 0.5 and 0.9
mechs = list(results[0.5].index)
r05 = [results[0.5].get(m, 0) for m in mechs]
r09 = [results[0.9].get(m, 0) for m in mechs]
from scipy.stats import spearmanr
rho, p = spearmanr(r05, r09)
print(f"\nMechanism rank Spearman rho (0.5 vs 0.9): {rho:.3f}, p={p:.4f}")
print("If rho > 0.95: mechanism hierarchy is robust to threshold choice")

print("\n[Done]")

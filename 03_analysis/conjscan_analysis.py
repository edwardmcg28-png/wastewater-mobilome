import pandas as pd
import numpy as np
import glob, os

RGI_CSV = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv"
GENOMAD_BASE = "/QRISdata/Q6636/data"

# load RGI (plasmid hits only)
rgi = pd.read_csv(RGI_CSV, low_memory=False)
rgi = rgi.drop_duplicates(subset=["sample_id","Contig","Best_Hit_ARO"])
rgi_plasmid = rgi[rgi["on_plasmid"].fillna(False)].copy()
print(f"Plasmid ARG hits: {len(rgi_plasmid)}")

# load all plasmid_genes.tsv and extract contig-level conjscan/amr flags
rows = []
files = glob.glob(os.path.join(GENOMAD_BASE,"ww_*","genomad","*","*_contigs_summary","*_contigs_plasmid_genes.tsv"))
print(f"Loading {len(files)} plasmid_genes files...")

for i, f in enumerate(files):
    if i % 20 == 0: print(f"  {i}/{len(files)}...")
    sid = f.split(os.sep)[f.split(os.sep).index("genomad")+1]
    try:
        t = pd.read_csv(f, sep="\t", usecols=["gene","annotation_conjscan","annotation_amr"])
        # extract contig name from gene id (format: contig_genenum)
        t["Contig"] = t["gene"].str.rsplit("_", n=1).str[0]
        t["sample_id"] = sid
        t["has_conj"] = t["annotation_conjscan"].notna() & (t["annotation_conjscan"] != "NA")
        t["has_amr_genomad"] = t["annotation_amr"].notna() & (t["annotation_amr"] != "NA")
        # collapse to contig level
        contig_flags = t.groupby(["sample_id","Contig"]).agg(
            has_conj=("has_conj","any"),
            has_amr_genomad=("has_amr_genomad","any"),
            n_conj_genes=("has_conj","sum")
        ).reset_index()
        rows.append(contig_flags)
    except Exception as e:
        pass

contig_annot = pd.concat(rows, ignore_index=True)
print(f"Contig annotations: {len(contig_annot)}")

# join with RGI plasmid hits
merged = rgi_plasmid.merge(contig_annot, on=["sample_id","Contig"], how="left")
matched = merged["has_conj"].notna().sum()
print(f"Matched: {matched}/{len(merged)}")

# key stats
n_total = len(merged)
n_with_conj = merged["has_conj"].fillna(False).sum()
n_without_conj = n_total - n_with_conj

print(f"\n=== ARG on plasmid contigs WITH conjugation genes ===")
print(f"With conj genes:    {n_with_conj}/{n_total} ({n_with_conj/n_total*100:.1f}%)")
print(f"Without conj genes: {n_without_conj}/{n_total} ({n_without_conj/n_total*100:.1f}%)")

# by drug class
print(f"\nConjugation gene co-occurrence by drug class (top 10):")
drug_col = "Drug Class"
if drug_col in merged.columns:
    for kw in ["tetracycline","aminoglycoside","sulfonamide","phenicol","macrolide",
               "beta-lactam","fluoroquinolone","glycopeptide"]:
        sub = merged[merged[drug_col].fillna("").str.lower().str.contains(kw)]
        if len(sub) > 0:
            pct = sub["has_conj"].fillna(False).mean()*100
            print(f"  {kw:<16}: {pct:.1f}% with conj genes (n={len(sub)})")

# by resistance mechanism
print(f"\nConjugation gene co-occurrence by mechanism:")
mech_col = "Resistance Mechanism"
if mech_col in merged.columns:
    for mech, sub in merged.groupby(mech_col):
        pct = sub["has_conj"].fillna(False).mean()*100
        print(f"  {mech[:50]:<50}: {pct:.1f}% (n={len(sub)})")

# by host phylum (if available)
if "phylum" in merged.columns:
    print(f"\nConjugation gene co-occurrence by host phylum (top 8):")
    phy = merged.groupby("phylum").agg(
        n=("has_conj","count"),
        pct_conj=("has_conj", lambda x: x.fillna(False).mean()*100)
    ).sort_values("n", ascending=False).head(8)
    print(phy.to_string())

print("\n[Done]")

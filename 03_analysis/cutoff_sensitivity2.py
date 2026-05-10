import pandas as pd
import numpy as np
import glob, os, re

RGI = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv"
df = pd.read_csv(RGI, low_memory=False)
df = df.drop_duplicates(subset=['sample_id','Contig','Best_Hit_ARO'])
df['is_mobile'] = (df['on_plasmid'].fillna(False) | df['on_virus'].fillna(False)).astype(int)

print("Loading contig lengths...")
fai_files = glob.glob('/QRISdata/Q6636/data/ww_*/assemblies/*.fai')
print(f"Found {len(fai_files)} .fai files")

length_rows = []
for f in fai_files:
    # 从文件名提取sample_id
    basename = os.path.basename(f)
    # 格式1: SRR12345_contigs.fa.fai
    # 格式2: 24_eff_E250051725_L01_562_contigs.fa.fai
    sid = basename.replace('_contigs.fa.fai','').replace('_contigs.fai','')
    try:
        t = pd.read_csv(f, sep='\t', header=None, usecols=[0,1],
                        names=['Contig','contig_length'])
        t['sample_id'] = sid
        length_rows.append(t)
    except:
        pass

lengths = pd.concat(length_rows, ignore_index=True)
print(f"Total contigs loaded: {len(lengths)}")
print(f"Unique sample_ids in lengths: {lengths['sample_id'].nunique()}")
print(f"Unique sample_ids in RGI: {df['sample_id'].nunique()}")

# 检查sample_id匹配情况
rgi_ids = set(df['sample_id'].unique())
fai_ids = set(lengths['sample_id'].unique())
overlap = rgi_ids & fai_ids
print(f"Overlapping sample_ids: {len(overlap)}")
print(f"Sample RGI IDs: {list(rgi_ids)[:3]}")
print(f"Sample FAI IDs: {list(fai_ids)[:3]}")

if len(overlap) == 0:
    print("\nNo overlap - checking ID formats...")
    print("RGI sample_ids (first 5):", list(rgi_ids)[:5])
    print("FAI sample_ids (first 5):", list(fai_ids)[:5])
else:
    merged = df.merge(lengths, on=['sample_id','Contig'], how='left')
    matched = merged['contig_length'].notna().sum()
    print(f"\nMatched: {matched}/{len(merged)}")

    cutoffs = [1000, 3000, 5000, 10000]
    print(f"\nAll ARG hits:")
    print(f"{'Cutoff':<10} {'N_ARG':<10} {'Escherichia%':<15} {'Entero%':<12} {'Mobility%':<12}")
    print("-"*60)
    for cutoff in cutoffs:
        sub = merged[merged['contig_length'] >= cutoff]
        sub_host = sub[sub['genus'].notna()]
        if len(sub_host) == 0:
            continue
        e_pct = (sub_host['genus'] == 'Escherichia').mean() * 100
        entero_pct = sub_host['family'].fillna('').str.contains('Enterobacteriaceae').mean() * 100
        mob_pct = sub['is_mobile'].mean() * 100
        print(f"{cutoff:<10} {len(sub):<10} {e_pct:<15.1f} {entero_pct:<12.1f} {mob_pct:<12.1f}")

    print(f"\nMobile ARG only:")
    print(f"{'Cutoff':<10} {'N_mobile':<12} {'Escherichia%':<15} {'Entero%':<12}")
    print("-"*50)
    for cutoff in cutoffs:
        sub = merged[(merged['contig_length'] >= cutoff) & (merged['is_mobile']==1)]
        sub_host = sub[sub['genus'].notna()]
        if len(sub_host) == 0:
            continue
        e_pct = (sub_host['genus'] == 'Escherichia').mean() * 100
        entero_pct = sub_host['family'].fillna('').str.contains('Enterobacteriaceae').mean() * 100
        print(f"{cutoff:<10} {len(sub):<12} {e_pct:<15.1f} {entero_pct:<12.1f}")


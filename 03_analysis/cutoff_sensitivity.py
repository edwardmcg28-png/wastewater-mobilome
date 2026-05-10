import pandas as pd
import numpy as np
import glob, os

RGI = "/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv"
df = pd.read_csv(RGI, low_memory=False)
df = df.drop_duplicates(subset=['sample_id','Contig','Best_Hit_ARO'])
df['is_mobile'] = (df['on_plasmid'].fillna(False) | df['on_virus'].fillna(False)).astype(int)

# 从.fai文件读取contig长度
print("Loading contig lengths from .fai files...")
fai_files = glob.glob('/QRISdata/Q6636/data/ww_*/assembly/*/*.fai') + \
            glob.glob('/QRISdata/Q6636/data/ww_*/megahit/*/*.fai')
print(f"Found {len(fai_files)} .fai files")

length_rows = []
for f in fai_files:
    sid = None
    # 从路径推断sample_id
    parts = f.split('/')
    for p in parts:
        if p.startswith('SRR') or p.startswith('ERR') or p.startswith('DRR'):
            sid = p
            break
    if sid is None:
        continue
    try:
        t = pd.read_csv(f, sep='\t', header=None, usecols=[0,1],
                        names=['Contig','contig_length'])
        t['sample_id'] = sid
        length_rows.append(t)
    except:
        pass

if len(length_rows) == 0:
    print("No .fai files found, trying alternative paths...")
    fai_files = glob.glob('/QRISdata/Q6636/data/*/*/*.fai')
    print(f"Alternative search found: {len(fai_files)} files")
    print("Sample paths:")
    for f in fai_files[:5]:
        print(f"  {f}")
else:
    lengths = pd.concat(length_rows, ignore_index=True)
    print(f"Contig lengths loaded: {len(lengths)}")
    
    # merge
    merged = df.merge(lengths, on=['sample_id','Contig'], how='left')
    matched = merged['contig_length'].notna().sum()
    print(f"Matched: {matched}/{len(merged)}")
    
    if matched > 0:
        print(f"\nContig length distribution:")
        print(merged['contig_length'].describe())
        
        # 敏感性分析
        cutoffs = [1000, 3000, 5000, 10000]
        print(f"\n{'Cutoff':<10} {'N_ARG':<10} {'Escherichia%':<15} {'Entero%':<12} {'Mobility%':<12}")
        print("-"*60)
        
        for cutoff in cutoffs:
            sub = merged[merged['contig_length'] >= cutoff]
            sub_host = sub[sub['genus'].notna()]
            n = len(sub_host)
            if n == 0:
                continue
            e_pct = (sub_host['genus'] == 'Escherichia').mean() * 100
            entero_pct = sub_host['family'].fillna('').str.contains('Enterobacteriaceae').mean() * 100
            mob_pct = sub['is_mobile'].mean() * 100
            print(f"{cutoff:<10} {len(sub):<10} {e_pct:<15.1f} {entero_pct:<12.1f} {mob_pct:<12.1f}")

        # mobile ARG only
        print(f"\nMobile ARG only:")
        print(f"{'Cutoff':<10} {'N_mobile':<12} {'Escherichia%':<15} {'Entero%':<12}")
        print("-"*50)
        for cutoff in cutoffs:
            sub = merged[(merged['contig_length'] >= cutoff) & (merged['is_mobile']==1)]
            sub_host = sub[sub['genus'].notna()]
            n = len(sub_host)
            if n == 0:
                continue
            e_pct = (sub_host['genus'] == 'Escherichia').mean() * 100
            entero_pct = sub_host['family'].fillna('').str.contains('Enterobacteriaceae').mean() * 100
            print(f"{cutoff:<10} {len(sub):<12} {e_pct:<15.1f} {entero_pct:<12.1f}")


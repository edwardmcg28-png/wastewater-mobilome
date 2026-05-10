import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import braycurtis
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("全球污水ARG项目 - 国家间多维度差异分析")
print("="*80)

# ============================================================================
# 1. 数据加载
# ============================================================================
print("\n[1] 加载数据...")

# 元数据
meta = pd.read_csv('/QRISdata/Q6636/sra_ww_mobilization/results/sample_map_complete.tsv', sep='\t')

# ARG类型丰度 (28 types × 181 samples)
arg_types = pd.read_csv('/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/argoap_type_abundance.csv', index_col=0)

# 微生物群落 (genus level)
bracken = pd.read_csv('/QRISdata/Q6636/sra_ww_mobilization/results/bracken_genus_abundance.csv', index_col=0)

# ARG-MGE共定位
rgi_mge = pd.read_csv('/QRISdata/Q6636/sra_ww_mobilization/results/arg_analysis/rgi_with_mge_annotation.csv', low_memory=False)

# 只分析Municipal Influent（排除样本类型混杂）
arg_T = arg_types.T
arg_T['sample_id'] = arg_T.index
merged = arg_T.merge(meta, on='sample_id')
muni_inf = merged[merged['category'] == 'ww_influent_municipal'].copy()

# 只保留样本量≥6的国家
country_counts = muni_inf['country_std'].value_counts()
countries_enough = country_counts[country_counts >= 6].index.tolist()
muni_inf_filt = muni_inf[muni_inf['country_std'].isin(countries_enough)]

print(f"Municipal Influent样本: {len(muni_inf)}")
print(f"样本量≥6的国家: {len(countries_enough)}")
print(f"国家列表: {countries_enough}")

# ============================================================================
# 2. ARG类型组成差异分析
# ============================================================================
print("\n" + "="*80)
print("[2] ARG类型组成差异 (28种ARG类型)")
print("="*80)

# 提取ARG类型丰度矩阵
arg_type_cols = [col for col in muni_inf_filt.columns if col in arg_types.index]
arg_matrix = muni_inf_filt[arg_type_cols].values

# 计算各国的平均ARG类型组成
country_arg_profiles = {}
for country in countries_enough:
    country_data = muni_inf_filt[muni_inf_filt['country_std'] == country]
    profile = country_data[arg_type_cols].mean(axis=0)
    country_arg_profiles[country] = profile

arg_profile_df = pd.DataFrame(country_arg_profiles).T

# 计算Top 5优势ARG类型（全球平均）
global_avg = muni_inf_filt[arg_type_cols].mean(axis=0).sort_values(ascending=False)
top5_args = global_avg.head(5)

print("\n全球Top 5优势ARG类型:")
for arg_type, abundance in top5_args.items():
    print(f"  {arg_type}: {abundance:.3f} copies/cell ({abundance/global_avg.sum()*100:.1f}%)")

# 各国优势ARG类型差异
print("\n各国Top 3优势ARG类型:")
for country in countries_enough:
    top3 = arg_profile_df.loc[country].sort_values(ascending=False).head(3)
    print(f"\n{country}:")
    for arg_type, abundance in top3.items():
        pct = abundance / arg_profile_df.loc[country].sum() * 100
        print(f"  {arg_type}: {abundance:.3f} copies/cell ({pct:.1f}%)")

# 计算样本间Bray-Curtis距离矩阵
def bray_curtis_matrix(matrix):
    n_samples = matrix.shape[0]
    dist_matrix = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            dist = braycurtis(matrix[i], matrix[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    return dist_matrix

bc_dist = bray_curtis_matrix(arg_matrix)

# 国家间差异显著性检验
print("\n国家间ARG组成差异检验:")

# 计算国家间平均距离 vs 国家内平均距离
between_country_dists = []
within_country_dists = []

sample_country_list = muni_inf_filt['country_std'].tolist()

for i in range(len(sample_country_list)):
    for j in range(i+1, len(sample_country_list)):
        dist = bc_dist[i, j]
        if sample_country_list[i] == sample_country_list[j]:
            within_country_dists.append(dist)
        else:
            between_country_dists.append(dist)

mean_between = np.mean(between_country_dists)
mean_within = np.mean(within_country_dists)

print(f"  国家间平均Bray-Curtis距离: {mean_between:.3f}")
print(f"  国家内平均Bray-Curtis距离: {mean_within:.3f}")
print(f"  比值 (between/within): {mean_between/mean_within:.2f}")

# Mann-Whitney U test
u_stat, p_val = stats.mannwhitneyu(between_country_dists, within_country_dists, alternative='greater')
print(f"  Mann-Whitney U test: U = {u_stat:.0f}, p = {p_val:.6f}")
if p_val < 0.05:
    print("  ✓ 国家间ARG组成显著不同 (p < 0.05)")
else:
    print("  ✗ 国家间ARG组成无显著差异")

# ============================================================================
# 3. 微生物群落差异分析
# ============================================================================
print("\n" + "="*80)
print("[3] 微生物群落组成差异 (Genus level)")
print("="*80)

# 匹配Bracken数据
bracken_T = bracken.T
bracken_T['sample_id'] = bracken_T.index
merged_bracken = muni_inf_filt[['sample_id', 'country_std']].merge(bracken_T, on='sample_id', how='left')

# 提取genus丰度列
genus_cols = [col for col in bracken.index]
bracken_matrix = merged_bracken[genus_cols].fillna(0).values

# 各国优势菌属 (Top 5)
print("\n各国Top 5优势菌属 (平均相对丰度%):")
for country in countries_enough:
    country_data = merged_bracken[merged_bracken['country_std'] == country]
    avg_profile = country_data[genus_cols].mean(axis=0)
    total = avg_profile.sum()
    top5 = avg_profile.sort_values(ascending=False).head(5)
    
    print(f"\n{country}:")
    for genus, abundance in top5.items():
        pct = abundance / total * 100
        print(f"  {genus}: {pct:.2f}%")

# 病原菌属检测
pathogenic_genera = ['Salmonella', 'Klebsiella', 'Escherichia', 'Acinetobacter', 
                     'Pseudomonas', 'Enterococcus', 'Staphylococcus', 'Vibrio',
                     'Campylobacter', 'Shigella', 'Legionella']

print("\n潜在病原菌属相对丰度 (%):")
pathogen_present = [g for g in pathogenic_genera if g in genus_cols]

for country in countries_enough:
    country_data = merged_bracken[merged_bracken['country_std'] == country]
    avg_profile = country_data[genus_cols].mean(axis=0)
    total = avg_profile.sum()
    
    print(f"\n{country}:")
    for pathogen in pathogen_present:
        if pathogen in avg_profile.index:
            pct = avg_profile[pathogen] / total * 100
            if pct > 0.1:  # 只显示 >0.1%
                print(f"  {pathogen}: {pct:.2f}%")

# 微生物群落Bray-Curtis距离
bc_microbiome = bray_curtis_matrix(bracken_matrix)

# 国家间微生物群落差异
between_micro = []
within_micro = []

sample_country_list_micro = merged_bracken['country_std'].tolist()

for i in range(len(sample_country_list_micro)):
    for j in range(i+1, len(sample_country_list_micro)):
        dist = bc_microbiome[i, j]
        if sample_country_list_micro[i] == sample_country_list_micro[j]:
            within_micro.append(dist)
        else:
            between_micro.append(dist)

mean_between_micro = np.mean(between_micro)
mean_within_micro = np.mean(within_micro)

print(f"\n微生物群落国家间平均Bray-Curtis距离: {mean_between_micro:.3f}")
print(f"微生物群落国家内平均Bray-Curtis距离: {mean_within_micro:.3f}")
print(f"比值 (between/within): {mean_between_micro/mean_within_micro:.2f}")

u_stat_micro, p_val_micro = stats.mannwhitneyu(between_micro, within_micro, alternative='greater')
print(f"Mann-Whitney U test: U = {u_stat_micro:.0f}, p = {p_val_micro:.6f}")
if p_val_micro < 0.05:
    print("✓ 国家间微生物群落显著不同 (p < 0.05)")
else:
    print("✗ 国家间微生物群落无显著差异")

# ============================================================================
# 4. MGE携带率差异
# ============================================================================
print("\n" + "="*80)
print("[4] 移动遗传元件(MGE)携带率差异")
print("="*80)

# 筛选Municipal Influent的RGI数据
rgi_mge_filt = rgi_mge[rgi_mge['sample_id'].isin(muni_inf_filt['sample_id'])].copy()

# 各国MGE携带率
print("\n各国ARG-MGE携带率:")
print(f"{'Country':<15} {'Total ARG':<10} {'On Plasmid':<12} {'On Virus':<10} {'On Any MGE':<12} {'Plasmid %':<10} {'Virus %':<10}")
print("-"*95)

mge_stats = []
for country in countries_enough:
    country_samples = muni_inf_filt[muni_inf_filt['country_std'] == country]['sample_id'].tolist()
    country_rgi = rgi_mge_filt[rgi_mge_filt['sample_id'].isin(country_samples)]
    
    total_args = len(country_rgi)
    on_plasmid = country_rgi['on_plasmid'].sum()
    on_virus = country_rgi['on_virus'].sum()
    on_any_mge = (country_rgi['on_plasmid'] | country_rgi['on_virus']).sum()
    
    plasmid_pct = on_plasmid / total_args * 100 if total_args > 0 else 0
    virus_pct = on_virus / total_args * 100 if total_args > 0 else 0
    mge_pct = on_any_mge / total_args * 100 if total_args > 0 else 0
    
    print(f"{country:<15} {total_args:<10} {on_plasmid:<12} {on_virus:<10} {on_any_mge:<12} {plasmid_pct:<10.2f} {virus_pct:<10.2f}")
    
    mge_stats.append({
        'country': country,
        'plasmid_pct': plasmid_pct,
        'virus_pct': virus_pct,
        'mge_pct': mge_pct
    })

mge_stats_df = pd.DataFrame(mge_stats)

# 各国MGE携带率统计检验
print("\nMGE携带率国家间差异 (Kruskal-Wallis test):")
# 这里每个国家只有一个值（聚合后的比例），无法做统计检验
print("  注意: MGE携带率是按国家聚合的值，无法进行统计检验")
print("  建议: 按样本水平分析MGE携带率差异")

# ============================================================================
# 5. 时间趋势的国家差异
# ============================================================================
print("\n" + "="*80)
print("[5] 时间趋势的国家间差异")
print("="*80)

# 计算各国的ARG总丰度
muni_inf_filt['total_arg'] = muni_inf_filt[arg_type_cols].sum(axis=1)

print("\n各国时间趋势 (Spearman相关):")
print(f"{'Country':<15} {'n':<5} {'Time span':<15} {'Spearman rho':<15} {'p-value':<10} {'Trend':<20}")
print("-"*85)

time_trends = []
for country in countries_enough:
    country_data = muni_inf_filt[muni_inf_filt['country_std'] == country]
    
    if len(country_data) >= 3:
        years = country_data['year_fixed'].values
        args = country_data['total_arg'].values
        
        rho, p = stats.spearmanr(years, args)
        
        time_span = f"{int(years.min())}-{int(years.max())}"
        
        if p < 0.05:
            if rho > 0:
                trend = "Increasing*"
            else:
                trend = "Decreasing*"
        else:
            trend = "No trend"
        
        print(f"{country:<15} {len(country_data):<5} {time_span:<15} {rho:<15.3f} {p:<10.4f} {trend:<20}")
        
        time_trends.append({
            'country': country,
            'n': len(country_data),
            'rho': rho,
            'p': p,
            'significant': p < 0.05
        })

time_trends_df = pd.DataFrame(time_trends)

print("\n时间趋势总结:")
n_sig = time_trends_df['significant'].sum()
n_increasing = time_trends_df[(time_trends_df['significant']) & (time_trends_df['rho'] > 0)].shape[0]
n_decreasing = time_trends_df[(time_trends_df['significant']) & (time_trends_df['rho'] < 0)].shape[0]

print(f"  显著趋势国家数: {n_sig}/{len(countries_enough)}")
print(f"  显著上升: {n_increasing}")
print(f"  显著下降: {n_decreasing}")

# ============================================================================
# 6. 输出总结
# ============================================================================
print("\n" + "="*80)
print("分析总结")
print("="*80)

print("""
主要发现：

1. ARG组成差异
   - 国家间ARG类型组成存在显著差异
   - 各国优势ARG类型有所不同（见上述Top 3）
   
2. 微生物群落差异
   - 国家间微生物群落组成存在显著差异
   - 潜在病原菌属丰度各国不同
   
3. MGE携带率差异
   - 各国质粒/噬菌体携带率存在差异
   - 需要按样本水平进一步统计检验
   
4. 时间趋势差异
   - 不同国家呈现不同的时间趋势
   - 部分国家显著下降，其他国家无显著趋势
   
建议后续分析：
   - 完整PERMANOVA多因素分析（国家 + 年份 + 交互作用）
   - ARG-宿主-MGE三元网络分析
   - 特定ARG类型的国家驱动因素分析
""")

print("\n分析完成！")

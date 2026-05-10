#!/usr/bin/env python3
"""
Beta多样性分析 + PCoA可视化
包括PERMANOVA统计检验
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
import os

print("=" * 80)
print("📊 Beta多样性分析 + PCoA")
print("=" * 80)

OUTPUT_DIR = "./beta_diversity_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 设置绘图样式
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# ============================================
# 1. 读取物种丰度数据
# ============================================
print("\n📊 [1/4] 读取物种丰度数据...")

species_df = pd.read_csv('results_2025_flood_analysis/05_species_abundance_all.csv')
print(f"  ✓ 读取了 {len(species_df)} 行物种数据")

# 创建物种丰度矩阵（样本 × 物种）
abundance_matrix = species_df.pivot_table(
    index='Sample',
    columns='name',
    values='fraction_total_reads',
    fill_value=0
)

print(f"  ✓ 丰度矩阵: {abundance_matrix.shape[0]} 样本 × {abundance_matrix.shape[1]} 物种")

# 样本信息
sample_info = species_df[['Sample', 'Label', 'Time', 'Matrix']].drop_duplicates()
sample_info = sample_info.set_index('Sample')

# 确保矩阵行顺序与样本信息一致
abundance_matrix = abundance_matrix.reindex(sample_info.index)

# ============================================
# 2. 计算距离矩阵（Bray-Curtis）
# ============================================
print("\n📊 [2/4] 计算Bray-Curtis距离...")

def bray_curtis_distance(matrix):
    """计算Bray-Curtis距离"""
    distances = []
    n_samples = matrix.shape[0]
    
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            num = np.sum(np.abs(matrix.iloc[i] - matrix.iloc[j]))
            denom = np.sum(matrix.iloc[i] + matrix.iloc[j])
            if denom > 0:
                bc = num / denom
            else:
                bc = 0
            distances.append(bc)
    
    return squareform(distances)

bc_distances = bray_curtis_distance(abundance_matrix)
print(f"  ✓ Bray-Curtis距离矩阵: {bc_distances.shape}")

# 保存距离矩阵
bc_df = pd.DataFrame(
    bc_distances,
    index=abundance_matrix.index,
    columns=abundance_matrix.index
)
bc_df.to_csv(f"{OUTPUT_DIR}/bray_curtis_distance_matrix.csv")

# ============================================
# 3. PCoA分析
# ============================================
print("\n📊 [3/4] 进行PCoA分析...")

# 使用PCA近似PCoA（对于欧氏距离等价）
# 对距离矩阵进行中心化
def pcoa(distance_matrix):
    """经典多维尺度分析（PCoA）"""
    n = distance_matrix.shape[0]
    
    # Double centering
    row_means = distance_matrix.mean(axis=1, keepdims=True)
    col_means = distance_matrix.mean(axis=0, keepdims=True)
    grand_mean = distance_matrix.mean()
    
    B = -0.5 * (distance_matrix - row_means - col_means + grand_mean)
    
    # 特征分解
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    
    # 按特征值降序排列
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # 只保留正特征值
    positive_idx = eigenvalues > 1e-10
    eigenvalues = eigenvalues[positive_idx]
    eigenvectors = eigenvectors[:, positive_idx]
    
    # 计算坐标
    coordinates = eigenvectors * np.sqrt(eigenvalues)
    
    # 计算解释方差比例
    explained_var = eigenvalues / eigenvalues.sum() * 100
    
    return coordinates, explained_var

pcoa_coords, explained_var = pcoa(bc_distances)

print(f"  ✓ PCoA完成")
print(f"  📊 PC1解释: {explained_var[0]:.2f}%")
print(f"  📊 PC2解释: {explained_var[1]:.2f}%")

# 创建PCoA结果DataFrame
pcoa_df = pd.DataFrame(
    pcoa_coords[:, :2],
    index=abundance_matrix.index,
    columns=['PC1', 'PC2']
)

# 添加样本信息
pcoa_df = pcoa_df.join(sample_info)

# 保存PCoA结果
pcoa_df.to_csv(f"{OUTPUT_DIR}/pcoa_coordinates.csv")

# ============================================
# 4. PERMANOVA简化版（组间vs组内距离比较）
# ============================================
print("\n📊 [4/4] 统计检验...")

# 简化的PERMANOVA：比较组间和组内距离
def simple_permanova(distance_matrix, groups):
    """简化的PERMANOVA分析"""
    within_distances = []
    between_distances = []
    
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            dist = distance_matrix[i, j]
            if groups.iloc[i] == groups.iloc[j]:
                within_distances.append(dist)
            else:
                between_distances.append(dist)
    
    if len(within_distances) > 0 and len(between_distances) > 0:
        stat, p_value = mannwhitneyu(between_distances, within_distances, alternative='greater')
        
        return {
            'within_mean': np.mean(within_distances),
            'between_mean': np.mean(between_distances),
            'statistic': stat,
            'p_value': p_value
        }
    return None

# 按时间分组
time_groups = sample_info['Time']
time_permanova = simple_permanova(bc_distances, time_groups)

# 按基质分组
matrix_groups = sample_info['Matrix']
matrix_permanova = simple_permanova(bc_distances, matrix_groups)

print("\nPERMANOVA结果（简化版）:")
print("\n按时间分组:")
if time_permanova:
    print(f"  组内平均距离: {time_permanova['within_mean']:.4f}")
    print(f"  组间平均距离: {time_permanova['between_mean']:.4f}")
    print(f"  P-value: {time_permanova['p_value']:.4f}")
    print(f"  显著性: {'是' if time_permanova['p_value'] < 0.05 else '否'}")

print("\n按基质分组:")
if matrix_permanova:
    print(f"  组内平均距离: {matrix_permanova['within_mean']:.4f}")
    print(f"  组间平均距离: {matrix_permanova['between_mean']:.4f}")
    print(f"  P-value: {matrix_permanova['p_value']:.4f}")
    print(f"  显著性: {'是' if matrix_permanova['p_value'] < 0.05 else '否'}")

# 保存统计结果
stats_results = []
if time_permanova:
    stats_results.append({
        'Grouping': 'Time',
        'Within_Distance': round(time_permanova['within_mean'], 4),
        'Between_Distance': round(time_permanova['between_mean'], 4),
        'P_value': round(time_permanova['p_value'], 4),
        'Significant': 'Yes' if time_permanova['p_value'] < 0.05 else 'No'
    })

if matrix_permanova:
    stats_results.append({
        'Grouping': 'Matrix',
        'Within_Distance': round(matrix_permanova['within_mean'], 4),
        'Between_Distance': round(matrix_permanova['between_mean'], 4),
        'P_value': round(matrix_permanova['p_value'], 4),
        'Significant': 'Yes' if matrix_permanova['p_value'] < 0.05 else 'No'
    })

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv(f"{OUTPUT_DIR}/permanova_results.csv", index=False)

# ============================================
# 5. 可视化
# ============================================
print("\n📊 生成可视化...")

fig = plt.figure(figsize=(16, 6))

# 子图1: 按时间着色
ax1 = plt.subplot(131)
time_colors = {'T0': '#2E7D32', 'T1': '#FFA726', 'T2': '#E53935', 'T3': '#5E35B1'}
time_order = ['T0', 'T1', 'T2', 'T3']

for time in time_order:
    mask = pcoa_df['Time'] == time
    ax1.scatter(pcoa_df.loc[mask, 'PC1'], 
               pcoa_df.loc[mask, 'PC2'],
               c=time_colors[time], 
               label=time, 
               s=200, 
               alpha=0.8, 
               edgecolors='black',
               linewidths=2)

ax1.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)', fontsize=13, fontweight='bold')
ax1.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)', fontsize=13, fontweight='bold')
ax1.set_title('A. PCoA by Time Point', fontsize=14, fontweight='bold')
ax1.legend(title='Time', fontsize=11, title_fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
ax1.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)

# 子图2: 按基质着色
ax2 = plt.subplot(132)
matrix_colors = {'Soil': '#8B4513', 'Water': '#1E90FF'}

for matrix in ['Soil', 'Water']:
    mask = pcoa_df['Matrix'] == matrix
    ax2.scatter(pcoa_df.loc[mask, 'PC1'],
               pcoa_df.loc[mask, 'PC2'],
               c=matrix_colors[matrix],
               label=matrix,
               s=200,
               alpha=0.8,
               edgecolors='black',
               linewidths=2)

ax2.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)', fontsize=13, fontweight='bold')
ax2.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)', fontsize=13, fontweight='bold')
ax2.set_title('B. PCoA by Matrix Type', fontsize=14, fontweight='bold')
ax2.legend(title='Matrix', fontsize=11, title_fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
ax2.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)

# 子图3: 组合（时间+基质）
ax3 = plt.subplot(133)

# Soil样本
for time in time_order:
    mask = (pcoa_df['Time'] == time) & (pcoa_df['Matrix'] == 'Soil')
    if mask.sum() > 0:
        ax3.scatter(pcoa_df.loc[mask, 'PC1'],
                   pcoa_df.loc[mask, 'PC2'],
                   c=time_colors[time],
                   marker='o',
                   s=200,
                   alpha=0.8,
                   edgecolors='black',
                   linewidths=2,
                   label=f'Soil-{time}')

# Water样本
for time in time_order:
    mask = (pcoa_df['Time'] == time) & (pcoa_df['Matrix'] == 'Water')
    if mask.sum() > 0:
        ax3.scatter(pcoa_df.loc[mask, 'PC1'],
                   pcoa_df.loc[mask, 'PC2'],
                   c=time_colors[time],
                   marker='s',
                   s=200,
                   alpha=0.8,
                   edgecolors='black',
                   linewidths=2,
                   label=f'Water-{time}')

ax3.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)', fontsize=13, fontweight='bold')
ax3.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)', fontsize=13, fontweight='bold')
ax3.set_title('C. PCoA Combined View', fontsize=14, fontweight='bold')
ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
ax3.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/pcoa_plot.png", bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/pcoa_plot.pdf", bbox_inches='tight')
plt.close()

print("  ✅ PCoA图已保存")

# ============================================
# 总结
# ============================================
print("\n" + "=" * 80)
print("✅ Beta多样性分析完成！")
print("=" * 80)

print(f"\n📁 输出文件位于: {OUTPUT_DIR}/")
print("\n生成的文件：")
files = [
    'bray_curtis_distance_matrix.csv',
    'pcoa_coordinates.csv',
    'permanova_results.csv',
    'pcoa_plot.png',
    'pcoa_plot.pdf'
]

for f in files:
    if os.path.exists(f"{OUTPUT_DIR}/{f}"):
        size = os.path.getsize(f"{OUTPUT_DIR}/{f}") / 1024
        print(f"  - {f:50s} ({size:8.1f} KB)")

print("\n🔬 主要发现：")
print(f"  • PC1解释: {explained_var[0]:.2f}%")
print(f"  • PC2解释: {explained_var[1]:.2f}%")
if time_permanova and time_permanova['p_value'] < 0.05:
    print(f"  • 时间点间群落结构显著不同 (p={time_permanova['p_value']:.4f})")
if matrix_permanova and matrix_permanova['p_value'] < 0.05:
    print(f"  • Soil与Water群落结构显著不同 (p={matrix_permanova['p_value']:.4f})")

print("\n📝 这张PCoA图可以作为Figure 6！")

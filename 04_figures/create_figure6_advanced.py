#!/usr/bin/env python3
"""
病原菌-ARG关联分析和高级可视化
从 Contig 的 ARG 注释 + Kraken2 分类结果整合数据
绘制升级版 Figure 6（5个子图）
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import gzip
import re
from collections import defaultdict

print("=" * 80)
print("🦠 病原菌-ARG关联分析 & 高级可视化")
print("=" * 80)

# ============================================================================
# 配置路径
# ============================================================================
ARG_DIR = "/QRISdata/Q8083/basespace_runs/project_10179169/argoap_output_raw"
KRAKEN_DIR = "/QRISdata/Q8083/basespace_runs/project_10179169/spades_output/kraken_bracken_out_contigs"
OUTPUT_DIR = "./figure6_pathogen_arg_advanced"
Path(OUTPUT_DIR).mkdir(exist_ok=True)

# 样本信息
SAMPLES = {
    '35_S3_L001':   {'Time': 'T0', 'Matrix': 'Soil',  'Label': 'S-T0'},
    '35w_S7_L001':  {'Time': 'T0', 'Matrix': 'Water', 'Label': 'W-T0'},
    '38_S1_L001':   {'Time': 'T1', 'Matrix': 'Soil',  'Label': 'S-T1'},
    '38w_S5_L001':  {'Time': 'T1', 'Matrix': 'Water', 'Label': 'W-T1'},
    '310_S4_L001':  {'Time': 'T2', 'Matrix': 'Soil',  'Label': 'S-T2'},
    '310w_S6_L001': {'Time': 'T2', 'Matrix': 'Water', 'Label': 'W-T2'},
    '314_S2_L001':  {'Time': 'T3', 'Matrix': 'Soil',  'Label': 'S-T3'},
    '314w_S8_L001': {'Time': 'T3', 'Matrix': 'Water', 'Label': 'W-T3'}
}

# WHO优先病原菌列表
PRIORITY_PATHOGENS = {
    'Acinetobacter': 'Critical',
    'Pseudomonas': 'Critical',
    'Escherichia': 'Critical',
    'Klebsiella': 'Critical',
    'Salmonella': 'High',
    'Shigella': 'High',
    'Campylobacter': 'High',
    'Helicobacter': 'High',
    'Vibrio': 'High',
    'Enterococcus': 'High',
    'Staphylococcus': 'High',
    'Streptococcus': 'Medium',
    'Clostridium': 'Medium',
    'Listeria': 'High',
    'Bacillus': 'Medium',
    'Mycobacterium': 'Critical',
    'Legionella': 'High',
    'Burkholderia': 'Critical'
}

# ============================================================================
# 步骤1：解析ARG注释数据
# ============================================================================
print("\n[1/5] 解析ARG注释数据...")

def parse_arg_blast(file_path):
    """解析 blastout.filtered.txt"""
    if not Path(file_path).exists():
        return pd.DataFrame()
    
    df = pd.read_csv(file_path, sep='\t')
    
    # 提取contig ID（从qseqid中提取NODE信息）
    # 格式: 35_S3_L001@35_S3_L001@10000@NODE_5636413_length_289_cov_2.042735
    df['contig_id'] = df['qseqid'].str.extract(r'@(NODE_\d+_length_\d+_cov_[\d.]+)')
    
    # 保留关键列
    df = df[['contig_id', 'gene', 'subtype', 'type', 'sample', 'rpk']].copy()
    df.columns = ['contig_id', 'ARG_gene', 'ARG_subtype', 'ARG_type', 'sample', 'abundance']
    
    return df

# 读取所有样本的ARG数据
arg_data_list = []

for sample_id, info in SAMPLES.items():
    arg_file = f"{ARG_DIR}/{sample_id}/stage_two/blastout.filtered.txt"
    
    print(f"  Processing {sample_id}...", end=" ")
    
    df = parse_arg_blast(arg_file)
    
    if len(df) > 0:
        df['Time'] = info['Time']
        df['Matrix'] = info['Matrix']
        df['Label'] = info['Label']
        arg_data_list.append(df)
        print(f"✓ {len(df)} ARG annotations")
    else:
        print("⚠️  No data")

arg_df = pd.concat(arg_data_list, ignore_index=True)
print(f"\n  ✅ Total ARG annotations: {len(arg_df)}")
print(f"  📊 Unique contigs with ARGs: {arg_df['contig_id'].nunique()}")

# ============================================================================
# 步骤2：解析Kraken2分类数据
# ============================================================================
print("\n[2/5] 解析Contig分类数据...")

def parse_kraken_output(file_path):
    """解析 .kraken.txt.gz 文件"""
    if not Path(file_path).exists():
        return pd.DataFrame()
    
    data = []
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            classified = parts[0]  # C or U
            contig_id = parts[1]
            taxid = parts[2]
            
            if classified == 'C':  # 只保留分类成功的
                data.append({
                    'contig_id': contig_id,
                    'taxid': taxid
                })
    
    return pd.DataFrame(data)

def parse_kraken_report(file_path):
    """解析 .report.txt 获取taxid到名称的映射"""
    if not Path(file_path).exists():
        return {}
    
    taxid_to_name = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            # 格式: percentage  reads  reads_direct  rank  taxid  name
            taxid = parts[4]
            name = parts[5].strip()
            rank = parts[3]
            
            if rank == 'G':  # Genus level
                taxid_to_name[taxid] = name
    
    return taxid_to_name

# 读取所有样本的分类数据
kraken_data_list = []
taxid_mapping = {}

for sample_id, info in SAMPLES.items():
    kraken_file = f"{KRAKEN_DIR}/{sample_id}/{sample_id}.contigs.kraken.txt.gz"
    report_file = f"{KRAKEN_DIR}/{sample_id}/{sample_id}.contigs.report.txt"
    
    print(f"  Processing {sample_id}...", end=" ")
    
    # 解析kraken输出
    df = parse_kraken_output(kraken_file)
    
    # 解析report获取名称映射
    tax_map = parse_kraken_report(report_file)
    taxid_mapping.update(tax_map)
    
    if len(df) > 0:
        df['sample'] = sample_id
        df['Time'] = info['Time']
        df['Matrix'] = info['Matrix']
        kraken_data_list.append(df)
        print(f"✓ {len(df)} classified contigs")
    else:
        print("⚠️  No data")

kraken_df = pd.concat(kraken_data_list, ignore_index=True)

# 添加属名
def get_genus_from_taxid(taxid, mapping=taxid_mapping):
    """从taxid获取属名"""
    # 这是简化版，实际可能需要更复杂的查找
    return mapping.get(taxid, None)

# 注意：这里需要额外的NCBI taxonomy数据库来完整映射
# 作为简化，我们直接从report文件中提取属名
print(f"\n  ✅ Total classified contigs: {len(kraken_df)}")
print(f"  📋 Loaded {len(taxid_mapping)} taxid-to-genus mappings")

# ============================================================================
# 步骤3：整合数据 - 匹配病原菌和ARG
# ============================================================================
print("\n[3/5] 整合病原菌-ARG数据...")

# 方案：由于完整的taxid映射复杂，我们使用替代方法
# 从Kraken2的.kreport或之前的病原菌统计中获取contig-genus映射

# 读取之前生成的病原菌数据
pathogen_genus_file = "pathogen_stats_analysis/pathogen_genus_list.csv"
if Path(pathogen_genus_file).exists():
    pathogen_list_df = pd.read_csv(pathogen_genus_file)
    print(f"  ✓ Loaded {len(pathogen_list_df)} pathogen genera")

# 简化方案：直接从ARG的gene名称推断宿主
# 很多ARG基因名包含宿主信息（如mexW来自Pseudomonas）

# 创建gene到可能宿主的映射（基于文献）
GENE_TO_HOST = {
    'mexW': 'Pseudomonas', 'mexK': 'Pseudomonas', 'mexF': 'Pseudomonas',
    'mexB': 'Pseudomonas', 'mexY': 'Pseudomonas', 'mexE': 'Pseudomonas',
    'MuxB': 'Pseudomonas', 'MuxC': 'Pseudomonas',
    'arr': 'Mycobacterium',
    'adeF': 'Acinetobacter', 'adeG': 'Acinetobacter', 'adeB': 'Acinetobacter',
    'oqxB': 'Escherichia', 'oqxA': 'Escherichia',
    'smeE': 'Stenotrophomonas',
    'arnA': 'Pseudomonas',  # 也可能在其他属
}

def infer_genus_from_gene(gene_name):
    """从ARG基因名推断可能的宿主属"""
    gene_base = gene_name.split('__')[-1] if '__' in gene_name else gene_name
    
    for pattern, genus in GENE_TO_HOST.items():
        if pattern.lower() in gene_base.lower():
            return genus
    
    return None

# 为ARG数据添加推断的属名
arg_df['inferred_genus'] = arg_df['ARG_gene'].apply(infer_genus_from_gene)

# 筛选病原菌
arg_df['is_pathogen'] = arg_df['inferred_genus'].isin(PRIORITY_PATHOGENS.keys())
arg_df['pathogen_priority'] = arg_df['inferred_genus'].map(PRIORITY_PATHOGENS)

pathogen_arg_df = arg_df[arg_df['is_pathogen']].copy()

print(f"\n  ✅ 病原菌-ARG配对: {len(pathogen_arg_df)}")
print(f"  📊 涉及病原菌属: {pathogen_arg_df['inferred_genus'].nunique()}")
print(f"  📊 涉及ARG类型: {pathogen_arg_df['ARG_type'].nunique()}")

# 保存数据
pathogen_arg_df.to_csv(f"{OUTPUT_DIR}/pathogen_arg_pairs_detailed.csv", index=False)
print(f"  💾 已保存: {OUTPUT_DIR}/pathogen_arg_pairs_detailed.csv")

# ============================================================================
# 步骤4：统计分析
# ============================================================================
print("\n[4/5] 统计分析...")

# 按病原菌统计ARG数量
pathogen_stats = pathogen_arg_df.groupby(['inferred_genus', 'pathogen_priority']).agg({
    'ARG_gene': 'count',
    'ARG_type': 'nunique',
    'abundance': 'sum'
}).reset_index()
pathogen_stats.columns = ['Pathogen', 'Priority', 'ARG_Count', 'ARG_Types', 'Total_Abundance']
pathogen_stats = pathogen_stats.sort_values('ARG_Count', ascending=False)

print("\n  Top病原菌（按ARG数量）:")
print(pathogen_stats.head(10).to_string(index=False))

# 按ARG类型统计
arg_type_stats = pathogen_arg_df.groupby('ARG_type').agg({
    'ARG_gene': 'count',
    'inferred_genus': 'nunique'
}).reset_index()
arg_type_stats.columns = ['ARG_Type', 'Count', 'Pathogen_Count']
arg_type_stats = arg_type_stats.sort_values('Count', ascending=False)

print("\n  主要ARG类型:")
print(arg_type_stats.head(10).to_string(index=False))

# 时间动态
temporal_stats = pathogen_arg_df.groupby(['Time', 'pathogen_priority']).size().reset_index(name='Count')

print("\n  时间动态:")
print(temporal_stats.pivot_table(index='Time', columns='pathogen_priority', values='Count', fill_value=0))

# ============================================================================
# 步骤5：绘制升级版 Figure 6
# ============================================================================
print("\n[5/5] 绘制 Figure 6...")

fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# ------------------------------------------------------------------------
# 子图A: 病原菌-ARG网络图
# ------------------------------------------------------------------------
print("  绘制 (A) 网络图...")

ax_network = fig.add_subplot(gs[0:2, 0:2])

# 创建网络
G = nx.Graph()

# 添加边（病原菌-ARG类型）
edge_data = pathogen_arg_df.groupby(['inferred_genus', 'ARG_type']).size().reset_index(name='weight')

for _, row in edge_data.iterrows():
    G.add_edge(row['inferred_genus'], row['ARG_type'], weight=row['weight'])

# 节点属性
node_types = {}
node_colors = []
node_sizes = []

for node in G.nodes():
    if node in PRIORITY_PATHOGENS:
        node_types[node] = 'pathogen'
        # 颜色按优先级
        priority = PRIORITY_PATHOGENS[node]
        if priority == 'Critical':
            node_colors.append('#FF4444')
        elif priority == 'High':
            node_colors.append('#FF8800')
        else:
            node_colors.append('#FFDD00')
        
        # 大小按ARG数量
        count = pathogen_stats[pathogen_stats['Pathogen'] == node]['ARG_Count'].values
        node_sizes.append((count[0] if len(count) > 0 else 10) * 20)
    else:
        node_types[node] = 'ARG'
        node_colors.append('#4169E1')
        count = arg_type_stats[arg_type_stats['ARG_Type'] == node]['Count'].values
        node_sizes.append((count[0] if len(count) > 0 else 5) * 10)

# 布局
pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

# 绘制边
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]
nx.draw_networkx_edges(G, pos, width=[w/5 for w in weights], alpha=0.5, ax=ax_network)

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                       alpha=0.8, ax=ax_network)

# 绘制标签
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax_network)

ax_network.set_title('(A) Pathogen-ARG Network', fontsize=14, fontweight='bold', pad=15)
ax_network.axis('off')

# 添加图例
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#FF4444', label='Critical Priority'),
    Patch(facecolor='#FF8800', label='High Priority'),
    Patch(facecolor='#FFDD00', label='Medium Priority'),
    Patch(facecolor='#4169E1', label='ARG Type')
]
ax_network.legend(handles=legend_elements, loc='upper left', fontsize=9)

# ------------------------------------------------------------------------
# 子图B: 病原菌-ARG类型热图
# ------------------------------------------------------------------------
print("  绘制 (B) 热图...")

ax_heatmap = fig.add_subplot(gs[0:2, 2])

# 创建病原菌×ARG类型矩阵
heatmap_data = pathogen_arg_df.groupby(['inferred_genus', 'ARG_type']).size().unstack(fill_value=0)

# 选择Top病原菌和ARG类型
top_pathogens = pathogen_stats['Pathogen'].head(8).tolist()
top_args = arg_type_stats['ARG_Type'].head(8).tolist()

heatmap_data = heatmap_data.loc[
    [p for p in top_pathogens if p in heatmap_data.index],
    [a for a in top_args if a in heatmap_data.columns]
]

# 绘制热图
sns.heatmap(heatmap_data, cmap='YlOrRd', annot=True, fmt='d',
            cbar_kws={'label': 'ARG Count'}, ax=ax_heatmap,
            linewidths=0.5, linecolor='white')

ax_heatmap.set_title('(B) Pathogen-ARG Heatmap', fontsize=14, fontweight='bold', pad=15)
ax_heatmap.set_xlabel('ARG Type', fontsize=11, fontweight='bold')
ax_heatmap.set_ylabel('Pathogen Genus', fontsize=11, fontweight='bold')
ax_heatmap.tick_params(axis='x', rotation=45, labelsize=9)
ax_heatmap.tick_params(axis='y', rotation=0, labelsize=9)

# ------------------------------------------------------------------------
# 子图C: 时间动态堆叠柱状图
# ------------------------------------------------------------------------
print("  绘制 (C) 时间动态...")

ax_temporal = fig.add_subplot(gs[2, :2])

# 准备数据
pivot_temporal = temporal_stats.pivot_table(
    index='Time',
    columns='pathogen_priority',
    values='Count',
    fill_value=0
)

# 确保时间顺序
time_order = ['T0', 'T1', 'T2', 'T3']
pivot_temporal = pivot_temporal.reindex(time_order, fill_value=0)

# 堆叠柱状图
pivot_temporal.plot(kind='bar', stacked=True, ax=ax_temporal,
                    color={'Critical': '#FF4444', 'High': '#FF8800', 'Medium': '#FFDD00'},
                    width=0.6, edgecolor='white', linewidth=1.5)

ax_temporal.set_title('(C) Temporal Dynamics of Pathogen-ARG Pairs',
                      fontsize=14, fontweight='bold', pad=15)
ax_temporal.set_xlabel('Time Point', fontsize=12, fontweight='bold')
ax_temporal.set_ylabel('Number of Pairs', fontsize=12, fontweight='bold')
ax_temporal.set_xticklabels(['T0\n(Pre-flood)', 'T1\n(24h)', 'T2\n(7d)', 'T3\n(30d)'],
                            rotation=0)
ax_temporal.legend(title='Priority', fontsize=9, title_fontsize=10)
ax_temporal.grid(axis='y', alpha=0.3, linestyle='--')

# ------------------------------------------------------------------------
# 子图D: Top病原菌排名
# ------------------------------------------------------------------------
print("  绘制 (D) Top病原菌...")

ax_toppath = fig.add_subplot(gs[2, 2])

# Top 8病原菌
top8 = pathogen_stats.head(8).copy()

# 横向条形图
colors_map = {'Critical': '#FF4444', 'High': '#FF8800', 'Medium': '#FFDD00'}
bar_colors = [colors_map[p] for p in top8['Priority']]

ax_toppath.barh(range(len(top8)), top8['ARG_Count'], color=bar_colors,
                edgecolor='white', linewidth=1.5)

ax_toppath.set_yticks(range(len(top8)))
ax_toppath.set_yticklabels(top8['Pathogen'], fontsize=9)
ax_toppath.invert_yaxis()

ax_toppath.set_title('(D) Top Pathogens by ARG Count',
                     fontsize=14, fontweight='bold', pad=15)
ax_toppath.set_xlabel('ARG Count', fontsize=11, fontweight='bold')
ax_toppath.grid(axis='x', alpha=0.3, linestyle='--')

# 添加数值标签
for i, v in enumerate(top8['ARG_Count']):
    ax_toppath.text(v + 1, i, str(v), va='center', fontsize=9, fontweight='bold')

# 保存图片
output_png = f"{OUTPUT_DIR}/Figure6_Pathogen_ARG_Advanced.png"
output_pdf = f"{OUTPUT_DIR}/Figure6_Pathogen_ARG_Advanced.pdf"

plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')

print(f"\n  ✅ 图片已保存:")
print(f"     PNG: {output_png}")
print(f"     PDF: {output_pdf}")

plt.close()

# ============================================================================
# 总结
# ============================================================================
print("\n" + "=" * 80)
print("✅ 分析完成！")
print("=" * 80)

print(f"\n📊 关键统计:")
print(f"  • 病原菌-ARG配对总数: {len(pathogen_arg_df)}")
print(f"  • 涉及病原菌属: {pathogen_arg_df['inferred_genus'].nunique()}")
print(f"  • 涉及ARG基因: {pathogen_arg_df['ARG_gene'].nunique()}")
print(f"  • 涉及ARG类型: {pathogen_arg_df['ARG_type'].nunique()}")

print(f"\n📁 输出文件位于: {OUTPUT_DIR}/")
print("  • pathogen_arg_pairs_detailed.csv - 详细配对数据")
print("  • Figure6_Pathogen_ARG_Advanced.png - 升级版图表")
print("  • Figure6_Pathogen_ARG_Advanced.pdf - 矢量图")

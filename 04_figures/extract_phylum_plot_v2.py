#!/usr/bin/env python3
"""
从 Kraken2/GTDB kreport 文件提取门级分类组成
并绘制 Figure 5C & 5D
改进版：添加详细日志和错误处理
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import re
import os

print("=" * 80)
print("🎨 Figure 5C & 5D - 细菌群落门级分类组成")
print("=" * 80)

# ============================================================================
# 配置
# ============================================================================
KREPORT_DIR = "/QRISdata/Q8083/basespace_runs/project_10179169/spades_output/kraken_bracken_out"
OUTPUT_DIR = os.path.expanduser("~/figure5_output")  # 明确使用home目录

print(f"\n📁 配置:")
print(f"  输入目录: {KREPORT_DIR}")
print(f"  输出目录: {OUTPUT_DIR}")

# 创建输出目录
Path(OUTPUT_DIR).mkdir(exist_ok=True)
print(f"  ✅ 输出目录已创建: {OUTPUT_DIR}")

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

# ============================================================================
# 第1步：从 kreport 文件提取门级数据
# ============================================================================
print("\n[1/4] 从 kreport 文件提取门级丰度...")

def parse_gtdb_kreport(file_path):
    """
    解析 GTDB 格式的 kreport 文件
    格式: k__Kingdom|p__Phylum|c__Class|...    count
    """
    phylum_counts = {}
    total_reads = 0
    kingdom_reads = 0
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    continue
                
                taxonomy = parts[0]
                try:
                    count = int(parts[1])
                except ValueError:
                    continue
                
                # 统计kingdom级别的reads（用于计算总数）
                if taxonomy.startswith('k__') and '|' not in taxonomy:
                    kingdom_reads += count
                
                # 提取门级信息（包含 |p__ 的行）
                if '|p__' in taxonomy:
                    # 提取门名
                    match = re.search(r'\|p__([^|]+)', taxonomy)
                    if match:
                        phylum = match.group(1)
                        
                        # 只统计恰好到门级的行（深度=1，即 k__|p__）
                        depth = taxonomy.count('|')
                        if depth == 1:  # kingdom|phylum
                            phylum_counts[phylum] = count
                            total_reads += count
        
        # 如果没有提取到门级数据，使用kingdom reads作为总数
        if total_reads == 0 and kingdom_reads > 0:
            total_reads = kingdom_reads
            print(f"    ⚠️  使用kingdom级reads: {kingdom_reads:,}")
        
        return phylum_counts, total_reads
    
    except Exception as e:
        print(f"    ❌ 错误: {e}")
        return {}, 0

# 提取所有样本的门级数据
all_data = []
missing_files = []

for sample_id, info in SAMPLES.items():
    kreport_file = f"{KREPORT_DIR}/{sample_id}/{sample_id}.kreport"
    
    print(f"\n  {sample_id} ({info['Label']}):", end=" ")
    
    if not Path(kreport_file).exists():
        print(f"❌ 文件不存在")
        missing_files.append(sample_id)
        continue
    
    phylum_counts, total_reads = parse_gtdb_kreport(kreport_file)
    
    if total_reads == 0:
        print(f"⚠️  没有读取到数据")
        continue
    
    print(f"✓ {len(phylum_counts)}门, {total_reads:,} reads")
    
    # 显示主要门类
    if phylum_counts:
        top3 = sorted(phylum_counts.items(), key=lambda x: x[1], reverse=True)[:3]
        for phylum, count in top3:
            pct = (count / total_reads) * 100
            print(f"      • {phylum:25s} {count:10,} ({pct:5.1f}%)")
    
    # 计算相对百分比
    for phylum, count in phylum_counts.items():
        relative_pct = (count / total_reads) * 100
        all_data.append({
            'Sample': sample_id,
            'Label': info['Label'],
            'Time': info['Time'],
            'Matrix': info['Matrix'],
            'Phylum': phylum,
            'Read_Count': count,
            'Relative_Percent': relative_pct
        })

# 检查是否有缺失文件
if missing_files:
    print(f"\n⚠️  警告：缺失 {len(missing_files)} 个样本的文件:")
    for sample in missing_files:
        print(f"    • {sample}")

# 转换为DataFrame
df = pd.DataFrame(all_data)

if len(df) == 0:
    print("\n❌ 错误：没有提取到任何数据！")
    print("请检查：")
    print("1. kreport文件是否存在")
    print("2. 文件格式是否正确")
    exit(1)

print(f"\n✅ 成功提取数据:")
print(f"  • 总记录数: {len(df)}")
print(f"  • 样本数: {df['Sample'].nunique()}")
print(f"  • 门数: {df['Phylum'].nunique()}")

# 保存原始数据
raw_csv = f"{OUTPUT_DIR}/phylum_abundance_raw.csv"
df.to_csv(raw_csv, index=False)
print(f"  💾 已保存: {raw_csv}")

# ============================================================================
# 第2步：数据汇总和Top门类选择
# ============================================================================
print("\n[2/4] 选择主要门类...")

# 计算每个门的平均相对丰度
phylum_mean = df.groupby('Phylum')['Relative_Percent'].mean().sort_values(ascending=False)

print(f"\n  检测到的所有门（共{len(phylum_mean)}个）:")
for i, (phylum, abundance) in enumerate(phylum_mean.items(), 1):
    print(f"    {i:2d}. {phylum:30s} {abundance:6.2f}%")

# 选择 Top 门类
TOP_N = 6
top_phyla = phylum_mean.head(TOP_N).index.tolist()

print(f"\n  ✅ 选择 Top {TOP_N} 门类用于绘图:")
for i, phylum in enumerate(top_phyla, 1):
    print(f"    {i}. {phylum:30s} (平均 {phylum_mean[phylum]:.2f}%)")

# 创建分组
df['Phylum_Group'] = df['Phylum'].apply(
    lambda x: x if x in top_phyla else 'Others'
)

# 重新汇总
df_grouped = df.groupby(['Sample', 'Label', 'Time', 'Matrix', 'Phylum_Group']).agg({
    'Read_Count': 'sum',
    'Relative_Percent': 'sum'
}).reset_index()

# 保存分组后数据
grouped_csv = f"{OUTPUT_DIR}/phylum_abundance_grouped.csv"
df_grouped.to_csv(grouped_csv, index=False)
print(f"  💾 已保存: {grouped_csv}")

# ============================================================================
# 第3步：数据透视
# ============================================================================
print("\n[3/4] 准备绘图数据...")

# 时间顺序
time_order = ['T0', 'T1', 'T2', 'T3']
df_grouped['Time'] = pd.Categorical(df_grouped['Time'], categories=time_order, ordered=True)
df_grouped = df_grouped.sort_values(['Matrix', 'Time'])

# 创建透视表
pivot_water = df_grouped[df_grouped['Matrix'] == 'Water'].pivot_table(
    index='Phylum_Group',
    columns='Time',
    values='Relative_Percent',
    aggfunc='sum',
    fill_value=0
)

pivot_soil = df_grouped[df_grouped['Matrix'] == 'Soil'].pivot_table(
    index='Phylum_Group',
    columns='Time',
    values='Relative_Percent',
    aggfunc='sum',
    fill_value=0
)

# 按平均丰度排序
pivot_water = pivot_water.loc[pivot_water.mean(axis=1).sort_values(ascending=False).index]
pivot_soil = pivot_soil.loc[pivot_soil.mean(axis=1).sort_values(ascending=False).index]

print("  ✓ 数据透视完成")

# ============================================================================
# 第4步：绘制堆叠柱状图
# ============================================================================
print("\n[4/4] 绘制图表...")

# 配色方案
color_palette = {
    'Proteobacteria': '#4169E1',
    'Pseudomonadota': '#4169E1',
    'Bacteroidota': '#FFA500',
    'Bacteroidetes': '#FFA500',
    'Actinomycetota': '#32CD32',
    'Actinobacteria': '#32CD32',
    'Bacillati': '#32CD32',
    'Firmicutes': '#FFD700',
    'Bacillota': '#FFD700',
    'Acidobacteriota': '#FF69B4',
    'Acidobacteria': '#FF69B4',
    'Verrucomicrobiota': '#8B008B',
    'Verrucomicrobia': '#8B008B',
    'Planctomycetota': '#00CED1',
    'Planctomycetes': '#00CED1',
    'Chloroflexota': '#FF6347',
    'Chloroflexi': '#FF6347',
    'Cyanobacteria': '#4682B4',
    'Cyanobacteriota': '#4682B4',
    'Myxococcota': '#DA70D6',
    'Myxococcus': '#DA70D6',
    'Others': '#D3D3D3'
}

# 创建图形
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for idx, (matrix, pivot, ax) in enumerate(zip(['Water', 'Soil'], 
                                                [pivot_water, pivot_soil], 
                                                axes)):
    # 获取颜色
    colors = [color_palette.get(phylum, '#CCCCCC') for phylum in pivot.index]
    
    # 绘制堆叠柱状图
    pivot.T.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=colors,
        width=0.7,
        edgecolor='white',
        linewidth=1.5,
        legend=False
    )
    
    # 标题
    subplot_letter = 'C' if matrix == 'Water' else 'D'
    ax.set_title(f"({subplot_letter}) {matrix} Samples", 
                 fontsize=16, fontweight='bold', pad=15)
    
    # 坐标轴
    ax.set_xlabel('Time Point', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=14, fontweight='bold')
    
    # X轴标签
    ax.set_xticklabels(['T0\n(Pre-flood)', 'T1\n(24h)', 'T2\n(7d)', 'T3\n(30d)'], 
                       rotation=0, ha='center', fontsize=12)
    
    # Y轴
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.tick_params(axis='y', labelsize=11)
    
    # 网格
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # 移除边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# 图例
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, 
           title='Phylum', 
           loc='center left',
           bbox_to_anchor=(1.0, 0.5),
           fontsize=11,
           title_fontsize=12,
           frameon=True,
           edgecolor='black')

plt.tight_layout(rect=[0, 0, 0.88, 1])

# 保存
output_png = f"{OUTPUT_DIR}/Figure5CD_Taxonomic_Composition.png"
output_pdf = f"{OUTPUT_DIR}/Figure5CD_Taxonomic_Composition.pdf"

try:
    plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✅ PNG已保存: {output_png}")
except Exception as e:
    print(f"  ❌ PNG保存失败: {e}")

try:
    plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
    print(f"  ✅ PDF已保存: {output_pdf}")
except Exception as e:
    print(f"  ❌ PDF保存失败: {e}")

plt.close()

# ============================================================================
# 数据摘要
# ============================================================================
print("\n" + "=" * 80)
print("📊 关键发现")
print("=" * 80)

for matrix in ['Water', 'Soil']:
    pivot = pivot_water if matrix == 'Water' else pivot_soil
    
    print(f"\n【{matrix}】")
    
    for phylum in pivot.index:
        if phylum == 'Others':
            continue
        t0_val = pivot.loc[phylum, 'T0']
        t1_val = pivot.loc[phylum, 'T1']
        
        if t0_val > 0:
            change_t1 = ((t1_val - t0_val) / t0_val) * 100
            if abs(change_t1) > 20:
                direction = "增加" if change_t1 > 0 else "减少"
                print(f"  • {phylum}: T0→T1 {direction} {abs(change_t1):.1f}% "
                      f"({t0_val:.1f}% → {t1_val:.1f}%)")

print("\n" + "=" * 80)
print("✅ 完成！")
print("=" * 80)

print(f"\n📁 所有文件保存在: {OUTPUT_DIR}")
print("\n生成的文件:")
for filename in ['phylum_abundance_raw.csv', 
                 'phylum_abundance_grouped.csv',
                 'Figure5CD_Taxonomic_Composition.png',
                 'Figure5CD_Taxonomic_Composition.pdf']:
    filepath = f"{OUTPUT_DIR}/{filename}"
    if os.path.exists(filepath):
        size = os.path.getsize(filepath) / 1024
        print(f"  ✅ {filename:45s} ({size:8.1f} KB)")
    else:
        print(f"  ❌ {filename:45s} (未生成)")

print(f"\n💡 下载命令（在本地电脑运行）:")
print(f"scp -r uqnzhai@bunya.rcc.uq.edu.au:~/figure5_output ~/Desktop/")

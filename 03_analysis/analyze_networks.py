#!/usr/bin/env python3
"""
网络分析：ARG-宿主和病原菌-ARG网络
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
import os

print("=" * 80)
print("🕸️  网络分析：ARG-宿主-病原菌")
print("=" * 80)

OUTPUT_DIR = "./network_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 设置绘图样式
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# ============================================
# 1. 读取数据
# ============================================
print("\n📊 [1/3] 读取数据...")

arg_host_df = pd.read_csv('arg_host_analysis/arg_host_complete.csv')
print(f"  ✓ ARG-宿主数据: {len(arg_host_df)} 条记录")

# 只看已分类的和病原菌
arg_host_classified = arg_host_df[arg_host_df['host_status'] == 'classified'].copy()
pathogen_args = arg_host_df[arg_host_df['is_pathogen']].copy()

print(f"  ✓ 已分类的ARG: {len(arg_host_classified)}")
print(f"  ✓ 病原菌上的ARG: {len(pathogen_args)}")

# ============================================
# 2. 构建病原菌-ARG网络
# ============================================
print("\n📊 [2/3] 构建病原菌-ARG网络...")

# 创建网络
G_pathogen = nx.Graph()

# 统计边的权重（相同配对出现的次数）
edge_weights = {}

for _, row in pathogen_args.iterrows():
    pathogen = row['host_genus']
    arg_type = row['arg_type']
    
    edge = (pathogen, arg_type)
    edge_weights[edge] = edge_weights.get(edge, 0) + 1

# 添加节点和边
for (pathogen, arg_type), weight in edge_weights.items():
    G_pathogen.add_node(pathogen, node_type='pathogen')
    G_pathogen.add_node(arg_type, node_type='arg')
    G_pathogen.add_edge(pathogen, arg_type, weight=weight)

print(f"  ✓ 网络节点: {G_pathogen.number_of_nodes()}")
print(f"  ✓ 网络边: {G_pathogen.number_of_edges()}")

# 保存网络统计
network_stats = {
    'Network': 'Pathogen-ARG',
    'Nodes': G_pathogen.number_of_nodes(),
    'Edges': G_pathogen.number_of_edges(),
    'Pathogen_Nodes': len([n for n, d in G_pathogen.nodes(data=True) if d['node_type'] == 'pathogen']),
    'ARG_Nodes': len([n for n, d in G_pathogen.nodes(data=True) if d['node_type'] == 'arg']),
    'Density': nx.density(G_pathogen)
}

# ============================================
# 3. 可视化网络
# ============================================
print("\n📊 [3/3] 生成网络图...")

fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# ===== 图1: 病原菌-ARG网络 =====
ax1 = axes[0]

# 使用spring layout
pos = nx.spring_layout(G_pathogen, k=2, iterations=50, seed=42)

# 分离病原菌和ARG节点
pathogen_nodes = [n for n, d in G_pathogen.nodes(data=True) if d['node_type'] == 'pathogen']
arg_nodes = [n for n, d in G_pathogen.nodes(data=True) if d['node_type'] == 'arg']

# 节点大小根据degree
pathogen_sizes = [G_pathogen.degree(n) * 500 for n in pathogen_nodes]
arg_sizes = [G_pathogen.degree(n) * 500 for n in arg_nodes]

# 绘制病原菌节点
nx.draw_networkx_nodes(G_pathogen, pos, 
                       nodelist=pathogen_nodes,
                       node_color='#E53935',
                       node_size=pathogen_sizes,
                       alpha=0.8,
                       ax=ax1,
                       label='Pathogen')

# 绘制ARG节点
nx.draw_networkx_nodes(G_pathogen, pos,
                       nodelist=arg_nodes,
                       node_color='#1E88E5',
                       node_size=arg_sizes,
                       alpha=0.8,
                       ax=ax1,
                       label='ARG Type')

# 绘制边
edge_widths = [G_pathogen[u][v]['weight'] * 0.5 for u, v in G_pathogen.edges()]
nx.draw_networkx_edges(G_pathogen, pos,
                       width=edge_widths,
                       alpha=0.3,
                       ax=ax1)

# 标签
labels = {}
for node in G_pathogen.nodes():
    if G_pathogen.degree(node) > 2:  # 只标注连接度>2的节点
        labels[node] = node

nx.draw_networkx_labels(G_pathogen, pos, 
                       labels=labels,
                       font_size=9,
                       font_weight='bold',
                       ax=ax1)

ax1.set_title('A. Pathogen-ARG Network', fontsize=16, fontweight='bold', pad=20)
ax1.legend(fontsize=12, loc='upper left')
ax1.axis('off')

# ===== 图2: 高连接度节点的子网络 =====
ax2 = axes[1]

# 提取高连接度节点（degree >= 3）
high_degree_nodes = [n for n in G_pathogen.nodes() if G_pathogen.degree(n) >= 3]
G_high = G_pathogen.subgraph(high_degree_nodes).copy()

if G_high.number_of_nodes() > 0:
    pos_high = nx.spring_layout(G_high, k=3, iterations=50, seed=42)
    
    # 分离节点类型
    pathogen_high = [n for n, d in G_high.nodes(data=True) if d['node_type'] == 'pathogen']
    arg_high = [n for n, d in G_high.nodes(data=True) if d['node_type'] == 'arg']
    
    # 节点大小
    pathogen_sizes_high = [G_high.degree(n) * 800 for n in pathogen_high]
    arg_sizes_high = [G_high.degree(n) * 800 for n in arg_high]
    
    # 绘制
    nx.draw_networkx_nodes(G_high, pos_high,
                          nodelist=pathogen_high,
                          node_color='#E53935',
                          node_size=pathogen_sizes_high,
                          alpha=0.9,
                          ax=ax2)
    
    nx.draw_networkx_nodes(G_high, pos_high,
                          nodelist=arg_high,
                          node_color='#1E88E5',
                          node_size=arg_sizes_high,
                          alpha=0.9,
                          ax=ax2)
    
    edge_widths_high = [G_high[u][v]['weight'] for u, v in G_high.edges()]
    nx.draw_networkx_edges(G_high, pos_high,
                          width=edge_widths_high,
                          alpha=0.5,
                          ax=ax2)
    
    # 所有节点都标注
    nx.draw_networkx_labels(G_high, pos_high,
                           font_size=11,
                           font_weight='bold',
                           ax=ax2)
    
    ax2.set_title('B. High-Connectivity Hub Network (degree ≥ 3)', 
                 fontsize=16, fontweight='bold', pad=20)
else:
    ax2.text(0.5, 0.5, 'No high-connectivity nodes\n(degree ≥ 3)',
            ha='center', va='center', fontsize=14)

ax2.axis('off')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/pathogen_arg_network.png", bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/pathogen_arg_network.pdf", bbox_inches='tight')
plt.close()

print("  ✅ 病原菌-ARG网络图已保存")

# ============================================
# 生成Sankey-style条形图
# ============================================
print("\n生成病原菌-ARG流向图...")

fig, ax = plt.subplots(figsize=(14, 10))

# 统计每个病原菌携带的ARG类型
pathogen_arg_counts = pathogen_args.groupby(['host_genus', 'arg_type']).size().reset_index(name='count')
pathogen_arg_counts = pathogen_arg_counts.sort_values(['host_genus', 'count'], ascending=[True, False])

# 获取病原菌列表（按总ARG数排序）
pathogen_totals = pathogen_args.groupby('host_genus').size().sort_values(ascending=True)
pathogens_ordered = pathogen_totals.index.tolist()

# 颜色映射
arg_colors = {
    'multidrug': '#E53935',
    'beta_lactam': '#FB8C00',
    'tetracycline': '#FDD835',
    'aminoglycoside': '#7CB342',
    'polymyxin': '#26C6DA',
    'bacitracin': '#AB47BC',
    'rifamycin': '#8D6E63',
    'macrolide-lincosamide-streptogramin': '#78909C'
}

# 绘制堆叠条形图
y_pos = np.arange(len(pathogens_ordered))
left = np.zeros(len(pathogens_ordered))

for arg_type in arg_colors.keys():
    widths = []
    for pathogen in pathogens_ordered:
        count = pathogen_arg_counts[
            (pathogen_arg_counts['host_genus'] == pathogen) &
            (pathogen_arg_counts['arg_type'] == arg_type)
        ]['count'].sum()
        widths.append(count)
    
    if sum(widths) > 0:
        ax.barh(y_pos, widths, left=left, 
               label=arg_type.replace('_', ' ').replace('-', '-\n'),
               color=arg_colors[arg_type],
               alpha=0.8,
               edgecolor='white',
               linewidth=1.5)
        left += widths

ax.set_yticks(y_pos)
ax.set_yticklabels(pathogens_ordered, fontsize=12, fontweight='bold')
ax.set_xlabel('Number of ARGs', fontsize=13, fontweight='bold')
ax.set_title('Pathogen-ARG Association Profile', fontsize=15, fontweight='bold', pad=20)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10, title='ARG Type', title_fontsize=11)
ax.grid(axis='x', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/pathogen_arg_barplot.png", bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/pathogen_arg_barplot.pdf", bbox_inches='tight')
plt.close()

print("  ✅ 病原菌-ARG条形图已保存")

# ============================================
# 节点中心性分析
# ============================================
print("\n计算网络中心性...")

# Degree centrality
degree_cent = nx.degree_centrality(G_pathogen)
betweenness_cent = nx.betweenness_centrality(G_pathogen)

centrality_df = pd.DataFrame({
    'Node': list(degree_cent.keys()),
    'Degree_Centrality': list(degree_cent.values()),
    'Betweenness_Centrality': list(betweenness_cent.values()),
    'Degree': [G_pathogen.degree(n) for n in degree_cent.keys()],
    'Node_Type': [G_pathogen.nodes[n]['node_type'] for n in degree_cent.keys()]
}).sort_values('Degree', ascending=False)

centrality_df.to_csv(f"{OUTPUT_DIR}/network_centrality.csv", index=False)

print("\n网络中心性Top 10:")
print(centrality_df.head(10).to_string(index=False))

# ============================================
# 保存边列表
# ============================================
edge_list = []
for u, v, data in G_pathogen.edges(data=True):
    edge_list.append({
        'Source': u,
        'Target': v,
        'Weight': data['weight'],
        'Source_Type': G_pathogen.nodes[u]['node_type'],
        'Target_Type': G_pathogen.nodes[v]['node_type']
    })

edge_df = pd.DataFrame(edge_list).sort_values('Weight', ascending=False)
edge_df.to_csv(f"{OUTPUT_DIR}/network_edges.csv", index=False)

# ============================================
# 总结
# ============================================
print("\n" + "=" * 80)
print("✅ 网络分析完成！")
print("=" * 80)

print(f"\n📁 输出文件位于: {OUTPUT_DIR}/")
print("\n生成的文件：")
files = [
    'pathogen_arg_network.png',
    'pathogen_arg_network.pdf',
    'pathogen_arg_barplot.png',
    'pathogen_arg_barplot.pdf',
    'network_centrality.csv',
    'network_edges.csv'
]

for f in files:
    if os.path.exists(f"{OUTPUT_DIR}/{f}"):
        size = os.path.getsize(f"{OUTPUT_DIR}/{f}") / 1024
        print(f"  - {f:50s} ({size:8.1f} KB)")

print("\n🕸️  网络统计：")
print(f"  • 节点总数: {G_pathogen.number_of_nodes()}")
print(f"  • 边总数: {G_pathogen.number_of_edges()}")
print(f"  • 病原菌节点: {len(pathogen_nodes)}")
print(f"  • ARG节点: {len(arg_nodes)}")
print(f"  • 网络密度: {nx.density(G_pathogen):.3f}")
print(f"  • 平均度: {sum(dict(G_pathogen.degree()).values()) / G_pathogen.number_of_nodes():.2f}")

print("\n🎯 Hub节点 (degree > 5):")
hubs = [(n, G_pathogen.degree(n)) for n in G_pathogen.nodes() if G_pathogen.degree(n) > 5]
for node, deg in sorted(hubs, key=lambda x: x[1], reverse=True):
    node_type = G_pathogen.nodes[node]['node_type']
    print(f"  • {node} ({node_type}): degree = {deg}")

print("\n📝 这些网络图可以作为Figure 7-8！")

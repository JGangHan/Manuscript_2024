########### example dataset: progenitor populations ####################
import os
import click
import sys
import igraph
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc

mpl.use('Agg')
sc.settings.verbosity = 3
results_file = os.path.join(outdir,f"{samplename}_paga.h5ad")
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(12, 12), facecolor='white')



### 1. 数据读取
# sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(4, 4), facecolor='white')
# sc.settings.reset()
adata = sc.read_h5ad('./progenitors.h5ad')
# 主成分
# 做不做这一步对结果影响不大
# sc.tl.pca(adata2, svd_solver='arpack')


### 2. 数据预处理
sc.pp.neighbors(adata, n_neighbors =5, n_pcs=30)


### 3. # 必须按照下边顺序
# draw_graph - paga - plot paga - draw_graph - plot
# 1.细胞空间分布
sc.tl.draw_graph(adata)
# 2.PAGA分析
sc.tl.paga(adata, groups='celltype_final')
sc.pl.paga(adata, color=['celltype_final'],
           #cmap = 'green',
           # save=f"_{samplename}.pdf",
           show=True,
           threshold=0.07,
           fontoutline =1,
           edge_width_scale=1)
# 3. 修改配色，PAGA plot
adata.uns['celltype_final_colors'] = ["#ffc682", "#fbc9c9", "#f57777", 
                                      "#b1ff91","#2e7d32", "#a9cf54", "#66bb6a", "#43a047", "#96ed89"]
sc.pl.paga(adata, color=['celltype_final'],
           #cmap = 'green',
           save= "paga_progenitors_line.pdf",
           show=True,
           threshold=0.07,
           fontoutline =1,
           edge_width_scale=1)
# 4.细胞 plot
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['celltype_final'], 
                 legend_loc='on data',
                 save= "paga_progenitors_cell.pdf",
                 show=True, size=30)

### 4. 相关性分析
ax = sc.pl.correlation_matrix(adata, 'celltype_final', show_correlation_numbers = False, 
                              save = 'progenitors_correlation.pdf', figsize=(4.5,3))

### 5. 保存
import dill
dill.dump(adata, open('paga_progenitors.pkl', 'wb'), protocol=-1)
adata = dill.load(open('paga_progenitors.pkl', 'rb'))

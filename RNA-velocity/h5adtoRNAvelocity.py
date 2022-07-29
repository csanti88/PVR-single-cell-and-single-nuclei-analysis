#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 22:55:02 2021

@author: megan
"""

############################################################################
#                create vectors for color palettes                                    
############################################################################

pip install -U ggplot
pip install palettable
import palettable
import ggplot
import random
random.seed(1234)
#color palette for MG
#import seaborn as sns; sns.set_theme()
#sns.palplot(sns.light_palette("purple"))
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm


colors = ("#F8766D", "#00BA38", "#619CFF")
colors = ("#F8766D", "#00BFC4")
colors = ("#2457B2", "#0DBCAF", "#F7D011") #treatment
colors = ("#461B47", "#DB71B8", "#F79B70") #state
colors= ("#2BA9E1", "#FBAF42", "#B7295C") #state
colors =pd.read_csv("/Volumes/MYG_5TB/Pvr/analysis/sc_v_sn_manuscript/colors_pseudotime.csv")
colors = plot_color_gradients('Sequential', ['Purples',])

import pandas as pd

mg_sc_genes = pd.read_csv("/Users/megan/Documents/Blackshaw_Lab/R/PVR/Objects/scVelo/MG_top500.csv")
mg_sn_genes = pd.read_csv("/Users/megan/Documents/Blackshaw_Lab/R/PVR/Objects/scVelo/snMG_top500.csv")

############################################################################
#                  create RNA velocity plots with scVelo 
############################################################################

import scvelo as scv 
#pip install scanpy==1.7.1
#!pip install scvelo --upgrade --quiet
import scanpy as sc

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo', fontsize=0)  # for beautified visualization

################################
#     following tutorial
################################

#comb
adata = scv.read("/Volumes/MYG_5TB/Pvr/analysis/sc_v_sn_manuscript/RNAVelocity/mg_sc.h5ad")
adata

#filter and normalize. Keep genes with minumum 20 counts and top 3000 genes significant genes.
    #####change genes to include to be top 3000 variable genes from the MG subsetted data & remove rod genes.
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=500)
#compute moments for velocity estimation using 5 PCs and 30 neighbors.
    ##### 5 PCs
scv.pp.moments(adata, n_pcs=5, n_neighbors=30) #5 and 30 say same story in comb
#estimate velocity gene-specifically.
scv.tl.velocity(adata)
#compute velocity graph based on cosine similarities
scv.tl.velocity_graph(adata)
#stream, grid, and arrow plot of velocities on embedding
stream = scv.pl.velocity_embedding_stream(adata, linewidth= 3, palette= colors, dpi= 300, basis="umap", color="State", size=100, alpha=0.75)
grid = scv.pl.velocity_embedding_grid(adata, basis="umap", palette= colors, color="State", arrow_color="black", arrow_length=4, arrow_size=3, dpi=300, size=100, alpha=0.75)
arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= colors, color="State", arrow_length=4, arrow_size=5, dpi=300, alpha=0.3)

#color_map = 'Blues'
################################
#          custom genes
################################

#sc_bm
adata = scv.read("/Users/megan/Documents/Blackshaw_Lab/R/PVR/Objects/mg_sc_bm_cc.h5ad")
adata

#filter and normalize. Keep genes with minumum 20 counts and top 3000 genes significant genes.
    #####change genes to include to be top 3000 variable genes from the MG subsetted data & remove rod genes.
scv.pp.filter_and_normalize(adata, min_shared_counts=20, retain_genes= mg_sc_genes, subset_highly_variable=(False))
#compute moments for velocity estimation using 5 PCs and 30 neighbors.
    ##### 5 PCs
scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
#estimate velocity gene-specifically.
scv.tl.velocity(adata)
#compute velocity graph based on cosine similarities
scv.tl.velocity_graph(adata)
#stream, grid, and arrow plot of velocities on embedding
sc_bm_stream = scv.pl.velocity_embedding_stream(adata, linewidth= 3, palette= sc_colors, dpi= 300, basis="umap", color="seurat_clusters")
sc_bm_grid = scv.pl.velocity_embedding_grid(adata, basis="umap", palette= sc_colors, color="seurat_clusters", arrow_color="black", arrow_length=4, arrow_size=3, dpi=500)
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="seurat_clusters", arrow_length=4, arrow_size=3, dpi=1000)

#color by timepoint, plot all clusters
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", arrow_length=4, arrow_size=3, dpi=1000)

#Keep overall shape, doesn't color by cluster
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", groups='0', arrow_length=4, arrow_size=3, dpi=1000)
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", groups='1', arrow_length=4, arrow_size=3, dpi=1000)
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", groups='2', arrow_length=4, arrow_size=3, dpi=1000)
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", groups='3', arrow_length=4, arrow_size=3, dpi=1000)
sc_bm_arrow = scv.pl.velocity_embedding(adata, basis="umap", palette= sc_colors, color="treatment_num", groups='4', arrow_length=4, arrow_size=3, dpi=1000)
import harmony
import palantir
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# Plotting and miscellaneous imports
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import mjc_functions as mjc

# Initialize random seed
import random
random.seed(101)






greatestPalette = ["#f44336","#265bcf","#36b3ba","#ffeb3b","#e91e63","#00cc92",
"#4caf50","#ffb65c","#9c27b0","#03a9f4","#43d61f","#ff9800","#673ab7","#cddc39",
"#81bdc5","#ff5722","#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb",
"#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff",
"#edf3ba","#ffccbd"]

## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
#feature_colors = [(230,230,230), (35,35,142), (255,127,0)]
feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
position=[0, 0.019999, 0.02, 0.55, 1]
my_feature_cmap = mjc.make_cmap(feature_colors, position=position, bit=True)
dot_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
my_dot_cmap = mjc.make_cmap(dot_colors, position=position, bit=True)

n_neighbors = 75
n_jobs = 8
n_components = 10
dot_size = 15 # size of each point in Harmony scatters

adata = sc.read_h5ad('./data/total.mesenchyme.anndata.h5ad')

tp = 'age'

timepoints = ['47', '59', '72', '80', '101', '122', '127', '132']
timepoint_connections = pd.DataFrame(np.array([timepoints[:-1], timepoints[1:]]).T)

data_df = adata.to_df()

# compute the augmented and non-augmented affinity matrices
aug_aff, aff = harmony.core.augmented_affinity_matrix(data_df=adata.to_df(), timepoints=adata.obs[tp], timepoint_connections=timepoint_connections, n_neighbors=n_neighbors, n_jobs=n_jobs, pc_components=n_components)



layout = harmony.plot.force_directed_layout(aug_aff, data_df.index)

#harmFig = harmony.plot.plot_timepoints(layout, adata.obs[tp])

adata.obsm["X_harmony"] = np.asarray(layout)
#adata.obsp["harmony_aff"] = aff
#adata.obsp["harmony_aff_aug"] = aug_aff
adata.uns["harmony_timepoint_var"] = tp
adata.uns["harmony_timepoint_connections"] = np.asarray(timepoint_connections)

figure_dir = './figures/total_mesenchyme'
sc.settings.figdir = figure_dir

sc.pl.embedding(adata, basis='harmony', color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', palette = greatestPalette, edges = False, size = dot_size, alpha = 0.95)

expressed_dict = dict()
	
for gene in adata.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

fig_2A_genes = ['ACTA2','TAGLN','DLL1','F3','NPY','GPX3','FRZB','PDGFRA','NRG1','EGF','EGFR','ERBB2','ERBB3']	
		
genes_to_plot = []	
	
for gene in fig_2A_genes:	
	if gene in expressed_dict:	
		genes_to_plot.append(gene)	
	else:	
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
	
print('Plotting genes:', ', '.join(genes_to_plot),'\n')	

sc.pl.embedding(adata, basis='harmony', color=genes_to_plot, save = '_featurePlots.pdf', show = False, cmap = my_feature_cmap, edges = False, size = dot_size*3, alpha = 0.95)

figure_1_mes_list = ['DCN','COL1A1','COL1A2','RGS5','PDGFRB','ANO1','KIT','HAND2','HAND1','ACTA2','TAGLN','PDGFRA','DLL1','F3','NPY','GPX3','RSPO1','RSPO2','RSPO3','WNT2B','EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','NOG','CHRD','FRZB','GREM1','GREM2','SOX6','CD34']

genes_to_plot = []	
	
for gene in figure_1_mes_list:	
	if gene in expressed_dict:	
		genes_to_plot.append(gene)	
	else:	
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
	
print('Plotting genes:', ', '.join(genes_to_plot),'\n')	

sc.pl.embedding(adata, basis='harmony', color=genes_to_plot, save = '_fig1_mesenchyme_featurePlots.pdf', show = False, cmap = my_feature_cmap, edges = False, size = dot_size*3, alpha = 0.95)

adata.write('./data/Harmony.mesenchyme.anndata.h5ad')















adata = sc.read_h5ad('./data/total.epithelium.anndata.h5ad')

tp = 'age'

timepoints = ['47', '59', '72', '80', '101', '122', '127', '132']
timepoint_connections = pd.DataFrame(np.array([timepoints[:-1], timepoints[1:]]).T)

data_df = adata.to_df()

# compute the augmented and non-augmented affinity matrices
aug_aff, aff = harmony.core.augmented_affinity_matrix(data_df=adata.to_df(), timepoints=adata.obs[tp], timepoint_connections=timepoint_connections, n_neighbors=n_neighbors, n_jobs=n_jobs, pc_components=n_components)



layout = harmony.plot.force_directed_layout(aug_aff, data_df.index)

#harmFig = harmony.plot.plot_timepoints(layout, adata.obs[tp])

adata.obsm["X_harmony"] = np.asarray(layout)
#adata.obsp["harmony_aff"] = aff
#adata.obsp["harmony_aff_aug"] = aug_aff
adata.uns["harmony_timepoint_var"] = tp
adata.uns["harmony_timepoint_connections"] = np.asarray(timepoint_connections)

figure_dir = './figures/total_epithelium'
sc.settings.figdir = figure_dir

sc.pl.embedding(adata, basis='harmony', color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', palette = greatestPalette, edges = False, size = dot_size, alpha = 0.95)

expressed_dict = dict()
	
for gene in adata.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

fig_2A_genes = ['LGR5','OLFM4','FABP2','SI','DPP4','F3','NPY','ACTA2','TAGLN','EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','EGFR','ERBB2','ERBB3','ERBB4','HAND1']	
		
genes_to_plot = []	
	
for gene in fig_2A_genes:	
	if gene in expressed_dict:	
		genes_to_plot.append(gene)	
	else:	
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
	
print('Plotting genes:', ', '.join(genes_to_plot),'\n')	

sc.pl.embedding(adata, basis='harmony', color=genes_to_plot, save = '_featurePlots.pdf', show = False, cmap = my_feature_cmap, edges = False, size = dot_size*3, alpha = 0.95)

adata.write('./data/Harmony.epithelium.anndata.h5ad')

















































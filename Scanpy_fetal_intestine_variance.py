import sys
import bbknn
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
from pathlib import Path
import mjc_functions as mjc
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 3			 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.settings.set_figure_params(dpi_save=1200, dpi=1200)
sns.set(style="white", color_codes=True)

greatestPalette = ["#f44336","#265bcf","#36b3ba","#ffeb3b","#e91e63","#00cc92",
"#4caf50","#ffb65c","#9c27b0","#03a9f4","#43d61f","#ff9800","#673ab7","#cddc39",
"#81bdc5","#ff5722","#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb",
"#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff",
"#edf3ba","#ffccbd"]

palestPalette = ["#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb","#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff","#edf3ba","#ffccbd"]



#############################################################################
##         Flag to toggle rerunning the analysis or just plotting          ##
#############################################################################

rerun = True

#############################################################################
##                          Plot section flags                             ##
#############################################################################

redraw_featureplots = True
redraw_umaps = True

run_marker_analysis = True

expression_cutoff = 0.01  # Sets the threshold for min expression to be shown in plots

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
mistorage_mount_point = '/home/mike/mistorage/'

#############################################################################
##        Adjust these parameters to get a nicer looking UMAP plot         ##
#############################################################################
# UMAP arguments
num_neighbors_use = 30
num_pcs_use = 9
umap_spread = 2
umap_min_dist = 0.5
maxiter = None
umap_gamma=1
random_state = 20120608

paga_init = True

dot_size = 15 # size of each point un UMAP scatters

# Louvain arguments
louv_res = 0.6

# PAGA arguments
size=20
paga_layout='fr'
threshold=0.01
node_size_scale=3

#############################################################################
## Change this to contain all genes you want to see expression for         ##
#############################################################################
genes_of_interest = ['EPCAM','VIM','CDX2','ACTA2','COL1A1','SOX17','T','TOP2A','WT1','CDH5','PECAM1','VWF','KDR','CD34','RGS5', 
'CSPG4','ITGAM','PTPRC','HBB','STMN2','S100B','TUBB3','SPDEF','CHGA','LYZ','MUC2','MUC5','VIL1','SHH','IHH','DHH', 
'HHIP','GLI1','GLI2','SMO','PTCH1','LGR5','OLFM4','PDGFRA','DLL1','F3','NPY','GPX3','BMP4','ANO1','KIT','HAND2','WT1','UPK3B','HES1', 
'HEY1','ID1','ID2','ID3','ID4','NOG','CHRD','GREM1','FOXF1','FOXF2','FOXL1','FOXL2','VEGFA','LOXL2','LAMC1','CYGB','FRZB',
'CTGF','CTSC','C1S','SYNPO2','EVA1B','ACKR3','XIST','DKK1','EGF','NRG1','EGFR']

mes_score_genes = ['TCF21','COL1A1','COL1A2','DCN','ACTA2','TAGLN']

#############################################################################
##                          End of global settings                         ##
#############################################################################

## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
#feature_colors = [(230,230,230), (35,35,142), (255,127,0)]
feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
position=[0, 0.019999, 0.02, 0.55, 1]
my_feature_cmap = mjc.make_cmap(feature_colors, position=position, bit=True)
dot_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
my_dot_cmap = mjc.make_cmap(dot_colors, position=position, bit=True)




if Path('./data/Raw.concatenated.anndata.h5ad').is_file():
	print('Found [./data/Raw.concatenated.anndata.h5ad] loading data from there')
	adata = sc.read_h5ad('./data/Raw.concatenated.anndata.h5ad')
	print('\nConcatenated samples contain...\n', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
else:
	mapping_genome = 'hg19'
	
	adata_47 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2598-31', mapping_genome)
	#sc.pp.downsample_counts(adata_47, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_47 = mjc.Filter_New_Anndata(adata_47, 'Day_47')
	
	adata_59 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2757-2', mapping_genome)
	#sc.pp.downsample_counts(adata_59, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_59 = mjc.Filter_New_Anndata(adata_59, 'Day_59')
	
	adata_72 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2856-1', mapping_genome)
	#sc.pp.downsample_counts(adata_72, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_72 = mjc.Filter_New_Anndata(adata_72, 'Day_72')
	
	adata_80 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2598-24', mapping_genome)
	#sc.pp.downsample_counts(adata_80, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_80 = mjc.Filter_New_Anndata(adata_80, 'Day_80')
	
	adata_101 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2511-2', mapping_genome)
	#sc.pp.downsample_counts(adata_101, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_101 = mjc.Filter_New_Anndata(adata_101, 'Day_101')
	
	adata_122 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2182-5', mapping_genome).concatenate(mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2182-6', mapping_genome))
	#sc.pp.downsample_counts(adata_122, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_122 = mjc.Filter_New_Anndata(adata_122, 'Day_122')
	
	adata_127 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2250-1', mapping_genome).concatenate(mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2250-2', mapping_genome))
	#sc.pp.downsample_counts(adata_127, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_127 = mjc.Filter_New_Anndata(adata_127, 'Day_127')
	
	adata_132 = mjc.Create_Scanpy_Anndata(mistorage_mount_point, '2598-28', mapping_genome)
	#sc.pp.downsample_counts(adata_132, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	adata_132 = mjc.Filter_New_Anndata(adata_132, 'Day_132')
	
	adata = adata_47.concatenate(adata_59, adata_72, adata_80, adata_101, adata_122, adata_127, adata_132)
	
	adata.write('./data/Raw.concatenated.anndata.h5ad')
	
	print('\nConcatenated samples contain...\n', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')


#sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.normalize_total(adata)

## Log transform the data.
sc.pp.log1p(adata)

## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
adata.write('./data/Filtered.concatenated.anndata.h5ad')
adata.raw = adata

sc.tl.score_genes(adata, mes_score_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='mes_score', random_state=0, copy=False, use_raw=False)

## Identify highly-variable genes based on dispersion relative to expression level.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)

## Filter the genes to remove non-variable genes since they are uninformative
adata = adata[:, adata.var['highly_variable']]

## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
sc.pp.regress_out(adata, ['n_counts', 'S_score', 'G2M_score'])

## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
sc.pp.scale(adata, max_value=10)

## Run PCA to compute the default number of components
sc.tl.pca(adata, svd_solver='arpack')

## Rank genes according to contributions to PCs.
sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')

## Draw the PCA elbow plot to determine which PCs to use
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)

## Compute nearest-neighbors
sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)

## fix batch differences based on XX/XY
bbknn.bbknn(adata, batch_key='age', n_pcs=num_pcs_use, neighbors_within_batch=3, copy=False)

## Calculate cell clusters via Louvain algorithm
sc.tl.louvain(adata, resolution = louv_res)

sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color='louvain', save=False, show=False, threshold=threshold, node_size_scale=node_size_scale, node_size_power=0.9, layout=paga_layout)

#sc.tl.umap(adata, init_pos='paga', min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma, random_state=random_state)
sc.tl.umap(adata, min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma, random_state=random_state)

#sc.pl.umap(adata, color='louvain', save = '_clusterIdentity.pdf', show = False, legend_loc = 'on data', edges = True, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_noEdge.pdf', show = False, legend_loc = 'on data', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
sc.pl.umap(adata, color=['louvain', 'age'], save = '_clusterIdentity_age.pdf', show = False, legend_loc = 'right margin', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
sc.pl.umap(adata, color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
sc.pl.umap(adata, color='sex', save = '_sex.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
sc.pl.umap(adata, color='sampleName', save = '_sample.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
sc.pl.umap(adata, color=['n_genes','n_counts','percent_mito'], save = '_stats.pdf', show = False, edges = False, cmap = my_feature_cmap, size = dot_size+10)
sc.pl.umap(adata, color='mes_score', save = '_mesenchyme_score.pdf', show = False, edges = False, cmap = my_feature_cmap, size = dot_size+10)



ageFigumap = plt.figure(dpi=80, figsize=(18,7))
ax1 = ageFigumap.add_subplot(1,3,1)
ax2 = ageFigumap.add_subplot(1,3,2)

sc.pl.umap(adata, color='louvain', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6, ax=ax1)
sc.pl.umap(adata, color='age', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6, ax=ax2)

ax1.set_title('Louvain clusters')
ax2.set_title('Age (days)')


ageFigumap.savefig('UMAP_louvain_age_panels.pdf')













expressed_dict = dict()
	
for gene in adata.raw.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

genes_to_plot = []

for gene in genes_of_interest:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.umap(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)

genes_to_plot = []
	
erbb_pathway = ['EGF','EGFR','ERBB2','ERBB3','ERBB4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','NRG1','NRG2','NRG3','NRG4']

for gene in erbb_pathway:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting ERBB pathway genes:', ' '.join(genes_to_plot),'\n')

sc.pl.umap(adata, color=genes_to_plot, save = '_ERBB_pathway_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)



genes_to_plot = []
	
erbb_pathway = ['COL1A1','ACTA2','TAGLN','DCN']

for gene in erbb_pathway:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting ERBB pathway genes:', ' '.join(genes_to_plot),'\n')

sc.pl.umap(adata, color=genes_to_plot, save = '_fig1_genes_featureplots.pdf', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)


fig_2A_genes = ['ACTA2','TAGLN','DLL1','F3','NPY','GPX3','FRZB']
	
genes_to_plot = []

for gene in fig_2A_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='age', mean_only_expressed=True, save = '_fig2A_DotPlot.pdf', show = False, color_map=my_dot_cmap, dendrogram=False)


fig_3B_umap_genes = ['EPCAM','ALPI','EGF','LGR5']
	
genes_to_plot = []

for gene in fig_3B_umap_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.umap(adata, color=genes_to_plot, save = '_fig3B_UMAP_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)

fig_3B_genes = ['LGR5','OLFM4','FABP2','SI','DPP4','F3','NPY','ACTA2','TAGLN','EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','EGFR','ERBB2','ERBB3','ERBB4']
	
genes_to_plot = []

for gene in fig_3B_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B_DotPlot_logScale.pdf', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.5, log=True)
sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B_DotPlot_linearScale.pdf', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.5, log=False)



fig_3B1_genes = ['LGR5','OLFM4']

genes_to_plot = []

for gene in fig_3B1_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B1_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)

fig_3B2_genes = ['FABP2','SI','DPP4']

genes_to_plot = []

for gene in fig_3B2_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B2_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)


fig_3B3_genes = ['F3','NPY','ACTA2','TAGLN']

genes_to_plot = []

for gene in fig_3B3_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B3_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)




fig_3B4_genes = ['EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG']

genes_to_plot = []

for gene in fig_3B4_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B4_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.15, log=True)



fig_3B5_genes = ['EGFR','ERBB2','ERBB3','ERBB4']

genes_to_plot = []

for gene in fig_3B5_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B5_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)

adata.write('./data/Processed.concatenated.anndata.h5ad')

print('Checking for expression of genes of interest\n')
expressed_dict = dict()

for gene in adata.raw.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

fig_1C_genes = ['EPCAM','CDH1','CDX2','CLDN4','PTPRC','HLA-DRA','ARHGDIB','CORO1A','S100B','PLP1','STMN2','ELAVL4','CDH5','KDR','ECSCR','CLDN5','COL1A2','COL1A2','DCN','ACTA2','TAGLN','ACTG2','MYLK']

genes_to_plot = []

for gene in fig_1C_genes:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')

sc.pl.dotplot(adata, genes_to_plot, color_map = my_feature_cmap, groupby='louvain', var_group_positions=[(0,3),(4,7),(8,11),(12,15),(16,23)], var_group_labels=['Epithelial','Immune','ENS','Endothelial','Mesenchymal'], var_group_rotation=45, use_raw=True, log=True, dendrogram=True, expression_cutoff=expression_cutoff, mean_only_expressed=True, standard_scale='var', show=False, save='_fig_1C.pdf')





hand_picked_markers = ['DLL1','F3','NPY','GPX3','ADAMDEC1','ACTA2','TAGLN','LGR5','OLFM4','SOX6','FRZB','WNT2B','WNT2','WNT3','EGF','ERBB2','ERBB3','ERBB4','EGFR','FOXL1','PDGFRA','DMRT1','HPRT1','MKI67','TOP2A','PCNA','GLI1','PTCH1','SHH','IHH','RSPO2','RSPO3','NRG1','BMP4','BMP2','BMP7']

genes_to_plot = []

for gene in hand_picked_markers:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ', '.join(genes_to_plot),'\n')


sc.tl.dendrogram(adata, 'louvain', n_pcs=num_pcs_use, use_raw=True, cor_method='pearson', linkage_method='complete', key_added='dendrogram_louvain')

sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', n_genes=50, use_raw=True)
	
sc.tl.filter_rank_genes_groups(adata, groupby='louvain', use_raw=True, log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=0.25, min_fold_change=1.25, max_out_group_fraction=0.25)
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered', groupby='louvain', mean_only_expressed=True,  n_genes=6, save = '_markerDotPlots.pdf', show = False, color_map=my_dot_cmap, dendrogram=True)

mjc.write_marker_file(adata)

adata.write('./data/Processed.concatenated.anndata.h5ad')

#mjc.__rank_genes(adata=adata, groupby='louvain')


def standard_plots(adata):
	sc.tl.score_genes(adata, mes_score_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='mes_score', random_state=0, copy=False, use_raw=False)
	## Identify highly-variable genes based on dispersion relative to expression level.
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)
	
	## Filter the genes to remove non-variable genes since they are uninformative
	adata = adata[:, adata.var['highly_variable']]
	
	## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
	sc.pp.regress_out(adata, ['n_counts', 'S_score', 'G2M_score'])
	
	## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
	sc.pp.scale(adata, max_value=10)
	
	## Run PCA to compute the default number of components
	sc.tl.pca(adata, svd_solver='arpack')
	
	## Rank genes according to contributions to PCs.
	sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')
	
	## Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)
	
	## Compute nearest-neighbors
	sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)
	
	## Calculate cell clusters via Louvain algorithm
	sc.tl.louvain(adata, resolution = louv_res)
	
	sc.tl.paga(adata, groups='louvain')
	sc.pl.paga(adata, color='louvain', save=False, show=False, threshold=threshold, node_size_scale=node_size_scale, node_size_power=0.9, layout=paga_layout)
	
	#sc.tl.umap(adata, init_pos='paga', min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma, random_state=random_state)
	sc.tl.umap(adata, min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma, random_state=random_state)
	
	#sc.pl.umap(adata, color='louvain', save = '_clusterIdentity.pdf', show = False, legend_loc = 'on data', edges = True, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_noEdge.pdf', show = False, legend_loc = 'on data', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color=['louvain', 'age'], save = '_clusterIdentity_age.pdf', show = False, legend_loc = 'right margin', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='sex', save = '_sex.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='sampleName', save = '_sample.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color=['n_genes','n_counts','percent_mito'], save = '_stats.pdf', show = False, edges = False, cmap = my_feature_cmap, size = dot_size+10)
	sc.pl.umap(adata, color='mes_score', save = '_mesenchyme_score.pdf', show = False, edges = False, cmap = my_feature_cmap, size = dot_size+10)
	
	expressed_dict = dict()
	
	for gene in adata.raw.var_names.values.tolist():
		if gene not in expressed_dict:
			expressed_dict[str(gene)] = 1
	
	genes_to_plot = []
	
	for gene in genes_of_interest:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)
	
	genes_to_plot = []	
	
	erbb_pathway = ['EGF','EGFR','ERBB2','ERBB3','ERBB4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','NRG1','NRG2','NRG3','NRG4']	
		
	for gene in erbb_pathway:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting ERBB pathway genes:', ' '.join(genes_to_plot),'\n')	
		
	sc.pl.umap(adata, color=genes_to_plot, save = '_ERBB_pathway_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)	
		
		
		
	genes_to_plot = []	
			
	erbb_pathway = ['COL1A1','ACTA2','TAGLN','DCN']	
		
	for gene in erbb_pathway:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting ERBB pathway genes:', ' '.join(genes_to_plot),'\n')	
		
	sc.pl.umap(adata, color=genes_to_plot, save = '_fig1_genes_featureplots.pdf', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)	
		
		
	fig_2A_genes = ['ACTA2','TAGLN','DLL1','F3','NPY','GPX3','FRZB']	
			
	genes_to_plot = []	
		
	for gene in fig_2A_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='age', mean_only_expressed=True, save = '_fig2A_DotPlot.pdf', show = False, color_map=my_dot_cmap, dendrogram=False)	
		
		
	fig_3B_umap_genes = ['EPCAM','ALPI','EGF','LGR5']	
			
	genes_to_plot = []	
		
	for gene in fig_3B_umap_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.umap(adata, color=genes_to_plot, save = '_fig3B_UMAP_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*2, use_raw = True)	
		
	fig_3B_genes = ['LGR5','OLFM4','FABP2','SI','DPP4','F3','NPY','ACTA2','TAGLN','EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG','EGFR','ERBB2','ERBB3','ERBB4']	
			
	genes_to_plot = []	
		
	for gene in fig_3B_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B_DotPlot_logScale.pdf', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.5, log=True)	
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B_DotPlot_linearScale.pdf', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.5, log=False)	
		
		
		
	fig_3B1_genes = ['LGR5','OLFM4']	
		
	genes_to_plot = []	
		
	for gene in fig_3B1_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B1_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)	
		
	fig_3B2_genes = ['FABP2','SI','DPP4']	
		
	genes_to_plot = []	
		
	for gene in fig_3B2_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B2_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)	
		
		
	fig_3B3_genes = ['F3','NPY','ACTA2','TAGLN']	
		
	genes_to_plot = []	
		
	for gene in fig_3B3_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B3_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)	
		
		
		
		
	fig_3B4_genes = ['EGF','NRG1','NRG2','NRG3','NRG4','TGFA','HBEGF','AREG','BTC','EPGN','EREG']	
		
	genes_to_plot = []	
		
	for gene in fig_3B4_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B4_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, dot_max=0.15, log=True)	
		
		
		
	fig_3B5_genes = ['EGFR','ERBB2','ERBB3','ERBB4']	
		
	genes_to_plot = []	
		
	for gene in fig_3B5_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_fig3B5_DotPlot_logScale.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)
	
	figure_1_mes_list = ['DCN','COL1A1','COL1A2','RGS5','PDGFRB','ANO1','KIT','HAND2','ACTA2','TAGLN','PDGFRA','DLL1','F3','NPY','GPX3']
	
	genes_to_plot = []	
		
	for gene in figure_1_mes_list:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_mesenchyme.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)
	
	figure_1_mes_list2 = ['DCN','COL1A1','COL1A2','RGS5','PDGFRB','ACTA2','TAGLN','PDGFRA','DLL1','F3','NPY','GPX3']
	
	genes_to_plot = []	
		
	for gene in figure_1_mes_list2:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_mesenchyme2.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)
	
	cell_type_genes = ['S100B','PLP1','STMN2','ELAVL4','CHD5','KDR','ECSCR','CLDN5','COL1A1','COL1A2','DCN','ACTA2','TAGLN','ACTG2','MYLK','EPCAM','CDH1','CDX2','CLDN4','PTPRC','HLA-DRA','ARHGDIB','CORO1A']
	
	genes_to_plot = []	
		
	for gene in cell_type_genes:	
		if gene in expressed_dict:	
			genes_to_plot.append(gene)	
		else:	
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')	
		
	print('Plotting genes:', ', '.join(genes_to_plot),'\n')	
		
	sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', mean_only_expressed=True, save = '_cell_type_genes.pdf', standard_scale='var', show = False, color_map=my_dot_cmap, dendrogram=True, log=True)
	
	sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', n_genes=50, use_raw=True)
	sc.tl.filter_rank_genes_groups(adata, groupby='louvain', use_raw=True, log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=0.5, min_fold_change=1.25, max_out_group_fraction=0.9)
	sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered', groupby='louvain', mean_only_expressed=True,  n_genes=6, save = '_markerDotPlots.pdf', show = False, color_map=my_dot_cmap, dendrogram=True)
	mjc.write_marker_file(adata)




init_adata = sc.read_h5ad('./data/Filtered.concatenated.anndata.h5ad')

init_adata.obs['louvain'] = adata.obs['louvain']

init_adata[init_adata.obs['louvain'].isin(['0','1','2','3'])].write('./data/total.mesenchyme.anndata.h5ad')

init_adata[init_adata.obs['louvain'].isin(['3'])].write('./data/sec.mesenchyme.anndata.h5ad')

mes_adata = sc.read_h5ad('./data/total.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['47'])].write('./data/day47.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['59'])].write('./data/day59.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['72'])].write('./data/day72.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['80'])].write('./data/day80.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['101'])].write('./data/day101.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['122'])].write('./data/day122.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['127'])].write('./data/day127.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['132'])].write('./data/day132.mesenchyme.anndata.h5ad')


init_adata[init_adata.obs['louvain'].isin(['4','5'])].write('./data/total.epithelium.anndata.h5ad')

epi_adata = sc.read_h5ad('./data/total.epithelium.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['47'])].write('./data/day47.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['59'])].write('./data/day59.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['72'])].write('./data/day72.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['80'])].write('./data/day80.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['101'])].write('./data/day101.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['122'])].write('./data/day122.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['127'])].write('./data/day127.epi.anndata.h5ad')
epi_adata[epi_adata.obs['age'].isin(['132'])].write('./data/day132.epi.anndata.h5ad')


init_adata[init_adata.obs['age'].isin(['47'])].write('./data/day47.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['59'])].write('./data/day59.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['72'])].write('./data/day72.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['80'])].write('./data/day80.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['101'])].write('./data/day101.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['122'])].write('./data/day122.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['127'])].write('./data/day127.total.anndata.h5ad')
init_adata[init_adata.obs['age'].isin(['132'])].write('./data/day132.total.anndata.h5ad')


outFile = open('cellcounts.csv','w')



figure_dir = './figures/total_mesenchyme'
sc.settings.figdir = figure_dir
adata = sc.read_h5ad('./data/total.mesenchyme.anndata.h5ad')
adata.raw = adata
standard_plots(adata)
outFile.write('total_mesenchyme\n')
outFile.write(str(len(adata.obs_names)) + '\n')

'''
sce.tl.harmony_timeseries(adata, tp="age", n_components=None)
harmfig = plt.figure(figsize=(10,4))
harmFig = sce.pl.harmony_timeseries(adata, return_fig=True)
harmFig.savefig(figure_dir + '/harmony.pdf')
'''


figure_dir = './figures/total_epithelium'
sc.settings.figdir = figure_dir
adata = sc.read_h5ad('./data/total.epithelium.anndata.h5ad')
adata.raw = adata
standard_plots(adata)
outFile.write('total_epithelium\n')
outFile.write(str(len(adata.obs_names)) + '\n')

'''
sce.tl.harmony_timeseries(adata, tp="age", n_components=None)
harmfig = plt.figure(figsize=(10,4))
harmFig = sce.pl.harmony_timeseries(adata, return_fig=True)
harmFig.savefig(figure_dir + '/harmony.pdf')
'''

stage_list = ['47','59','72','80','101','122','127','132']

for stage in stage_list:
	figure_dir = './figures/' + stage + '_mesenchyme'
	sc.settings.figdir = figure_dir
	adata = sc.read_h5ad('./data/day' + stage + '.mesenchyme.anndata.h5ad')
	adata.raw = adata
	adata.T.to_df().to_csv('./harmony_in/day' + stage + '_matrix.csv')
	standard_plots(adata)
	outFile.write(stage + ' mesenchyme\n')
	outFile.write(str(len(adata.obs_names)) + '\n')
	

for stage in stage_list:
	figure_dir = './figures/' + stage + '_epithelium'
	sc.settings.figdir = figure_dir
	adata = sc.read_h5ad('./data/day' + stage + '.epi.anndata.h5ad')
	adata.raw = adata
	standard_plots(adata)
	outFile.write(stage + ' epithelium\n')
	outFile.write(str(len(adata.obs_names)) + '\n')

for stage in stage_list:
	figure_dir = './figures/' + stage + '_total'
	sc.settings.figdir = figure_dir
	adata = sc.read_h5ad('./data/day' + stage + '.total.anndata.h5ad')
	adata.raw = adata
	standard_plots(adata)
	outFile.write(stage + ' total\n')
	outFile.write(str(len(adata.obs_names)) + '\n')


outFile.close()


def sec_analysis(adata):
	## Identify highly-variable genes based on dispersion relative to expression level.
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)
	
	## Filter the genes to remove non-variable genes since they are uninformative
	adata = adata[:, adata.var['highly_variable']]
	
	## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
	sc.pp.regress_out(adata, ['n_counts', 'S_score', 'G2M_score'])
	
	## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
	sc.pp.scale(adata, max_value=10)
	
	## Run PCA to compute the default number of components
	sc.tl.pca(adata, svd_solver='arpack')
	
	## Rank genes according to contributions to PCs.
	sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')
	
	## Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)
	
	## Compute nearest-neighbors
	sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)
	
	## Calculate cell clusters via Louvain algorithm
	adata.obs['louvain'] = adata.obs['age']
	
	sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', n_genes=50, use_raw=True)
	sc.tl.filter_rank_genes_groups(adata, groupby='louvain', use_raw=True, log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=0.25, min_fold_change=1.25, max_out_group_fraction=0.9)
	sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered', groupby='louvain', mean_only_expressed=True,  n_genes=6, save = '_markerDotPlots.pdf', show = False, color_map=my_dot_cmap, dendrogram=True)
	mjc.write_marker_file(adata)



init_adata = sc.read_h5ad('./data/sec.mesenchyme.anndata.h5ad')
init_adata.obs['louvain'] = adata.obs['age']

mes_adata[mes_adata.obs['age'].isin(['47','59'])].write('./data/day47_59.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['59','72'])].write('./data/day59_72.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['72','80'])].write('./data/day72_80.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['80','101'])].write('./data/day80_101.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['101','122'])].write('./data/day101_122.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['122','127'])].write('./data/day122_127.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['127','132'])].write('./data/day127_132.mesenchyme.anndata.h5ad')

'''
mes_adata[mes_adata.obs['age'].isin(['49','132'])].write('./data/day132_47.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['59','132'])].write('./data/day132_59.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['72','132'])].write('./data/day132_72.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['80','132'])].write('./data/day132_80.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['101','132'])].write('./data/day132_101.mesenchyme.anndata.h5ad')
mes_adata[mes_adata.obs['age'].isin(['122','132'])].write('./data/day132_122.mesenchyme.anndata.h5ad')

sec_list = ['47_59','59_72','72_80','80_101','101_122','122_127','127_132','132_47','132_59','132_72','132_80','132_101','132_122']
'''


sec_list = ['47_59','59_72','72_80','80_101','101_122','122_127','127_132']

for comp in sec_list:
	figure_dir = './figures/' + comp + '_mesenchyme'
	sc.settings.figdir = figure_dir
	adata = sc.read_h5ad('./data/day' + comp + '.mesenchyme.anndata.h5ad')
	adata.raw = adata
	sec_analysis(adata)





















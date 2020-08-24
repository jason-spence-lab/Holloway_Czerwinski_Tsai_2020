import sys
import bbknn
import copy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

greatestPalette = ["#f44336","#265bcf","#36b3ba","#ffeb3b","#e91e63","#00cc92",
"#4caf50","#ffb65c","#9c27b0","#03a9f4","#43d61f","#ff9800","#673ab7","#cddc39",
"#81bdc5","#ff5722","#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb",
"#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff",
"#edf3ba","#ffccbd"]

palestPalette = ["#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb","#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff","#edf3ba","#ffccbd"]

## Define function to generate a colormap from rgb colors and positions 0.0-1.0
def make_cmap(colors, position=None, bit=False):
	'''
	make_cmap takes a list of tuples which contain RGB values. The RGB
	values may either be in 8-bit [0 to 255] (in which bit must be set to
	True when called) or arithmetic [0 to 1] (default). make_cmap returns
	a cmap with equally spaced colors.
	Arrange your tuples so that the first color is the lowest value for the
	colorbar and the last is the highest.
	position contains values from 0 to 1 to dictate the location of each color.
	'''
	import matplotlib as mpl
	import numpy as np
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit("position must start with 0 and end with 1")
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],
						 bit_rgb[colors[i][1]],
						 bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))

	cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap

## Helper function to make a anndata object from a sample ID
'''
def Create_Scanpy_Anndata(mistorage_mount_point, sampleID, mappingGenome):
	annotation_dict = dict()
	for line in open(''.join([mistorage_mount_point, '01_RNAseq_RAW_Data/single_cell_meta_data_table.tsv']), 'r'):
		#print(line)
		elem = str.split(line.rstrip())
		#print(elem)
		if elem:
			if elem[0] not in annotation_dict:
				annotation_dict[elem[0]] = elem[1:]
	metadata_list = annotation_dict[sampleID][1:]
	newAdata = sc.read_10x_h5(''.join([mistorage_mount_point, annotation_dict[sampleID][0]]), genome=mappingGenome)
	## Set gene names to be unique since there seem to be duplicate names from Cellranger
	newAdata.var_names_make_unique()
	## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
	print('\nAdding Metadata for sample',sampleID,'\n')
	for field in metadata_list:
		field_list = str.split(field, ':')
		meta_name = field_list[0]
		meta_value = field_list[1]
		newAdata.obs[meta_name] = meta_value
	return(newAdata)
'''

def Create_Scanpy_Anndata(mistorage_mount_point, sampleID, mappingGenome):
	annotation_dict = dict()
	for line in open('/mnt/black/scRNA-seq/single_cell_meta_data_table_cleaned.tsv', 'r'):
		#print(line)
		elem = str.split(line.rstrip())
		#print(elem)
		if elem:
			if elem[0] not in annotation_dict:
				annotation_dict[elem[0]] = elem[1:]
	metadata_list = annotation_dict[sampleID][1:]
	newAdata = sc.read_10x_h5(''.join([mistorage_mount_point, annotation_dict[sampleID][0]]))
	## Set gene names to be unique since there seem to be duplicate names from Cellranger
	newAdata.var_names_make_unique()
	## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
	print('\nAdding Metadata for sample',sampleID,'\n')
	for field in metadata_list:
		field_list = str.split(field, ':')
		meta_name = field_list[0]
		meta_value = field_list[1]
		newAdata.obs[meta_name] = meta_value
	return(newAdata)

## Helper function to get unique values from a list with duplicates 
def unique(list1): 
	# intilize a null list 
	unique_list = [] 
	# traverse for all elements 
	for x in list1: 
		# check if exists in unique_list or not 
		if x not in unique_list: 
			unique_list.append(x) 
	return(unique_list)


def Filter_New_Anndata(adata, prefix):
	sc.pp.filter_cells(adata, min_genes=1000)
	#sc.pp.filter_genes(adata, min_cells=3)
	sc.pp.filter_cells(adata, min_counts=3500)
	
	print('\nDoing initial filtering...\nSample has', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	if len(adata.obs_names) >=3200:
		sc.pp.subsample(adata, n_obs=3200, copy=False)
	
	sc.pp.downsample_counts(adata, counts_per_cell=5000, random_state=20120608, replace=False, copy=False)
	
	cell_cycle_genes = [x.strip() for x in open('/mnt/black/scRNA-seq/regev_lab_cell_cycle_genes.txt')]
	
	s_genes = cell_cycle_genes[:43]
	g2m_genes = cell_cycle_genes[43:]
	cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
	
	sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
	
	mito_genes = adata.var_names.str.startswith('MT-')
	# Calculate the percent of genes derived from mito vs genome
	# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
	adata.obs['percent_mito'] = np.sum(
		adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	# add the total counts per cell as observations-annotation to adata
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save = ''.join(['_', prefix, '_PRE_Filter_plot.pdf']), show = False)
	
	## Actually do the filtering.
	adata = adata[adata.obs['n_genes'] > 1000, :]   # Keep cells with more than 1000 genes
	adata = adata[adata.obs['n_genes'] < 9000, :]   # Keep cells with less than 5000 genes to remove most doublets
	adata = adata[adata.obs['n_counts'] < 25000, :] # Keep cells with less than 15000 UMIs to catch a few remaining doublets
	adata = adata[adata.obs['n_counts'] > 3500, :] # Keep cells with less than 15000 UMIs to catch a few remaining doublets
	#adata = adata[adata.obs['percent_mito'] < 0.15, :]   # Keep cells with less than 0.1 mito/genomic gene ratio
	
	print('\nAfter downsampling and final filtering...\nSample has', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save = ''.join(['_', prefix, '_POST_Filter_plot.pdf']), show = False)
	
	return(adata)


## Pasring function to get marker gene data similar to Seurat's marker output
def write_marker_file(adata, file_out='markers_output.csv', n_genes=50):
	print('Parsing markers...')
	marker_dict = adata.uns['rank_genes_groups']
	unique_list = [] 
	for x in adata.obs['louvain'].values.tolist(): 
		if x not in unique_list: 
			unique_list.append(str(x))
	
	outFile = open(sc.settings.figdir + file_out, 'w')
	outFile.write('logFC,gene,score,pval,padj,cluster\n')
	
	#i = 0
	
	parsed_dict = dict()
	
	for item in marker_dict:
		if type(marker_dict[item]) is not dict:
			cluster_data = []
			for subitem in marker_dict[item]:
				cluster_data.append(subitem.tolist())
			
			if str(item) not in parsed_dict:
				parsed_dict[str(item)] = cluster_data
	for cluster in range(0, len(unique_list)):
		for i in range(0, n_genes):
			line_out = []
			for marker_value in marker_dict:
				if type(marker_dict[marker_value]) is not dict:
					line_out.append(str(marker_dict[marker_value][i].tolist()[cluster]))
			line_out.append(str(cluster))
			outFile.write(','.join(line_out) + '\n')
	
	print('Saving marker data to:', file_out)
	outFile.close()

## Writes results of rank genes analysis to multiple csv files, each representing a Louvain cluster
def rank_genes(adata, groupby, clusters2_compare=None, figdir='./figures/'):
	'''
	groupby: Adata observation metadata categories to compare
	clusters2_compare: Selection of either 2 clusters to compare - if none, then do 1vAll comparison
	Need to add file clearing
	'''
	if clusters2_compare == 'All': # Does 1 to 1 comparison between of all of the clusters
		print("Functionality not developed yet")
		return 0 # Functionality not available yet
	elif clusters2_compare == None: # Do default 1vAll comparison
		print("No clusters selected for comparison. Doing default 1vAll comparison")
		sc.tl.rank_genes_groups(adata,groupby ,method='wilcoxon', rankby_abs=False, n_genes=200)
		__write_rank_genes(adata, groupby, clusters2_compare, figdir)
	else: # Compare 1v1
		adata_temp = adata[adata.obs['louvain'].isin(clusters2_compare)]
		sc.tl.rank_genes_groups(adata_temp, groupby, method='wilcoxon', n_genes=200)
		__write_rank_genes(adata_temp, groupby, clusters2_compare, figdir)
	return 0
## Actually does the writing to csv files of the rank genes analysis
def __write_rank_genes(adata, groupby, clusters2_compare, figdir='./figures/'):
	rank_genes_data = copy.deepcopy(adata.uns['rank_genes_groups']) # create copy of data to manipulate
	rank_genes_data.pop('params')
	if clusters2_compare == None:
		clusters2_compare=['all']
	for cluster in adata.obs[groupby].cat.categories:
		csv_fileName = '/'.join([figdir,'csv_files','_'.join([groupby]+clusters2_compare),
			'_'.join([cluster,'compare.csv'])])
		os.makedirs(os.path.dirname(csv_fileName), exist_ok=True) # Make file if it doesn't exist already
		with open(csv_fileName,'w',newline='') as f:
			wr = csv.writer(f)
			wr.writerow(__ele_swap(list(rank_genes_data.keys()),0,1))
			wr.writerows(zip(*__ele_swap([params[cluster] for params in rank_genes_data.values()],0,1)))

## Swaps the elements at the proposed indices in an applicable data structure
def __ele_swap(structure, index1, index2):
	structure[index1], structure[index2] = structure[index2], structure[index1]
	return structure













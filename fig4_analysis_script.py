'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu
4 June, 2019

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

Contains analysis of HIO samples for Emily Holloway

In progress ---
Moving to encapsulate parameters and relevant functions using class sca_set()
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca_run/sca_run')

from sca_run_ml import *
import csv
#from tools.pipelines import *
def run_analysis():
	figdir = './figures_081920/'
	an_run = sca_run_ml()
	#############################################################################
	## Change this to point toward your mount location for our MiStorage share ##
	#############################################################################
	an_run.storage_mount_point = 'Z:/'

	# ## List of interesting genes
	an_run.add_gene_list(markers = ['OLFM4','LGR5'],
						 label='olfm4_lgr5')

	an_run.add_gene_list(markers = ['MKI67','TOP2A'], label='mki67_top2a')

	an_run.add_gene_list(markers = ['NRG1','EGF','CDX2','SOX9','ELF3','VIL1','MUC2','CHGA','DEFA5',
									'LYZ','DPP4','MKI67','TPI1','SPDEF','MGAM','PDGFA',
									'NOTCH1','NOTCH2','NOTCH3','EPCAM',
									'NOTCH4','DLL1','DLL3','DLL4','JAG1','JAG2',
									'SMOC2','IGFBP4','BMI1','LRIG1','TERT','MSI1','PHLDA1','PROM1',
									'TNFRSF19','EPHB2','POFUT1','SOX9','ASCL2','NEUROG3'], label='basic_list')

	an_run.add_gene_list(markers = ['CHGA','CPE','SCGN','NEUROD1','PCSK1N','TUBA1A','TM4SF4','SCG2',
									'FEV','BEX1','CHGB','CRYBA2','PCSK1','CADPS','INSM1','SSR4','KCTD12',
									'GC','SMIM6','NGFRAP1','NPDC1','RGS2','C9orf16','PTMS','SCG3',
									'PAM','HEPACAM2','LY6H','HOPX','AES','EID1','SCG5','UCP2','BEX2',
									'SEC11C','H3F3B','DDX5','GABARAPL2','GNAS','MAP1B','GAD2','KIAA1244',
									'RAB3B','S100A6','QPCT','VWA5B2','PROX1','STMN2','NKX2-2','RAB26',
									'MIR7-3HG','ALDH1A1','CBX6','RFX6','GTF2I','NLRP1','STMN1','RIMBP2',
									'DNAJC12','PKM','TUSC3','DSP','MAP1LC3A','AK1','VAMP2','QDPR','C1QL1',
									'CNIH2','LRRFIP1','CACNA1A','PPDPF','SCT','APLP1','NOVA1','TTR',
									'SYP','MYL6','LINC00261','CACNA2D1','MARCKSL1','FXYD3','ABCC8',
									'MARCKS','ZKSCAN1','PCBP4','KIAA1324','PFN2','TBCB','IDS','SYT13',
									'SSTR2','GHRL','SNAP25','VTN','TUBA4A','MDK','EGR1','FOXA2','GCH1',
									'C10orf10'], label='enteroendocrine')

	an_run.add_gene_list(markers = ['S100A6','LGALS4','GSN','SH3BGRL3','TFF3','HEPACAM2','HSPB1','MDK',
									 'LYZ','PPDPF','S100A11','FXYD3','SELM','FCGBP','TMSB4X','CLCA1',
									 'KRT19','RNASE1','GUCA2A','MYL6','HPCAL1','NEURL1','TPM1','SPINK1',
									 'TPD52','HES6','GUCA2B','SPIB','SERF2','ITLN1','FRZB','KRT8',
									 'LRRC26','MUC2','BCAS1','H3F3B','KLF4','PLP2','REP15','CD9','JUN',
									 'STARD10','ST6GALNAC1','YWHAZ','QSOX1','AGR3','MT2A','IER2','CALM2',
									 'GPX2','AGR2','KRT20','TMSB10','MGLL','FOS','KLK1','CDC42EP5','DLL1',
									 'BTG2','CTSE','HMGN2','TXNIP','WFDC2','MARCKSL1','LINC00261','FOXA3',
									 'CEACAM6','UCP2','ALDH2','XBP1','PRSS1','TCEA3','EGR1','ELF3','H3F3A',
									 'CA7','DYNLL1','RBPJ','RAB11FIP1','SPATS2L','TSPAN1','TSTA3',
									 'SLC38A2','C12orf57','GMDS','FOSB','RASSF6','SSR4','CEACAM5','DYRK4',
									 'WSB1','KIAA1324','SELK','ENHO','LPCAT4','BEST4','ABCA4','FOXP1',
									 'CFTR','SYTL2'], label='secretory')

	an_run.add_gene_list(markers = ['TSPAN8','HSP90AB1','GNB2L1','AGR2','SLC12A2','ELF3','IMPDH2','NPM1',
									'JUN','PTMA','EEF1A1','PABPC1','HNRNPA1','BTF3','MYC','BEX1','TMSB10',
									'ZFP36L2','ASCL2','EEF2','HMGCS2','NACA','NAP1L1','IER2','MARCKSL1',
									'ZFAS1','TXN','HES1','CLU','EDN1','OLFM4','LGR5'], label='stem_cell')

	an_run.add_gene_list(markers = ['DCN','TCF21','COL1A1','VIM'], label='mesenchyme')

	an_run.add_gene_list(markers = ['APOC3','SEPP1','APOA4','APOA1','PRAP1','SERPINA1','B2M','CAPN3',
									'FTH1','FABP2','APOB','SLC7A7','IL32','SLC40A1','ITM2B','MUC13',
									'C19orf77','NEAT1','VAMP8','LGMN','SAT1','MAMDC4','AMN','CTSA',
									'CDHR5','RHOC','SAT2','FAM3C','GRN','C8G','ID2','CTSH','GCHFR',
									'TDP2','ASAH1','DAB2','PCSK1N','CLTB','MYL9','PHGR1','ATOX1','NEU1',
									'FBP1','ETHE1','RTN4','ACSL5','CIDEB','NPC2','HLA-C','SLC39A4',
									'FOLH1','LPGAT1','HLA-B','G0S2','CTSB','CD63','MAX','MFI2','KHK',
									'GNS','DHRS7','PEPD','RARRES1','SLC5A12','SLC46A3','MME','AGPAT2',
									'TMEM176A','SLC25A5','FN3K','FAM132A','ALDOB','NR0B2','CD74','HSD17B11',
									'BTNL3','ALPI','CST3','REEP3','SLC15A1','ADA','RBP2','MAF','GPX4',
									'TMEM92','ANPEP','SPG21','FBXO2','PLAC8','OCIAD2','TNFSF10','MXRA8',
									'COX7A2L','SLC6A19','C19orf69','NR1H4','GLRX','RBP2','FTL','SLC7A7',
									'FTH1','CTSA','AMN','MAMDC4','SLC26A3','RARRES1','S100A14','SERPINA1',
									'APOA4','NPC2','FABP2','MYL9','C19orf77','MUC13','ALDOB','SEPP1',
									'APOA1','APOC3','OAT','SAT2','GUK1','PRAP1','PEPD','B2M','ITM2B',
									'C19orf33','ASS1','FCGRT','ASAH1','CDHR5','CD68','G0S2','FAM3C',
									'SLC25A5','GCHFR','ID2','FBP1','SCT','FABP1','KHK','DHRS7','SMIM1',
									'TTR','FOLH1','APOB','MS4A8','SLC5A12','CD63','PCSK1N','RHOC','ATP1B3',
									'ESPN','ENPP7','MAF','COX7A2L','COX7C','ETHE1','TM4SF20','IL18','CIDEB',
									'MTTP','IL32','SMLR1','FBXO2','F10','CBR1','SLC9A3R1','GRN','NR0B2',
									'GOLT1A','CTSH','TDP2','ANPEP','SLC6A19','VAMP8','ANXA4','CST3',
									'HSD17B11','PRR13','TMEM92','LGALS2','ACSL5','ADA','SLC40A1','PHPT1',
									'SLC2A2','NEAT1','ATP5E','SCPEP1','SLC15A1','HADH','LGMN','TMEM176A',
									'PLAC8','ATP6V1F'], label='absoprtive')

	an_run.add_gene_list(markers = ['EGF','EGFR','ERBB2','ERBB3','ERBB4','TGFA','HBEGF','AREG','BTC',
									'EREG','NRG1','NRG2','NRG3','NRG4'], label='erbb_pathway')#'EPGN',

	an_run.add_gene_list(markers = ['EPCAM','CDX2','LGR5','OLFM4','CHGA','LYZ','MUC2',
									'VIL1','DPP4','FABP1','ALPI','SI','MGAM','SPDEF','SHH','IHH',
									'PDGFA','SOX9','AXIN2','LGR5','ASCL2','CD44',
									'LEF1','PPARD','MMP7','MSI1','LRP5',
									'LRP6','DEFA6','DEFA5','DEFA4','REG3A','REG3G'], 
									label='epi_cell_type_genes') #'CYCLIND','CMYC','BHLHA15''MUC5','TCF1','REG3B','CJUN'

	an_run.add_gene_list(markers = ['LGR5','OLFM4','FABP2','ALPL','RBP2','BEST4','SPIB','MUC2','SPDEF','DLL1','TRPM5',
									'TAS1R3','CHGA','NEUROD1','PAX6','ARX','REG4'],
									label='epithelial annotation')

	## Parameters used for initial clustering analysis
	an_run.set_analysis_params(n_neighbors = 20, # Size of the local neighborhood used for manifold approximation
							   n_pcs = 15, # Number of principle components to use in construction of neighborhood graph
							   spread = 1, # In combination with min_dist determines how clumped embedded points are
							   min_dist = 0.4, # Minimum distance between points on the umap graph
							   resolution = 0.4) # High resolution attempts to increases # of clusters identified

	## Parameters used to filter the data - Mainly used to get rid of bad cells
	an_run.set_filter_params(min_cells = 0, # Filter out cells 
							 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
							 max_genes = 3000, # Filter out cells with more genes to remove most doublets
							 max_counts = 10000, # Filter out cells with more UMIs to catch a few remaining doublets
							 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

	an_run.set_plot_params(umap_obs = ['louvain','sampleName'],
						   exp_grouping = ['louvain'],
						   size=5)#,
						   # final_quality=True)

	an_run.sample_list=['2182-5','2182-6','2250-1','2250-2','2511-2','2598-28']
	# an_run.pipe_basic(''.join([figdir,'intestine/']))#,load_save='adata_save.p')#, save_transform=True)
	
	## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
	# New analysis parameters for the subset of parameters
	analysis_params_ext = dict(n_neighbors = 20,
							n_pcs = 15,
							spread = 1,
							min_dist = 0.4,
							resolution = 0.4)

	an_run.size=20
	an_out = an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'intestine/']), label='epithelium',
				extracted=['8','17'], load_save='adata_save.p', save_transform=True)

	[an_run, highly_variable, transform_X, transform_PCA] = an_out
	adata_int = an_run.adata.copy()


	## Parameters used to filter the data - Mainly used to get rid of bad cells
	an_run.set_filter_params(min_cells = 0, # Filter out cells 
							 min_genes = 250, # Filter out cells with fewer genes to remove dead cells
							 max_genes = 8000, # Filter out cells with more genes to remove most doublets
							 max_counts = 100000, # Filter out cells with more UMIs to catch a few remaining doublets
							 max_mito = 0.2) # Filter out cells with high mitochondrial gene content

	## Parameters used for initial clustering analysis
	an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
							   n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
							   spread = 1, # In combination with min_dist determines how clumped embedded points are
							   min_dist = 0.4, # Minimum distance between points on the umap graph
							   resolution = 0.4) # High resolution attempts to increases # of clusters identified
	an_run.size=5

	analysis_dict = {'combined':['150-2','150-5','150-8'],
					 'NRG1':['150-2'],
					 'EGF':['150-5'],
					 'EGF_NRG1':['150-8']}

	# Combined samples
	for key in analysis_dict:
		an_run.sample_list = analysis_dict[key]
		an_run.pipe_basic(''.join([figdir,key,'/']), highly_variable = highly_variable)
		adata_tomap = an_run.adata.copy()

		sc.tl.ingest(adata_tomap, adata_int, obs='louvain')
		an_run.plot_sca(adata_tomap, figdir=''.join([figdir,key,'_on_int/']))
		an_run.adata = adata_tomap.copy()
		an_run.write_summary(figdir=''.join([figdir,key,'_on_int/']))

	an_run.umap_cluster_color = 'default'
	# an_run.gene_dict['new_list']['groupby_positions'] = None									

if __name__ == "__main__":
	run_analysis()
# AS-BandMEF-reprogramming
Analysis of the alternative splicing (AS) changes in B cell reprogramming, correlation of expression of RNA-binding proteins and comparison to MEF reprogramming.
Published in [Vivori et al., 2020](https://www.biorxiv.org/content/10.1101/2020.09.17.299867v1).

### 00. Run vast-tools and edgeR and apply initial filtering
- Datasets:
  - B cell reprogramming dataset: Stadhouders et al., 2018.
  - MEF reprogramming dataset: Cieply et al., 2016.
  - List of RNA-binding proteins (RBPs) and splicing-associated RBPs: generated from Uniprot ([RBP-SPL_Lists](https://github.com/cvivori/AS-BandMEF-reprogramming/tree/master/RBP-SPL_Lists)).
  
- Alternative splicing analysis ([00_Run_VTS.sh](00_Run_VTS.sh)):
  - `Vast-tools` v2.2.2 align (mm10), combine, compare. See [vast-tools webpage](https://github.com/vastgroup/vast-tools) for details.
  - Filtering of vast-tools output for reads coverage with [VTS_INCL_filtering.R](https://github.com/cvivori/useful-cluster-scripts/blob/master/README.md#vts_incl_filteringr), to consider only AS events with a minimum of 10 actual reads per sample (0 “N” values allowed for each dataset).
  - Extension of the filtered INCLUSION tables, including all the AS events differentially spliced in at least one comparison of each dataset (and their dPSI in all the comparisons). See [VTS_add_dPSI_toINCL.R](https://github.com/cvivori/useful-cluster-scripts/blob/master/README.md#vts_add_dpsi_toinclr).

- Gene Expression analysis ([00_Run_edgeR.R](00_Run_edgeR.R)):
  - Import of gene counts (from STAR mapping) 
  - Filtering for minimum 5 cpms in at least in 33% of samples (5 for B cell reprogramming, 6 for MEF reprogramming)
  - Calculation of cpm values and differential expression analysis.

### 01. Import and processing of AS and gene expression analyses outputs
- Import vast-tools outputs ([01_Import_VTS.R](01_Import_VTS.R)).
- Import edgeR outputs ([01_Import_edgeRcpm.R](01_Import_edgeRcpm.R)).

### 02. Analysis of AS events and gene expression of RNA-binding proteins
- Extract gene expression profiles from RNA-binding proteins and splicing factors ([02_Extract_RBPscpm.R](02_Extract_RBPscpm.R)).
- Count and extract different types of AS events ([02_VTS_TypeEvents.R](02_VTS_TypeEvents.R)).
- Melt dataframes for plots ([02_mdfs_toplot.R](02_mdfs_toplot.R)).

### 03. Correlation of reprogramming stages, PCA analysis and other plots
- Correlation of B cell reprogramming stages according to AS or gene expression ([03_Plot_CorrelationStages.R](03_Plot_CorrelationStages.R)).
- Calculate and plot numbers of each type of AS events (stacked bar plot, heatmap) ([03_Plot_StackedBarplot_VTSevents.R](03_Plot_StackedBarplot_VTSevents.R)).
- Calculate overlap between AS events in B cell and MEF reprogramming & plot Venn diagram ([03_Plot_VennDiagrams_BandMEFstoiPS.R](03_Plot_VennDiagrams_BandMEFstoiPS.R)).
- Principal component analysis (PCA) of B cell and MEF reprogramming on most variable genes & generate plot ([03_PCA_BandMEFstoiPS.R](03_PCA_BandMEFstoiPS.R)).
- Plot heatmaps of scaled PSI values of differential cassette exons ([03_Plot_Heatmap_VTS_CEx_BandMEFstoiPS.R](03_Plot_Heatmap_VTS_CEx_BandMEFstoiPS.R)).

- Plot heatmap of scaled expression of RNA-binding proteins involved in pluripotency (data from Han et al., 2013), reprogramming or development ([03_Plot_Heatmap_RBPs_HanDev_BandMEFstoiPS.R](03_Plot_Heatmap_RBPs_HanDev_BandMEFstoiPS.R)).
- Plot expression of genes in B cell and MEF reprogramming according to the PCA-derived 'reprogramming pseudotime' ([03_Plot_Expression_BandMEFstoiPS.R](03_Plot_Expression_BandMEFstoiPS.R)).

### 04. Clustering of AS exons in B cell reprogramming and correlation of RBPs expression
- _Mfuzz_ clustering of differentially spliced cassette exons ([04_Clustering_VTS_CEx.R](04_Clustering_VTS_CEx.R)).
- Correlation of gene expression profiles of RBPs to each AS cluster ([04_Clusters_Assign_RBPs_GEx.R](04_Clusters_Assign_RBPs_GEx.R)).

### 05. Analyses on AS clusters
- Expression of genes containing the AS exons belonging to each AS cluster ([05_Clusters_Assign_GEx_VTS.R](05_Clusters_Assign_GEx_VTS.R)).
- Prediction of ORF disruption effects for exons in AS clusters ([05_ORFdisruption_Clusters_CEx.R](05_ORFdisruption_Clusters_CEx.R)).
- Calculate overlap between AS clusters in Bcell reprogramming and MEF reprogramming ([05_Overlap_VTS_BandMEFstoiPS.R](05_Overlap_VTS_BandMEFstoiPS.R)).

# AS-reprogramming

**00. Run vast-tools and edgeR**
- B cell reprogramming dataset: Stadhouders et al., 2018
- MEFs reprogramming dataset: Cieply et al., 2016
- Generate list of RNA-binding proteins (RBPs) and splicing-associated RBPs (extended) from Uniprot

**01. Import vast-tools and edgeR outputs in R**

**02. Processing of vast-tools and edgeR tables**
- Extract gene expression profiles from RNA-binding proteins and splicing factors (Uniprot)
- Count and extract different types of AS events
- Melt dataframes for plots

**03. Plots**
- Plot heatmap with correlation of B cell reprogramming stages
- Plot types of AS events (area plot, barplot and heatmap)
- Plot heatmaps of scaled PSI values of differential cassette exons
- Plot PCA of B cell and MEFs reprogramming on most variable genes/exons
- Plot Venn Diagrams with overlaps between AS events in B cell and MEFs reprogramming
- Plot heatmap of scaled expression of RNA-binding proteins involved in pluripotency (from Han et al., 2013), reprogramming or development.

**04. Clustering**
- Mfuzz (Kumar et al., 2007) clustering of differentially spliced cassette exons
- Correlation of gene expression profiles of genes containing the exons to each AS cluster

**05. Analyses on AS clusters**
- Prediction of ORF disruption effects for exons in AS clusters
- Calculate overlap between AS events in B cell and MEFs reprogramming

require(org.Mm.eg.db)   #keytypes(org.Mm.eg.db)
require("AnnotationDbi")
require("GenomicFeatures")
require(GO.db)

### UNIPROT TABLEs DESCRIPTION ###
## RBPgenes_mouse_table.txt  (nrows= 834)
  #Uniprot: (keyword:"RNA-binding [KW-0694]" OR keyword:"mRNA splicing [KW-0508]" OR keyword:"mRNA processing [KW-0507]" OR keyword:"Spliceosome [KW-0747]") AND organism:"Mus musculus (Mouse) [10090]"
  #Columns: Uniprot_Entry  Gene_names_(primary)  Gene_names
## SPLgenes_mouse_table.txt   (nrows= 242)
  #Uniprot: (keyword:"mRNA splicing [KW-0508]" OR keyword:"Spliceosome [KW-0747]") AND organism:"Mus musculus (Mouse) [10090]"
  #Columns: Uniprot_Entry  Gene_names_(primary)  Gene_names

### Open Uniprot-derived tables
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/RBPs_lists/FINAL/")
# RBPgenes_mouse_table=read.table("RBPgenes_mouse_table.txt",sep = "\t",quote="",header = TRUE)
# SPLgenes_mouse_table=read.table("SPLgenes_mouse_table.txt",sep = "\t",quote="",header = TRUE)
RBPgenes_mouse_table <- read.table("Uniprot_RBPGenes.tab",sep="\t",quote="",header = T,fill=T)
SPLgenes_mouse_table <- read.table("Uniprot_SplicingGenes.tab",sep="\t",quote="",header = T,fill=T)


### To generate a VECTOR with ALL possible RBP/SPL GENE NAMES
RBPgenes=as.vector(RBPgenes_mouse_table[,"Gene.names"])
SPLgenes=as.vector(SPLgenes_mouse_table[,"Gene.names"])
RBPgenes_mouse_vector=unlist(strsplit(RBPgenes," "))
SPLgenes_mouse_vector=unlist(strsplit(SPLgenes," "))
# write.csv(RBPgenes_mouse_vector,file = "RBPgenes_mouse_vector.txt")
# write.csv(SPLgenes_mouse_vector,file = "SPLgenes_mouse_vector.txt") #Use _ext!

### To expand list of SPLs (take all possible names of the gene from RBPgenes_mouse_table!)
SPL_ext=c("Rbm10","Kiaa0122","Snrnp70","Snrp70","Ddx18","mCG_1040626","Ewsr1","Ews","Ewsh","Matr3","Cpsf2","Cpsf100","Mcpsf","Dek","Tia1")
SPLgenes_mouse_ext=c(SPLgenes_mouse_vector,SPL_ext)
# write.csv(SPLgenes_mouse_ext,file = "SPLgenes_mouse_ext.txt")

#### CONVERT GENE NAMES
mRBP=mapIds(org.Mm.eg.db, keys=(RBPgenes_mouse_vector), column='SYMBOL', keytype='ALIAS',multiVals="first")		
mSPL=mapIds(org.Mm.eg.db, keys=(SPLgenes_mouse_ext), column='SYMBOL', keytype='ALIAS',multiVals="first")		# SPL is a vector containing Splicing Factors Gene Names


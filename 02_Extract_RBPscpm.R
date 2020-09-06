require(tidyverse)
require(org.Mm.eg.db)

#MOUSE vectors with ALL possible RBP/SPL GENE NAMES  from Uniprot database
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/RBPs_lists/")
RBPs=read.table("RBPgenes_mouse_vector.txt",sep=",",header = TRUE,row.names = 1) 
mRBP=mapIds(org.Mm.eg.db, keys=as.vector(RBPs[,1]), column='SYMBOL', keytype='ALIAS',multiVals="first")		# RBP is a vector containing RBP Gene Names 
SPLs=read.table("SPLgenes_mouse_ext.txt",sep=",",header = TRUE,row.names = 1)
SPLs <- c(as.vector(SPLs$x), "Rbm47","Cpsf2","Cpsf3","Tia1","Zrsr1","Zrsr2","U2surp")
mSPL=mapIds(org.Mm.eg.db, keys=as.vector(SPLs), column='SYMBOL', keytype='ALIAS',multiVals="first")		# SPL is a vector containing Splicing Factors Gene Names


## EXTRACT ENSGIDs from RBPs/SPL
ENS2Gene_RBPs <- ENS2Gene %>%
  filter(GeneName %in% mRBP)
length(mRBP); dim(ENS2Gene_RBPs)
ENS2Gene_SPLs <- ENS2Gene %>%
  filter(GeneName %in% mSPL)
length(mSPL); dim(ENS2Gene_SPLs)


## EXTRACT CPMs of RBPs
RBPs_B_CPMs_av_fil <- B_CPMs_av_fil[intersect(ENS2Gene_RBPs$ENSGID,rownames(B_CPMs_av_fil)),]
RBPs_M_CPMs_av_fil <- M_CPMs_av_fil[intersect(ENS2Gene_RBPs$ENSGID,rownames(M_CPMs_av_fil)),]
RBPs_B_CPMs_av_fil_scaled <- t(scale(t(RBPs_B_CPMs_av_fil)))
RBPs_M_CPMs_av_fil_scaled <- t(scale(t(RBPs_M_CPMs_av_fil)))
RBPs_BM_CPMs_av_fil <- BM_CPMs_av_fil[intersect(ENS2Gene_RBPs$ENSGID,rownames(BM_CPMs_av_fil)),]
RBPs_BM_CPMs_av_fil_scaled <- t(scale(t(RBPs_BM_CPMs_av_fil)))

head(RBPs_B_CPMs_av_fil)
dim(RBPs_B_CPMs_av_fil); dim(RBPs_M_CPMs_av_fil)
dim(RBPs_BM_CPMs_av_fil); dim(RBPs_BM_CPMs_av_fil_scaled)


## EXTRACT CPMs of SPLs
SPLs_B_CPMs_av_fil <- B_CPMs_av_fil[intersect(ENS2Gene_SPLs$ENSGID,rownames(B_CPMs_av_fil)),]
SPLs_M_CPMs_av_fil <- M_CPMs_av_fil[intersect(ENS2Gene_SPLs$ENSGID,rownames(M_CPMs_av_fil)),]
SPLs_B_CPMs_av_fil_scaled <- t(scale(t(SPLs_B_CPMs_av_fil)))
SPLs_M_CPMs_av_fil_scaled <- t(scale(t(SPLs_M_CPMs_av_fil)))
SPLs_BM_CPMs_av_fil <- BM_CPMs_av_fil[intersect(ENS2Gene_SPLs$ENSGID,rownames(BM_CPMs_av_fil)),]
SPLs_BM_CPMs_av_fil_scaled <- t(scale(t(SPLs_BM_CPMs_av_fil)))

dim(SPLs_B_CPMs_av_fil); dim(SPLs_M_CPMs_av_fil)
dim(SPLs_BM_CPMs_av_fil); dim(SPLs_BM_CPMs_av_fil_scaled)


#----- Outputs
message("Outputs for B and MEFs:
        RBPs_ and & SPLs_ expression (CPMs) dataframes")
#



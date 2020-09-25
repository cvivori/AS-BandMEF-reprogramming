require(edgeR)
require(org.Mm.eg.db)

## Mapping to mm10 genome was performed with STAR v2.7.1a with the following parameters:
# --limitGenomeGenerateRAM 31000000000 --limitBAMsortRAM 31000000000 --runThreadN 10 --readFilesCommand zcat 
# --outSJfilterReads Unique --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 6 --alignSJDBoverhangMin 3 
# --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated 
# --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --seedSearchStartLmax 50 --quantMode GeneCounts

## READ COUNT FILES FROM STAR
setwd("~/CEBPa_NEW/EdgeR_GEx/")
filesB <- list.files (path= "BtoiPS/", pattern = "_*ReadsPerGene.out.tab")
filesM <- list.files (path= "MEFstoiPS/", pattern = "_*ReadsPerGene.out.tab")
  datasetB <- sapply(filesB, function(x) unlist(strsplit(x,split="_"))[1])
  datasetM <- sapply(filesM, function(x) unlist(strsplit(x,split="_"))[1])
  samplesB <- sapply(filesB, function(x) unlist(strsplit(x,split="_"))[2])
  samplesM <- sapply(filesM, function(x) unlist(strsplit(x,split="_"))[2])
  replicatesB <- sapply(filesB, function(x) unlist(strsplit(x,split="_"))[3])
  replicatesB <- gsub(pattern = "rep", replacement = "", x = replicatesB)
  replicatesM <- sapply(filesM, function(x) unlist(strsplit(x,split="_"))[3])
  replicatesM <- gsub(pattern = "rep", replacement = "", x = replicatesM)
labelsB <- paste(datasetB,samplesB,replicatesB, sep="_")
labelsM <- paste(datasetM,samplesM,replicatesM, sep="_")


## CREATE DGELIST OBJECT FROM MANY FILES CONTAINING COUNTS OF SAMPLES
dgListB <- readDGE(skip=4, ## skip first 4 lines cause they contain the summary
                  files = filesB, path = "~/CEBPa_NEW/EdgeR_GEx/BtoiPS/", 
                  columns=c(1,2), ## the 2nd column contains the unstranded reads!!!!!!!  
                  group = samplesB, labels = labelsB)
dgListB
head(dgListB$counts)
summary(dgListB$counts)
#total number of readcounts per sample       colSums(dgListB$counts)

dgListM <- readDGE(skip=4, ## skip first 4 lines cause they contain the summary
                   files = filesM, path = "~/CEBPa_NEW/EdgeR_GEx/MEFstoiPS/", 
                   columns=c(1,2), ## the 2nd column contains the unstranded reads!!!!!!!  
                   group = samplesM, labels = labelsM)
dgListM
head(dgListM$counts)
summary(dgListM$counts)

## CONVERT ENSGIDs in human-readable GENE NAMES
IDsB <- rownames(dgListB$counts); IDsM <- rownames(dgListM$counts)
ENSGIDsB <- sapply(IDsB, function(x) unlist(strsplit(x,split="[.]"))[1]); ENSGIDsM <- sapply(IDsM, function(x) unlist(strsplit(x,split="[.]"))[1])
SymbolB <- mapIds(org.Mm.eg.db, keys=ENSGIDsB, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
SymbolM <- mapIds(org.Mm.eg.db, keys=ENSGIDsM, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")

## ADD IT TO THE LIST
dgListB$genes <- data.frame (Symbol = SymbolB)
dgListM$genes <- data.frame (Symbol = SymbolM)
dgListB


# SAVE ORIGINAL dgList
# dgListB_orig <- dgListB
# dgListB_orig <- calcNormFactors(dgListB_orig)
# dgListM_orig <- dgListM
# dgListM_orig <- calcNormFactors(dgListM_orig)


#### SET FILTERS HERE
threshold_cpm <- 5
threshold_nsamplesB <- 5    # 5 cpms at least in 33% of samples? 5 for BtoiPS 
threshold_nsamplesM <- 6    # 5 cpms at least in 33% of samples? 6 for MEFstoiPS

## Filtering for minimal gene expression
keepB <- rowSums(cpm(dgListB)>threshold_cpm) >= threshold_nsamplesB    
keepM <- rowSums(cpm(dgListM)>threshold_cpm) >= threshold_nsamplesM    
dgListB <- dgListB[keepB,]
dgListM <- dgListM[keepM,]
dim(dgListB_orig); dim(dgListB)
dim(dgListM_orig); dim(dgListM)
## Recalculate the library size accordingly
dgListB$samples$lib.size <- colSums(dgListB$counts)
dgListM$samples$lib.size <- colSums(dgListM$counts)



## Calculate Normalization Factor
dgListB <- calcNormFactors(dgListB); #ps_dgListB <- calcNormFactors(ps_dgListB)
dgListB$samples
dgListM <- calcNormFactors(dgListM); #ps_dgListM <- calcNormFactors(ps_dgListM)
dgListM$samples


## MD PLOT
filter_titleB = paste("Filter: cpm >", threshold_cpm, "in", threshold_nsamplesB, "samples", sep=" ")
  for (x in c(1:ncol(dgListB$counts))) {
    sample_name = colnames(dgListB$counts)[x]
    pdf(file = paste("MD_cpm",threshold_cpm,"_nsamples",threshold_nsamplesB,"_",sample_name,".pdf", sep=""), width =10, height=7)
      par(mfrow=c(1,2))
      plotMD(cpm(dgListB_orig, log=T),column=x, main = "Before filtering"); abline(h=0,col="red",lty=2,lwd=1)
      plotMD(cpm(dgListB, log=T),column=x, main = filter_titleB ); abline(h=0,col="red",lty=2,lwd=1)
      title(sample_name, outer = TRUE, cex = 1.5, line = -1)
    dev.off()
  }

  filter_titleM = paste("Filter: cpm >", threshold_cpm, "in", threshold_nsamplesM, "samples", sep=" ")
  for (x in c(1:ncol(dgListM$counts))) {
    sample_name = colnames(dgListM$counts)[x]
    pdf(file = paste("MD_cpm",threshold_cpm,"_nsamples",threshold_nsamplesM,"_",sample_name,".pdf", sep=""), width =10, height=7)
    par(mfrow=c(1,2))
    plotMD(cpm(dgListM_orig, log=T),column=x, main = "Before filtering"); abline(h=0,col="red",lty=2,lwd=1)
    plotMD(cpm(dgListM, log=T),column=x, main = filter_titleM ); abline(h=0,col="red",lty=2,lwd=1)
    title(sample_name, outer = TRUE, cex = 1.5, line = -1)
    dev.off()
  }


### MDS PLOTS
  pdf(file = paste("MDS_cpm",threshold_cpm,"_nsamples",threshold_nsamplesB,"_BtoiPS.pdf", sep=""),width=10,height=7)
  par(mfrow=c(1,2))
    plotMDS(dgListB_orig, method="bcv",col=replicatesB, main = "Before filtering")
    plotMDS(dgListB, method="bcv",col=replicatesB, main = filter_titleB )
    title("BtoiPS", outer = TRUE, cex = 1.5, line = -1)
  dev.off()
  pdf(file = paste("MDS_cpm",threshold_cpm,"_nsamples",threshold_nsamplesM,"_MEFstoiPS.pdf", sep=""),width=10,height=7)
  par(mfrow=c(1,2))
    plotMDS(dgListM_orig, method="bcv",col=replicatesM, main = "Before filtering")
    plotMDS(dgListM, method="bcv",col=replicatesM, main = filter_titleM )
    title("MEFstoiPS", outer = TRUE, cex = 1.5, line = -1)
  dev.off()


## OUTPUT TABLES WITH CPMs 
cpm_B <- cpm(dgListB)
ENSGID_B <- sapply(rownames(cpm_B), function(x) unlist(strsplit(x,split="[.]"))[1])
GeneName_B <- mapIds(org.Mm.eg.db, keys=ENSGID_B, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
df_B <- data.frame(ENSGID=ENSGID_B, GeneName=GeneName_B, cpm_B)
head(df_B)   
write.table(x= df_B,file = "FilteredGeneCpm_BtoiPS.txt", row.names = F, sep="\t",quote =F)

cpm_M <- cpm(dgListM)
ENSGID_M <- sapply(rownames(cpm_M), function(x) unlist(strsplit(x,split="[.]"))[1])
GeneName_M <- mapIds(org.Mm.eg.db, keys=ENSGID_M, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
df_M <- data.frame(ENSGID=ENSGID_M,GeneName= GeneName_M, cpm_M)
head(df_M)   
write.table(x= df_M,file = "FilteredGeneCpm_MEFstoiPS.txt", row.names = F, sep="\t",quote =F)




### ESTIMATE DISPERSION ACCORDING TO EDGER
d1_B <- estimateCommonDisp(dgListB, verbose = T)
d1_M <- estimateCommonDisp(dgListM, verbose = T)
d1_B <- estimateTagwiseDisp(d1_B)
d1_M <- estimateTagwiseDisp(d1_M)

plotBCV(d1_B)
plotBCV(d1_M)

design.matB <- model.matrix(~ 0+ samplesB) 
rownames(design.matB) <- rownames(dgListB$samples)
d2_B <- estimateGLMCommonDisp(dgListB,design.matB)
d2_B <- estimateGLMTrendedDisp(d2_B,design.matB, method="auto")
d2_B <- estimateGLMTagwiseDisp(d2_B,design.matB)
plotBCV(d2_B)

design.matM <- model.matrix(~  0+ samplesM) 
rownames(design.matM) <- rownames(dgListM$samples)
d2_M <- estimateGLMCommonDisp(dgListM,design.matM)
d2_M <- estimateGLMTrendedDisp(d2_M,design.matM, method="auto")
d2_M <- estimateGLMTagwiseDisp(d2_M,design.matM)
plotBCV(d2_M)


## CALCULATE DIFFERENTIALLY EXPRESSED GENES between B cells and iPS
et1B_iPS <- exactTest(d1_B, pair = c("B","iPS"))
topTags(et1B_iPS)
de1 <- decideTestsDGE(et1B_iPS, adjust.method="BH", p.value=0.001)
summary(de1)

de1tagsB_iPS <- rownames(d1_B)[as.logical(de1)] 
plotSmear(et1B_iPS, de.tags=de1tagsB_iPS)
abline(h = c(-2, 2), col = "blue")


et2B_iPS <- exactTest(d2_B, pair = c("B","iPS"))
topTags(et2B_iPS)
de2 <- decideTestsDGE(et2B_iPS, adjust.method="BH", p.value=0.001)
summary(de2)

de2tagsB_iPS <- rownames(d2_B)[as.logical(de1)] 
plotSmear(et2B_iPS, de.tags=de1tagsB_iPS)
abline(h = c(-2, 2), col = "blue")

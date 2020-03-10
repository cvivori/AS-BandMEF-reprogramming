require(tidyverse)
require(org.Mm.eg.db)
require(factoextra)
#require(princurve)
#require(corrplot)
require(RColorBrewer)
#require(scales)

##########################################
## PCA ON MOST VARIABLE GENES (EXPRESSION)

## CALCULATE COEFFVAR IN BOTH DATASETS
BM_CPMs_av_tofil <- na.omit(BM_CPMs_av)
BM_CPMs_av_tofil$CoeffVar <- (apply (BM_CPMs_av_tofil, 1, function(x) sd(x)/mean(x)))
head(BM_CPMs_av_tofil)
summary(BM_CPMs_av_tofil$CoeffVar)

#3rd quartile
mostvar <- rownames(subset(BM_CPMs_av_tofil, CoeffVar > 0.66))
# mostvar <- top_n(rownames_to_column(BM_CPMs_av_tofil,var="ENSGID"), n=500,wt = CoeffVar)$ENSGID
# mostvar <- rownames(subset(BM_CPMs_av_tofil, CoeffVar > 0))

##########################################
## EXTRACT SCALED CPMs of MOST VARIABLE GENES
vmat=subset(BM_CPMs_av_scaled, rownames(BM_CPMs_av_scaled) %in% mostvar)
vmat=na.omit(vmat)
dim(vmat)
head(vmat)

## PCA on SAMPLES
pca=prcomp((vmat),scale. = F,center=F)
plot(pca, type="l")
summary(pca)
#head(pca$rotation)

 
## PLOT PCA ON SCALED CPMs (COEFFVAR > 3rd Quartile)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_CoeffVar_066/")
col=rep("maroon",nrow(pca$rotation))
col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
pdf("PCA_ScaledCPMs_Samples_CoeffVar_066.pdf")
par(pty="s")
plot(pca$rotation[,1:2], col=col, asp=0.5)
title(main = paste("PCA on samples with most variable GEx profiles",sep = ""), cex.main=0.9)
title(xlab = "41.42% of variance explained",
      ylab = "18.06% of variance explained", cex.main=0.8, line=4)
text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
dev.off()

# ## PLOT PCA ON SCALED CPMs (top 500 COEFFVAR)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_Top500_CoeffVar/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_ScaledCPMs_Samples_Top500_CoeffVar.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with most variable GEx profiles",sep = ""), cex.main=0.9)
# title(xlab = "49.24% of variance explained",
#       ylab = "20.30% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()

## PLOT PCA ON SCALED CPMs (No Filter)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_CoeffVar_066/")
col=rep("maroon",nrow(pca$rotation))
col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
pdf("PCA_ScaledCPMs_Samples_NoFilter.pdf")
par(pty="s")
plot(pca$rotation[,1:2], col=col, asp=0.5)
title(main = paste("PCA on samples with all GEx profiles",sep = ""), cex.main=0.9)
title(xlab = "37.86% of variance explained",
      ylab = "14.99% of variance explained", cex.main=0.8, line=4)
text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
dev.off()


# Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
res.var$cor           # Correlations between variables and dimensions 
# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

## PERFORM K MEANS CLUSTERING on SCALED CPMs
k <- kmeans(res.var$coord, centers = 4, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-1.8,0.8),ylim=c(-0.8,0.8))
title(main = paste("PCA on samples with most variable GEx profiles - kmeans",sep = ""), cex.main=0.9)
title(xlab = "41.42% of variance explained",
      ylab = "18.06% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# NoFilter
k <- kmeans(res.var$coord, centers = 4, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-1.2,1.5),ylim=c(-0.6,0.7))
title(main = paste("PCA on samples with all GEx profiles",sep = ""), cex.main=0.9)
title(xlab = "37.86% of variance explained",
      ylab = "14.99% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)


##########################################
## EXTRACT RAW CPMs of MOST VARIABLE GENES
vmat=subset(BM_CPMs_av, rownames(BM_CPMs_av) %in% mostvar)
vmat=na.omit(vmat)
dim(vmat)
head(vmat)


## PCA on SAMPLES
pca=prcomp((vmat),scale. = F,center=F)
plot(pca, type="l")
summary(pca)
#head(pca$rotation)

## PLOT PCA ON RAW CPMs (COEFFVAR > 3rd Quartile)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_CoeffVar_066/")
col=rep("maroon",nrow(pca$rotation))
col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
pdf("PCA_RawCPMs_Samples_CoeffVar_066.pdf")
par(pty="s")
plot(pca$rotation[,1:2], col=col, asp=0.5)
title(main = paste("PCA on samples with most variable GEx profiles",sep = ""), cex.main=0.9)
title(xlab = "78.93% of variance explained",
      ylab = "08.615% of variance explained", cex.main=0.8, line=4)
text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
dev.off()

# ## PLOT PCA ON RAW CPMs (top 500 COEFFVAR)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_Top500_CoeffVar/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_RawCPMs_Samples_Top500_CoeffVar.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with most variable GEx profiles",sep = ""), cex.main=0.9)
# title(xlab = "89.13% of variance explained",
#       ylab = "05.17% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()

# ## PLOT PCA ON RAW CPMs (NoFilter)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/PCA_analysis/Most_variable_genes_CoeffVar_066/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_RawCPMs_Samples_NoFilter.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with all GEx profiles",sep = ""), cex.main=0.9)
# title(xlab = "80.24% of variance explained",
#       ylab = "13.38% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()


# Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
res.var$cor           # Correlations between variables and dimensions 
# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

## PERFORM K MEANS CLUSTERING on RAW CPMs
k <- kmeans(res.var$coord, centers = 5, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-650,50),ylim=c(-100,150))
title(main = paste("PCA on samples with most variable GEx profiles - kmeans",sep = ""), cex.main=0.9)
title(xlab = "78.93% of variance explained",
      ylab = "08.615% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# NoFilter
k <- kmeans(res.var$coord, centers = 4, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-350,-200),ylim=c(-100,270))
title(main = paste("PCA on samples with all GEx profiles",sep = ""), cex.main=0.9)
title(xlab = "80.24% of variance explained",
      ylab = "13.38% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)









##############################
### PCA ON MOST VARIABLE EXONS
## CALCULATE COEFFVAR IN BOTH DATASETS
BM_PSIs_VTS10_av_CEx_tofil <- BM_PSIs_VTS10_av_CEx
BM_PSIs_VTS10_av_CEx_tofil$CoeffVar <- (apply (BM_PSIs_VTS10_av_CEx_tofil, 1, function(x) sd(x)/mean(x)))
head(BM_PSIs_VTS10_av_CEx_tofil)
summary(BM_PSIs_VTS10_av_CEx_tofil$CoeffVar)

#3rd quartile
mostvar <- rownames(subset(BM_PSIs_VTS10_av_CEx_tofil, CoeffVar > 0.46525))
# mostvar <- top_n(rownames_to_column(BM_PSIs_VTS10_av_CEx_tofil,var="EVENT"), n=500,wt = CoeffVar)$EVENT
# mostvar <- rownames(subset(BM_PSIs_VTS10_av_CEx_tofil, CoeffVar > 0))

###############################################
## EXTRACT SCALED VALUES of MOST VARIABLE EXONS
vmat=subset(BM_PSIs_VTS10_av_CEx_scaled, rownames(BM_PSIs_VTS10_av_CEx_scaled) %in% mostvar)
vmat=na.omit(vmat)
dim(vmat)
head(vmat)

## PCA on SAMPLES
pca=prcomp((vmat),scale. = F,center=F)
plot(pca, type="l")
summary(pca)
#head(pca$rotation)


## PLOT PCA ON SCALED PSIs (COEFFVAR > 3rd Quartile)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
col=rep("maroon",nrow(pca$rotation))
col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
pdf("PCA_ScaledPSIs_Samples_CoeffVar_046.pdf")
par(pty="s")
plot(pca$rotation[,1:2], col=col, asp=0.5)
title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "45.35% of variance explained",
      ylab = "13.81% of variance explained", cex.main=0.8, line=4)
text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
dev.off()

## PLOT PCA ON SCALED PSIs (top 500 COEFFVAR)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_ScaledPSIs_Samples_Top500_CoeffVar.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
# title(xlab = "36.11% of variance explained",
#       ylab = "15.20% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()

# ## PLOT PCA ON SCALED PSIs (No Filter)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_ScaledPSIs_Samples_NoFilter.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with all PSIs profiles",sep = ""), cex.main=0.9)
# title(xlab = "33.93% of variance explained",
#       ylab = "14.67% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()
 
# Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
res.var$cor           # Correlations between variables and dimensions 
# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

## PERFORM K MEANS CLUSTERING on SCALED PSI VALUES
k <- kmeans(res.var$coord, centers = 4, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-1.5,1.5),ylim=c(-1.5,0.8))
title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "45.35% of variance explained",
      ylab = "13.81% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
#Nofiltering
k <- kmeans(res.var$coord, centers = 6, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
#plot(res.var$coord, col=k$clust, pch=16)
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-1,1),ylim=c(-1.2,0.8))
title(main = paste("PCA on samples with all PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "33.93% of variance explained",
      ylab = "14.67% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)



###############################################
## EXTRACT RAW PSI VALUES of MOST VARIABLE EXONS
vmat=subset(BM_PSIs_VTS10_av_CEx, rownames(BM_PSIs_VTS10_av_CEx) %in% mostvar)
vmat=na.omit(vmat)
dim(vmat)
head(vmat)

## PCA on SAMPLES
pca=prcomp((vmat),scale. = F,center=F)
plot(pca, type="l")
summary(pca)
#head(pca$rotation)


## PLOT PCA ON RAW PSIs (COEFFVAR > 3rd Quartile)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
col=rep("maroon",nrow(pca$rotation))
col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
pdf("PCA_RawPSIs_Samples_CoeffVar_046.pdf")
par(pty="s")
plot(pca$rotation[,1:2], col=col, asp=0.5, ylim=c(-0.4,0.4))
title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "73.97% of variance explained",
      ylab = "15.53% of variance explained", cex.main=0.8, line=4)
text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
dev.off()

## PLOT PCA ON RAW PSIs (top 500 COEFFVAR)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_RawPSIs_Samples_Top500_CoeffVar.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
# title(xlab = "93.25% of variance explained",
#       ylab = "02.37% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()

# ## PLOT PCA ON RAW PSIs (NoFilter)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/PCA analysis/")
# col=rep("maroon",nrow(pca$rotation))
# col[ grep("B_", rownames(pca$rotation)) ] = "darkslategray"
# pdf("PCA_RawPSIs_Samples_NoFilter.pdf")
# par(pty="s")
# plot(pca$rotation[,1:2], col=col, asp=0.5)
# title(main = paste("PCA on samples with all PSIs profiles",sep = ""), cex.main=0.9)
# title(xlab = "96.43% of variance explained",
#       ylab = "01.689% of variance explained", cex.main=0.8, line=4)
# text(pca$rotation[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
# dev.off()


# Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
res.var$cor           # Correlations between variables and dimensions 
# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

## PERFORM K MEANS CLUSTERING on SCALED RAW VALUES
k <- kmeans(res.var$coord, centers = 5, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-35,-17),ylim=c(-20,20))
title(main = paste("PCA on samples with most variable PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "73.97% of variance explained",
      ylab = "15.53% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)
#NoFilter
k <- kmeans(res.var$coord, centers = 4, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
#plot(res.var$coord, col=k$clust, pch=16)
plot(res.var$coord, col=k$clust, pch=16,xlim=c(-63,-56),ylim=c(-15,15))
title(main = paste("PCA on samples with all PSIs profiles",sep = ""), cex.main=0.9)
title(xlab = "96.43% of variance explained",
      ylab = "01.689% of variance explained", cex.main=0.8, line=4)
text(res.var$coord[,1:2], labels=rownames(pca$rotation),pos=3,cex=0.7,col = col)


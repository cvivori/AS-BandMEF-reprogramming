require(tidyverse)
require(org.Mm.eg.db)
require(reshape2)


setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/")

## MAKE MELTED TABLES TO PLOT GEx VALUES
# melted_BM_EdgeR=melt(merge(B_EdgeR,M_EdgeR,by=c("ENSGID","GeneName")), id.vars = c("ENSGID","GeneName"), variable.name = "Sample",value.name = "cpm")
# melted_BM_EdgeR$Dataset=sapply(as.character(melted_BM_EdgeR$Sample), function(x) unlist(strsplit(x,split="_"))[1])
# melted_BM_EdgeR$Condition=sapply(as.character(melted_BM_EdgeR$Sample), function(x) paste(unlist(strsplit(x,split="_"))[1],unlist(strsplit(x,split="_"))[2],sep = "_"))
# melted_BM_EdgeR$Replicate=sapply(as.character(melted_BM_EdgeR$Sample), function(x) unlist(strsplit(x,split="_"))[3])
# head(melted_BM_EdgeR)
# unique(melted_BM_EdgeR$Replicate)
# write.table(melted_BM_EdgeR,file="melted_BM_EdgeR.txt",row.names = F, quote=F,sep = "\t")
melted_BM_EdgeR <- read.table(file="melted_BM_EdgeR.txt",header = T,sep = "\t")

# melted_BM_EdgeR_av=aggregate(cpm ~ ENSGID+GeneName+Dataset+Condition, melted_BM_EdgeR, mean)
# melted_BM_EdgeR_sd=aggregate(cpm ~ ENSGID+GeneName+Dataset+Condition, melted_BM_EdgeR, sd)
# melted_BM_EdgeR_avsd=merge(x=melted_BM_EdgeR_av,y=melted_BM_EdgeR_sd,by=c("ENSGID","GeneName","Dataset","Condition"),suffixes=c("_av","_sd"))
# melted_BM_EdgeR_avsd$GeneName_Dataset=apply(melted_BM_EdgeR_avsd, 1, function(x) paste(x["GeneName"], x["Dataset"],sep="_"))
# head(melted_BM_EdgeR_avsd)
# write.table(melted_BM_EdgeR_avsd,file="melted_BM_EdgeR_avsd.txt",row.names = F, quote=F,sep = "\t")
melted_BM_EdgeR_avsd <- read.table(file="melted_BM_EdgeR_avsd.txt", header=T,sep = "\t")

# melted_BM_EdgeR_avsd_srep <- merge(x=melted_BM_EdgeR_avsd,y=melted_BM_EdgeR,by=c("ENSGID","GeneName","Dataset","Condition"))
# head(melted_BM_EdgeR_avsd_srep)
# write.table(melted_BM_EdgeR_avsd_srep,file="melted_BM_EdgeR_avsd_srep.txt",row.names = F, quote=F,sep = "\t")
melted_BM_EdgeR_avsd_srep <- read.table(file="melted_BM_EdgeR_avsd_srep.txt", header=T,sep = "\t")


## MAKE MELTED TABLES TO PLOT SCALED GEx VALUES
# melted_BM_EdgeR_scaled=melt(merge(x=ENS2Gene,y=BM_CPMs_av_scaled,by.x="ENSGID",by.y="row.names"), id.vars = c("ENSGID","GeneName","Gene"), variable.name = "Sample",value.name = "Scaled_cpm")
# melted_BM_EdgeR_scaled$Dataset=sapply(as.character(melted_BM_EdgeR_scaled$Sample), function(x) unlist(strsplit(x,split="_"))[1])
# melted_BM_EdgeR_scaled$Condition=sapply(as.character(melted_BM_EdgeR_scaled$Sample), function(x) unlist(strsplit(x,split="_"))[2])
# melted_BM_EdgeR_scaled$GeneName_Dataset=apply(melted_BM_EdgeR_scaled, 1, function(x) paste(x["GeneName"], x["Dataset"],sep="_"))
# head(melted_BM_EdgeR_scaled)
# unique(melted_BM_EdgeR_scaled$Condition)
# write.table(melted_BM_EdgeR_scaled,file="melted_BM_EdgeR_scaled.txt",row.names = F, quote=F,sep = "\t")
melted_BM_EdgeR_scaled <- read.table(file="melted_BM_EdgeR_scaled.txt", header=T,sep = "\t")



### BIG TABLE WITH PSIs
# melted_BM_PSIs=melt(merge(B_INCL_dPSI10[,c("EVENT","GENE",PSIcols_B)],M_INCL_dPSI10[,c("EVENT","GENE",PSIcols_M)],by=c("EVENT","GENE")), id.vars = c("EVENT","GENE"), variable.name = "Sample",value.name = "PSI")
# melted_BM_PSIs$Dataset=sapply(as.character(melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[1])
# melted_BM_PSIs$Condition=sapply(as.character(melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[2])
# melted_BM_PSIs$Condition <- paste(melted_BM_PSIs$Dataset,melted_BM_PSIs$Condition,sep = "_")
# melted_BM_PSIs$Replicate=sapply(as.character(melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[3])
# head(melted_BM_PSIs)
# write.table(melted_BM_PSIs,file="melted_BM_PSIs.txt",row.names = F, quote=F,sep = "\t")
melted_BM_PSIs <- read.table(file="melted_BM_PSIs.txt",header=T,sep = "\t")

# scB <- B_INCL_dPSI10[,c("EVENT","GENE",PSIcols_B)];   scM <- M_INCL_dPSI10[,c("EVENT","GENE",PSIcols_M)]
# scB[,PSIcols_B] <- t(scale(t(scB[,PSIcols_B])));   scM[,PSIcols_M] <- t(scale(t(scM[,PSIcols_M])))
# melted_BM_PSIs_scaled <- melt(merge(scB,scM,by=c("EVENT","GENE")), id.vars = c("EVENT","GENE"), variable.name = "Sample",value.name = "PSI")
# melted_BM_PSIs_scaled$Dataset=sapply(as.character(melted_BM_PSIs_scaled$Sample), function(x) unlist(strsplit(x,split="_"))[1])
# melted_BM_PSIs_scaled$Condition=sapply(as.character(melted_BM_PSIs_scaled$Sample), function(x) unlist(strsplit(x,split="_"))[2])
# melted_BM_PSIs_scaled$Condition <- paste(melted_BM_PSIs_scaled$Dataset,melted_BM_PSIs_scaled$Condition,sep = "_")
# melted_BM_PSIs_scaled$Replicate=sapply(as.character(melted_BM_PSIs_scaled$Sample), function(x) unlist(strsplit(x,split="_"))[3])
# head(melted_BM_PSIs_scaled)
# write.table(melted_BM_PSIs_scaled,file="melted_BM_PSIs_scaled.txt",row.names = F, quote=F,sep = "\t")
melted_BM_PSIs_scaled <- read.table(file="melted_BM_PSIs_scaled.txt",header=T,sep = "\t")


# melted_BM_PSIs_av=aggregate(PSI ~ EVENT+GENE+Dataset+Condition, melted_BM_PSIs, mean)
# melted_BM_PSIs_sd=aggregate(PSI ~ EVENT+GENE+Dataset+Condition, melted_BM_PSIs, sd)
# melted_BM_PSIs_avsd=merge(x=melted_BM_PSIs_av,y=melted_BM_PSIs_sd,by=c("EVENT","GENE","Dataset","Condition"),suffixes=c("_av","_sd"))
# melted_BM_PSIs_avsd$EVENT_Dataset=apply(melted_BM_PSIs_avsd, 1, function(x) paste(x["EVENT"], x["Dataset"],sep="_"))
# head(melted_BM_PSIs_avsd)
# write.table(melted_BM_PSIs_avsd,file="melted_BM_PSIs_avsd.txt",row.names = F, quote=F,sep = "\t")
melted_BM_PSIs_avsd <- read.table(file="melted_BM_PSIs_avsd.txt", header=T,sep = "\t")


### BIG TABLE WITH PSIs of NON_FILTERED INCL TABLE
# NonFil_melted_BM_PSIs=melt(merge(NonFil_B_INCL_dPSI10[,c("EVENT","GENE",PSIcols_B)],NonFil_M_INCL_dPSI10[,c("EVENT","GENE",PSIcols_M)],by=c("EVENT","GENE")), id.vars = c("EVENT","GENE"), variable.name = "Sample",value.name = "PSI")
# NonFil_melted_BM_PSIs$Dataset=sapply(as.character(NonFil_melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[1])
# NonFil_melted_BM_PSIs$Condition <- sapply(as.character(NonFil_melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[2])
# NonFil_melted_BM_PSIs$Condition <- paste(NonFil_melted_BM_PSIs$Dataset,NonFil_melted_BM_PSIs$Condition,sep = "_")
# NonFil_melted_BM_PSIs$Replicate=sapply(as.character(NonFil_melted_BM_PSIs$Sample), function(x) unlist(strsplit(x,split="_"))[3])
# head(NonFil_melted_BM_PSIs)
# write.table(NonFil_melted_BM_PSIs,file="NonFil_melted_BM_PSIs.txt",row.names = F, quote=F,sep = "\t")
NonFil_melted_BM_PSIs <- read.table(file="NonFil_melted_BM_PSIs.txt", header=T,sep = "\t")

# NonFil_melted_BM_PSIs_av=aggregate(PSI ~ EVENT+GENE+Dataset+Condition, NonFil_melted_BM_PSIs, mean)
# NonFil_melted_BM_PSIs_sd=aggregate(PSI ~ EVENT+GENE+Dataset+Condition, NonFil_melted_BM_PSIs, sd)
# NonFil_melted_BM_PSIs_avsd=merge(x=NonFil_melted_BM_PSIs_av,y=NonFil_melted_BM_PSIs_sd,by=c("EVENT","GENE","Dataset","Condition"),suffixes=c("_av","_sd"))
# NonFil_melted_BM_PSIs_avsd$EVENT_Dataset=apply(NonFil_melted_BM_PSIs_avsd, 1, function(x) paste(x["EVENT"], x["Dataset"],sep="_"))
# head(NonFil_melted_BM_PSIs_avsd)
# write.table(NonFil_melted_BM_PSIs_avsd,file="NonFil_melted_BM_PSIs_avsd.txt",row.names = F, quote=F,sep = "\t")
NonFil_melted_BM_PSIs_avsd <- read.table(file="NonFil_melted_BM_PSIs_avsd.txt",header=T,sep = "\t")

#----- Outputs
message("Outputs:
        melted_BM_EdgeR & melted_BM_EdgeR_avsd are melted dataframes (no cv filter) 
            to plot expression values respectively of single or average/SD of replicates
        melted_BM_EdgeR_scaled is a melted dataframe (no cv filter) to plot scaled expression values respectively of average of replicates")
#

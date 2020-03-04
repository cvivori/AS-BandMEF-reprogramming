require(tidyverse)
require(org.Mm.eg.db)


## MAKE VECTORS WITH NAMES OF CONDITIONS SO THAT THEY'RE UNIFORM!
conds_B_iPSES <- c("B_Bcells_1","B_Bcells_2","B_Bpulse_1","B_Bpulse_2","B_Day2_1","B_Day2_2","B_Day4_1","B_Day4_2","B_Day6_1","B_Day6_2","B_Day8_1","B_Day8_2","B_iPS_1","B_iPS_2","B_ES_1","B_ES_2")
conds_B_ESiPS <- c("B_Bcells_1","B_Bcells_2","B_Bpulse_1","B_Bpulse_2","B_Day2_1","B_Day2_2","B_Day4_1","B_Day4_2","B_Day6_1","B_Day6_2","B_Day8_1","B_Day8_2","B_ES_1","B_ES_2","B_iPS_1","B_iPS_2")
conds_av_B_iPSES <- c("B_Bcells","B_Bpulse","B_Day2","B_Day4","B_Day6","B_Day8","B_iPS","B_ES")
conds_av_B_ESiPS <- c("B_Bcells","B_Bpulse","B_Day2","B_Day4","B_Day6","B_Day8","B_ES","B_iPS")
conds_M_CloneFirst <- c("M_Clone_1","M_Clone_2","M_Clone_3","M_Day0_1","M_Day0_2","M_Day0_3","M_Day4_1","M_Day4_2","M_Day4_3","M_Day7_1","M_Day7_2","M_Day7_3","M_Day10_1","M_Day10_2","M_Day10_3","M_Day15_1","M_Day15_2","M_Day15_3","M_Day20_1","M_Day20_2","M_Day20_3") 
conds_M_CloneLast <- c("M_Day0_1","M_Day0_2","M_Day0_3","M_Day4_1","M_Day4_2","M_Day4_3","M_Day7_1","M_Day7_2","M_Day7_3","M_Day10_1","M_Day10_2","M_Day10_3","M_Day15_1","M_Day15_2","M_Day15_3","M_Day20_1","M_Day20_2","M_Day20_3","M_Clone_1","M_Clone_2","M_Clone_3") 
conds_av_M_CloneFirst <- c("M_Clone","M_Day0","M_Day4","M_Day7","M_Day10","M_Day15","M_Day20","M_Clone") 
conds_av_M_CloneLast <- c("M_Day0","M_Day4","M_Day7","M_Day10","M_Day15","M_Day20","M_Clone") 


## IMPORT
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/")
B_EdgeR <- read.table(file="FilteredGeneCpm_BtoiPS.txt",header = T,sep="\t") #filter min cpm 5
colnames(B_EdgeR) <- c("ENSGID","GeneName",conds_B_ESiPS)
rownames(B_EdgeR) <- B_EdgeR$ENSGID
head(B_EdgeR)
dim(B_EdgeR)

M_EdgeR <- read.table(file="FilteredGeneCpm_MEFstoiPS.txt",header = T,sep="\t") #filter min cpm 5
colnames(M_EdgeR) <- c("ENSGID","GeneName",conds_M_CloneFirst)
rownames(M_EdgeR) <- M_EdgeR$ENSGID
head(M_EdgeR)
dim(M_EdgeR)

## MAKE CONVERSION TABLE ENSGIDs-GENENAMES
ENS2Gene_B <- B_EdgeR[,c("ENSGID","GeneName")]
ENS2Gene_M <- M_EdgeR[,c("ENSGID","GeneName")]
ENS2Gene <- na.omit(rbind(ENS2Gene_B,ENS2Gene_M) %>%  distinct())
ENS2Gene$Gene <- mapIds(org.Mm.eg.db, keys=as.character(ENS2Gene$GeneName), column='SYMBOL', keytype='ALIAS',multiVals="first")



## CALCULATE AVERAGE CPMs and ORDER CONDITIONS
B_CPMs <- B_EdgeR %>% 
  dplyr::select(conds_B_iPSES)
head(B_CPMs)
B_CPMs_av <- B_EdgeR %>% 
  transmute(B_Bcells = (B_Bcells_1+B_Bcells_2)/2,
            B_Bpulse = (B_Bpulse_1+B_Bpulse_2)/2,
            B_Day2 = (B_Day2_1+B_Day2_2)/2,
            B_Day4 = (B_Day4_1+B_Day4_2)/2,
            B_Day6 = (B_Day6_1+B_Day6_2)/2,
            B_Day8 = (B_Day8_1+B_Day8_2)/2,
            B_iPS = (B_iPS_1+B_iPS_2)/2,
            B_ES = (B_ES_1+B_ES_2)/2 )
rownames(B_CPMs_av) <- rownames(B_EdgeR)
head(B_CPMs_av)

M_CPMs <- M_EdgeR %>% 
  dplyr::select(conds_M_CloneLast)
head(M_CPMs)
M_CPMs_av <- M_EdgeR %>% 
  transmute(M_Day0 = (M_Day0_1+M_Day0_2+M_Day0_3)/3,
            M_Day4 = (M_Day4_1+M_Day4_2+M_Day4_3)/3,
            M_Day7 = (M_Day7_1+M_Day7_2+M_Day7_3)/3,
            M_Day10 = (M_Day10_1+M_Day10_2+M_Day10_3)/3,
            M_Day15 = (M_Day15_1+M_Day15_2+M_Day15_3)/3,
            M_Day20 = (M_Day20_1+M_Day20_2+M_Day20_3)/3,
            M_Clone = (M_Clone_1+M_Clone_2+M_Clone_3)/3
  )
rownames(M_CPMs_av) <- rownames(M_EdgeR)
head(M_CPMs_av)





### FILTERING
## CALCULATE COEFFICIENT OF VARIATION *ON AVERAGE VALUES!*
filter_coeffvar <- 0.2

B_CPMs_av_tofil <- B_CPMs_av
B_CPMs_av_tofil$CoeffVar <- apply (B_CPMs_av, 1, function(x) sd(x)/mean(x))
B_CPMs_av_tofil <- subset(B_CPMs_av_tofil, CoeffVar >= filter_coeffvar)
B_CPMs_av_fil <- B_CPMs_av_tofil %>% dplyr::select(-CoeffVar)
dim(B_CPMs_av); dim(B_CPMs_av_fil); 
message("Filtered out: "); 1-nrow(B_CPMs_av_fil)/nrow(B_CPMs_av)
head(B_CPMs_av_fil)

M_CPMs_av_tofil <- M_CPMs_av
M_CPMs_av_tofil$CoeffVar <- apply (M_CPMs_av, 1, function(x) sd(x)/mean(x))
M_CPMs_av_tofil <- subset(M_CPMs_av_tofil, CoeffVar >= filter_coeffvar)
M_CPMs_av_fil <- M_CPMs_av_tofil %>% dplyr::select(-CoeffVar)
dim(M_CPMs_av); dim(M_CPMs_av_fil); 
message("Filtered out: "); 1-nrow(M_CPMs_av_fil)/nrow(M_CPMs_av)
head(M_CPMs_av_fil)

## OBTAIN SINGLE REPLICATE VALUES OF FILTERED GENES
B_CPMs_fil <- B_CPMs[rownames(B_CPMs_av_fil),]
M_CPMs_fil <- M_CPMs[rownames(M_CPMs_av_fil),]

## MERGE IN BIG TABLEs WITH BOTH DATASETS
BM_CPMs <- (merge(x=B_CPMs, y=M_CPMs, by="row.names",all = T))
rownames(BM_CPMs) <- BM_CPMs$Row.names
BM_CPMs <- BM_CPMs %>% dplyr::select(-Row.names)
head(BM_CPMs)

BM_CPMs_fil <- (merge(x=B_CPMs_fil, y=M_CPMs_fil, by="row.names",all = T))
rownames(BM_CPMs_fil) <- BM_CPMs_fil$Row.names
BM_CPMs_fil <- BM_CPMs_fil %>% dplyr::select(-Row.names)
head(BM_CPMs_fil)


BM_CPMs_av <- (merge(x=B_CPMs_av, y=M_CPMs_av, by="row.names",all = T))
rownames(BM_CPMs_av) <- BM_CPMs_av$Row.names
BM_CPMs_av <- BM_CPMs_av %>% dplyr::select(-Row.names)
head(BM_CPMs_av)


BM_CPMs_av_fil <- (merge(x=B_CPMs_av_fil, y=M_CPMs_av_fil, by="row.names",all=T))
rownames(BM_CPMs_av_fil) <- BM_CPMs_av_fil$Row.names
BM_CPMs_av_fil <- BM_CPMs_av_fil %>% dplyr::select(-Row.names)
head(BM_CPMs_av_fil)


dim(BM_CPMs); dim(BM_CPMs_av)
dim(BM_CPMs_fil); dim(BM_CPMs_av_fil)





#### SCALE VALUES
## SCALE SINGLE REPLICATE VALUES
B_CPMs_scaled <- t(scale(t(B_CPMs)))
head(B_CPMs_scaled)
M_CPMs_scaled <- t(scale(t(M_CPMs)))
head(M_CPMs_scaled)

sum(is.na(B_CPMs_scaled)); sum(is.na(M_CPMs_scaled))
dim(B_CPMs_scaled); dim(M_CPMs_scaled)


## SCALE SINGLE REPLICATE VALUES AFTER FILTERING
B_CPMs_fil_scaled <- t(scale(t(B_CPMs_fil)))
head(B_CPMs_fil_scaled)
M_CPMs_fil_scaled <- t(scale(t(M_CPMs_fil)))
head(M_CPMs_fil_scaled)

sum(is.na(B_CPMs_fil_scaled));sum(is.na(M_CPMs_fil_scaled))
dim(B_CPMs_fil_scaled); dim(M_CPMs_fil_scaled)


## SCALE AVERAGE VALUES
B_CPMs_av_scaled <- t(scale(t(B_CPMs_av)))
head(B_CPMs_av_scaled)
M_CPMs_av_scaled <- t(scale(t(M_CPMs_av)))
head(M_CPMs_av_scaled)

sum(is.na(B_CPMs_av_scaled)); sum(is.na(M_CPMs_av_scaled))
dim(B_CPMs_av_scaled); dim(M_CPMs_av_scaled)


## SCALE AVERAGE VALUES AFTER FILTERING
B_CPMs_av_fil_scaled <- t(scale(t(B_CPMs_av_fil)))
head(B_CPMs_av_fil_scaled)
M_CPMs_av_fil_scaled <- t(scale(t(M_CPMs_av_fil)))
head(M_CPMs_av_fil_scaled)

sum(is.na(B_CPMs_av_fil_scaled)); sum(is.na(M_CPMs_av_fil_scaled))
dim(B_CPMs_av_fil_scaled); dim(M_CPMs_av_fil_scaled)






## MERGE IN BIG TABLEs WITH BOTH DATASETS
BM_CPMs_scaled <- (merge(x=B_CPMs_scaled, y=M_CPMs_scaled, by="row.names",all=T))
rownames(BM_CPMs_scaled) <- BM_CPMs_scaled$Row.names
BM_CPMs_scaled <- BM_CPMs_scaled %>% dplyr::select(-Row.names)
head(BM_CPMs_scaled)
dim(BM_CPMs_scaled)

BM_CPMs_av_scaled <- (merge(x=B_CPMs_av_scaled, y=M_CPMs_av_scaled, by="row.names",all=T))
rownames(BM_CPMs_av_scaled) <- BM_CPMs_av_scaled$Row.names
BM_CPMs_av_scaled <- BM_CPMs_av_scaled %>% dplyr::select(-Row.names)
head(BM_CPMs_av_scaled)
dim(BM_CPMs_av_scaled)

BM_fil_CPMs_av_scaled <- (merge(x=B_CPMs_av_fil_scaled, y=M_CPMs_av_fil_scaled, by="row.names",all=T))
rownames(BM_fil_CPMs_av_scaled) <- BM_fil_CPMs_av_scaled$Row.names
BM_fil_CPMs_av_scaled <- BM_fil_CPMs_av_scaled %>% dplyr::select(-Row.names)
head(BM_fil_CPMs_av_scaled)
dim(BM_fil_CPMs_av_scaled)


#----- Outputs
message("Outputs for B and MEFs:
        * all tables are filtered by min 5 cpm (in 5/6 samples = 33%) 
        * _av tables contain average values among replicates
        * _fil tables are also filtered for coeff.var >= 0.2
        fil_B_CPMs_scaled & fil_M_CPMs_scaled = scaled cpm values, filtered 
        BM_fil_CPMs = cpm values of genes expressed in both dataset
        BM_fil_CPMs_scaled = scaled cpm values of genes expressed in both dataset
        conds(_av)_B/M_iPSES/CloneLast = vectors containing correct orders of conditions!")
#




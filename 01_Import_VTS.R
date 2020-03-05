require(tidyverse)

#### B2iPS
## IMPORT TABLE FROM VTS_add_dPSI_toINCL.R (Cluster)
  #dPSI10
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/")
B_INCL_dPSI10=read.table(file="Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noM_withdPSI_28DiffAS.txt",header = T, sep="\t",fill=T)
rownames(B_INCL_dPSI10) <- B_INCL_dPSI10$EVENT
colnames(B_INCL_dPSI10) <- gsub("_ALT_use25_","-",colnames(B_INCL_dPSI10))
summary(B_INCL_dPSI10$freq_N_sr)   #check that filtering 
  #dPSI15
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI15/")
B_INCL_dPSI15=read.table(file="Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noM_withdPSI_28DiffAS.txt",header = T, sep="\t",fill=T)
rownames(B_INCL_dPSI15) <- B_INCL_dPSI15$EVENT
colnames(B_INCL_dPSI15) <- gsub("_ALT_use25_","-",colnames(B_INCL_dPSI15))
summary(B_INCL_dPSI15$freq_N_sr)   #check that filtering 
  #dPSI10 with no N filtering
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/")
NonFil_B_INCL_dPSI10=read.table(file="INCLUSION_LEVELS_FULL-Mm237_withdPSI_28DiffAS.txt",header = T, sep="\t",fill=T)
rownames(NonFil_B_INCL_dPSI10) <- NonFil_B_INCL_dPSI10$EVENT
colnames(NonFil_B_INCL_dPSI10) <- gsub("_ALT_use25_","-",colnames(NonFil_B_INCL_dPSI10))
summary(NonFil_B_INCL_dPSI10$freq_N_sr)   #check that filtering 
head(NonFil_B_INCL_dPSI10)

## CATEGORIES OF COLUMNS
cols <- colnames(B_INCL_dPSI10)
STDcols <- cols[1:6]
Qcols <- cols[grep(".Q",cols)]
dPSIcols_B <- cols[grep("dPSI",cols)]
PSIcols <- c("B_Bcells_1","B_Bcells_2","B_Bpulse_1","B_Bpulse_2","B_Day2_1","B_Day2_2","B_Day4_1","B_Day4_2","B_Day6_1","B_Day6_2","B_Day8_1","B_Day8_2","B_iPS_1","B_iPS_2","B_ES_1","B_ES_2",
             "M_Day0_1","M_Day0_2","M_Day0_3","M_Day4_1","M_Day4_2","M_Day4_3","M_Day7_1","M_Day7_2","M_Day7_3","M_Day10_1","M_Day10_2","M_Day10_3","M_Day15_1","M_Day15_2","M_Day15_3","M_Day20_1","M_Day20_2","M_Day20_3","M_Clone_1","M_Clone_2","M_Clone_3") 
PSIcols_B <- PSIcols[grep("^B_",PSIcols)]

## CALCULATE MAX_dPSI
  #dPSI10
B_INCL_dPSI10$Max_dPSI <- apply(B_INCL_dPSI10[,dPSIcols_B], 1, function(x) max(x,na.rm = T))
B_INCL_dPSI10 <- B_INCL_dPSI10 %>%
  mutate(Max_dPSI = replace(Max_dPSI, Max_dPSI == -Inf, NA))
dim(B_INCL_dPSI10)
head(B_INCL_dPSI10)
summary(B_INCL_dPSI10)
# tibble(B_INCL_dPSI10)
summary(B_INCL_dPSI10$Max_dPSI)
  #dPSI15
B_INCL_dPSI15$Max_dPSI <- apply(B_INCL_dPSI15[,dPSIcols_B], 1, function(x) max(x,na.rm = T))
B_INCL_dPSI15 <- B_INCL_dPSI15 %>%
  mutate(Max_dPSI = replace(Max_dPSI, Max_dPSI == -Inf, NA))
dim(B_INCL_dPSI15)
head(B_INCL_dPSI15)
summary(B_INCL_dPSI15)
# tibble(B_INCL_dPSI15)
summary(B_INCL_dPSI15$Max_dPSI)

  #dPSI10
## NUMBER OF MAPPED EVENTS
nrow(B_INCL_dPSI10)
## NUMBER OF NON-CHANGING EVENTS
sum(is.na(B_INCL_dPSI10$Max_dPSI))
## NUMBER OF EVENTS DIFFERENTIALLY SPLICED AT LEAST IN ONE COMPARISON
nrow(B_INCL_dPSI10)-sum(is.na(B_INCL_dPSI10$Max_dPSI))
  #dPSI15
## NUMBER OF MAPPED EVENTS
nrow(B_INCL_dPSI15)
## NUMBER OF NON-CHANGING EVENTS
sum(is.na(B_INCL_dPSI15$Max_dPSI))
## NUMBER OF EVENTS DIFFERENTIALLY SPLICED AT LEAST IN ONE COMPARISON
nrow(B_INCL_dPSI15)-sum(is.na(B_INCL_dPSI15$Max_dPSI))


## CREATE LISTS with dPSI10 and 15
## EXTRACT ONLY DIFFERENTIAL EVENTS
B_diffEV <- list(dPSI10 = B_INCL_dPSI10[!is.na(B_INCL_dPSI10$Max_dPSI),], 
                 dPSI15 = B_INCL_dPSI15[!is.na(B_INCL_dPSI15$Max_dPSI),])
B_diffEV <- lapply(B_diffEV, function(x) {rownames(x) = x$EVENT; x})
lapply(B_diffEV, function(x) dim(x))

## EXTRACT ONLY PSIs of DIFFERENTIAL EVENTS
B_PSIs_VTS <- list(dPSI10 = B_diffEV$dPSI10[,PSIcols_B],
                   dPSI15 = B_diffEV$dPSI15[,PSIcols_B])
head(B_PSIs_VTS$dPSI10)
lapply(B_PSIs_VTS, function(x) dim(x))
## CHECK THAT THERE ARE NO NAs 
which(rowSums(is.na(B_PSIs_VTS$dPSI10))>0)

## CALCULATE AVERAGES OF PSIs
B_PSIs_VTS_av <-lapply(B_PSIs_VTS, function(y) { df = data.frame (
  B_Bcells = apply(y, 1, function(x) mean(as.numeric(x[c("B_Bcells_1","B_Bcells_2")])) ),
  B_Bpulse=apply(y, 1,function(x) mean(as.numeric(x[c("B_Bpulse_1","B_Bpulse_2")])) ),
  B_Day2=apply(y, 1,function(x) mean(as.numeric(x[c("B_Day2_1","B_Day2_2")])) ),
  B_Day4=apply(y, 1,function(x) mean(as.numeric(x[c("B_Day4_1","B_Day4_2")])) ),
  B_Day6=apply(y, 1,function(x) mean(as.numeric(x[c("B_Day6_1","B_Day6_2")])) ),
  B_Day8=apply(y, 1,function(x) mean(as.numeric(x[c("B_Day8_1","B_Day8_2")])) ),
  B_iPS=apply(y, 1,function(x) mean(as.numeric(x[c("B_iPS_1","B_iPS_2")])) ),
  B_ES=apply(y, 1,function(x) mean(as.numeric(x[c("B_ES_1","B_ES_2")])) ) );
  return(df)})
  



#### MEFs2iPS

## IMPORT TABLE FROM VTS_add_dPSI_toINCL.R (Cluster)
  #dPSI10
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/MEFs2iPS/dPSI10/")
M_INCL_dPSI10=read.table(file="Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noB_withdPSI_21DiffAS.txt",header = T, sep="\t",fill=T)
rownames(M_INCL_dPSI10) <- M_INCL_dPSI10$EVENT
colnames(M_INCL_dPSI10) <- gsub("_ALT_use25_","-",colnames(M_INCL_dPSI10))
summary(M_INCL_dPSI10$freq_N_sr)   #check that filtering 
  #dPSI15
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/MEFs2iPS/dPSI15/")
M_INCL_dPSI15=read.table(file="Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noB_withdPSI_21DiffAS.txt",header = T, sep="\t",fill=T)
rownames(M_INCL_dPSI15) <- M_INCL_dPSI15$EVENT
colnames(M_INCL_dPSI15) <- gsub("_ALT_use25_","-",colnames(M_INCL_dPSI15))
summary(M_INCL_dPSI15$freq_N_sr)   #check that filtering 
  #dPSI with no N filtering
setwd ("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/MEFs2iPS/dPSI10/")
NonFil_M_INCL_dPSI10=read.table(file="INCLUSION_LEVELS_FULL-Mm237_withdPSI_21DiffAS.txt",header = T, sep="\t",fill=T)
rownames(NonFil_M_INCL_dPSI10) <- NonFil_M_INCL_dPSI10$EVENT
colnames(NonFil_M_INCL_dPSI10) <- gsub("_ALT_use25_","-",colnames(NonFIl_M_INCL_dPSI10))
summary(NonFil_M_INCL_dPSI10$freq_N_sr)   #check that filtering 
head(NonFil_M_INCL_dPSI10)


## CATEGORIES OF COLUMNS
cols <- colnames(M_INCL_dPSI10)
STDcols <- cols[1:6]
Qcols <- cols[grep(".Q",cols)]
dPSIcols_M <- cols[grep("dPSI",cols)]
#PSIcols <- setdiff(cols,c(STDcols,Qcols,dPSIcols,"freq_N_sr"))
PSIcols <- c("B_Bcells_1","B_Bcells_2","B_Bpulse_1","B_Bpulse_2","B_Day2_1","B_Day2_2","B_Day4_1","B_Day4_2","B_Day6_1","B_Day6_2","B_Day8_1","B_Day8_2","B_iPS_1","B_iPS_2","B_ES_1","B_ES_2",
             "M_Day0_1","M_Day0_2","M_Day0_3","M_Day4_1","M_Day4_2","M_Day4_3","M_Day7_1","M_Day7_2","M_Day7_3","M_Day10_1","M_Day10_2","M_Day10_3","M_Day15_1","M_Day15_2","M_Day15_3","M_Day20_1","M_Day20_2","M_Day20_3","M_Clone_1","M_Clone_2","M_Clone_3") 
PSIcols_M <- PSIcols[grep("^M_",PSIcols)]


## CALCULATE MAX_dPSI
#dPSI10
M_INCL_dPSI10$Max_dPSI <- apply(M_INCL_dPSI10[,dPSIcols_M], 1, function(x) max(x,na.rm = T))
M_INCL_dPSI10 <- M_INCL_dPSI10 %>%
  mutate(Max_dPSI = replace(Max_dPSI, Max_dPSI == -Inf, NA))
dim(M_INCL_dPSI10)
head(M_INCL_dPSI10)
summary(M_INCL_dPSI10)
#dPSI15
M_INCL_dPSI15$Max_dPSI <- apply(M_INCL_dPSI15[,dPSIcols_M], 1, function(x) max(x,na.rm = T))
M_INCL_dPSI15 <- M_INCL_dPSI15 %>%
  mutate(Max_dPSI = replace(Max_dPSI, Max_dPSI == -Inf, NA))
dim(M_INCL_dPSI15)
head(M_INCL_dPSI15)
summary(M_INCL_dPSI15)

#dPSI10
## NUMBER OF MAPPED EVENTS
nrow(M_INCL_dPSI10)
## NUMBER OF NON-CHANGING EVENTS
sum(is.na(M_INCL_dPSI10$Max_dPSI))
## NUMBER OF EVENTS DIFFERENTIALLY SPLICED AT LEAST IN ONE COMPARISON
nrow(M_INCL_dPSI10)-sum(is.na(M_INCL_dPSI10$Max_dPSI))
#dPSI15
## NUMBER OF MAPPED EVENTS
nrow(M_INCL_dPSI15)
## NUMBER OF NON-CHANGING EVENTS
sum(is.na(M_INCL_dPSI15$Max_dPSI))
## NUMBER OF EVENTS DIFFERENTIALLY SPLICED AT LEAST IN ONE COMPARISON
nrow(M_INCL_dPSI15)-sum(is.na(M_INCL_dPSI15$Max_dPSI))


## CREATE LISTS with dPSI10 and 15
## EXTRACT ONLY DIFFERENTIAL EVENTS
M_diffEV <- list(dPSI10 = M_INCL_dPSI10[!is.na(M_INCL_dPSI10$Max_dPSI),], 
                 dPSI15 = M_INCL_dPSI15[!is.na(M_INCL_dPSI15$Max_dPSI),])
M_diffEV <- lapply(M_diffEV, function(x) {rownames(x) = x$EVENT; x})
lapply(M_diffEV, function(x) dim(x))

## EXTRACT ONLY PSIs of DIFFERENTIAL EVENTS
M_PSIs_VTS <- list(dPSI10 = M_diffEV$dPSI10[,PSIcols_M],
                   dPSI15 = M_diffEV$dPSI15[,PSIcols_M])
head(M_PSIs_VTS$dPSI10)
lapply(M_PSIs_VTS, function(x) dim(x))
## CHECK THAT THERE ARE NO NAs 
which(rowSums(is.na(M_PSIs_VTS$dPSI10))>0)

## CALCULATE AVERAGES OF PSIs
M_PSIs_VTS_av <-lapply(M_PSIs_VTS, function(y) { df = data.frame (
  M_Day0=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day0_1","M_Day0_2","M_Day0_3")]))),
  M_Day4=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day4_1","M_Day4_2","M_Day4_3")])) ),
  M_Day7=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day7_1","M_Day7_2","M_Day7_3")])) ),
  M_Day10=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day10_1","M_Day10_2","M_Day10_3")])) ),
  M_Day15=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day15_1","M_Day15_2","M_Day15_3")])) ),
  M_Day20=apply(y, 1,function(x) mean(as.numeric(x[c("M_Day20_1","M_Day20_2","M_Day20_3")])) ),
  M_Clone=apply(y, 1,function(x) mean(as.numeric(x[c("M_Clone_1","M_Clone_2","M_Clone_3")])) ));
return(df)})


  


#----- Outputs
message("Outputs for B and MEFs:\n\tB_diffEV & M_diffEV = only events differentially spliced in at least one comparison of each dataset,\tbut all values in both
        B_PSIs_VTS & M_PSIs_VTS = only PSI columns of diffEV tables in the corresponding dataset
        B_PSIs_VTS_av & M_PSIs_VTS_av = only average PSI columns of diffEV in the corresponding dataset
        BM_PSIs_VTS = only PSI columns of diffEV tables in both datasets
        \t*all are lists with dPSI10/15*")
#


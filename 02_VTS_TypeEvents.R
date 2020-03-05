## FUNCTIONs
Wrap_VTS_Events = function (input, microexons = c("in_CEx","out_CEx")) {
  
  all_counts <- table(input$COMPLEX)
  
  if (sum(grepl("IR-S",unique(input$COMPLEX)))>0) {
    RIs <- c("IR-C","IR-S") } else { RIs <- "IR" }
  
  CEx_inclMIC <- c("S","C1","C2","C3","MIC","ANN") 
  CEx_all <- subset(input, COMPLEX %in% CEx_inclMIC)
  MICs <- subset(input, LENGTH <= 27 & COMPLEX %in% CEx_inclMIC)
  CEx_long <- subset(CEx_all, !EVENT %in% MICs$EVENT)
  
  if (microexons == "in_CEx") {
    out <- c(all_counts["Alt3"],all_counts["Alt5"],
             IR=sum(all_counts[RIs]),
             CEx=nrow(CEx_all) )
  } else if (microexons == "out_CEx") {
    out <- c(all_counts["Alt3"],all_counts["Alt5"],
             IR=sum(all_counts[RIs]),
             CEx=nrow(CEx_long), µEx=nrow(MICs) )
  }
  return(out) 
}



###### NUMBER OF EVENTS IN THE WHOLE DATASETS ######

## EXTRACT EVENT IDS FROM ALL COMPARISONS TOGETHER
B_diffEV_list <- lapply(B_diffEV, function(x) { y = list(CEx = x %>% filter(COMPLEX %in% c("S","C1","C2","C3","ANN","MIC")),
                                                     #MIC = x %>% filter(COMPLEX == "MIC"),
                                                     RI = x %>% filter(COMPLEX %in% c("IR-S","IR-C","IR")),
                                                     Alt3 = x %>% filter(COMPLEX == "Alt3"),
                                                     Alt5 = x %>% filter(COMPLEX == "Alt5"));
                                            return(y) })
lapply(B_diffEV_list$dPSI10, length)

M_diffEV_list <- lapply(M_diffEV, function(x) { y = list(CEx = x %>% filter(COMPLEX %in% c("S","C1","C2","C3","ANN","MIC")),
                                                     #MIC = x %>% filter(COMPLEX == "MIC"),
                                                     RI = x %>% filter(COMPLEX %in% c("IR-S","IR-C","IR")),
                                                     Alt3 = x %>% filter(COMPLEX == "Alt3"),
                                                     Alt5 = x %>% filter(COMPLEX == "Alt5"));
                                            return(y) })
lapply(M_diffEV_list$dPSI10, length)



## CALCULATE NUMBER OF EVENTS
# Calculate how many events in each table
B_Numbers=lapply(B_diffEV,  function(x) table(x[,"COMPLEX"]))
B_Numbers
M_Numbers=lapply(M_diffEV, function(x) table(x[,"COMPLEX"]))
M_Numbers

# Summarize events in CEx and IR + Alt3 and Alt5
B_nEV=lapply(B_diffEV,function(x) Wrap_VTS_Events(x,microexons = "in_CEx") )
B_nEV
M_nEV=lapply(M_diffEV,function(x) Wrap_VTS_Events(x,microexons = "in_CEx") )
M_nEV




###### NUMBER OF EVENTS IN EACH COMPARISON ######

## CREATE TABLE WITH ALL POSSIBLE COMPARISONS
B_couples <- data.frame(str_split_fixed(str_split_fixed(dPSIcols_B,pattern = "-",n=2)[,2], pattern="[.]", n=3)[,c(1,3)])
    rownames(B_couples) <- dPSIcols_B
    B_couples <- rownames_to_column(B_couples)
    colnames(B_couples) <- c("Comparison","Condition1","Condition2")
    B_couples$Condition1 <- factor(B_couples$Condition1,levels=conds_av_B_iPSES)
    B_couples$Condition2 <- factor(B_couples$Condition2,levels=conds_av_B_iPSES)
  head(B_couples)
M_couples <- data.frame(str_split_fixed(str_split_fixed(dPSIcols_M,pattern = "_ALT_use25_",n=2)[,2], pattern="[.]", n=3)[,c(1,3)])
    rownames(M_couples) <- dPSIcols_M
    M_couples <- rownames_to_column(M_couples)
    colnames(M_couples) <- c("Comparison","Condition1","Condition2")
    M_couples$Condition1 <- factor(M_couples$Condition1,levels=conds_av_M_CloneLast)
    M_couples$Condition2 <- factor(M_couples$Condition2,levels=conds_av_M_CloneLast)
  head(M_couples)

## CREATE TABLE WITH EVENTS FROM single COMPARISONS
B_diffEV_comparisons <- list()
    for (c in dPSIcols_B) { B_diffEV_comparisons$dPSI10[[c]] <- B_diffEV$dPSI10[!is.na(B_diffEV$dPSI10[,c]),]}
    for (c in dPSIcols_B) { B_diffEV_comparisons$dPSI15[[c]] <- B_diffEV$dPSI15[!is.na(B_diffEV$dPSI15[,c]),]}
M_diffEV_comparisons <- list()
    for (c in dPSIcols_M) { M_diffEV_comparisons$dPSI10[[c]] <- M_diffEV$dPSI10[!is.na(M_diffEV$dPSI10[,c]),]}
    for (c in dPSIcols_M) { M_diffEV_comparisons$dPSI15[[c]] <- M_diffEV$dPSI15[!is.na(M_diffEV$dPSI15[,c]),]}

## CALCULATE NUMBER OF EVENTS in each COMPARISON
tmpB <-list( dPSI10 = data.frame( do.call(rbind,lapply(B_diffEV_comparisons$dPSI10, function(x) table(x[,"COMPLEX"]))) ),
             dPSI15 = data.frame( do.call(rbind,lapply(B_diffEV_comparisons$dPSI15, function(x) table(x[,"COMPLEX"]))) ))
  tmpB <- lapply(tmpB, function(x) rownames_to_column(x))
B_diffEV_comparisons_Numbers <- lapply(tmpB, function(t) merge(B_couples,t, by.x="Comparison", by.y="rowname",all=T))

tmpB <-list( dPSI10 = data.frame( do.call(rbind,lapply(B_diffEV_comparisons$dPSI10, function(x) Wrap_inclMIC(table(x[,"COMPLEX"])))) ),
             dPSI15 = data.frame( do.call(rbind,lapply(B_diffEV_comparisons$dPSI15, function(x) Wrap_inclMIC(table(x[,"COMPLEX"])))) ))
  tmpB <- lapply(tmpB, function(x) rownames_to_column(x))
B_diffEV_comparisons_EV <- lapply(tmpB, function(t) merge(B_couples,t, by.x="Comparison", by.y="rowname",all=T))
      B_diffEV_comparisons_EV           

tmpM <-list( dPSI10 = data.frame( do.call(rbind,lapply(M_diffEV_comparisons$dPSI10, function(x) table(x[,"COMPLEX"]))) ),
             dPSI15 = data.frame( do.call(rbind,lapply(M_diffEV_comparisons$dPSI15, function(x) table(x[,"COMPLEX"]))) ))
    tmpM <- lapply(tmpM, function(x) rownames_to_column(x))
M_diffEV_comparisons_Numbers <- lapply(tmpM, function(t) merge(M_couples,t, by.x="Comparison", by.y="rowname",all=T))

tmpM <-list( dPSI10 = data.frame( do.call(rbind,lapply(M_diffEV_comparisons$dPSI10, function(x) Wrap_inclMIC(table(x[,"COMPLEX"])))) ),
             dPSI15 = data.frame( do.call(rbind,lapply(M_diffEV_comparisons$dPSI15, function(x) Wrap_inclMIC(table(x[,"COMPLEX"])))) ))
    tmpM <- lapply(tmpM, function(x) rownames_to_column(x))
M_diffEV_comparisons_EV <- lapply(tmpM, function(t) merge(M_couples,t, by.x="Comparison", by.y="rowname",all=T))
  M_diffEV_comparisons_EV           

rm(tmpB); rm(tmpM)

## WRITE OUTPUTS  
# for (d in names(B_diffEV)) {
# setwd(paste("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS",d,"NumberEvents/",sep="/"))
#   write.table(B_Numbers[[d]],file = "NumberEvents_WholeDataset_BtoiPS.txt",sep="\t",quote=F,row.names = F)
#   write.table(B_diffEV_comparisons_Numbers[[d]],file = "NumberEvents_EachComparis_BtoiPS.txt",sep="\t",quote=F,row.names = F)
# }
# for (d in names(M_diffEV)) {
# setwd(paste("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/MEFs2iPS",d,"NumberEvents/",sep="/"))
#   write.table(M_Numbers[[d]],file = "NumberEvents_WholeDataset_MEFstoiPS.txt",sep="\t",quote=F,row.names = F)
#   write.table(M_diffEV_comparisons_Numbers[[d]],file = "NumberEvents_EachComparis_MEFstoiPS.txt",sep="\t",quote=F,row.names = F)
# }  
  
  
  
  ## EXTRACT PSIs of EACH SET OF EXONS in B/M_PSIs_VTS - including GENEname
  B_PSIs_VTS_list = M_PSIs_VTS_list = list()
  for (d in names(B_PSIs_VTS)) { 
    B_PSIs_VTS_list[[d]] <- list(CEx = subset(B_PSIs_VTS[[d]], row.names(B_PSIs_VTS[[d]]) %in% B_diffEV_list[[d]]$CEx$EVENT),
                                 RI = subset(B_PSIs_VTS[[d]], row.names(B_PSIs_VTS[[d]]) %in% B_diffEV_list[[d]]$RI$EVENT),
                                 Alt3 = subset(B_PSIs_VTS[[d]], row.names(B_PSIs_VTS[[d]]) %in% B_diffEV_list[[d]]$Alt3$EVENT),
                                 Alt5 = subset(B_PSIs_VTS[[d]], row.names(B_PSIs_VTS[[d]]) %in% B_diffEV_list[[d]]$Alt5$EVENT)) }
  for (d in names(M_PSIs_VTS)) { 
    M_PSIs_VTS_list[[d]] <- list(CEx = subset(M_PSIs_VTS[[d]], row.names(M_PSIs_VTS[[d]]) %in% M_diffEV_list[[d]]$CEx$EVENT),
                                 RI = subset(M_PSIs_VTS[[d]], row.names(M_PSIs_VTS[[d]]) %in% M_diffEV_list[[d]]$RI$EVENT),
                                 Alt3 = subset(M_PSIs_VTS[[d]], row.names(M_PSIs_VTS[[d]]) %in% M_diffEV_list[[d]]$Alt3$EVENT),
                                 Alt5 = subset(M_PSIs_VTS[[d]], row.names(M_PSIs_VTS[[d]]) %in% M_diffEV_list[[d]]$Alt5$EVENT)) }
  
  lapply(B_PSIs_VTS_list$dPSI10, dim)
  lapply(M_PSIs_VTS_list$dPSI10, dim)
  head(M_PSIs_VTS_list$dPSI10$CEx)
  
  ## EXTRACT AVERAGE PSIs EACH SET OF EXONS in B/M_PSIs_VTS_av_av
  B_PSIs_VTS_av_list = M_PSIs_VTS_av_list = list()
  for (d in names(B_PSIs_VTS_av)) { 
    B_PSIs_VTS_av_list[[d]] <- list(CEx = subset(B_PSIs_VTS_av[[d]], row.names(B_PSIs_VTS_av[[d]]) %in% B_diffEV_list[[d]]$CEx$EVENT),
                                    RI = subset(B_PSIs_VTS_av[[d]], row.names(B_PSIs_VTS_av[[d]]) %in% B_diffEV_list[[d]]$RI$EVENT),
                                    Alt3 = subset(B_PSIs_VTS_av[[d]], row.names(B_PSIs_VTS_av[[d]]) %in% B_diffEV_list[[d]]$Alt3$EVENT),
                                    Alt5 = subset(B_PSIs_VTS_av[[d]], row.names(B_PSIs_VTS_av[[d]]) %in% B_diffEV_list[[d]]$Alt5$EVENT)) }
  for (d in names(M_PSIs_VTS_av)) { 
    M_PSIs_VTS_av_list[[d]] <- list(CEx = subset(M_PSIs_VTS_av[[d]], row.names(M_PSIs_VTS_av[[d]]) %in% M_diffEV_list[[d]]$CEx$EVENT),
                                    RI = subset(M_PSIs_VTS_av[[d]], row.names(M_PSIs_VTS_av[[d]]) %in% M_diffEV_list[[d]]$RI$EVENT),
                                    Alt3 = subset(M_PSIs_VTS_av[[d]], row.names(M_PSIs_VTS_av[[d]]) %in% M_diffEV_list[[d]]$Alt3$EVENT),
                                    Alt5 = subset(M_PSIs_VTS_av[[d]], row.names(M_PSIs_VTS_av[[d]]) %in% M_diffEV_list[[d]]$Alt5$EVENT)) }
  lapply(B_PSIs_VTS_av_list$dPSI10, dim)
  lapply(M_PSIs_VTS_av_list$dPSI10, dim)
  head(M_PSIs_VTS_av_list$dPSI10$CEx)  

#----- Outputs
message("Outputs for B and MEFs:\n\tB_Numbers & M_Numbers = number of events in each category,
        B_EV & M_EV = Wrapped in *Events (CEx, IR, Alt3 and Alt5)
        \nlist of tables including all events/CEx/RI/Alt3/Alt5:
        \n\tB_diffEV_list & M_diffEV_list = only events differentially spliced in at least one comparison of each dataset,\tbut all values in both
        B_PSIs_VTS_list & M_PSIs_VTS_list = only PSI columns of diffEV tables in the corresponding dataset
        B_PSIs_VTS_av_list & M_PSIs_VTS_av_list = only average PSI columns of diffEV in the corresponding dataset
        *all are lists with dPSI10/15*")
#


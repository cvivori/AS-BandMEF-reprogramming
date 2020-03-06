require(tidyverse)

## FUNCTIONs
Plot_ONTO <- function (df_superonto, df_supersuperonto, conditions_to_plot, palette_superonto) {
  
  ll=list(SO=df_superonto,SSO=df_supersuperonto)
  pp=list(SO=palette_superonto,SSO=palette_superonto[c(1,3,6)])
  cc=list(SO=SO,SSO=SSO)
  
  for (n in names(ll)) {
    df = ll[[n]]
    df_perc = as.data.frame(apply(df, 1, function(x) (x/sum(x))))
    df_cum = as.data.frame(t(cumsum(df_perc)))
    df_perc= as.data.frame(t(df_perc))
    #colnames(df_cum)=paste(colnames(df_cum),"_cum",sep="")
    #df_perc=cbind.data.frame(t(df_perc),df_cum)
    df$Conditions=rownames(df)
    df_perc$Conditions=rownames(df_perc)
    df_cum$Conditions=rownames(df_cum)
    
    mdf=melt(df,
             id.vars = "Conditions", 
             measure.vars = cc[[n]], 
             variable.name = "ORF_disruption", 
             value.name = "Number_of_events")
    if (sum(grepl("variable", colnames(mdf))) || sum(grepl("value", colnames(mdf))) >0) {
      colnames(mdf)[2:3]=c("ORF_disruption","Number_of_events")
    } else {}
    mdfp=melt(df_perc, 
              id.vars = "Conditions", 
              variable.name = "ORF_disruption", 
              measure.vars = 1:(ncol(df_perc)-1),
              value.name = "Fraction_of_events_perc")
    if (sum(grepl("variable", colnames(mdfp))) || sum(grepl("value", colnames(mdfp))) >0) {
      colnames(mdfp)[2:3]=c("ORF_disruption","Fraction_of_events_perc")
    } else {}
    mdfc=melt(df_cum,
              id.vars = "Conditions", 
              variable.name = "ORF_disruption", 
              measure.vars = 1:(ncol(df_cum)-1), 
              value.name = "Fraction_of_events_cum")
    if (sum(grepl("variable", colnames(mdfc))) || sum(grepl("value", colnames(mdfc))) >0) {
      colnames(mdfc)[2:3]=c("ORF_disruption","Fraction_of_events_cum")
    } else {}
    mdfpc=merge.data.frame(mdfp,mdfc,by=c("Conditions","ORF_disruption"))
    
    tab = subset(mdf, Conditions %in% conditions_to_plot)
    longestlabel=max(sapply(tab$Conditions, function(x) nchar(x)))
    if (longestlabel > 10 ) { ang = 20; h = 1; v = 1 }    else { ang = 0; h = 0.5; v = 0.5 }
    
    p <- ggplot(tab) + 
      geom_bar(aes(x=Conditions, y=Number_of_events, fill=ORF_disruption), 
               stat="identity", position = "fill") 
    
    plot_tab = ggplot_build(p)$data[[1]]
    plot_tab = plot_tab[order(plot_tab$group),]
    tabp = subset(mdfpc, Conditions %in% conditions_to_plot)
    tabp$coord=apply(plot_tab, 1, function(x) mean(as.numeric(x[c("ymax","ymin")])))
    
    
    p +
      scale_fill_viridis(option = "D",discrete = T) +
      scale_x_discrete(limits=conditions_to_plot) +
      scale_y_continuous(name = c("Percentage of events"), 
                         breaks = c(0,0.25,0.5,0.75,1), 
                         labels=c(0,25,50,75,100)) +
      geom_text(data=tabp, aes(x= Conditions, y=coord,
                               label=paste(100*round(Fraction_of_events_perc,2),"%",sep="") ,
                               hjust=0.5,vjust=0.5,size=7),
                position = position_dodge(width=1),size=3) + 
      theme_classic() +
      theme (text = element_text(color="grey20",size=10),
             axis.title = element_text(face="bold"), 
             axis.text.x = element_text(angle = ang, hjust =h, vjust =v),
             legend.position = "bottom", legend.title = element_blank(),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
             strip.background = element_blank(),strip.text = element_text(face="bold"))
    
    if (length(conditions_to_plot) <= 4) { ggsave(paste("ONTO_",paste(conditions_to_plot,collapse="-"),"_",n,".pdf",sep = ""),width=6,height=6,device = cairo_pdf) 
      } else { ggsave(paste("ONTO_",paste(c(conditions_to_plot[1],conditions_to_plot[length(conditions_to_plot)]),collapse="--"),"_",n,".pdf",sep = ""),width=10,height=6,device = cairo_pdf) }
  }
}
FisherTest_ORF <- function(SUPERONTO_table, which_cluster, which_prediction) {
  
  SUPERONTO_table$total=rowSums(SUPERONTO_table)
  totest <- SUPERONTO_table %>%
    rownames_to_column(var = "cluster") %>%
    filter(cluster %in% c("all",which_cluster)) %>%
    dplyr::select(c("cluster",which_prediction,"total")) 
  totest <- totest %>%
    mutate(rest = total - totest[,which_prediction]) %>%
    column_to_rownames(var = "cluster")
  contingency <- totest %>%
    rownames_to_column(var = "cluster") %>%
    dplyr::select(-total) %>%
    column_to_rownames(var = "cluster")
  p <- fisher.test(x = contingency, alternative = "two.sided")$p.value
  return(list(proportions = contingency/totest$total,p_value = p))
}



########## IMPORT ONTO TABLE ##########
source("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/VTS_ToolsAndPlots/VTS_Import_ONTO_Mm2_mm10_.R")
head(ONTO)

ONTO$ONTOGENY= factor(ONTO$ONTOGENY, levels=O)
ONTO$SUPERONTO= factor(ONTO$SUPERONTO, levels=SO)
ONTO$SUPERSUPERONTO= factor(ONTO$SUPERSUPERONTO, levels=SSO)
ONTO$PlotDisruption= factor(ONTO$PlotDisruption, levels = D)
ONTO$PlotNonCoding= factor(ONTO$PlotNonCoding,levels= N)
head(ONTO)

########## IMPORT DATA ##########
i=12
SETS_ev=cl_events_VTS
SETS_genes=cl_genes_VTS
SETS_ev_EX=SETS_ev
all <- rownames(rc2)


# CALCULATE numbers of each (CASSETTE EXONS)
SETS_EX_ONTO=lapply(SETS_ev_EX, function(x) subset(ONTO, EVENT %in% x))
SETS_EX_O=lapply(SETS_EX_ONTO, function(x) table(x[,"ONTOGENY"])[O])
SETS_EX_SO=lapply(SETS_EX_ONTO, function(x) table(x[,"SUPERONTO"])[SO])
SETS_EX_SSO=lapply(SETS_EX_ONTO, function(x) table(x[,"SUPERSUPERONTO"])[SSO])
all_EX_ONTO <-  subset(ONTO, EVENT %in% all)

# Merge tables in one for plotting
ONTOGENY_EX=as.data.frame(rbind(do.call(rbind,SETS_EX_O),table(all_EX_ONTO[,"ONTOGENY"])[O]))
    ONTOGENY_EX <- ONTOGENY_EX[c(13,in_clusters,out_clusters),]
    rownames(ONTOGENY_EX) <- sapply(rownames(ONTOGENY_EX), function(x) paste("clVTS",x,sep = "_"))
    rownames(ONTOGENY_EX) <- str_replace(rownames(ONTOGENY_EX),pattern="clVTS_13",replacement="all")
SUPERONTO_EX=as.data.frame(rbind(do.call(rbind,SETS_EX_SO),table(all_EX_ONTO[,"SUPERONTO"])[SO]))
    SUPERONTO_EX <- SUPERONTO_EX[c(13,in_clusters,out_clusters),]
    rownames(SUPERONTO_EX) <- sapply(rownames(SUPERONTO_EX), function(x) paste("clVTS",x,sep = "_"))
    rownames(SUPERONTO_EX) <- str_replace(rownames(SUPERONTO_EX),pattern="clVTS_13",replacement="all")
SSUPERONTO_EX=as.data.frame(rbind(do.call(rbind,SETS_EX_SSO),table(all_EX_ONTO[,"SUPERSUPERONTO"])[SSO]))
    SSUPERONTO_EX <- SSUPERONTO_EX[c(13,in_clusters,out_clusters),]
    rownames(SSUPERONTO_EX)=sapply(rownames(SSUPERONTO_EX), function(x) paste("clVTS",x,sep = "_"))
    rownames(SSUPERONTO_EX) <- str_replace(rownames(SSUPERONTO_EX),pattern="clVTS_13",replacement="all")
    
cl_Data_VTS_ORF=list()
for (x in c(1:i)) {
  cl_Data_VTS_ORF[[x]] = merge.data.frame(cl_Data_VTS[[x]],SETS_EX_ONTO[[x]],by = c("GENE","EVENT","COORD","LENGTH","FullCO","COMPLEX"))
}
head(cl_Data_VTS_ORF[[x]])

## PLOT CASSETTE EXONS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/ORF_Disruption/")
Plot_ONTO(df_superonto = SUPERONTO_EX, df_supersuperonto = SSUPERONTO_EX,
          conditions_to_plot = rownames(SUPERONTO_EX), palette_superonto = viridis(6))

## OUTPUT TABLES CASSETTE EXONS
# write.table(ONTOGENY_EX, file= "NumbersONTO_ONTOGENY.txt", sep="\t", quote =F)
# write.table(SUPERONTO_EX, file= "NumbersONTO_SUPERONTO.txt", sep="\t", quote =F)
# for (j in 1:i) {
#   write.csv(cl_Data_VTS_ORF[[j]], file=paste("cl_",j,"_Data_ORF.csv",sep = ""),row.names = F)
# }


## STATS ON SUPERONTO

SUPERONTO_table <- SUPERONTO_EX
which_cluster <- "clVTS_4"
which_prediction <- "CDS_DISR_uEXC"

# PULSE  
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
                which_cluster = "clVTS_4", which_prediction = "CDS_DISR_uEXC")
# EARLY
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
               which_cluster = "clVTS_10", which_prediction = "CDS_DISR_uEXC")
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
               which_cluster = "clVTS_3", which_prediction = "CDS_DISR_uINC")
# MIDDLE
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
               which_cluster = "clVTS_9", which_prediction = "CDS_PROT")
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
               which_cluster = "clVTS_1", which_prediction = "CDS_PROT")
FisherTest_ORF(SUPERONTO_table = SSUPERONTO_EX, 
               which_cluster = "clVTS_1", which_prediction = "CDS_DISR")
# LATE
FisherTest_ORF(SUPERONTO_table = SUPERONTO_EX, 
               which_cluster = "clVTS_12", which_prediction = "CDS_PROT")
FisherTest_ORF(SUPERONTO_table = SSUPERONTO_EX, 
               which_cluster = "clVTS_12", which_prediction = "CDS_DISR")

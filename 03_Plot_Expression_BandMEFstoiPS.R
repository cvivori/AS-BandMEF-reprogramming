require(tidyverse)

## FUNCTIONs
Plot_GEx_inBM = function (mdf_avsd, toplot, x_labels) {
  
    sapply(toplot, function(tp) {
  
  #POSITIVE REGULATORS
  mdf_toplot <-  mdf_avsd %>%
    filter(Condition != "ES") %>%
    filter(GeneName == tp) 
  mdf_toplotB=subset( mdf_toplot,Dataset == "B")
  mdf_toplotM=subset( mdf_toplot,Dataset == "M")
  
  min=min(mdf_toplot$cpm_av)
  max=max(mdf_toplot$cpm_av)
  if (nrow(mdf_toplot)== 0) {
    message("No plot generated for ",tp,".. most probably not expressed in MEFs2iPS?")
  } else {
  
  ggplot(mdf_toplot,
         aes(x=Time,y=cpm_av, color=Dataset, group=GeneName_Dataset)) +
      geom_vline(xintercept =c(2.5,4.5,7.5) , color = "gray80", size= 0.4) +
      annotate("text", x=1.5, y=-1.6, label="Diff.", hjust = 0.5, vjust=1, color = "grey20", size=3, fontface="bold") +
      annotate("text", x=3.5, y=-1.6, label="Early", hjust = 0.5, vjust=1, color = "grey20", size=3, fontface="bold") +
      annotate("text", x=6, y=-1.6, label="Late", hjust = 0.5, vjust=1, color = "grey20", size=3, fontface="bold") +
      annotate("text", x=8, y=-1.6, label="iPS", hjust = 0.5, vjust=1, color = "grey20", size=3, fontface="bold") +
    geom_point(data=mdf_toplotB,
               aes(x=Time,y=cpm),size=3) +
    geom_line(data=mdf_toplotB[!is.na(mdf_toplotB$cpm_av),],
              aes(x=Time,y=cpm_av),size=1.5) +
    geom_point(data=mdf_toplotM[!is.na(mdf_toplotM$cpm_av),],
               size=3) +
    geom_line(data=mdf_toplotM[!is.na(mdf_toplotM$cpm_av),],
              aes(x=Time,y=cpm_av),size=1.5) +
    geom_errorbar(data=mdf_toplotM[!is.na(mdf_toplotM$cpm_av),],
                  aes(ymin=cpm_av-cpm_sd, ymax=cpm_av+cpm_sd), width=0.3, position = position_dodge(0)) +
    scale_x_discrete(labels=x_labels) +
    scale_y_continuous( limits = c(-2,max)) +
    scale_color_manual(values=c("#2e4e4e","#af305f"),
                       labels=c("Bcell_reprogramming","MEFs_reprogramming")) +
    theme_classic() +
    theme (text = element_text(color="grey20",size=11),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(angle=30,hjust=1,vjust=1),
           legend.position = "bottom", legend.title = element_text(face="bold"),
           panel.grid.minor.x=element_blank(),
           panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
           strip.background = element_blank(),strip.text = element_text(face="bold"))
    ggsave(paste("GEx_CPMavsd_BandMEFs_",tp,".pdf",sep = ""), width=7,height=6,device = cairo_pdf)
  }
 })
}



## ADD PCA TIME TO mdf of CPM values
time_PCA=c("B_Bcells/M_Day0","B_Bpulse","B_Day2/M_Day4","B_Day4/M_Day7","B_Day6/M_Day10","M_Day15","B_Day8/M_Day20","iPS")

head(melted_BM_EdgeR_avsd_srep)
mdf_toplot <- melted_BM_EdgeR_avsd_srep
    mdf_toplot$Condition <- factor(mdf_toplot$Condition,levels = c(conds_av_B_iPSES[-8],conds_av_M_CloneLast))
    mdf_toplot$Sample <- factor(mdf_toplot$Sample,levels = c(conds_B_iPSES[-c(15:16)],conds_av_M_CloneLast))
    mdf_toplot$Time=rep(NA,nrow(melted_BM_EdgeR_avsd_srep))
mdf_toplot <- mdf_toplot %>%
  mutate(Time = replace(Time, Condition=="B_Bcells", "time01")) %>%
  mutate(Time = replace(Time, Condition=="B_Bpulse", "time02")) %>%
  mutate(Time = replace(Time, Condition=="B_Day2", "time03")) %>%
  mutate(Time = replace(Time, Condition=="B_Day4", "time04")) %>%
  mutate(Time = replace(Time, Condition=="B_Day6", "time05")) %>%
  mutate(Time = replace(Time, Condition=="B_Day8", "time07")) %>%
  mutate(Time = replace(Time, Condition=="B_iPS", "time08")) %>%
  mutate(Time = replace(Time, Condition=="M_Day0", "time01")) %>%
  mutate(Time = replace(Time, Condition=="M_Day4", "time03")) %>%
  mutate(Time = replace(Time, Condition=="M_Day7", "time04")) %>%
  mutate(Time = replace(Time, Condition=="M_Day10", "time05")) %>%
  mutate(Time = replace(Time, Condition=="M_Day15", "time06")) %>%
  mutate(Time = replace(Time, Condition=="M_Day20", "time07")) %>%
  mutate(Time = replace(Time, Condition=="M_Clone", "time08"))
  head(mdf_toplot) 


## PLOT SCALED EXPRESSION OF pos&neg REGULATORS of each cluster
setwd("~/Desktop/")
Plot_GEx_inBM(mdf_avsd = mdf_toplot, toplot= c("Cpsf3, Hnrnpul1"), x_labels = time_PCA)




#
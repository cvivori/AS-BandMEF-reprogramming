#----- Packages required
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("org.Mm.eg.db","org.Hs.eg.db", "AnnotationDbi","GenomicFeatures","GO.db"))
# BiocManager::install(c("tidyverse","reshape2", "Mfuzz","RColorBrewer","viridis"))
#install.packages("plotrix")
require(tidyverse)
require(reshape2)
require(viridis)
#
#----- Dependencies
message("Dependencies:
        01_NumberEvents_CEBPaNEW_Cieply -> B_Numbers & M_Numbers,
                                        -> B_EV & M_EV")
#
#
#----- Functions to set
Plot_AreaPlot <- function(table_numberevents,condition_compared,events=c("all","wrap"),suffix,pdf_output=c("yes","no")) {
  data <- table_numberevents
  
  p <- ggplot(data, aes(x=Condition2_numbers, y=Number, fill=Event_type)) +
    scale_x_continuous(breaks=data$Condition2_numbers,labels=data$Condition2,
                       name=paste("Condition (compared to ",condition_compared,")",sep="")) +
    scale_y_continuous(breaks=seq(0,10000,1000),name="Number of events (cumulative)") +
    geom_area() +
    theme_minimal() +
    theme (text = element_text(color="grey20",size=11),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(angle=45,hjust=1,vjust=1),
           legend.position = "bottom", legend.title = element_text(face="bold"),
           panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
           panel.grid.major.y=element_line(linetype = "dashed",color="gray80")) +
    if (events=="wrap") { scale_fill_viridis(option="D",discrete = T,direction = 1,begin=0,end=1, name = "Event type") 
    }  else if (events=="all")  {   
      scale_fill_manual(values = c(magma(2,begin =.55,end =.7),"gray70",  #plasma(1,begin=0.95),
                                   viridis(6,begin=0,end=0.8)), name = "Event type") }
  if (pdf_output == "yes") { p;
    ggsave(plot = p, paste("AreaPlot_",events,"_",suffix,".pdf",sep = ""), width=7,height=6,device = cairo_pdf) 
  } else { p }
}

Plot_StackedBarPlot <- function(table_numberevents,condition_compared,events=c("all","wrap"),suffix,pdf_output=c("yes","no")) {
  data <- table_numberevents
  
  p <- ggplot(data, aes(x=Condition2_numbers, y=Number, fill=Event_type)) +
    scale_x_continuous(breaks=data$Condition2_numbers,labels=data$Condition2,
                       name=paste("Condition (compared to ",condition_compared,")",sep="")) +
    scale_y_continuous(breaks=seq(0,10000,1000),name="Number of events (cumulative)") +
    geom_bar(position="stack",stat = "identity") +
    theme_minimal() +
    theme (text = element_text(color="grey20",size=11),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(angle=45,hjust=1,vjust=1),
           legend.position = "bottom", legend.title = element_text(face="bold"),
           panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
           panel.grid.major.y=element_line(linetype = "dashed",color="gray80")) +
    if (events=="wrap") { scale_fill_viridis(option="D",discrete = T,direction = 1,begin=0,end=1, name = "Event type") 
    }  else if (events=="all")  {   
      scale_fill_manual(values = c(magma(2,begin =.55,end =.7),"gray70",  #plasma(1,begin=0.95),
                                   viridis(6,begin=0,end=0.8)), name = "Event type") }
  if (pdf_output == "yes") { p
  ggsave(paste("StackedBarPlot_",events,"_",suffix,".pdf",sep = ""), width=7,height=6,device = cairo_pdf)
  } else { p }
}

Plot_Heatmap_NumberEvents <- function(table_numberevents,suffix){
  for (e in unique(table_numberevents$Event_type)) {
    data <- table_numberevents %>%
      filter(Event_type == e)
    p <-  ggplot(data, aes(x=Condition2, y=Condition1, fill=Number)) +
      facet_grid(~Event_type,scales="free")+
      geom_tile() +
      scale_fill_viridis(option = "D",limits=c(0,4100),guide=guide_colorbar(title = "Number of events",raster = F,nbin=40)) +
      theme_minimal()  + 
      theme(axis.text.x = element_text(size = 8,angle = 30,hjust=1,vjust=1),
            axis.text.y = element_text(size = 8,angle = 30,hjust=1,vjust=1),
            axis.title= element_blank(),axis.ticks=element_blank(),
            legend.title = element_text(size = 7))
    
    plot_tab = ggplot_build(p)$data[[1]]
    #plot_tab = plot_tab[order(plot_tab$group),]
    tabp = data
    tabp$xcoord=apply(plot_tab, 1, function(x) mean(as.numeric(x[c("xmax","xmin")])))
    tabp$ycoord=apply(plot_tab, 1, function(x) mean(as.numeric(x[c("ymax","ymin")])))
    plt <- p+ geom_text(data=tabp, aes(x= xcoord, y=ycoord,
                             label=Number, 
                             hjust=0.5,vjust=0.5,size=7),
              position = position_dodge(width=1),size=3,color="white")
    ggsave(plot = plt, paste("Heatmap_Number",e,"_",suffix,".pdf",sep = ""), width=7,height=6,device = cairo_pdf)
    }
}
#


## dPSI10
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/NumberEvents/")

## NUMBER OF EVENTS dPSI10 WRAPPED
data_Bcells_dPSI10_wrap <- melt(B_diffEV_comparisons_EV$dPSI10, variable.name = "Event_type", value.name = "Number") %>% 
  filter(Condition1 == "B_Bcells")
data_Bcells_dPSI10_wrap$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_wrap$Condition2, "B_Day","")
data_Bcells_dPSI10_wrap$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_wrap$Condition2_numbers, "B_Bpulse","0")
data_Bcells_dPSI10_wrap$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_wrap$Condition2_numbers, "B_iPS","10")
data_Bcells_dPSI10_wrap$Condition2_numbers <- as.numeric(str_replace_all(data_Bcells_dPSI10_wrap$Condition2_numbers, "B_ES","12"))
data_Bcells_dPSI10_wrap$Event_type <- factor(data_Bcells_dPSI10_wrap$Event_type, levels= c("Alt3","Alt5","RI","MIC","CEx"))
head(data_Bcells_dPSI10_wrap)
## STACKED AREA & BAR CHARTs
Plot_AreaPlot(table_numberevents = data_Bcells_dPSI10_wrap, condition_compared = "Bcells",events = "wrap",suffix = "dPSI10", pdf_output = "yes")
Plot_StackedBarPlot(table_numberevents = data_Bcells_dPSI10_wrap, condition_compared = "Bcells",events = "wrap",suffix = "dPSI10", pdf_output = "yes")


## NUMBER OF EVENTS dPSI10 ALL
data_Bcells_dPSI10_all <- melt(B_diffEV_comparisons_Numbers$dPSI10, variable.name = "Event_type", value.name = "Number") %>% 
  filter(Condition1 == "B_Bcells") 
data_Bcells_dPSI10_all$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_all$Condition2, "B_Day","")
data_Bcells_dPSI10_all$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_all$Condition2_numbers, "B_Bpulse","0")
data_Bcells_dPSI10_all$Condition2_numbers <- str_replace_all(data_Bcells_dPSI10_all$Condition2_numbers, "B_iPS","10")
data_Bcells_dPSI10_all$Condition2_numbers <- as.numeric(str_replace_all(data_Bcells_dPSI10_all$Condition2_numbers, "B_ES","12"))
data_Bcells_dPSI10_all$Event_type <- factor(data_Bcells_dPSI10_all$Event_type, levels= c("Alt3","Alt5","IR","S","C1","C2","C3","MIC","ANN"))
head(data_Bcells_dPSI10_all)
## STACKED AREA & BAR CHARTs
Plot_AreaPlot(table_numberevents = data_Bcells_dPSI10_all, condition_compared = "Bcells",events = "all",suffix = "dPSI10", pdf_output = "yes")
Plot_StackedBarPlot(table_numberevents = data_Bcells_dPSI10_all, condition_compared = "Bcells",events = "all",suffix = "dPSI10", pdf_output = "yes")

## HEATMAPs
data_B_dPSI10_wrap <- melt(B_diffEV_comparisons_EV$dPSI10, variable.name = "Event_type", value.name = "Number")
Plot_Heatmap_NumberEvents(table_numberevents = data_B_dPSI10_wrap,suffix = "dPSI10")




require(tidyverse)
require(VennDiagram)


## PSI tables divided by type of events
lapply(B_diffEV_list$dPSI10, dim)
  head(B_diffEV_list$dPSI10$CEx)
lapply(M_diffEV_list$dPSI10, dim)
  head(M_diffEV_list$dPSI10$CEx)

## OVERLAP EVENTS FROM B cell AND MEFs REPROGRAMMING
  # a1 = elements unique of the 1st set, a2 = elements unique of the 2nd set, a3 = intersection of the 2 sets
  overlaps_BM <- list()
for (t in names(B_diffEV_list$dPSI10)) {
  overlaps_BM[[t]] <- calculate.overlap(list(B_diffEV_list$dPSI10[[t]]$EVENT,M_diffEV_list$dPSI10[[t]]$EVENT))
}
  lapply(overlaps_BM, function(x) summary(x))

## PERCENTAGE OF MEFs EVENTS overlapping
perc_BM <- lapply(overlaps_BM, function(x) round(100*length(x$a3)/length(x$a2),digits = 2))

## DRAW VENNS
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/Overlap_MEFs_Venns/")
# for (t in names(B_EV_list$dPSI10)) {
# pdf(file=paste("VennDiagram_",t,"_BandMEFstoiPS_",perc_BM[[t]],"perc.pdf",sep = ""),width = 8, height=8)
#   draw.pairwise.venn( area1 = length(overlaps_BM[[t]]$a1),
#                       area2 = length(overlaps_BM[[t]]$a2),
#                       cross.area = length(overlaps_BM[[t]]$a3),
#                       category = c("B cell reprogramming","MEFs reprogramming"),
#                       lwd=0,
#                       fill=c("darkslategray","maroon"),
#                       cat.fontface = "bold",
#                       cat.fontfamily = "Helvetica",
#                       cat.pos =0,
#                       #print.mode=c('raw','percent'),
#                       fontfamily = "Helvetica"
#   )
#   dev.off()
# }




## OVERLAP OF B cell CLUSTERS with MEFs AS EVENTS 
# from 04_Clustering_VTS_CEx.R
summary(cl_memb_VTS)
summary(cl_events_VTS)
summary(cl_genes_VTS) 
summary(cl_Data_VTS)

num_events <- unlist(lapply(cl_events_VTS, function(x) length(x)))
ov_events <- lapply(cl_events_VTS, function(x) sapply(x, function(y) is.element(y,rownames(M_PSIs_VTS_av_list$dPSI10$CEx))))
sum_ov_events <- unlist(lapply(ov_events, function(x) sum(x)))

overlap <- data.frame(cluster = c(1:i),
                      num_events = num_events,
                      overlap_events = sum_ov_events)
overlap$overlap_percentage <- overlap$overlap_events/overlap$num_events
overlap$cluster <- factor(overlap$cluster, levels=c(4,10,3,9,1,12,11,2,5,6,7,8))
overlap$name_cluster <- c("MidDOWN","iPSDO","EarlyDOWN","PulseUP","D8UP","ESDO","D8UP","D8DO","MidUP","EarlyUP","TransEarly","Late")

moverlap <- melt(overlap)
head(moverlap)


### CALCULATE AVERAGE OF ALL CLUSTERS
av_all <- mean(overlap$overlap_percentage)

## PLOT BARPLOT
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/Overlap_MEFs2iPS/")
ggplot(subset(moverlap, variable == "overlap_percentage"), aes(x=cluster,y=value*100)) +
  geom_bar(stat="identity") +
  scale_x_discrete() +
  scale_y_continuous(limits=c(0,50)) +
  ylab(label = "Percentage of overlap with MEFs reprogramming") +
  geom_hline(yintercept = av_all*100, linetype=2, color="maroon",show.legend = T) +
  geom_text(data = subset(moverlap, variable == "num_events"),
            aes(x=cluster,y=1, label=value, 
                hjust=0.5,vjust=0.5),
            position = position_dodge(width=1), size=2,color="white") +
  geom_text(data =  subset(moverlap, variable == "num_events"),
            aes(x=cluster,y=2, label=name_cluster, 
                hjust=0.5,vjust=0.5),
            position = position_dodge(width=1), size=2,color="white") +
  theme_classic() +
  theme (text = element_text(color="grey20",size=18),
         title = element_text(color="grey20",size=11,face="bold"),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave(paste("OverlapMEFs_Percentage.pdf",sep = ""), width=7,height=6,device = cairo_pdf)





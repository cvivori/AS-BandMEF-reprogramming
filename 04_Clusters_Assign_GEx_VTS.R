require(tidyverse)

## FUNCTIONs
## ClustersMembers - Calculates mamberships scores for elements of lists
ClustersMembers <- function (Zscores_matrixes_list, clusters_number, clusters_centroids)   # Zscores_matrixes_list = list of matrices containing GEx values 
{                                                                                      # clusters_number = number of clusters (e.g. i)
  mem_list=list()                                                                      # clusters_centroids = e.g. cl2$centers or -cl2$centers
  for (j in 1:clusters_number) {
    if (nrow(Zscores_matrixes_list[[j]])==0) next
    mem_temp=membership(Zscores_matrixes_list[[j]], clusters=clusters_centroids, m=m)
    mem_list[[j]]=mem_temp
  }
  return(mem_list)
}
## MembersOnly_M - Returns a list containing members above a minimum membership threshold to each cluster:   
MembersOnly <- function (clusters_members_list, clusters_number, minmem)
{ 
  only_list=list()
  for (j in 1:clusters_number) {
    if (is.null(clusters_members_list[[j]])>0) next
    only_temp=which(clusters_members_list[[j]][,j]>minmem)
    only_list[[j]]=clusters_members_list[[j]][only_temp,j]
  }
  return(only_list)
}




## EXTRACT EXPRESSION OF GENES FROM EACH CLUSTER
head(B_CPMs_av)
head(cl_genes_VTS)

cl_ENSGIDs_VTS_CEx <- lapply(cl_genes_VTS, function(x) ENS2Gene[which(ENS2Gene$GeneName %in% x),"ENSGID"])
cl_leftout_genes_CEx <- lapply(cl_genes_VTS, function(x) setdiff(x, ENS2Gene$GeneName))

cl_GEx_CEx <- lapply(cl_ENSGIDs_VTS_CEx, function(x) B_CPMs_av[intersect(x,rownames(B_CPMs_av)),])
  
# Standardise to Z SCOREs
cl_GEx_CEx_Z= lapply(cl_GEx_CEx, function(x) t(scale(t(x))))


## ATTRIBUTE MEMBERSHIP TO SCALED GEx of GENES (to CEx clusters dPSI10)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/Overlap_GEx_CEx/")
# Assess membership of ZScores of CLUSTER Gene Expression to each cluster (CONCORDANT & DISCORDANT gene exp)
Cl_GE_memberships_conc=ClustersMembers(cl_GEx_CEx_Z, clusters_number = i, clusters_centroids = cl2_VTS$centers)
Cl_GE_memberships_disc=ClustersMembers(cl_GEx_CEx_Z, clusters_number = i, clusters_centroids = -cl2_VTS$centers)
# Select the Genes that show a minmem0.3 in their own cluster
minmemGEx=0.3   # Set the minimum membership threshold
Cl_GE_membersonly_conc=MembersOnly(Cl_GE_memberships_conc, clusters_number = i, minmem =minmemGEx)
Cl_GE_membersonly_disc=MembersOnly(Cl_GE_memberships_disc, clusters_number = i, minmem =minmemGEx)
# Calculate number of genes passing the threshold (should be lower than 30%!)
N_conc=unlist(lapply(Cl_GE_membersonly_conc,length))
N_disc=unlist(lapply(Cl_GE_membersonly_disc,length))
# Calculate percentage of genes passing the threshold (should be lower than 30%!)
Percentage_conc=unlist(lapply(Cl_GE_membersonly_conc,length))/unlist(lapply(Cl_GE_memberships_conc,length))
Percentage_conc
Percentage_others=rep(1,12)  
Percentage_others=Percentage_others-Percentage_conc
Percentage_others

Percentage_disc=unlist(lapply(Cl_GE_membersonly_disc,length))/unlist(lapply(Cl_GE_memberships_disc,length))
Percentage_disc
Percentage_others=rep(1,12)  
Percentage_others=Percentage_others-Percentage_disc
Percentage_others


## MAKE DATAFRAME TO PLOT
overlap_GEx <- data.frame(cluster = c(1:12),
                      num_genes_considered = unlist(lapply(cl_ENSGIDs_VTS_CEx, function(x) length(unique(x)))),
                      num_concordant_GEx = N_conc,
                      num_discordant_GEx = N_disc,
                      perc_concordant_GEx = Percentage_conc,
                      perc_discordant_GEx = Percentage_disc,
                      perc_others = 1-Percentage_conc-Percentage_disc)
overlap_GEx$cluster <- factor(overlap_GEx$cluster, levels=c(4,10,3,9,1,12,11,2,5,6,7,8))
overlap_GEx$name_cluster <- c("MidDOWN","iPSDO","EarlyDOWN","PulseUP","D8UP","ESDO","D8UP","D8DO","MidUP","EarlyUP","TransEarly","Late")
overlap_GEx

moverlap_GEx <- melt(overlap_GEx)
moverlap_GEx$variable <- factor(moverlap_GEx$variable, levels =c("num_genes_considered","num_discordant_GEx","num_concordant_GEx","perc_others","perc_discordant_GEx","perc_concordant_GEx"))
head(moverlap_GEx)

### select only relevant clusters and calculate average of all clusters
in_clusters <- c(4,10,3,9,1,12)
out_clusters <- c(11,2,5,6,7,8)
order_clusters_3x4 <- c(10,3,5,7,9,1,8,6,4,12,11,2)
av_all_conc <- mean(overlap_GEx$perc_concordant_GEx)
av_all_disc <- mean(overlap_GEx$perc_discordant_GEx)
av_others_conc <- mean(overlap_GEx[out_clusters,"perc_concordant_GEx"])
av_others_disc <- mean(overlap_GEx[out_clusters,"perc_discordant_GEx"])



## PLOT STACKED BARPLOT
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/Overlap_GEx_CEx/")
p <- ggplot(subset(moverlap_GEx, grepl("perc_", variable)), aes(x=cluster,y=value*100,fill=variable)) +
  geom_bar(stat="identity") +
  scale_x_discrete() 
  #scale_fill_manual(values = c("grey80",viridis(4)[2:3])) 
  
plot_tab = ggplot_build(p)$data[[1]]
plot_tab = plot_tab[order(plot_tab$group),]
moverlap_GEx$coord=apply(plot_tab, 1, function(x) mean(as.numeric(x[c("ymax","ymin")])))
# lines <- subset(mdfc, Conditions == control_cond)$Fraction_of_events_cum
# lines <- c(1-lines,1)

p +  ylab(label = "Percentage of events with conc/discordant gene expression") +
     geom_hline(yintercept = av_others_conc*100, linetype=1, color="maroon",show.legend = T) +
     geom_hline(yintercept = (av_others_conc*100 + av_others_disc*100), linetype=1, color="maroon",show.legend = T) +
    geom_hline(yintercept = av_all_conc*100, linetype=2, color="maroon",show.legend = T) +
    geom_hline(yintercept = (av_all_conc*100 + av_all_disc*100), linetype=2, color="maroon",show.legend = T) +
       geom_text(data = subset(moverlap_GEx, variable == "num_genes_considered"),
            aes(x=cluster,y=95,label=value,hjust=0.5,vjust=0.5),
            position = position_dodge(width=1), size=3,color="grey20",show.legend = F) +
     # geom_text(data =  subset(moverlap_GEx, variable == "num_genes_considered"),
     #        aes(x=cluster,y=5, label=name_cluster, hjust=0.5,vjust=0.5),
     #        position = position_dodge(width=1), size=2,color="white",show.legend = F) +
  scale_fill_manual(values = c("black",viridis(4)[2:3],"grey80")) +
  guides(fill = guide_legend(nrow = 2), color=F) +
  theme_classic() +
  theme (text = element_text(color="grey20",size=18),
         title = element_text(color="grey20",size=11,face="bold"),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
         legend.position = "bottom",legend.title = element_text(face="bold"),legend.text = element_text(color="grey20",size=10),
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
#ggsave(paste("OverlapGEx_Percentage.pdf",sep = ""), width=7,height=6,device = cairo_pdf)








## MAKE A BIG TABLE WITH EXPRESSION VALUES AND MEMBERSHIP VALUES OF GEx PROFILES
membersonly_conc <- melt(lapply(Cl_GE_membersonly_conc, function(x) data.frame(ENSGID=names(x),
                                                      membership=x)), level=1)
membersonly_conc <- membersonly_conc %>% dplyr::select(-variable)
colnames(membersonly_conc) <- c("ENSGID", "membership_conc", "cluster")
membersonly_disc <- melt(lapply(Cl_GE_membersonly_disc, function(x) data.frame(ENSGID=names(x),
                                                      membership=x)), level=1)
membersonly_disc <- membersonly_disc %>% dplyr::select(-variable)
colnames(membersonly_disc) <- c("ENSGID", "membership_disc", "cluster")
membersonly <- melt(cl_GEx_CEx_Z,level = 1)
colnames(membersonly) <- c("ENSGID","Condition","Scaled_CPM","cluster")

membersonly_merged <- merge(membersonly,membersonly_conc,all=T)
membersonly_merged <- merge(membersonly_merged,membersonly_disc,all=T)
membersonly_merged$number_genes <- sapply(membersonly_merged$cluster, function(x) nrow(Cl_GE_memberships_conc[[x]])) 
membersonly_merged$Facetcluster <- paste0("Cluster ",formatC(membersonly_merged$cluster, width = 2, format = "d", flag = "0"),
                                          " - n=",membersonly_merged$number_genes)
head(membersonly_merged)



ggplot(subset(membersonly_merged, is.na(membership_conc) & is.na(membership_disc)), 
       aes(x=Condition,y=Scaled_CPM, group=ENSGID)) + 
  facet_wrap(.~Facetcluster) +
  #geom_point() + 
  # scale_color_viridis(discrete = T,option="D",end = .9) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  geom_line(color="grey80") +
  geom_line(data=subset(membersonly_merged, membership_disc>=0.3), 
            aes(x=Condition,y=Scaled_CPM, group=ENSGID), color=viridis(4)[3]) +
    geom_line(data=subset(membersonly_merged, membership_conc>=0.3), 
              aes(x=Condition,y=Scaled_CPM, group=ENSGID), color=viridis(4)[2]) +
# scale_color_viridis(discrete = T,option="B",end=.9) + 
  #guides(color=F) +
  theme_classic() +
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=45,hjust=1,vjust=1),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
# ggsave("OverlapGEx_Profiles.pdf", width=7,height=6,device = cairo_pdf)






## PIE CHART OF ALL EVENTS
df <- data.frame(value=c(1-av_all_conc-av_all_disc,av_all_disc,av_all_conc),
                 group=c("Others","Concordant","Discordant"))
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0, direction=-1) +
  scale_fill_manual(values = c(viridis(4)[2:3],"grey80")) +
  theme_minimal()


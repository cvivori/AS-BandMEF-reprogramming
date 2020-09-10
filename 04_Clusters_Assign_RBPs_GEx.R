require(tidyverse)
require(ggrepel)

## FUNCTIONs
## MembersOnly_M - Returns a list containing members above a minimum membership threshold to each cluster, including the membership value   
MembersOnly_M <- function (memberships_matrix, clusters_number, minmem)   # memberships_matrix = matrix containing membership values
{ 																	                                  # clusters_number = number of clusters (e.g. i)
  only_temp=list()													                          # minmem = minimum membership threshold
  only_list=list()
  for (j in 1:clusters_number) {
    only_temp[[j]]=which(memberships_matrix[,j]>minmem)
    only_list[[j]]=noquote(names(only_temp[[j]]))
  }
  return(only_list)
}

Plot_posneg_B_regulators = function (mdf_scaled, highlight_genes="", prefix="") {
  
  
  #POSITIVE REGULATORS
  mdf_toplot <- mdf_scaled %>%
    filter(!is.na(cluster_posRBP)) %>%
    filter(Dataset=="B")
  nclusters <- length(unique(mdf_toplot$cluster_posRBP))
  ncolors <- length(unique(mdf_toplot$GeneName))
  pal <- colorRampPalette(brewer.pal(9, "Greens")[3:9])(ncolors+1)
    ggplot(mdf_toplot,
         aes(x=Sample,y=Scaled_cpm, color=GeneName, group=GeneName_Dataset)) +
    facet_wrap(~FacetposRBP) +
    geom_line(size=0.5) +
    geom_point(size=0.5) +
      geom_line(data= subset(mdf_toplot, GeneName %in% highlight_genes),
              size =1,show.legend =F) +
      geom_point(data= subset(mdf_toplot, GeneName %in% highlight_genes),
                 size =1,show.legend =F) +
      geom_label_repel(data = subset(mdf_toplot, GeneName %in% highlight_genes & Condition == "Day8"),
                       color="gray20",label = highlight_genes,nudge_x=1,nudge_y=2) +
    scale_color_manual(values = pal,guide=F) +
    theme_classic() +
    theme (text = element_text(color="grey20",size=11),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
           legend.position = "bottom", legend.title = element_text(face="bold"),
           panel.grid.minor.x=element_blank(),
           panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
           strip.background = element_blank(),strip.text = element_text(face="bold"))
  ggsave(paste(prefix,"ScaledCPM_posRBP_minmem",minmemRBP,"_cl",nclusters,".pdf",sep = ""), width=5.88,height=3.6,device = cairo_pdf)
  
  #NEGATIVE REGULATORS
  mdf_toplot <- mdf_scaled %>%
    filter(!is.na(cluster_negRBP)) %>%
    filter(Dataset=="B") 
  nclusters <- length(unique(mdf_toplot$cluster_negRBP))
  ncolors <- length(unique(mdf_toplot$GeneName))
  pal <- colorRampPalette(brewer.pal(9, "OrRd")[3:9])(ncolors+1)
  ggplot(mdf_toplot,
         aes(x=Sample,y=Scaled_cpm, color=GeneName, group=GeneName_Dataset)) +
    facet_wrap(~FacetnegRBP) +
    geom_line(size=0.5) +
    geom_point(size=0.5) +
    geom_line(data= subset(mdf_toplot, GeneName %in% highlight_genes),
              size =1,show.legend =F) +
    geom_point(data= subset(mdf_toplot, GeneName %in% highlight_genes),
               size =1,show.legend =F) +
    geom_label_repel(data = subset(mdf_toplot, GeneName %in% highlight_genes & Condition == "Day8"),
                     color="gray20",label = highlight_genes,nudge_x=1,nudge_y=2) +
    scale_color_manual(values = pal,guide=F) +
    theme_classic() +
    theme (text = element_text(color="grey20",size=11),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
           legend.position = "bottom", legend.title = element_text(face="bold"),
           panel.grid.minor.x=element_blank(),
           panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
           strip.background = element_blank(),strip.text = element_text(face="bold"))
  ggsave(paste(prefix,"ScaledCPM_negRBP_minmem",minmemRBP,"_cl",nclusters,".pdf",sep = ""), width=5.88,height=3.6,device = cairo_pdf)
  
}


## EXTRACT EXPRESSION OF GENES CODING FOR RBPs
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/RBP_membership/")
gRBPs_B_CPMs_av_fil_scaled <- RBPs_B_CPMs_av_fil_scaled 
rownames(gRBPs_B_CPMs_av_fil_scaled) <- ENS2Gene[which(ENS2Gene$ENSGID %in% rownames(RBPs_B_CPMs_av_fil_scaled)),"GeneName"]
gRBPs_B_CPMs_av_fil <- RBPs_B_CPMs_av_fil
rownames(gRBPs_B_CPMs_av_fil_scaled) <- ENS2Gene[which(ENS2Gene$ENSGID %in% rownames(RBPs_B_CPMs_av_fil)),"GeneName"]
gSPLs_B_CPMs_av_fil_scaled <- SPLs_B_CPMs_av_fil_scaled
rownames(gSPLs_B_CPMs_av_fil_scaled) <- ENS2Gene[which(ENS2Gene$ENSGID %in% rownames(SPLs_B_CPMs_av_fil_scaled)),"GeneName"]

head(gRBPs_B_CPMs_av_fil_scaled)
head(gSPLs_B_CPMs_av_fil_scaled)


## ATTRIBUTE MEMBERSHIP TO SCALED GEx of RBPs (to CEx clusters dPSI10)
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/RBP_membership/")
RBPs_memberships_pos=membership(na.omit(gRBPs_B_CPMs_av_fil_scaled), clusters=cl2_VTS$centers, m=m)    		# Attributes memberships of ZScores of RBPs Gene Expression 
  # write.csv(RBPs_memberships_pos, file="RBP_memberships_pos.csv")               # in the AS events clusters generated previously
RBPs_memberships_neg=membership(na.omit(gRBPs_B_CPMs_av_fil_scaled), clusters=-cl2_VTS$centers, m=m) 	    # in both directions (POSITIVE and NEGATIVE regulators) 
  # write.csv(RBPs_memberships_neg, file="RBP_memberships_neg.csv")
SPLs_memberships_pos=membership(na.omit(gSPLs_B_CPMs_av_fil_scaled), clusters=cl2_VTS$centers, m=m)    		# Attributes memberships of ZScores of SPLs Gene Expression 
  # write.csv(SPLs_memberships_pos, file="SPL_memberships_pos.csv")         	    # in the AS events clusters generated previously
SPLs_memberships_neg=membership(na.omit(gSPLs_B_CPMs_av_fil_scaled), clusters=-cl2_VTS$centers, m=m) 		# in both directions (POSITIVE and NEGATIVE regulators)
  # write.csv(SPLs_memberships_neg, file="SPL_memberships_neg.csv")

## MAKE LISTS of POSITIVE and NEGATIVE REGULATORS of each CLUSTER
minmemRBP=0.3   #Set minimum membership threshold for members assignment and plotting
RBP_members_pos_cl=MembersOnly_M(RBPs_memberships_pos, clusters_number = i, minmem = minmemRBP)   	# RBP_members_pos_cl is a LIST containing RBPs that positively correlate with each cluster (>minmem)
RBP_members_neg_cl=MembersOnly_M(RBPs_memberships_neg, clusters_number = i, minmem = minmemRBP) 		# RBP_members_neg_cl is a LIST containing RBPs that negatively correlate with each cluster (>minmem) 
SPL_members_pos_cl=MembersOnly_M(SPLs_memberships_pos, clusters_number = i, minmem = minmemRBP) 		# SPL_members_pos_cl is a LIST containing SPLs that positively correlate with each cluster (>minmem) 
SPL_members_neg_cl=MembersOnly_M(SPLs_memberships_neg, clusters_number = i, minmem = minmemRBP) 		# SPL_members_neg_cl is a LIST containing SPLs that negatively correlate with each cluster (>minmem)

## CHECK NUMBERS
numbers_regulators <- data.frame(cluster=c("MiddleDOWN","iPS__","EarlyDOWN","PulseUP","Day8UP__","ESDOWN__","D82ESUP__","Day8DOWN__","MiddleUP","EarlyUP","TransientDOWN","Late"),
                                 posRBP=unlist(lapply(RBP_members_pos_cl, function(x) length(x))),
                                 negRBP=unlist(lapply(RBP_members_neg_cl, function(x) length(x))))
numbers_regulators

## MAKE 1 BIG TABLE WITH ALL REGULATORS
tmp_Rp <- melt(RBP_members_pos_cl,level=1);  colnames(tmp_Rp) <- c("GeneName","cluster_posRBP")
tmp_Rn <- melt(RBP_members_neg_cl,level=1);  colnames(tmp_Rn) <- c("GeneName","cluster_negRBP")
tmp_Sp <- melt(SPL_members_pos_cl,level=1);  colnames(tmp_Sp) <- c("GeneName","cluster_posSPL")
tmp_Sn <- melt(SPL_members_neg_cl,level=1);  colnames(tmp_Sn) <- c("GeneName","cluster_negSPL")
Regulators <- merge(tmp_Rp,tmp_Rn,all=T,by="GeneName")
  Regulators <- merge(Regulators,tmp_Sp,all=T,by="GeneName")
  Regulators <- merge(Regulators,tmp_Sn,all=T,by="GeneName")
  head(Regulators) 
  dim(Regulators)
Reg_BM_CPMs <- merge(Regulators,melted_BM_EdgeR_scaled,by="GeneName",all=T)
Reg_BM_CPMs$Sample <- factor(Reg_BM_CPMs$Sample,levels = c(conds_av_B_iPSES,conds_av_M_CloneLast))
Reg_BM_CPMs$FacetposRBP <- paste0("Cluster ",formatC(Reg_BM_CPMs$cluster_posRBP, width = 2, format = "d", flag = "0"),
                                  " - n=",numbers_regulators[Reg_BM_CPMs$cluster_posRBP,"posRBP"])
Reg_BM_CPMs$FacetnegRBP <- paste0("Cluster ",formatC(Reg_BM_CPMs$cluster_negRBP, width = 2, format = "d", flag = "0"),
                                  " - n=",numbers_regulators[Reg_BM_CPMs$cluster_negRBP,"negRBP"])
  head(subset(Reg_BM_CPMs,cluster_negRBP ==6)) 

## PLOT SCALED EXPRESSION OF pos&neg REGULATORS of each cluster
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/RBP_membership/")
Plot_posneg_B_regulators(Reg_BM_CPMs,prefix = "", highlight_genes = c("Cpsf3","Tia1","Hnrnpul1"))



## MAKE TABLES FOR WORDCLOUDS (Tagul/Wordart)
# setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/RBP_membership/WordClouds/")
# for (x in c(1:i)) {
#   memb <-  data.frame(Genes=RBP_members_pos_cl[[x]],
#                       membership= RBPs_memberships_pos[RBP_members_pos_cl[[x]],x])
#                       #ProtvsGEx = ProtvsGEx[RBP_members_pos_cl[[x]]])
#   memb_sort <- memb %>% arrange(desc(membership))
#   memb_sort$vcol <- rep("DimGray",nrow(memb_sort))
#   #  memb_sort[which(memb_sort$ProtvsGEx =="discordant" ),"vcol"] <- "LightGray"
#     memb_sort[which(memb_sort$Genes %in% intersect(memb_sort$Genes, SPL_members_pos_cl[[x]]) ),"vcol"] <- "ForestGreen"
#     memb_sort[which(memb_sort$Genes %in% intersect(memb_sort$Genes, SPL_members_pos_cl[[x]])
#                     & memb_sort$ProtvsGEx =="discordant" ),"vcol"] <- "LightGreen"
#   write.csv(memb_sort,file=paste("RBPSPL_members_pos_cl",x,"_mimmemRBP",minmemRBP,"_WordCloud.csv"),quote = F,row.names = F)
# 
#   memb <-  data.frame(Genes=RBP_members_neg_cl[[x]],
#                       membership= RBPs_memberships_pos[RBP_members_neg_cl[[x]],x])
#                       #ProtvsGEx = ProtvsGEx[RBP_members_neg_cl[[x]]])
#   memb_sort <- memb %>% arrange(desc(membership))
#   memb_sort$vcol <- rep("DimGray",nrow(memb_sort))
#  # memb_sort[which(memb_sort$ProtvsGEx =="discordant" ),"vcol"] <- "LightGray"
#   memb_sort[which(memb_sort$Genes %in% intersect(memb_sort$Genes, SPL_members_neg_cl[[x]]) ),"vcol"] <- "OrangeRed"
#   memb_sort[which(memb_sort$Genes %in% intersect(memb_sort$Genes, SPL_members_neg_cl[[x]])
#                   & memb_sort$ProtvsGEx =="discordant" ),"vcol"] <- "Orange"
#   write.csv(memb_sort,file=paste("RBPSPL_members_neg_cl",x,"_mimmemRBP",minmemRBP,"_WordCloud.csv"),quote = F,row.names = F)
# }


# ## MAKE TABLEs FOR SUPPLEMENTARIES
# # RBP_members_pos_toprint=lapply(RBP_members_pos_cl, function(x) as.data.frame(gRBPs_fil_B_CPMs_av[x,]))
#   # RBP_members_neg_toprint=lapply(RBP_members_neg_cl, function(x) as.data.frame(gRBPs_fil_B_CPMs_av[x,]))
#   RBP_members_pos_toprint=RBP_members_neg_toprint=list()
#   for (c in c(1:i)) { 
#     if (length(RBP_members_pos_cl[[c]])>1) {
#       RBP_members_pos_toprint[[c]] = as.data.frame(gRBPs_fil_B_CPMs_av[RBP_members_pos_cl[[c]],])
#     } else if (length(RBP_members_pos_cl[[c]])==1) {
#       RBP_members_pos_toprint[[c]] = as.data.frame(t(gRBPs_fil_B_CPMs_av[RBP_members_pos_cl[[c]],])) } 
#     RBP_members_pos_toprint[[c]]$Cluster = rep(c,nrow(RBP_members_pos_toprint[[c]]))
#     RBP_members_pos_toprint[[c]]$Membership = RBPs_memberships_pos[RBP_members_pos_cl[[c]],c]
#     if (length(RBP_members_neg_cl[[c]])>1) {
#       RBP_members_neg_toprint[[c]] = as.data.frame(gRBPs_fil_B_CPMs_av[RBP_members_neg_cl[[c]],])
#     } else if (length(RBP_members_neg_cl[[c]])==1) {
#       RBP_members_neg_toprint[[c]] = as.data.frame(t(gRBPs_fil_B_CPMs_av[RBP_members_neg_cl[[c]],])) } 
#     RBP_members_neg_toprint[[c]]$Cluster = rep(c,nrow(RBP_members_neg_toprint[[c]]))
#     RBP_members_neg_toprint[[c]]$Membership = RBPs_memberships_neg[RBP_members_neg_cl[[c]],c]
#   }
# # setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/RBP_membership/")
# # write.table(x = do.call(rbind,RBP_members_pos_toprint), file = "RBPs_pos_CPMs_ClusterMembership.txt", quote=F, row.names = T,col.names=NA,sep = "\t")
# # write.table(x = do.call(rbind,RBP_members_neg_toprint), file = "RBPs_neg_CPMs_ClusterMembership.txt", quote=F, row.names = T,col.names=NA,sep = "\t")
  

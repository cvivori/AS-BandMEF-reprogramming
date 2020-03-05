require(ggdendro)
require(ggplot2)
require(viridis)
require(grid)


### MAKE HEATMAP ON AS ON BOTH DATASETS (RAW + SCALED)
## TWO WITH HCLUST ON BOTH EVENTS AND SAMPLES
## TWO WITH ONLY EVENTS CLUSTERED

### Merge in a big table
tmpB <- B_PSIs_VTS_list$dPSI10$CEx
colnames(tmpB) <- conds_B_iPSES
tmpM <- M_PSIs_VTS_list$dPSI10$CEx
colnames(tmpM) <- conds_M_CloneLast

BM_PSIs_VTS10_CEx <- merge(tmpB,tmpM,by="row.names")
rownames(BM_PSIs_VTS10_CEx) <- BM_PSIs_VTS10_CEx$Row.names
BM_PSIs_VTS10_CEx <- BM_PSIs_VTS10_CEx %>% dplyr::select(-Row.names)
mat=t(scale(t((BM_PSIs_VTS10_CEx)))); values = "scaledPSI"
dim(mat)
head(mat)


### Merge in a big table (Average values)
tmpB <- B_PSIs_VTS_av_list$dPSI10$CEx
colnames(tmpB) <- conds_av_B_iPSES
tmpM <- M_PSIs_VTS_av_list$dPSI10$CEx
colnames(tmpM) <- conds_av_M_CloneLast

BM_PSIs_VTS10_av_CEx <- merge(tmpB,tmpM,by="row.names")
rownames(BM_PSIs_VTS10_av_CEx) <- BM_PSIs_VTS10_av_CEx$Row.names
BM_PSIs_VTS10_av_CEx <- BM_PSIs_VTS10_av_CEx %>% dplyr::select(-Row.names)
mat=t(scale(t((BM_PSIs_VTS10_av_CEx)))); values = "scaledPSI"
dim(mat)
head(mat)


## Create dendrograms and sort accordingly
col_clust=hclust( dist(t(mat), method = "euclidean"), method = "ward.D2" )
col_dendro=as.dendrogram(col_clust)
plot(col_clust,cex=0.5)
col_ord=order.dendrogram(col_dendro)
col_sorted=mat[,col_ord]
col_ddata=dendro_data(col_dendro)

row_clust=hclust( dist((mat), method = "euclidean"), method = "ward.D2" )
row_dendro=as.dendrogram(row_clust)
#plot(row_clust,cex=0.5)
row_ord=order.dendrogram(row_dendro)
row_sorted=mat[row_ord,]
double_sorted=col_sorted[row_ord,]
row_ddata=dendro_data(row_dendro)


## Create mdf to plot heatmaps with row dendrogram only
MBtoiPS=as.data.frame(row_sorted)
MBtoiPS$Gene=rownames(row_sorted)
MBtoiPS$Gene=with(MBtoiPS,factor(Gene,levels=Gene,ordered=T))
MBtoiPS_m=melt(MBtoiPS, id.vars = "Gene", variable.name = "Conditions")   #measure.vars = c("Day0","Day4","Day7","Day10","Day15","Day20","Clone"),
head(MBtoiPS_m)


## Plot heatmap
p <- ggplot(MBtoiPS_m, aes(Conditions, Gene,fill = value)) + 
  geom_tile() +
  geom_vline(xintercept = 16.5, color = "white",size=2.5) +
  theme_minimal()  + 
  theme(axis.text.x = element_text(size = 8,angle = 30,hjust=1,vjust=1),
        axis.text.y = element_blank(),  #element_text(size = 6,angle = 0,hjust=1,vjust=0.5),
        axis.title= element_blank(),axis.ticks=element_blank(),
        legend.title = element_text(size = 7))

if (values== "rawPSI") {plt <- p + scale_fill_viridis(option = "D",
                                                      guide=guide_colorbar(title = "PSI",raster = F,nbin=40))
} else if (values == "scaledPSI") {plt <- p + 
  scale_fill_gradient2(high = viridis(3)[1],low = viridis(3)[2], midpoint = 0,limits=c(-4,4),
                       guide=guide_colorbar(title = "Scaled\nPSI",raster = F,nbin=40))   
} else if (values== "rawFPKM") {plt <- p + scale_fill_viridis(option = "D",
                                                              guide=guide_colorbar(title = "FPKM",raster = F,nbin=40))
} else if (values == "scaledFPKM") {plt <- p + 
  scale_fill_gradient2(high = "#3E9941",low = "#8C1A6A", midpoint = 0,limits=c(-4,4),
                       guide=guide_colorbar(title = "Scaled\nFPKM",raster = F,nbin=40))   }
plt


drows <- ggplot(segment(row_ddata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)
) + coord_flip(
) +scale_y_reverse(
) + theme_minimal(
) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/")
pdf(file=paste("Heatmap_",values,"_CEx_dPSI10_BandMEFstoiPS_nocoldendro.pdf",sep = ""),width = 16, height=12.5)
grid.newpage()
print(drows, vp=viewport(0.2, 0.925, x=0.16, y=0.448))
print(plt, vp=viewport(0.6, 0.88, x=0.55, y=0.43))
dev.off() 


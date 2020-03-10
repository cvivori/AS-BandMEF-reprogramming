require(ggdendro)
require(org.Mm.eg.db)
require(viridis)
require(grid)

## FUNCTIONS
## Distance calculation
##Euclidean_dist
dist_fun=dist
##Pearson_dist
dist2 <- function(x, ...)
  as.dist(1-(cor(t(x), method="pearson"))) 
dist_fun=dist2

## Clustering method
##Complete
clust_method= hclust
##WardD2
hclust2 <- function(x, method="ward.D2", ...)
  hclust(x, method=method, ...)    
clust_method= hclust2



## IMPORT TABLE FROM HAN 2013
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Papers/2013 Han-Blencowe MBNLs repress AS reprogr SUPP/")
HanFull=read.csv("nature12270-s4 Mouse.csv",sep=";",header=T)
rownames(HanFull) <- HanFull$NAME
#Han=read.table("Han_S4_Mouse_pvalbelow005.txt",sep="\t",header=T)
Han=subset(HanFull, C_p.value.Wilcox < 0.05)

## EXTRACT ENSGIDs
ENSMUSG_H=Han[,"NAME"]
names(ENSMUSG_H)=mapIds(org.Mm.eg.db, keys=c(gsub(",.+","",ENSMUSG_H)), column='SYMBOL', keytype='ALIAS',multiVals="first")
names(ENSMUSG_H)[which(is.na(names(ENSMUSG_H)))]=Han[which(is.na(names(ENSMUSG_H))),"NAME"]
rownames(Han)=names(ENSMUSG_H)

## ADD RBPs INVOLVED IN DEVELOPMENT / REPROGRAMMING
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/RBPs_lists/")
ENSMUSG_DEV <- read.csv("RBPs_in_Development_Reprogramming.csv",header=T,row.names = 1)[,1]
names(ENSMUSG_DEV)=mapIds(org.Mm.eg.db, keys=c(gsub(",.+","",ENSMUSG_DEV)), column='SYMBOL', keytype='ALIAS',multiVals="first")
ENSMUSG_DEV

## UNION
ENSMUSG=unique(c(names(ENSMUSG_H),ENSMUSG_DEV))
head(ENSMUSG)




### MAKE TABLE OF GENES EXPRESSED and VARIABLE ONLY IN BOTH DATASETS
gBM_CPMs_scaled <- BM_CPMs_scaled
tmp <- mapIds(org.Mm.eg.db, keys=rownames(BM_CPMs_scaled), column='SYMBOL', keytype='ENSEMBL',multiVals="first")
tmp[which(is.na(tmp))] <- names(tmp[which(is.na(tmp))])
tmp[which(duplicated(tmp))] <- names(tmp[which(duplicated(tmp))])
rownames(gBM_CPMs_scaled) <- tmp

mat <- na.omit(gBM_CPMs_scaled) ; values = "scaledCPM"
mat <- mat %>%
  rownames_to_column(var="Gene") %>%
  filter(Gene %in% ENSMUSG) %>%
  column_to_rownames(Gene)
head(mat)
dim(mat)

clust=hclust( dist(mat, method = "euclidean"), method = "ward.D2" )
dendro=as.dendrogram(clust)
# plot(clust,cex=0.5)
ord=order.dendrogram(dendro)
sorted=mat[ord,]
#group_col=group_col[ord]
ddata=dendro_data(dendro)

MBtoiPS=as.data.frame(sorted)
MBtoiPS$Gene=rownames(sorted)
MBtoiPS$Gene=with(MBtoiPS,factor(Gene,levels=Gene,ordered=T))
MBtoiPS_m=melt(MBtoiPS, id.vars = "Gene", variable.name = "Conditions")   #measure.vars = c("Day0","Day4","Day7","Day10","Day15","Day20","Clone"),
names(MBtoiPS$Conditions)=NULL
head(MBtoiPS_m)



## MAKE HEATMAP
p <- ggplot(MBtoiPS_m, aes(Conditions, Gene,fill = value)) + geom_tile()

plt <- p + scale_fill_gradient2(low = "#3E9941",mid="white",high = "#8C1A6A", midpoint = 0,limits=c(-4,4),
                                guide=guide_colorbar(title = "Scaled\nFPKM",raster = F,nbin=40)
)  + theme_minimal(
)  + theme(axis.text.x = element_text(size = 10,angle = 90,hjust=1,vjust=0.5),
           axis.text.y = element_text(size = 6,angle = 0,hjust=1,vjust=0.5),
           axis.title= element_blank(),axis.ticks=element_blank(),
           legend.title = element_text(size = 7)
) + geom_vline(xintercept = 16.5, col="white",size=2
)
plt


d <- ggplot(segment(ddata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)
) + coord_flip(
) +scale_y_reverse(
) + theme_minimal(
) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

df=data.frame(Gene=rownames(sorted), RANK=rep(0,nrow(sorted)), ESCvsDIFF=rep(0,nrow(sorted)))
rownames(df)=df$Gene
df[intersect(rownames(df),HanFull$NAME),"RANK"]=HanFull[intersect(rownames(df),HanFull$NAME),"RANK"]
df[intersect(rownames(df),HanFull$NAME),"ESCvsDIFF"]=HanFull[intersect(rownames(df),HanFull$NAME),"ESCvsDIFF"]
df


names_heatmap=rownames(sorted)

b <- ggplot (df, aes(x=Gene)) + 
  geom_bar(data = (df), 
           aes(y=ESCvsDIFF, fill=RANK), width = 0.8,
           position=position_stack(), 
           stat="identity") +
  scale_fill_viridis(option = "C",direction = -1) +
  scale_x_discrete(limits=names_heatmap) +
  scale_y_log10() +
  ylab(label = "Fold Change ES vs. Diff (log10)") +
  coord_flip() +
  #geom_hline(yintercept = 1, colour = "white",size=0.6) +
  # geom_text(data = subset(df, ESCvsDIFF <=40 & ESCvsDIFF >0),
  #           aes(y=ESCvsDIFF+0.7, label=ESCvsDIFF,
  #               hjust="left",vjust=0.5),
  #           position = position_dodge(width=1),size=3) +
  # geom_text(data = subset(df, ESCvsDIFF >40),
  #           aes(y=ESCvsDIFF, label=ESCvsDIFF,
  #               hjust="right",vjust=0.5),color="white",
  #           position = position_dodge(width=1),size=3) +
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"), #axis.text.y = element_blank(), axis.ticks.y = element_blank(),
         legend.position = "bottom", legend.title = element_text(color="grey20",size=11), legend.direction = "horizontal",
         panel.background = element_blank(), panel.grid.major = element_line(color="grey90")) 
b
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/") 
#ggsave(filename = "Han_ESvsDiff_BarPlot.pdf",plot = b,width=7,height=6,device = cairo_pdf)
pdf(file=paste("RBPs_DevHan_inBOTH_EdgeR_BARPLOT.pdf",sep = ""),width = 16, height=12.5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(60, 5)))
print(d, vp=viewport(layout.pos.row = 1:57, layout.pos.col = 1))
print(b, vp=viewport(layout.pos.row = 3:60 , layout.pos.col = 2))
print(plt, vp=viewport(layout.pos.row = 3:59, layout.pos.col = 3:5))
dev.off() 





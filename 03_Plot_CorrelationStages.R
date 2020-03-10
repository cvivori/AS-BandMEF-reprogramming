require(tidyverse)
require(ggdendro)
require(ggplot2)
require(gplots)
require(viridis)
require(grid)
require(pheatmap)

######Functions for Heatmaps: clustering and pearson index for distance correlation
## dist_correlation= c("euclidean_dist","pearson_dist")
##Euclidean_dist
dist_fun=dist
##Pearson_dist
dist2 <- function(x, ...)
  as.dist(1-(cor(t(x), method="pearson"))) 
dist_fun=dist2

## clustering_method=c("complete_method","ward_D2_method"),
##Complete
clust_method= hclust
##WardD2
hclust2 <- function(x, method="ward.D2", ...)
  hclust(x, method=method, ...)    
clust_method= hclust2
##################
#


## CORRELATION OF STAGES ACCORDING TO GEx
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/Gene_Expression_EdgeR/CorrelationStages/")
B_CPMs_av_tofil <- na.omit(B_CPMs_av)
B_CPMs_av_tofil$CoeffVar <- (apply (B_CPMs_av_tofil, 1, function(x) sd(x)/mean(x)))
head(B_CPMs_av_tofil)
summary(B_CPMs_av_tofil$CoeffVar)

## 3rd quartile
v <- subset(B_CPMs_av_tofil,CoeffVar > 0.73864)  %>% 
  dplyr::select(-CoeffVar)
CorrGEx <- cor(v)

x <- CorrGEx
quantile.range=quantile(x, probs = seq(0, 1, 0.01))     
palette.breaks=seq(quantile.range["5%"], quantile.range["95%"], 0.01)
color.palette=colorRampPalette(c("white", "#2F5F6F"))(length(palette.breaks) - 1) 

heatmap.2 (x, #numeric matrix of the values to be plotted
           distfun = dist, #"euclidean_dist"
           hclustfun = hclust, #"ward_D2_method"
           reorderfun = function (d,w) reorder (d,w=c(10,0,6,0,6,0,6,0)),
           
           # dendrogram control
           Rowv = T,    #determines if and how the row dendrogram should be reordered. 
           #By default, it is TRUE, which implies dendrogram is computed and reordered based on row means. If NULL or FALSE, then no dendrogram is computed and no reordering is done. 
           #If a dendrogram, then it is used "as-is", ie without any reordering. If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
           Colv=T,  #if(symm)"Rowv" else TRUE,  #determines if and how the column dendrogram should be reordered. Options as Rowv above + when x is a square matrix, Colv="Rowv" means that columns should be treated identically to the rows.
           dendrogram = "both", #one among c("both","row","column","none")
           
           # colors
           col= color.palette,
           
           # block sepration (optional)
           #colsep=6,   # vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor
           #sepcolor="white",
           #sepwidth=c(0.08,0.05),   #Vector of length 2 giving the width (colsep) or height (rowsep) the separator box drawn by colsep and rowsep as a function of the width (colsep) or height (rowsep) of a cell
           
           # level trace
           trace="none", #c("column","row","both","none") character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'column'.
           tracecol="gainsboro",
           
           # Row/Column Labeling
           margins = c(5, 7),  #numeric vector of length 2 containing the margins for column (lower margin) and row names (right margin), respectively
           #RowSideColors=rc,  #character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x
           #colRow = rc,   #color of row labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
           #ColSideColors,  #character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x
           #colCol = cc,  #color of column labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
           cexRow = .8,  #positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively
           cexCol = .8,
           labRow = rownames(x),  #character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
           labCol = colnames(x),
           #srtRow = 0,  #angle of row/column labels, in degrees from horizontal
           srtCol = 0,
           #adjRow = c(0,NA),  #2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation)
           adjCol = c(0.5,NA),
           offsetRow = 0.005,  #Number of character-width spaces to place between row/column labels and the edge of the plotting region.
           offsetCol = 0.4,
           
           # color key + density info
           key = TRUE,  #logical indicating whether a color-key should be shown.
           keysize = 1.1,   #numeric value indicating the size of the key
           density.info="density",  #c("histogram","density","none") character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key
           denscol="gainsboro",  #color for the density display specified by density.info, defaults to the same value as tracecol
           #symkey = any(x < 0, na.rm=TRUE) || symbreaks,  #Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise
           #densadj = 0.5,  #Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key
           key.title = NA,   # titles of the color key
           key.xlab = NA,
           key.ylab = NA,
           #key.xtickfun = NULL, #function computing tick location and labels for the xaxis of the color key. Returns a named list containing parameters that can be passed to axis
           #key.ytickfun = NULL,
           #key.par=list(),  #graphical parameters for the color key. Named list that can be passed to par
           
           # plot labels
           main = NULL,  #main, x- and y-axis titles
           xlab = NULL,
           ylab = NULL
)



## CORRELATION OF STAGES ACCORDING TO PSIs
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/CorrelationStages/")
CEx_dPSI10 <- B_PSIs_VTS_av_list$dPSI10$CEx
CEx_dPSI10$CoeffVar <- apply (CEx_dPSI10, 1, function(x) sd(x)/mean(x))
summary(CEx_dPSI10$CoeffVar)

## 3rd quartile
fil <- subset(CEx_dPSI10,CoeffVar > 0.4107)  %>% 
  dplyr::select(-CoeffVar)
CorrAS <- cor(fil)

x <- CorrAS
quantile.range=quantile(x, probs = seq(0, 1, 0.01))     
palette.breaks=seq(quantile.range["5%"], quantile.range["95%"], 0.01)
color.palette=colorRampPalette(c("white","#CD5C5C"))(length(palette.breaks) - 1) 

heatmap.2 (x, #numeric matrix of the values to be plotted
           distfun = dist, #"euclidean_dist"
           hclustfun = hclust, #"ward_D2_method"
           reorderfun = function (d,w) reorder (d,w=c(0,0,0,0,0,0,0,1)),
           
           # dendrogram control
           Rowv = T,    #determines if and how the row dendrogram should be reordered. 
           #By default, it is TRUE, which implies dendrogram is computed and reordered based on row means. If NULL or FALSE, then no dendrogram is computed and no reordering is done. 
           #If a dendrogram, then it is used "as-is", ie without any reordering. If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
           Colv=T,  #if(symm)"Rowv" else TRUE,  #determines if and how the column dendrogram should be reordered. Options as Rowv above + when x is a square matrix, Colv="Rowv" means that columns should be treated identically to the rows.
           dendrogram = "both", #one among c("both","row","column","none")
           
           # colors
           col= color.palette,
           
           # block sepration (optional)
           #colsep=6,   # vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor
           #sepcolor="white",
           #sepwidth=c(0.08,0.05),   #Vector of length 2 giving the width (colsep) or height (rowsep) the separator box drawn by colsep and rowsep as a function of the width (colsep) or height (rowsep) of a cell
           
           # level trace
           trace="none", #c("column","row","both","none") character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'column'.
           tracecol="gainsboro",
           
           # Row/Column Labeling
           margins = c(5, 7),  #numeric vector of length 2 containing the margins for column (lower margin) and row names (right margin), respectively
           #RowSideColors=rc,  #character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x
           #colRow = rc,   #color of row labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
           #ColSideColors,  #character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x
           #colCol = cc,  #color of column labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
           cexRow = 1,  #positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively
           cexCol = 1,
           labRow = rownames(x),  #character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
           labCol = colnames(x),
           #srtRow = 0,  #angle of row/column labels, in degrees from horizontal
           srtCol = 0,
           #adjRow = c(0,NA),  #2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation)
           adjCol = c(0.5,NA),
           offsetRow = 0.005,  #Number of character-width spaces to place between row/column labels and the edge of the plotting region.
           offsetCol = 0.4,
           
           # color key + density info
           key = TRUE,  #logical indicating whether a color-key should be shown.
           keysize = 1.1,   #numeric value indicating the size of the key
           density.info="density",  #c("histogram","density","none") character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key
           denscol="gainsboro",  #color for the density display specified by density.info, defaults to the same value as tracecol
           #symkey = any(x < 0, na.rm=TRUE) || symbreaks,  #Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise
           #densadj = 0.5,  #Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key
           key.title = NA,   # titles of the color key
           key.xlab = NA,
           key.ylab = NA,
           #key.xtickfun = NULL, #function computing tick location and labels for the xaxis of the color key. Returns a named list containing parameters that can be passed to axis
           #key.ytickfun = NULL,
           #key.par=list(),  #graphical parameters for the color key. Named list that can be passed to par
           
           # plot labels
           main = NULL,  #main, x- and y-axis titles
           xlab = NULL,
           ylab = NULL
)




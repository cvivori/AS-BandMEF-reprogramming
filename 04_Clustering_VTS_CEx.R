require(tidyverse)
require(reshape2)
require(Mfuzz)
require(RColorBrewer)
require(plotrix)
require(org.Mm.eg.db)   #keytypes(org.Mm.eg.db)

## FUNCTIONS
# mfuzz.plot function (adapted)
mfuzz.plot3 <- function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0.5, time.labels, 
                         time.points, ylim = c(0, 0), xlab = "", ylab = "Scaled PSI",        ### change y axis label if needed!
                         x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", 
                         col.lab = "black", col.main = "black", col.sub = "black", 
                         col = "black", Xwidth = 5, Xheight = 5, single = FALSE, ...) 
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  if (missing(colo)) {
    colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", 
              "#FF6000", "#FF7800", "#FF8F00", "#FFA700", "#FFBF00", 
              "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", 
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", 
              "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28", 
              "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", 
              "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", 
              "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", 
              "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF", 
              "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF", 
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", 
              "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", 
              "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }
  else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), 
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in 1:dim(cl[[1]])[[1]]) {
    if (single) 
      j <- single
    tmp <- exprs(eset)[clusterindex == j, ]
    tmpmem <- memship[clusterindex == j, j]
    tmpmem[ tmpmem >= min.mem ]
    sizeCluster = length(tmpmem[ tmpmem >= min.mem ])
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
      if (x11) 
        X11(width = Xwidth, height = Xheight)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      #            if (sum(ylim == c(0, 0)) == 2) 
      ylim <- c(ymin, ymax)
      if (!is.na(sum(mfrow))) {
        par(mfrow = mfrow, bg = bg, col.axis = col.axis, 
            col.lab = col.lab, col.main = col.main, col.sub = col.sub, 
            col = col)
      }
      else {
        par(bg = bg, col.axis = col.axis, col.lab = col.lab, 
            col.main = col.main, col.sub = col.sub, col = col)
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points))) 
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, las=3,
                   xlab = xlab, ylab = ylab, main = paste("Cluster ",j,", n=",sizeCluster), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points, 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col, las=3, 
             ...)
        axis(2, col = ax.col, las=1, ...)
      }
      abline(h=c(-4,-3,-2,-1,0,1,2,3,4,5,6), lty=2, col="gray")
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      #            if (sum(ylim == c(0, 0)) == 2) 
      ylim <- c(ymin, ymax)
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points))) 
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, las=3, 
                   xlab = xlab, ylab = ylab, main = paste("Cluster ",j,", n=",sizeCluster), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points, 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
             col = ax.col, las=3, ...)
        axis(2, col = ax.col, las=1, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col, las=3, 
             ...)
        axis(2, col = ax.col, las=1, ...)
      }
      abline(h=c(-4,-3,-2,-1,0,1,2,3,4,5,6), lty=2, col="gray")
    }
    if (length(tmpmem) > 0) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                     colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            if (missing(time.points)) {
              lines(tmp[tmpind[k], ], col = colo[jj])
            }
            else lines(time.points, tmp[tmpind[k], ], 
                       col = colo[jj])
          }
        }
      }
    }
    if (single) 
      return()
  }
}





## SET UP
## STARTING MATRIX with AVERAGE PSIs (scaled)
rc2 <- t(scale(t(B_PSIs_VTS_av_list$dPSI10$CEx)))
acore_list_VTS = acore_sort_VTS = cl_Data_VTS = list()
cl_memb_VTS = cl_events_VTS = cl_genes_VTS = list()



      setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/Clustering/Clustering_CEx_dPSI10/")
      i <- 12
      
      ## CLUSTERING PART A
      ## Set up the Expression set object  (no changes needed!)
      exprs=new("ExpressionSet",exprs=as.matrix(rc2))   ## defines ExpressionSet
      #exprs.r = filter.NA(exprs, thres=0.25)    ## excludes ROWS with more than 25% of the measurements missing - Useless, since we did the same earlier
      #exprs.f = fill.NA(exprs.r,mode="knnw")       ## replaces remaining missing values by the knn algorithm, in which values averaged are weighted by the distance to the corresponding neighbour - Useless, if we already filtered out the NAs
      #exprs.s = filter.std(exprs.f,min.std=0)      ## excludes genes with low standard deviation - leave min.std=0 if no filter needed (i.e. when already considering PSIs significantly changing in at least one comparison)
      exprs.0 = exprs
      #exprs.0 <- standardise(exprs.s)            ## Z-score standardisation - No! If needed, start with the standardised matrix!
      m1 <- mestimate(exprs.0)              ## Estimates fuzzyfier parameter

      
      ## CLUSTERING PART B: CLUSTER NUMBER CHOICE
      ## Calculation of minimum centroid distance for a range of cluster numbers to estimate optimised number of clusters #
      ## (Run only the first time, choose the desired number of clusters and set "i" in PART A) #

            # Dmin <-function (eset, m, crange = seq(4, 40, 4), repeats = 3, visu = TRUE)
            #    {
            #     DminM <- matrix(0, nrow = length(crange), ncol = repeats)
            #      for (ii in 1:repeats) {
            #       j <- 0
            #       for (i in crange) {
            #         cl <- mfuzz(eset, c = i, m = m)
            #          DminM[j <- j + 1, ii] <- min(dist(cl[[1]]))
            #        }
            #      }
            #      DminMav <- apply(DminM, 1, mean)
            #      if (visu)
            #        plot(crange, DminMav, xlab = "Cluster number", ylab = "Min. centroid distance", type="l")
            #      return(DminMav)
            #    }
            # 
            # # Prints the pdf with graph of minimum centroid distances
            # pdf("Minimum_distance_to_cluster_centroid.pdf")
            #   tmp  <- Dmin( exprs.0,        ## ExpressionSet
            #                  m=m1,         ## Fuzzyfier parameter
            #                  crange=seq(3,15,1),   ## Vector of cluster results to test
            #                  repeats=2,       ## Repeats
            #                 visu=TRUE       ## Plot results
            #    )
            #   dev.off()


        ## CLUSTERING PART C: RUN IT!
        # Remember to set "i" to the number of clusters wanted!

        # Soft Clustering code
        set.seed(5678)      # Random setting of the starting centroids in the clusters -> alters the order of appearance of the clusters
        m = m1              # Fuzzifier parameter
        clusters = i        # Number of clusters (set in Part A)
        cl2 <- mfuzz(exprs.0,c=clusters,m=m, iter.max=1000)       # Clustering
        
        # Extracts rows forming the "alpha cores" of soft clusters - This code will extract a unique set of rows (junction coordinates) for each cluster
        acore_list_VTS <- acore(exprs.0, cl=cl2, min.acore=0)
        # str(acore_list_VTS) 

        
        ## CLUSTERING PART E: Summary and Print-Outs 
        # Plots Clusters
        minmemplot=0
        # pdf(paste("resolve_to_",i,"_c.mEst",round(m1,2),".zscorePlot_minmem",minmemplot,".pdf",sep=""),width=13,height=10)
        #   mfuzz.plot3(exprs.0, cl=cl2, mfrow=c(3,3), time.labels = colnames(rc2), x11=FALSE, min.mem = minmemplot)  #change minmem to print "stricter" graphs
        #   dev.off()

        # ColorBar for the graphical result
        # png("ColorBar.png", width=1.5*300, height=10*300, res=300)
        #  mfuzzColorBar()
        #  dev.off()

        # PCA plot showing how clusters are close to each other in a 2 variable space
        # pdf(paste("resolve_to_",i,"c.mEst",round(m1,2),".pcaPlot.pdf",sep=""))
        #   O <- overlap(cl2)
        #   Ptmp <- overlap.plot(cl2,over=O,thres=0.05) ## copied from the protocol
        #   dev.off()

        # For each row, extract all the memberships to the different clusters
        #write.csv(cl2[4], file="membership_emissions.csv")
        
        # Print a unique set of rows (junction coordinates) for each cluster with their membership score    # You can otherwise print the sorted version, acore_sort!
        # for(j in 1:length(acore_list_VTS)) {
           # table <-as.data.frame(cbind(as.character(acore_list_VTS[[j]][,1]),as.numeric(acore_list_VTS[[j]][,2])))
           # write.table(table, paste("resolve_to_",i,"clusters.c",j,".txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
           # }
           # 
        
        # SORT Events in clusters by membership
        acore_sort_VTS <- as.list(sapply (c(1:i), function(x) acore_sort_VTS[[x]]=as.matrix(acore_list_VTS[[x]][(order(factor(acore_list_VTS[[x]][,2]), decreasing = TRUE)),])))   # acore_sort is a list containing SKIPPING coordinates and membership values (sorted by membership)
        # for(j in 1:length(acore_list_VTS)) {
        #     table <-as.data.frame(cbind(as.character(acore_sort_VTS[[j]][,1]),as.numeric(acore_sort_VTS[[j]][,2])))
        #     write.table(table, paste("resolve_to_",i,"clusters.c",j,"_sorted.txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
        # }


        # SELECT EVENTS/GENES OF EACH CLUSTER by min membership
        minmemlist=0 # Set minimum membership threshold
        cl_memb_VTS <- noquote(sapply (c(1:i), function(x) cl_memb_VTS[[x]]=acore_sort_VTS[[x]][which(acore_sort_VTS[[x]][,2]>minmemlist),]))    
        cl_events_VTS <- noquote(sapply (c(1:i), function(x) cl_events_VTS[[x]]=as.vector(acore_sort_VTS[[x]][which(acore_sort_VTS[[x]][,2]>minmemlist),1])))    # cl_events is a list containing the SKIPPING coordinates of the events belonging to each cluster (>minmem) 
        cl_genes_VTS <- lapply (cl_events_VTS, function(x) subset(B_diffEV$dPSI10, EVENT %in% x)[,"GENE"])                       # cl_genes is a list containing the gene names of the events belonging to each cluster (>minmem) 
        
        # Create & PRINT a SORTED summary matrix for every cluster (Gene, Skipping coord, Inclusion coord, PSI/JEFF values)
        cl_Data_VTS=lapply(cl_events_VTS, function(x) subset(B_diffEV$dPSI10, EVENT %in% x))
        #cl_Data_VTS_withdPSI=lapply(cl_events_VTS, function(x) subset(B_EdgeR, EVENT %in% x))
        # for (j in 1:i) {
        #   write.table(cl_Data_VTS[[j]], file=paste("cl_",j,"_Data_","minmem",minmemlist,".txt",sep = ""),row.names = F,quote = F,sep="\t")
        # }
       
## PREPARE SUPP.TABLE1 
        cl_Data_VTS_SUPP1 <- sapply(c(1:i), function(x) {y = subset(B_diffEV$dPSI10, EVENT %in% cl_event[[x]]),
                                                  y$Cluster = x},
                                                  return(y))
  


cl2_VTS <- cl2

#----- Outputs
message("
Outputs: acore_sort_VTS = list with events belonging to each cluster, sorted my membership
         cl_Data_VTS = list of events belonging to each cluster, sorted my membership, including all details
         * All tables and plot are saved in the selected Clustering folder")
#
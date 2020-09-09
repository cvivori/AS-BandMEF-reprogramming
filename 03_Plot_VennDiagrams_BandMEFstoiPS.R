require(tidyverse)
require(VennDiagram)

## IMPORT
summary(B_diffEV_list$dPSI10)
  head(B_diffEV_list$dPSI10$CEx)
summary(M_diffEV_list$dPSI10)
  head(M_diffEV_list$dPSI10$CEx)

## OVERLAP EACH CATEGORY OF EVENTS IN B CELL AND MEFs
  overlaps_BM <- list()
for (t in names(B_diffEV_list$dPSI10)) {
  overlaps_BM[[t]] <- calculate.overlap(list(B_diffEV_list$dPSI10[[t]]$EVENT,M_diffEV_list$dPSI10[[t]]$EVENT))
}
  overlaps_BM$ALL <- calculate.overlap(list(B_diffEV$dPSI10$EVENT,M_diffEV$dPSI10$EVENT))
  lapply(overlaps_BM, function(x) summary(x))

## PERCENTAGE
perc_BM <- lapply(overlaps_BM, function(x) round(100*length(x$a3)/length(x$a2),digits = 2))

## DRAW VENNS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1601 CEBPa_NEW/VASTTOOLS_v2.2_FINAL_Mm10/B2iPS/dPSI10/Overlap_MEFs_Venns/")
for (t in names(overlaps_BM$dPSI10)) {
pdf(file=paste("VennDiagram_",t,"_BandMEFstoiPS_",perc_BM[[t]],"perc.pdf",sep = ""),width = 8, height=8)
  draw.pairwise.venn( area1 = length(overlaps_BM[[t]]$a1),
                      area2 = length(overlaps_BM[[t]]$a2),
                      cross.area = length(overlaps_BM[[t]]$a3),
                      category = c("B cell reprogramming","MEFs reprogramming"),
                      lwd=0,
                      fill=c("darkslategray","maroon"),
                      cat.fontface = "bold",
                      cat.fontfamily = "Helvetica",
                      cat.pos =0,
                      #print.mode=c('raw','percent'),
                      fontfamily = "Helvetica"
  )
  dev.off()
}


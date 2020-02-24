require(tidyverse)
require(ggrepel)
require(reshape2)
require(viridis)
require(VennDiagram)

## FUNCTIONs
## Find pairs of duplicates (in the FullCO column, for example)
FindPairs = function (table, column) {
  Dups=which(duplicated(table[,column]))
  Pairs=list()
  for (i in rownames(table)[Dups]) {
    n=unlist(strsplit(i, split="-"))[1]
    lines=grep(n, rownames(table))
    Pairs[[n]]=table[lines,]
  }
  return(Pairs)
}  
## Check if all the duplicates found are only ALTA / ALTD events
CheckAlt = function (listofpairs) {
  n=names(listofpairs)
  if (is.null(n)) {return("EMPTY")
  } else { test=sapply(n, function(x) grepl("ALT",x))
  if (sum(test)==length(test)) {return("OK")
  } else {return("NO!")}
  }
}
## Count events (Wrap_VTS_Events function)
source("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/useful-R-scripts/VTS_Count_EventTypes.R/")
## Separate events in a list
SplitList_inclMIC = function (input_list) {
  
  out <- list()
  CEx_inclMIC <- c("S","C1","C2","C3","MIC","ANN") 
  RIs <- c("IR-C","IR-S","IR")
  
  for(n in names(input_list)) {
    out[[n]] <- list(
      CEx = subset(input_list[[n]], COMPLEX %in% CEx_inclMIC),
      RI = subset(input_list[[n]], COMPLEX %in% RIs),
      Alt3 = subset(input_list[[n]], COMPLEX == "Alt3"),
      Alt5 = subset(input_list[[n]], COMPLEX == "Alt5")
    )}  
  return(out) }
## Plot scatterplots ddPSI (Plot_ScatterPlots_dPSI function)
source("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/useful-R-scripts/Plot_Scatterplot_dPSI_values.R/")

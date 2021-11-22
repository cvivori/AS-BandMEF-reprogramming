require("optparse")
require(dplyr)
require(stringr)
require(purrr)



#### INPUT ARGUMENTS
option_list = list(
  make_option(c("-p", "--path_vast_out"), type="character",  
              help="path to vast_out folder", metavar="character"),
  make_option(c("-i", "--inclusion_file"), type="character", 
              help="inclusion file name", metavar="character"),
  make_option(c("-s", "--single_replicates"), type="logical", default=FALSE,
              help="TRUE/FALSE: Should you consider Q columns of merge/single replicates samples (not only _1, _2 or _3)? [default %default]"),
  make_option(c("-d", "--path_DiffAS"), type="character", 
              help="path to folder containing all DiffAS files to be considered", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (length(opt) < length(option_list)){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n\nOUTPUT of THIS SCRIPT:
            
       Depending of the input INCLUSION table...     
       -  If you input the FILTERED INCLUSION table (see VTS_INCL_filtering.R) in the --inclusion_file argument,
          'Filtered_INCLUSION[...]_withdPSI_[...]' is a table with ALL CORRECTLY MAPPED events in the dataset (with dPSI from DiffAS files)
       -  If you input the original INCLUSION table in the --inclusion_file argument
          'Filtered_INCLUSION[...]_withdPSI_[...]' is a table with ALL events of the inclusion table (with dPSI from DiffAS files)

** Example run in the cluster: 
qsub_job 1 4 1  -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0.txt -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_with_dPSI/ -s FALSE
\n")
}
# -  'All_AS_events_withdPSI_...' is table with ALL DIFFERENTIALLY SPLICED events in the dataset (with dPSI from DiffAS files)
#     * Pay attention that if you ran vast-tools compare on the original INCLUSION table, there might be N values in some columns,
#     because the filter (min VLOW) is only applied to the columns of given group a and b!


#### IMPORT INCLUSION TABLE ARGUMENTS
message(">> Importing INCLUSION table...")
setwd(opt$path_vast_out)
INCL=read.table(file = opt$inclusion_file, header=T,sep="\t")
rownames(INCL)=INCL$EVENT

message(">> Extracting columns")
cols <- colnames(INCL)
STDcols <- cols[1:6]
Qcols <- cols[grep(".Q",cols)]
if (opt$single_replicates == TRUE) { Qcols_sr <- Qcols } else { Qcols_sr <- Qcols[grep("_[123].Q",Qcols)] }
PSIcols <- setdiff(cols,c(STDcols,Qcols))
if (opt$single_replicates == TRUE) { PSIcols_sr <- PSIcols } else { PSIcols_sr <- PSIcols[grep("_[123]$",PSIcols)] }
message(">> The columns that were considered to filter for N values are:")
print(Qcols_sr)


#### IMPORT DIFFAS FILES
message(">> Opening DiffAS files")
setwd(opt$path_DiffAS)
DiffAS_files <- data.frame(filename=list.files (path= opt$path_DiffAS, pattern = "^DiffAS-*"))
DiffAS_files$name <- gsub("-with_dPSI.tab","",DiffAS_files$filename)
DiffAS_files$comparison <- sapply(DiffAS_files$name, function(x) paste(unlist(strsplit(x,split = "_"))[-1], collapse="_"))
message(paste("... ", DiffAS_files$filename,"\n", sep=""))

message(">> Extracting dPSIs from all DiffAS files")
DiffAS <- list()
for (c in c(1:nrow(DiffAS_files))) {
  DiffAS[[DiffAS_files[c,"comparison"]]] <- read.table(file = as.character(DiffAS_files[c,"filename"]), sep="\t",header=T);
  colnames(DiffAS[[DiffAS_files[c,"comparison"]]])[which(colnames(DiffAS[[DiffAS_files[c,"comparison"]]])=="dPSI")] <- paste("dPSI",DiffAS_files[c,"comparison"],sep = "_")
} 
DiffAS <- lapply(DiffAS, function(x) {rownames(x) <- x$EVENT; return(x)})
#summary(DiffAS)
#head(DiffAS[[1]])


## JOIN TABLES 
# all DiffAS files together
message(">> Joining all DiffAS files")
allDiffAS <- reduce(DiffAS, full_join)
rownames(allDiffAS) <- allDiffAS$EVENT 
lengthDiffAS <- length(DiffAS)

# plus provided inclusion table
message(">> Joining INCLUSION table & all DiffAS files")
INCL_withDiffAS <- INCL %>%
  left_join(allDiffAS)
rownames(INCL_withDiffAS) <- INCL_withDiffAS$EVENT 


prefix <- str_split_fixed(opt$inclusion_file,pattern="\\.",n=2)[,1]
# out_D=allDiffAS
# out_name_D=paste("All_AS_events_withdPSI_",lengthDiffAS,"DiffAS.txt",sep="")
out_ID=INCL_withDiffAS
out_name_ID=paste(prefix,"_withdPSI_",lengthDiffAS,"DiffAS.txt",sep="")

setwd(opt$path_vast_out)
message(">> Printing (in vast_out folder):")
# message(paste(out_name_D,"\nNumber of rows: ",nrow(out_D),sep="\t"))
message(paste(out_name_ID,"\nNumber of rows: ",nrow(out_ID),sep="\t"))
# write.table(out_D, file = out_name_D,sep="\t",row.names = FALSE)
write.table(out_ID, file = out_name_ID,sep="\t",row.names = FALSE)




# if (opt$print_only_unionDiffAS == TRUE) {
#   message(">> Joining all DiffAS files")
#   out=allDiffAS
#   out_name=paste("All_AS_events_withdPSI_",lengthDiffAS,"DiffAS.txt",sep="")
# } else {
#   message(">> Joining INCLUSION table & all DiffAS files")
#   out=INCL_withDiffAS
#   out_name=paste("All_INCL_events_withdPSI_",lengthDiffAS,"DiffAS.txt",sep="")
# }
# message(paste(">> Printing (in vast_out folder)", out_name, sep=" "))
# message(paste0("Number of rows: ",nrow(out))) 
# setwd(opt$path_vast_out)
# write.table(out, file = out_name,sep="\t")


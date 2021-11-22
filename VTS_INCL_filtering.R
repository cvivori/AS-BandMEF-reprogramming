require("optparse")
require(stringr)


#### INPUT ARGUMENTS
option_list = list(
  make_option(c("-p", "--path_vast_out"), type="character", 
              help="Path to vast_out folder", metavar="character"),
  make_option(c("-i", "--inclusion_file"), type="character",  
              help="Inclusion file name", metavar="character"),
  make_option(c("-e", "--exclusion_pattern"), type="character",  
              help="Pattern included in columns to be excluded (e.g. 'M_') **can be 'nothing'!", metavar="character"),
  # make_option(c("-o", "--print_INCL_Qcols"), type="logical", default=FALSE,
  #             help="TRUE/FALSE: Should the program print INCLUSION table with Q columns? [default %default]"),
  make_option(c("-s", "--single_replicates"), type="logical", default=FALSE,
              help="TRUE/FALSE: Should you consider Q columns of merge/single replicates samples (not only _1, _2 or _3)? [default %default]"),
  make_option(c("-m", "--max_N_sr"), type="numeric", default="0", 
              help="Max number of N values allowed in read coverage of single replicates [default %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (length(opt) < length(option_list)){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}


#### IMPORT INCLUSION TABLE ARGUMENTS
message(">> Importing INCLUSION table...")
setwd(opt$path_vast_out)
INCL=read.table(file = opt$inclusion_file, header=T,sep="\t")
rownames(INCL)=INCL$EVENT
head(INCL)

#### EXTRACT COLUMNS TO FILTER FOR N VALUES
message(">> Extracting columns")
cols <- colnames(INCL)
STDcols <- cols[1:6]
Qcols <- cols[grep(".Q",cols)]
if  (opt$exclusion_pattern =="nothing" ) {Qcols_fil <- Qcols} else { Qcols_fil <- Qcols[-grep(opt$exclusion_pattern,Qcols)] }
if (opt$single_replicates == TRUE) { Qcols_sr <- Qcols_fil } else { Qcols_sr <- Qcols_fil[grep("_[123].Q",Qcols_fil)] }
PSIcols <- setdiff(cols,c(STDcols,Qcols))
if  (opt$exclusion_pattern =="nothing" ) {PSIcols_fil <- PSIcols} else { PSIcols_fil <- PSIcols[-grep(opt$exclusion_pattern,PSIcols)] }
if (opt$single_replicates == TRUE) { PSIcols_sr <- PSIcols_fil } else { PSIcols_sr <- PSIcols_fil[grep("_[123]$",PSIcols_fil)] }
message(">> The columns that were considered to filter for N values are:")
print(Qcols_sr)

#### COUNTING N VALUES IN READ COVERAGE (based on actual reads)
message(">> Counting number of Ns")
INCLQ <- INCL[,Qcols_sr]
Q_ReadCov_sr <- apply(INCLQ,2, function(x) str_split_fixed(x,pattern=",",n=6)[,1])  # Score 1: Read coverage, based on actual reads 
    rownames(Q_ReadCov_sr) <- rownames(INCL)
    # head(Q_ReadCov_sr)
freq_N_sr <- apply(Q_ReadCov_sr, 1, function(x) sum(grepl("^N",x)))
# head(freq_N_sr)
    # freq_VLOW <- apply(Q_ReadCov, 1, function(x) sum(grepl("^VLOW",x)))
    # freq_LOW <- apply(Q_ReadCov, 1, function(x) sum(grepl("^LOW",x)))
    # freq_OK <- apply(Q_ReadCov, 1, function(x) sum(grepl("^OK",x)))
    # freq_SOK <- apply(Q_ReadCov, 1, function(x) sum(grepl("^SOK",x)))
# Q_CorrReadCov <- apply(INCLQ,2, function(x) str_split_fixed(x,pattern=",",n=6)[,2])   # Score 2: Read coverage, based on corrected reads 
# Q_Reads <- apply(INCLQ,2, function(x) str_split_fixed(x,pattern=",",n=6)[,3])     # Score 3: see https://github.com/vastgroup/vast-tools#combine-output-format
# Q_Imbalance <- apply(INCLQ,2, function(x) str_split_fixed(x,pattern=",",n=6)[,4])     # Score 4: see https://github.com/vastgroup/vast-tools#combine-output-format
# Q_Complexity <- apply(INCLQ,2, function(x) str_split_fixed(x,pattern=",",n=6)[,5])    # Score 5: see https://github.com/vastgroup/vast-tools#combine-output-format

#### FILTER: NO MORE N VALUES THAN WHAT IS SET IN THE FILTER
message(">> Filtering for max number of Ns (in single replicate samples)")
to_filter<- INCL
to_filter$freq_N_sr <- freq_N_sr
filtered <- to_filter[which(to_filter$freq_N_sr <= opt$max_N_sr),]
#head(filtered)
# summary(INCL$freq_N_sr)
# summary(filtered$freq_N_sr)

#### REMOVE 'Q' COLUMNS CONTAINING QUALITY SCORES
message(">> Removing Q columns")
prefix <- str_replace_all(opt$inclusion_file,".tab","")
out_noQ=filtered[,c(STDcols,PSIcols,"freq_N_sr")]
out_name_noQ=paste("Filtered_",prefix,"_noQ_maxN_",opt$max_N_sr,".txt",sep="")
out=filtered
out_name=paste("Filtered_",prefix,"_maxN_",opt$max_N_sr,".txt",sep="")

#### OUTPUT FILTERED INCLUSION TABLES BOTH WITH AND WITHOUT Q COLUMNS
message(paste(">> Printing tables:", out_name, out_name_noQ, sep="\n"))
write.table(out, file = out_name,sep="\t", quote = FALSE,row.names = FALSE)
write.table(out_noQ, file = out_name_noQ,sep="\t", quote = FALSE,row.names = FALSE)



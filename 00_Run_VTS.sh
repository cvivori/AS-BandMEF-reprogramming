#!/bin/bash
# See https://github.com/vastgroup/vast-tools

#  VAST_TOOLS ALIGN, MM10

# B2iPS
# qsub_job 48 20 15 -N vtA_B1 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/B_1_12474_GCCAAT_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/B_1_12474_GCCAAT_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_B2 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/B_2_12475_ACAGTG_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/B_2_12475_ACAGTG_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_Ba1 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/Ba_1_12476_CTTGTA_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/Ba_1_12476_CTTGTA_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_Ba2 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/Ba_2_12477_GTGAAA_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/Ba_2_12477_GTGAAA_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D21 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D2_1_12478_TGACCA_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D2_1_12478_TGACCA_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D22 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D2_2_12479_CAGATC_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D2_2_12479_CAGATC_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D41 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D4_1_12480_ACTTGA_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D4_1_12480_ACTTGA_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D42 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D4_2_12481_GGCTAC_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D4_2_12481_GGCTAC_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D61 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D6_1_12482_ACTGAT_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D6_1_12482_ACTGAT_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D62 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D6_2_12483_AGTCAA_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D6_2_12483_AGTCAA_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D81 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D8_1_12484_GTCCGC_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D8_1_12484_GTCCGC_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_D82 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D8_2_12485_GAGTGG_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/D8_2_12485_GAGTGG_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_ES1 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/ES_1_12488_GATCAG_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/ES_1_12488_GATCAG_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_ES2 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/ES_2_12489_TAGCTT_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/ES_2_12489_TAGCTT_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_iPS1 vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/iPS_1_12486_CGATGT_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/iPS_1_12486_CGATGT_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
# qsub_job 48 20 15 -N vtA_iPS2 -m ea -M claudia.vivori@crg.eu vast-tools align /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/iPS_2_12487_ATCACG_read1.fastq.gz /users/jvalcarcel/cvivori/CEBPa_NEW/FASTQ_Links/iPS_2_12487_ATCACG_read2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
 
# MEFs2iPS
#  qsub_job 168 30 15 -N vtA_D01 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep1_SRR2070947_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep1_SRR2070947_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D02 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep2_SRR2070948_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep2_SRR2070948_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D03 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep3_SRR2070949_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day0_rep3_SRR2070949_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D41 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep1_SRR2070950_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep1_SRR2070950_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D42 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep2_SRR2070951_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep2_SRR2070951_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D43 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep3_SRR2070952_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day4_rep3_SRR2070952_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D71 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep1_SRR2070953_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep1_SRR2070953_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D72 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep2_SRR2070954_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep2_SRR2070954_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D73 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep3_SRR2070955_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day7_rep3_SRR2070955_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D101 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep1_SRR2070956_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep1_SRR2070956_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D102 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep2_SRR2070957_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep2_SRR2070957_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D103 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep3_SRR2070958_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day10_rep3_SRR2070958_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D151 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep1_SRR2070959_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep1_SRR2070959_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D152 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep2_SRR2070960_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep2_SRR2070960_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D153 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep3_SRR2070961_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day15_rep3_SRR2070961_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D201 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep1_SRR2070962_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep1_SRR2070962_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D202 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep2_SRR2070963_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep2_SRR2070963_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_D203 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep3_SRR2070964_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Day20_rep3_SRR2070964_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_CL1 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep1_SRR2070965_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep1_SRR2070965_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_CL2 vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep2_SRR2070966_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep2_SRR2070966_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 
#  qsub_job 168 30 15 -N vtA_CL3 -m ea -M claudia.vivori@crg.eu vast-tools align /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep3_SRR2070967_1.fastq.gz  /no_backup/jvalcarcel/cvivori/Cieply_MEFs2iPS_Apr2016/FASTQ_files/iPS_PairedEnd/Clone_rep3_SRR2070967_2.fastq.gz -sp Mm2 -c 15 -o vast_out --IR_version 2 

# COPY MEFs FILES
#cp ~/Cieply_MEFstoiPS_2016/VASTTOOLS_FINAL/vast_out/to_combine/* vast_out/to_combine/

## RENAMING OUTPUTS TO COMBINE
#source 00_Rename_VTSintermediates.sh

## no VASTTOOLS MERGE

## VASTTOOLS COMBINE
#qsub_job 12 8 2 -N vtC_BandMEFs vast-tools combine -sp Mm2




#### VASTTOOLS COMPARE ON ALL POSSIBLE COUPLES
# # B cells ∆PSI>15 ∆range>5
# perl ~/Scripts/VTS_run_compare_all_couples.pl /users/jvalcarcel/cvivori/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/INCLUSION_LEVELS_FULL-Mm237.tab 15 5 Mm2 DUPS B_Bcells B_Bpulse B_Day2 B_Day4 B_Day6 B_Day8 B_iPS B_ES
# # B cells ∆PSI>10 ∆range>5
# perl ~/Scripts/VTS_run_compare_all_couples.pl /users/jvalcarcel/cvivori/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/INCLUSION_LEVELS_FULL-Mm237.tab 10 5 Mm2 DUPS B_Bcells B_Bpulse B_Day2 B_Day4 B_Day6 B_Day8 B_iPS B_ES
# # MEFs ∆PSI>15 ∆range>5
# perl ~/Scripts/VTS_run_compare_all_couples.pl /users/jvalcarcel/cvivori/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/INCLUSION_LEVELS_FULL-Mm237.tab 15 5 Mm2 TRIPS M_Day0 M_Day4 M_Day7 M_Day10 M_Day15 M_Day20 M_Clone
# # MEFs ∆PSI>10 ∆range>5
# perl ~/Scripts/VTS_run_compare_all_couples.pl /users/jvalcarcel/cvivori/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/INCLUSION_LEVELS_FULL-Mm237.tab 10 5 Mm2 TRIPS M_Day0 M_Day4 M_Day7 M_Day10 M_Day15 M_Day20 M_Clone


## TIDY UP
#source ~/Scripts/VTS_tidyup_vastout.sh



## FILTER INCLUSION TABLE (no N values)
# qsub_job 4 24 8 -N R_VIf Rscript --vanilla ~/Scripts/VTS_INCL_filtering.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -e M_ -s FALSE -m 0
# qsub_job 4 24 8 -N R_VIf Rscript --vanilla ~/Scripts/VTS_INCL_filtering.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -e B_ -s FALSE -m 0

# ## CREATE UNION OF DIFFAS FILES (dPSI15) SEPARATE B / M
# qsub_job 1 4 1  -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI15_with_dPSI_B/ -s FALSE
# qsub_job 1 4 1  -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI15_with_dPSI_M/ -s FALSE
# qsub_job 1 4 1  -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noM.txt -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI15_with_dPSI_B/ -s FALSE
# qsub_job 1 4 1  -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noB.txt -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI15_with_dPSI_M/ -s FALSE

# ## CREATE UNION OF DIFFAS FILES (dPSI10) SEPARATE B / M
# qsub_job 1 4 1 -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI10_with_dPSI_B/ -s FALSE
# qsub_job 1 4 1 -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i INCLUSION_LEVELS_FULL-Mm237.tab -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI10_with_dPSI_M/ -s FALSE
# qsub_job 1 4 1 -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noM.txt -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI10_with_dPSI_B/ -s FALSE
# qsub_job 1 4 1 -N R_Vdpsi Rscript --vanilla ~/Scripts/VTS_add_dPSI_toINCL.R -p ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/ -i Filtered_INCLUSION_LEVELS_FULL-Mm237_maxN_0_noB.txt -d ~/CEBPa_NEW/VASTTOOLS_FINAL/vast_out/DiffAS_files_dPSI10_with_dPSI_M/ -s FALSE



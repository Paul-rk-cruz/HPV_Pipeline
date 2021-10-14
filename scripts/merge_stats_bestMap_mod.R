#!/usr/bin/env Rscript
# HPV Pipeline v1.0 NGS Analysis R Script - Modified for Nextflow compatibility

# Import libraries
library(tidyverse)
#NOTE: After a first pass, realized that ambiguous calls on some genomes are causing misses. Removing duplicate genomes (the ones less used) for 16,18,35,45,51,39,59. May end up doing same for non-high risk

# SETUP COMMAND LINE ARGUMENTS
#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
	print("No arguments supplied.")
}else{
    output_directory=args[[1]]
	path=args[[2]]
    run_name=args[[3]]
    }

threshold_filter <- 0.05 #percentage above which calls are considered real
ext <- "_hpvAll_scafstats.txt" #scafstats file name extension
glob <- paste0("*",ext)

setwd(path) 

#read in mapping stats and merge:
files <- fs::dir_ls(path, glob=glob)   
scaf_stats <- files %>% purrr::map_dfr(read_tsv, .id = "Sample")

#tidy sample names:
scaf_stats$Sample <- gsub(path, "",scaf_stats$Sample)
scaf_stats$Sample <- gsub(ext, "",scaf_stats$Sample)

#Rename Columns:
colnames(scaf_stats) <- c("Sample", "Reference", "Percent Unambiguous Reads", "x", "Percent Ambiguous Reads", "x", "Unambiguous Reads", "Ambiguous Reads", "Assigned Reads", "x")

ref_type <- str_split_fixed(scaf_stats$Reference, "_", 2)
colnames(ref_type) <- c("Reference Accession", "HPV Type")

scaf_stats <- cbind(scaf_stats, ref_type)


#Filter Relevant Columns 
scaf_stats <- scaf_stats[,c(1,11,12,3,5,7,8,9)]

filtered_scaf_stats <- filter(scaf_stats, `Percent Unambiguous Reads` > threshold_filter) 

top_hit <- scaf_stats %>% group_by(Sample) %>% filter(row_number()==1)

write.csv(filtered_scaf_stats, paste0("~${params.outdir}filtered_scafstats_", run_name, ".csv"))
write.csv(scaf_stats, paste0("~${params.outdir}all_scafstats_", run_name, ".csv"))
write.csv(top_hit, paste0("~${params.outdir}topHit_scafstats_", run_name, ".csv"))
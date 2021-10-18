    #!/usr/bin/env Rscript

    # IMPORT LIBRARIES
    library(tidyverse)

    #Get args from command line 
    args<-(commandArgs(TRUE));
    if(length(args)==0){
        print("No arguments supplied.")
    }else{
    run_name <- args[1]
    path <- args[2]
    outdir <- args[3]
    }

    # SETUP VARIABLES
    threshold_filter <- 0.05
    ext <- "_hpvAll_scafstats.txt" 
    glob <- paste0("*",ext)

    # SET PATH
    setwd(path)

    # READ BBMAP OUTPUT - SCAFFOLD STATS
    files <- fs::dir_ls(path, glob=glob)   
    scaf_stats <- files %>% purrr::map_dfr(read_tsv, .id = "Sample")

    # SAMPLE NAMES:
    scaf_stats$Sample <- gsub(path, "",scaf_stats$Sample)
    scaf_stats$Sample <- gsub(ext, "",scaf_stats$Sample)

    # RENAME COLUMNS:
    colnames(scaf_stats) <- c("Sample", "Reference", "Percent Unambiguous Reads", "x", "Percent Ambiguous Reads", "x", "Unambiguous Reads", "Ambiguous Reads", "Assigned Reads", "x")

    ref_type <- str_split_fixed(scaf_stats$Reference, "_", 2)
    colnames(ref_type) <- c("Reference Accession", "HPV Type")

    scaf_stats <- cbind(scaf_stats, ref_type)

    # FILTER RELEVANT COLUMNS 
    scaf_stats <- scaf_stats[,c(1,11,12,3,5,7,8,9)]

    filtered_scaf_stats <- filter(scaf_stats, `Percent Unambiguous Reads` > threshold_filter) 

    top_hit <- scaf_stats %>% group_by(Sample) %>% filter(row_number()==1)

    # WRITE CSV
    write.csv(filtered_scaf_stats, paste0(outdir, "analysis/filtered_scafstats_", run_name, ".csv"))
    write.csv(scaf_stats, paste0(outdir, "analysis/all_scafstats_", run_name, ".csv"))
    write.csv(top_hit, paste0(outdir, "analysis/topHit_scafstats_", run_name, ".csv"))
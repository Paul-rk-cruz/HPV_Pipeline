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

    # DIRECTORY CONTAINING SCAF STATS
    scaf_path <- paste0(path,"bbmap_scaf_stats/") #folder containing scaf stats
    scaf_ext <- "_R1_hpvAll_scafstats.txt" #scafstats file name extension
    glob_scaf <- paste0("*",scaf_ext)

    # DIRECTORY CONTAINING COV STATS
    cov_path <- paste0(path,"bbmap_cov_stats/") #folder containing coverage stats
    cov_ext <- "_R1_hpvAll_covstats.txt" #covstats file name extension
    glob_cov <- paste0("*",cov_ext)

    # READ BBMAP OUTPUT/SCAFFOLD STATS AND MERGE
    scaf_files <- fs::dir_ls(scaf_path, glob=glob_scaf)   
    scaf_stats <- scaf_files %>% purrr::map_dfr(read_tsv, .id = "Sample")

    # TIDY SAMPLE NAMES:
    scaf_stats$Sample <- gsub(scaf_path, "",scaf_stats$Sample)
    scaf_stats$Sample <- gsub(scaf_ext, "",scaf_stats$Sample)

    # RENAME COLUMNS:
    colnames(scaf_stats) <- c("Sample", "Reference", "Percent Unambiguous Reads", "x", "Percent Ambiguous Reads", "y", "Unambiguous Reads", "Ambiguous Reads", "Assigned Reads", "z")

    # READ BBMAP OUTPUT/COVERAGE STATS AND MERGE
    cov_files <- fs::dir_ls(cov_path, glob=glob_cov)   
    cov_stats <- cov_files %>% purrr::map_dfr(read_tsv, .id = "Sample")

    # TIDY SAMPLE NAMES:
    cov_stats$Sample <- gsub(cov_path, "",cov_stats$Sample)
    cov_stats$Sample <- gsub(cov_ext, "",cov_stats$Sample)
    colnames(cov_stats)[2] <- "Reference"
    colnames(cov_stats)[3] <- "Average Fold Coverage"
    colnames(cov_stats)[6] <- "Covered Percent of Reference"
    colnames(cov_stats)[11] <- "Median Fold Coverage"

    # MERGE SCAF & COV STATS
    # first paste sample name and ref to make unique key on each:

    scaf_stats$key <- paste0(scaf_stats$Sample, "_", scaf_stats$Reference)
    cov_stats$key <- paste0(cov_stats$Sample, "_", cov_stats$Reference)

    # inner join from scaf to cov by key:
    scaf_stats <- inner_join(scaf_stats, cov_stats, by="key")
    colnames(scaf_stats)[1:2] <- c("Sample", "Reference")

    # CLEAN DATA
    ref_type <- str_split_fixed(scaf_stats$Reference, "_", 2)
    colnames(ref_type) <- c("Reference Accession", "HPV Type")

    scaf_stats <- cbind(scaf_stats, ref_type)

    # FILTER RELEVANT COLUMNS
    scaf_stats <- scaf_stats[,c(1,24,25,3,5,7,8,9,14,22,17)]

    filtered_scaf_stats <- filter(scaf_stats, `Percent Unambiguous Reads` > threshold_filter) 

    top_hit <- scaf_stats %>% group_by(Sample) %>% filter(row_number()==1)

    # WRITE CSV FILES TO ANALYSIS DIRECTORY
    write.csv(filtered_scaf_stats, paste0(outdir, "analysis/filtered_scafstats_", run_name, ".csv"))
    write.csv(scaf_stats, paste0(outdir, "analysis/all_scafstats_", run_name, ".csv"))
    write.csv(top_hit, paste0(outdir, "analysis/topHit_scafstats_", run_name, ".csv"))
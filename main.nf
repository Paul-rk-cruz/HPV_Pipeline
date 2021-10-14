#!/usr/bin/env nextflow

/*
========================================================================================
                 Human Papilloma Virus Pipeline v1.0
========================================================================================
 Github Repo:
 Greninger Lab
 
 Author:
 Paul RK Cruz <kurtisc@uw.edu>
 Nicole A Lieberman <naliebe@uw.edu>
 Alex L Greninger <agrening@uw.edu>
 UW Medicine | Virology
 Department of Laboratory Medicine and Pathology
 University of Washington
 Created: October 12, 2021
 Updated: October 12, 2021
 LICENSE: GNU
----------------------------------------------------------------------------------------
*/
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                DISPLAY HELP MSG                    */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Pipeline version
version = '1.0'
def helpMsg() {
    log.info"""
	 _______________________________________________________________________________
     Human Rhinovirus Genome Mapping Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/HPV_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR --SingleEnd
    Valid CLI Arguments:
    REQUIRED:
      --reads                       Path to input fastq.gz folder).
      --outdir                      The output directory where the results will be saved
    OPTIONAL:
      --skipTrimming                Skips the fastq trimmming process   
    """.stripIndent()
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*              CONFIGURATION VARIABLES               */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Make sure outdir path ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Make sure reads path ends with trailing slash
if (!params.reads.endsWith("/")){
   params.reads = "${params.outdir}/"
}
// Trimmomatic Paths and variables
params.ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
params.ADAPTERS_EE = file("${baseDir}/adapters/TruSeq2-PE.fa")
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
ADAPTERS_PE = file("${baseDir}/adapters/TruSeq2-PE.fa")
params.SETTING = "2:30:10:1:true"
SETTING = "2:30:10:1:true"
params.LEADING = "3"
LEADING = "3"
params.TRAILING = "3"
TRAILING = "3"
params.SWINDOW = "4:20"
SWINDOW = "4:20"
params.MINLEN = "35"
MINLEN = "35"
// Setup Parameters to default values
params.skipTrimming = false
params.singleEnd = false
params.helpMsg = false
params.runName = false
params.ADAPTERS_PE = false
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*            HPV MULTI-FASTA REFERENCES              */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
REF_HPV_ALL = file("${baseDir}/ref_fasta/hpvAll.fasta")
REF_HPV_HIGHRISK = file("${baseDir}/ref_fasta/hpvHighRisk.fasta")
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                  R SCRIPT PATH                     */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
MERGE_STATS_R = file("${baseDir}/scripts/merge_stats_bestMap.R")
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                  SET UP CHANNELS                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}
if (! params.reads ) exit 1, "> Error: Fastq files not found. Please specify a valid path with --reads"
// Create channel for input reads.
// Import reads depending on single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}*_R{1,2}*.fastq.gz")
        // .ifEmpty { error "> Cannot locate paired-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // input: *.gz
    input_read_ch = Channel
        .fromPath("${params.reads}*.gz")
        .ifEmpty { error "> Cannot locate single-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> file(it)}
}
// Setup Run Name
if(params.runName != false) {
    RUN_NAME = params.runName
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*              WORKFLOW DISPLAY HEADER               */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
def hrvheader() {
    
    return """
    """.stripIndent()
}
// log files header
// log.info hrvheader()
log.info "_______________________________________________________________________________"
log.info " Human Papilloma Virus Pipeline :  v${version}"
log.info "_______________________________________________________________________________"
def summary = [:]
summary['Configuration Profile:'] = workflow.profile
summary['Current directory path:']        = "$PWD"
summary['HRV Pipeline directory path:']          = workflow.projectDir
summary['Input directory path:']               = params.reads
summary['Output directory path:']          = params.outdir
summary['Work directory path:']         = workflow.workDir
summary['Sequence type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
if (params.singleEnd) {
summary['Trimmomatic adapters:'] = params.ADAPTERS_SE
} else {
summary['Trimmomatic adapters:'] = params.ADAPTERS_PE
}
summary["Trimmomatic read length (minimum):"] = params.MINLEN
summary["Trimmomatic Setting:"] = params.SETTING
summary["Trimmomatic Sliding Window:"] = params.SWINDOW
summary["Trimmomatic Leading:"] = params.LEADING
summary["Trimmomatic Trailing:"] = params.TRAILING
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "_______________________________________________________________________________"
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                WORKFLOW PROCESSES                  */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*
 * STEP 1: Trim_Reads
 * Trimming of low quality and short NGS sequences.
 */
if (params.singleEnd) {
process Trimming {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3
    // echo true

    input:
        file R1 from input_read_ch
        file ADAPTERS_SE
        val MINLEN

    output:
        tuple env(base),file("*.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv") into Trim_out_SE, Trim_out_SE_FQC

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    
    base=`basename ${R1} ".fastq.gz"`
    echo \$base
    
    /usr/local/miniconda/bin/trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_SE}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    printf "\$num_trimmed" >> ${R1}_num_trimmed.txt
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed> \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr

    """
}
process Mapping {
    container "docker.io/paulrkcruz/hrv-pipeline:latest" 
    errorStrategy 'retry'
    maxRetries 3
    // echo true

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_SE
        file REF_HPV_ALL_FASTA from REF_HPV_ALL
        file REF_HPV_HIGHRISK_FASTA from REF_HPV_HIGHRISK
        RUN_NAME

    output:
        tuple val(base), file("${base}_summary2.csv"), file("${base}_hpvAll.sorted.bam"), file("${base}_hpvAll_covstats.txt"), file("${base}_hpvAll_scafstats.txt") into Mapping_files_ch
        tuple val (base), file("*") into Dump_files_ch

    publishDir "${params.outdir}bam_sorted", mode: 'copy', pattern:'*_hpvAll.sorted.bam*'
    publishDir "${params.outdir}bbmap_stats", mode: 'copy', pattern:'*.txt*'

    script:

    """
    #!/bin/bash

    /usr/local/miniconda/bin/bbmap.sh in=${base}.trimmed.fastq.gz ref=${REF_HPV_ALL_FASTA} outm=${base}_hpvAll.sam outu=nope.sam threads=${task.cpus} maxindel=9 ambiguous=best covstats=${base}_hpvAll_covstats.txt scafstats=${base}_hpvAll_scafstats.txt
    
    /usr/local/miniconda/bin/samtools view -S -b ${base}_hpvAll.sam > ${base}_hpvAll.bam
    /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}_hpvAll.bam > ${base}_hpvAll.sorted.bam

    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'filtered_scafstats_${RUN_NAME}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'all_scafstats_${RUN_NAME}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'topHit_scafstats_${RUN_NAME}.csv'

    mkdir ${params.outdir}analysis
    cp filtered_scafstats_${RUN_NAME}.csv ${params.outdir}/analysis/
    cp all_scafstats_${RUN_NAME}.csv ${params.outdir}/analysis/
    cp topHit_scafstats_${RUN_NAME}.csv ${params.outdir}/analysis/

    cp ${base}_summary.csv ${base}_summary2.csv

    """
}
process Analysis {
    container "docker.io/rocker/tidyverse:latest"
    errorStrategy 'ignore'
    // maxRetries 3
    // echo true

    input: 
    tuple val(base), file("${base}_summary2.csv"), file("${base}_hpvAll.sorted.bam"), file("${base}_hpvAll_covstats.txt"), file("${base}_hpvAll_scafstats.txt") from Mapping_files_ch
        file MERGE_STATS_R
        RUN_NAME

    output:
        tuple val(base), file("${base}_trim_stats.csv"), file("${base}_hpvAll.sorted.bam"), file("${base}_hpvAll_covstats.txt"), file("${base}_hpvAll_scafstats.txt") into Analysis_ch

    publishDir "${params.outdir}analysis", mode: 'copy', pattern:'*_trim_stats.csv*'
    
    script:

    """
    #!/usr/bin/env Rscript

    # IMPORT LIBRARIES
    library(tidyverse)

    # SETUP VARIABLES
    output_directory <- "${params.outdir}"
    run_name <- "${RUN_NAME}"
    threshold_filter <- 0.05
    path <- "${params.outdir}/bbmap_stats/"
    ext <- "_hpvAll_scafstats.txt" 
    glob <- paste0("*",ext)

    # SET PATH
    setwd(path)

    # READ BBMAP OUTPUT - SCAFFOLD STATS
    files <- fs::dir_ls(path, glob=glob)   
    scaf_stats <- files %>% purrr::map_dfr(read_tsv, .id = "Sample")

    # SAMPLE NAMES:
    scaf_stats\$Sample <- gsub(path, "",scaf_stats\$Sample)
    scaf_stats\$Sample <- gsub(ext, "",scaf_stats\$Sample)

    # RENAME COLUMNS:
    colnames(scaf_stats) <- c("Sample", "Reference", "Percent Unambiguous Reads", "x", "Percent Ambiguous Reads", "x", "Unambiguous Reads", "Ambiguous Reads", "Assigned Reads", "x")
    ref_type <- str_split_fixed(scaf_stats\$Reference, "_", 2)
    colnames(ref_type) <- c("Reference Accession", "HPV Type")
    scaf_stats <- cbind(scaf_stats, ref_type)

    # FILTER RELEVANT COLUMNS 
    scaf_stats <- scaf_stats[,c(1,11,12,3,5,7,8,9)]
    filtered_scaf_stats <- filter(scaf_stats, 'Percent Unambiguous Reads' > threshold_filter) 
    top_hit <- scaf_stats %>% group_by(Sample) %>% filter(row_number()==1)

    # WRITE CSV
    write.csv(filtered_scaf_stats, paste0("${params.outdir}analysis/filtered_scafstats_", run_name, ".csv"))
    write.csv(scaf_stats, paste0("${params.outdir}analysis/all_scafstats_", run_name, ".csv"))
    write.csv(top_hit, paste0("${params.outdir}analysis/topHit_scafstats_", run_name, ".csv"))

    # OUTPUT TRIM_STATS SUMMARY
    cp ${base}_summary2.csv ${base}_trim_stats.csv
    
    """
}
}
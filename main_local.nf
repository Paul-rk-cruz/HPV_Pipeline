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
// Set pipeline version
version = '1.0'
params.helpMsg = false
def helpMsg() {
    log.info"""
	 _______________________________________________________________________________
     Human Papilloma Virus Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/HPV_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR --SingleEnd
    Valid CLI Arguments:
    REQUIRED:
      --input                       Path to input fastq.gz folder
      --outdir                      The output directory where the results will be saved
      --runName                     Specifies the run name
      --singleEnd                   Specifies single-end illumina sequences
      --ref                         Specifies the HPV reference multifasta type (all OR highRisk)
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
// Setup Parameters to default values
params.input = false
params.outdir = false
params.ref = false
params.skipTrimming = false
params.singleEnd = false
params.runName = false
params.ADAPTERS_PE = false
// Make sure outdir path ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Make sure reads path ends with trailing slash
if (!params.input.endsWith("/")){
   params.input = "${params.outdir}/"
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
if (! params.input ) exit 1, "> Error: Fastq files not found. Please specify a valid path with --input"
// Create channel for input reads.
// Import reads depending on single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.input}*_R{1,2}*.fastq.gz")
        // .ifEmpty { error "> Cannot locate paired-end reads in: ${params.input}.\n> Please enter a valid file path." }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // input: *.gz
    input_read_ch = Channel
        .fromPath("${params.input}*.gz")
        .ifEmpty { error "> Cannot locate single-end reads in: ${params.input}.\n> Please enter a valid file path." }
        .map { it -> file(it)}
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*        SETUP CONFIGURATION VARIABLES               */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
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
params.MINLEN = "75"
MINLEN = "75"
// Reference multifastas
if(params.ref != false) {
    REF = params.ref
}
REF_HPV_ALL = file("${baseDir}/ref_fasta/hpvAll.fasta")
REF_HPV_HIGHRISK = file("${baseDir}/ref_fasta/hpvHighRisk.fasta")
// Setup Run Name
if(params.runName != false) {
    RUN_NAME = params.runName
}
// R script path
MERGE_STATS_R = file("${baseDir}/scripts/analysis.R")
// Workflow display header
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
summary['Input directory path:']               = params.input
summary['Output directory path:']          = params.outdir
summary['Work directory path:']         = workflow.workDir
summary['Sequence type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Reference type:']           	  = params.ref
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

if (params.singleEnd) {
/*
 * STEP 1: Trim_Reads
 * Trimming of low quality and short NGS sequences.
 */
if (params.singleEnd) {
process Trimming {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        file R1 from input_read_ch
        file ADAPTERS_SE
        val MINLEN

    output:
        tuple env(base),file("*.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv") into Trimming_ch

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
/*
 * STEP 2: Alignment
 * Align NGS Sequence reads to HPV-ALL multifasta
 */
process Aligning {
    // container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trimming_ch
        file REF_HPV_ALL_FASTA from REF_HPV_ALL
        file REF_HPV_HIGHRISK_FASTA from REF_HPV_HIGHRISK
        RUN_NAME
        REF

    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_stats1.csv"), file("${base}_hpvAll.sam") into Aligning_ch
        tuple val(base), file("${base}_hpvAll_scafstats.txt") into Bbmap_scaf_stats_ch
        tuple val(base), file("${base}_hpvAll_covstats.txt") into Bbmap_cov_stats_ch

    publishDir "${params.outdir}bbmap_scaf_stats", mode: 'copy', pattern:'*_hpvAll_scafstats.txt*'
    publishDir "${params.outdir}bbmap_cov_stats", mode: 'copy', pattern:'*_hpvAll_covstats.txt*'    

    script:
    if( REF == 'all' )
    """
    #!/bin/bash

    bbmap.sh -Xmx2g in=${base}.trimmed.fastq.gz ref=${REF_HPV_ALL_FASTA} outm=${base}_hpvAll.sam outu=${base}_nope.sam scafstats=${base}_hpvAll_scafstats.txt covstats=${base}_hpvAll_covstats.txt maxindel=9 ambiguous=best threads=${task.cpus} 

    cp ${base}_summary.csv ${base}_stats1.csv

    """
    else if( REF == 'highRisk' )
    """
    #!/bin/bash

    bbmap.sh -Xmx2g in=${base}.trimmed.fastq.gz ref=${REF_HPV_HIGHRISK_FASTA} outm=${base}_hpvAll.sam outu=${base}_nope.sam scafstats=${base}_hpvAll_scafstats.txt covstats=${base}_hpvAll_covstats.txt maxindel=9 ambiguous=best threads=${task.cpus}

    cp ${base}_summary.csv ${base}_stats1.csv

    """
   else
    """
    #!/bin/bash

    /usr/local/bin/bbmap.sh -Xmx2g in=${base}.trimmed.fastq.gz ref=${REF_HPV_ALL_FASTA} outm=${base}_hpvAll.sam outu=${base}_nope.sam scafstats=${base}_hpvAll_scafstats.txt covstats=${base}_hpvAll_covstats.txt maxindel=9 ambiguous=best threads=${task.cpus}

    cp ${base}_summary.csv ${base}_stats1.csv
    
    """
}
/*
 * STEP 2: Bam_Sorting
 * Sort bam file and collect summary statistics.
 */
process Bam_Sorting { 
    container "quay.io/greninger-lab/swift-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
      tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_stats1.csv"), file("${base}_hpvAll.sam") from Aligning_ch
    output:
      tuple val(base), file("${base}_hpvAll.sam"), file("${base}_hpvAll.sorted.bam"), file("${base}_trim_stats.csv")  into Analysis_ch   

    publishDir "${params.outdir}trim_stats", mode: 'copy', pattern:'*_trim_stats.csv*'
    publishDir "${params.outdir}bam_sorted", mode: 'copy', pattern:'*_hpvAll.sorted.bam*'

    script:
    """
    #!/bin/bash

    /usr/local/miniconda/bin/samtools view -S -b ${base}_hpvAll.sam > ${base}_hpvAll.bam
    /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}_hpvAll.bam > ${base}_hpvAll.sorted.bam

    cp ${base}_stats1.csv ${base}_trim_stats.csv

    """
}
/*
 * STEP 4: Analysis
 * Analysis summary creation utilizing R script.
 */
process Analysis {
    container "docker.io/rocker/tidyverse:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
    file("${base}_hpvAll_scafstats.txt") from Bbmap_scaf_stats_ch.collect()     
    file MERGE_STATS_R
    RUN_NAME

    script:
    """

    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'filtered_scafstats_${RUN_NAME}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'all_scafstats_${RUN_NAME}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'topHit_scafstats_${RUN_NAME}.csv'

    if [ ! -d ${params.outdir}analysis ]; then
    mkdir -p ${params.outdir}analysis;
    fi;

    cp filtered_scafstats_${RUN_NAME}.csv ${params.outdir}analysis/
    cp all_scafstats_${RUN_NAME}.csv ${params.outdir}analysis/
    cp topHit_scafstats_${RUN_NAME}.csv ${params.outdir}analysis/

    ls -latr
    Rscript --vanilla ${MERGE_STATS_R} \'${RUN_NAME}' \'${params.outdir}\' \'${params.outdir}\'
    """
}
}
}
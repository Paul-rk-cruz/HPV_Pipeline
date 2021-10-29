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
// Nextflow dsl2
nextflow.enable.dsl=2
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
        nextflow run FILE_PATH/HPV_Pipeline/main.nf --input PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR --SingleEnd
    Valid CLI Arguments:
    REQUIRED:
      --input                       Path to input fastq.gz folder
      --outdir                      The output directory where the results will be saved
      --runName                     Specifies the run name
      --ref                         Specifies the HPV reference multifasta type (all OR highRisk)
      --singleEnd                   Specifies single-end fastq.gz sequence files
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
params.runName = false
params.singleEnd = false
params.SETTING = "2:30:10:1:true"
params.LEADING = "3"
params.TRAILING = "3"
params.SWINDOW = "4:20"
params.MINLEN = "75"
// Check if input is set
if (params.input == false) {
    println( "Must provide an input directory with --input") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.input.endsWith("/")){
   params.input = "${params.input}/"
}
// if OUTDIR not set
if (params.outdir == false) {
    println( "Must provide an output directory with --outdir") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}
// Setup file paths
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
// Setup Run Name
if(params.runName) {
    runName = params.runName
}
// Setup reference multifasta
if (params.ref == "all") {
    REF_HPV = file("${baseDir}/ref_fasta/hpvAll.fasta")
}
else if (params.ref == "highRisk") {
    REF_HPV = file("${baseDir}/ref_fasta/hpvHighRisk.fasta")
}
// R script path
MERGE_STATS_R = file("${baseDir}/scripts/analysis.R")
// Workflow display header
def hpvheader() {
    return """
    """.stripIndent()
}
// log files header
// log.info hpvheader()
log.info "_______________________________________________________________________________"
log.info " Human Papilloma Virus Pipeline :  v${version}"
log.info "_______________________________________________________________________________"
def summary = [:]
summary['Configuration Profile:'] = workflow.profile
summary['Run Name:']           	  = params.runName
summary['Current directory path:']        = "$PWD"
summary['HRV Pipeline directory path:']          = workflow.projectDir
summary['Input directory path:']               = params.input
summary['Output directory path:']          = params.outdir
summary['Work directory path:']         = workflow.workDir
summary['Sequence type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Reference type:']           	  = params.ref
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
if (params.singleEnd) {
summary['Trimmomatic adapters:'] = ADAPTERS_SE
} else {
summary['Trimmomatic adapters:'] = 'Paired-end'
}
summary["Trimmomatic read length (minimum):"] = params.MINLEN
summary["Trimmomatic Setting:"] = params.SETTING
summary["Trimmomatic Sliding Window:"] = params.SWINDOW
summary["Trimmomatic Leading:"] = params.LEADING
summary["Trimmomatic Trailing:"] = params.TRAILING
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "_______________________________________________________________________________"


//
// Import processes
// 

include { Trimming } from './modules.nf'
include { Aligning } from './modules.nf'
include { Bam_Sorting } from './modules.nf'
include { Analysis } from './modules.nf'

// Create channel for input reads: single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.input}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.input} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.input}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    if(params.singleEnd == true) {
    Trimming (
        input_read_ch, 
        ADAPTERS_SE,
        params.MINLEN,
        params.SETTING, 
        params.LEADING,
        params.TRAILING,
        params.SWINDOW,
    )
    Aligning (
        Trimming.out[0],
        REF_HPV
    )
    Bam_Sorting (
        Aligning.out[0]
    )
    Analysis (
        Aligning.out[0].collect(),
        Aligning.out[1].collect(),
        MERGE_STATS_R,
        params.runName
    )
    }
}
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
        nextflow run FILE_PATH/HPV_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR
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
vapid_rhinovirus_sbt = file("${baseDir}/vapid/rhinovirus.sbt")
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
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*            VIRAL MULTI-FASTA REFERENCES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Setup MULTIFASTA Reference parameters.
params.skipTrimming = false
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
summary['Fastq type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
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
process Trim_Reads {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

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
    // container "docker.io/paulrkcruz/hrv-pipeline:latest" 
    errorStrategy 'retry'
    maxRetries 3
    // echo true

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_final_summary.csv") from Trim_out_SE

    output:
        tuple val(base), file("${base}_summary2.csv") into Everything_ch
        tuple val (base), file("*") into Dump_map1_ch

    publishDir "${params.outdir}Summary", mode: 'copy', pattern:'*.csv*'


    script:

    """
    #!/bin/bash


    cp ${base}_summary.csv ${base}_final_summary.csv

    """
}
}

// HPV Pipeline Workflow Processes
/*
 * STEP 1: Trim_Reads
 * Trimming of low quality and short NGS sequences.
 */
process Trimming {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        file R1 //from input_read_ch
        file ADAPTERS_SE
        val MINLEN
        val SETTING
        val LEADING
        val TRAILING
        val SWINDOW
    output:
        tuple env(base),file("*.trimmed.fastq.gz"), file("*summary.csv")// into Trimming_ch

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
    echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary.csv")// from Trimming_ch
        file REF_HPV

    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_stats1.csv"), file("${base}_hpvAll.sam"), file("${base}_hpvAll_scafstats.txt"), file("${base}_hpvAll_covstats.txt")// into Align_ch
    
    publishDir "${params.outdir}bbmap_scaf_stats", mode: 'copy', pattern:'*_hpvAll_scafstats.txt*'
    publishDir "${params.outdir}bbmap_cov_stats", mode: 'copy', pattern:'*_hpvAll_covstats.txt*'    

    script:

    """
    #!/bin/bash

    bbmap.sh in=${base}.trimmed.fastq.gz ref=${REF_HPV} outm=${base}_hpvAll.sam outu=${base}_nope.sam maxindel=9 ambiguous=best threads=${task.cpus} scafstats=${base}_hpvAll_scafstats.txt covstats=${base}_hpvAll_covstats.txt -Xmx6g > bbmap_out.txt 2>&1 

    cp ${base}_summary.csv ${base}_stats1.csv

    """
}
// /usr/local/bin/bbmap.sh
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
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_stats1.csv"), file("${base}_hpvAll.sam"), file("${base}_hpvAll_scafstats.txt"), file("${base}_hpvAll_covstats.txt")// from Align_ch

    output:
        tuple val(base), file("${base}_hpvAll.sorted.bam"), file("${base}_trim_stats.csv")//  into Analysis_ch   

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
    container "docker.io/cgrlab/tidyverse:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
    file("${base}_hpvAll_scafstats.txt")// from Bbmap_scaf_stats_ch.collect()     
    file("${base}_hpvAll_covstats.txt")// from Bbmap_cov_stats_ch.collect()   
    file MERGE_STATS_R
    val runName

    script:
    """
    #!/bin/bash

    ls -latr

    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'filtered_scafstats_${runName}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'all_scafstats_${runName}.csv'
    echo Sample, Reference, Percent_Unambiguous_Reads, x, Percent_Ambiguous_Reads, x, Unambiguous_Reads, Ambiguous Reads, Assigned Reads, x> 'topHit_scafstats_${runName}.csv'

    if [ ! -d ${params.outdir}analysis ]; then
    mkdir -p ${params.outdir}analysis;
    fi;

    cp filtered_scafstats_${runName}.csv ${params.outdir}analysis/
    cp all_scafstats_${runName}.csv ${params.outdir}analysis/
    cp topHit_scafstats_${runName}.csv ${params.outdir}analysis/

    Rscript --vanilla ${MERGE_STATS_R} \'${runName}' \'${params.outdir}\'
    """
}
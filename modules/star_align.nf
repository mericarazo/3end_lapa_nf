process STAR_ALIGN {
    tag "$sample"
    cpus 4
    memory '40 GB'
    time '24h'
    conda 'bioconda::star=2.7.11b'
    publishDir "${params.outdir}/03_star", mode: params.publish_mode, pattern: '*.Log.final.out'

    input:
    tuple val(sample), path(fq)
    path index

    output:
    tuple val(sample), path("${sample}.bam"), emit: bam
    path "${sample}.Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${fq} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM 10000000000 \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outFilterMatchNminOverLread 0.33 \\
        --alignEndsType Local \\
        --outSAMattributes NH HI AS nM NM MD

    mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.bam
    """
}

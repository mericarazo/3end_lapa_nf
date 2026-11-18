process INTERNAL_PRIMING {
    tag "$sample"
    cpus 2
    memory '12 GB'
    conda 'bioconda::samtools=1.20 bioconda::bedtools=2.31.1 conda-forge::gawk'
    publishDir "${params.outdir}/04_internal_priming", mode: params.publish_mode, pattern: '*.{bam}*'

    input:
    tuple val(sample), path(bam)
    path fasta
    path fai
    path sizes

    output:
    tuple val(sample), path("${sample}.ip.bam"), path("${sample}.ip.bam.bai"), emit: bam
    path "${sample}.ip.flagstat", emit: flagstat
    path "${sample}_remove_reads.txt"

    script:
    """
    internal_priming_filter.sh ${bam} ${sample} ${params.ip_window} ${params.ip_max_a} ${params.ip_max_run} ${fasta} ${sizes}
    """
}

process DEDUP {
    tag "$sample"
    cpus 2
    memory '24 GB'
    conda 'bioconda::umi_tools=1.1.5 bioconda::samtools=1.20'
    publishDir "${params.outdir}/05_dedup", mode: params.publish_mode, pattern: '*.{bam}*'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.dedup.bam"), path("${sample}.dedup.bam.bai"), emit: bam
    path "${sample}.dedup.flagstat", emit: flagstat
    path "${sample}.dedup.log"

    script:
    """
    umi_tools dedup -I ${bam} -S ${sample}.dedup.bam -L ${sample}.dedup.log
    samtools index ${sample}.dedup.bam
    samtools flagstat ${sample}.dedup.bam > ${sample}.dedup.flagstat
    """
}

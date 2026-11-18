process UMI_EXTRACT {
    tag "$sample"
    cpus 1
    memory '8 GB'
    conda 'bioconda::umi_tools=1.1.5'

    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("${sample}.umi.fq.gz"), emit: reads
    path "${sample}.umi_extract.log", emit: log

    // Solo se usa R1 aguas abajo (STAR solo con R1), asi que no se procesa R2.
    script:
    """
    umi_tools extract \\
        --bc-pattern=${params.umi_pattern} \\
        -I ${fq} \\
        -S ${sample}.umi.fq.gz \\
        -L ${sample}.umi_extract.log
    """
}

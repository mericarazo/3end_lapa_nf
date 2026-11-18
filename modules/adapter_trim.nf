process ADAPTER_TRIM {
    tag "$id"
    cpus 4
    memory '12 GB'
    conda 'bioconda::cutadapt=4.9'
    publishDir "${params.outdir}/01_adapter_trim", mode: params.publish_mode, pattern: '*.json'

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_1.trim.fq.gz"), path("${id}_2.trim.fq.gz"), emit: reads
    path "${id}.cutadapt.json", emit: json

    script:
    """
    cutadapt \\
        -j ${task.cpus} \\
        -g ${params.adapter_r1} \\
        -G ${params.adapter_r2} \\
        -e 0.1 \\
        -m 25 \\
        -o ${id}_1.trim.fq.gz \\
        -p ${id}_2.trim.fq.gz \\
        --json ${id}.cutadapt.json \\
        ${r1} ${r2}
    """
}

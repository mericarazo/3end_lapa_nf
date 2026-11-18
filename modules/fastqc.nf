process FASTQC {
    tag "$sample"
    cpus 2
    memory '4 GB'
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_mode

    input:
    tuple val(sample), path(fq)

    output:
    path '*_fastqc.zip', emit: zip
    path '*_fastqc.html'

    script:
    """
    fastqc -t ${task.cpus} ${fq}
    """
}

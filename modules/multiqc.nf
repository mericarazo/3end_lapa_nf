process MULTIQC {
    cpus 1
    memory '4 GB'
    conda 'bioconda::multiqc=1.25'
    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    path files

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc . -f
    """
}

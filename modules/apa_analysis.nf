process APA_ANALYSIS {
    cpus 1
    memory '8 GB'
    conda 'conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy conda-forge::scipy conda-forge::matplotlib'
    publishDir "${params.outdir}/07_apa_analysis", mode: params.publish_mode

    input:
    path clusters
    path counts
    path usage
    path gtf

    output:
    path "${params.pool_name}.*"

    script:
    def extra = params.require_detected_all ? '--require-detected-all' : ''
    """
    apa_analysis.py \\
        --clusters ${clusters} \\
        --counts ${counts} \\
        --usage ${usage} \\
        --gtf ${gtf} \\
        --case ${params.treatment} \\
        --control ${params.control} \\
        --min-mean ${params.min_mean} \\
        --min-usage ${params.min_usage} \\
        --padj ${params.padj} \\
        --delta ${params.delta} \\
        --highlight '${params.highlight}' \\
        --outprefix ${params.pool_name} ${extra}
    """
}

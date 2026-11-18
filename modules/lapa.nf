process LAPA {
    cpus 2
    memory '40 GB'
    time '48h'
    // el entorno con 'lapa' (pip install lapa) se indica en el perfil / config: process.withName:LAPA.conda
    publishDir "${params.outdir}/06_lapa", mode: params.publish_mode

    input:
    tuple val(samples), val(datasets), path(bams), path(bais)
    path fasta
    path fai
    path gtf
    path sizes

    output:
    path "lapa_out/${params.lapa_clusters_file}", emit: clusters
    path "lapa_out/${params.lapa_counts_file}", emit: counts
    path "lapa_out/${params.lapa_usage_file}", emit: usage
    path "alignment.csv"
    path "lapa_out"

    script:
    def names = bams.collect { it.name }
    def rows = [samples, datasets, names].transpose().collect { it.join(',') }
    def lines = rows.collect { "echo '${it}' >> alignment.csv" }.join('\n')
    """
    echo 'sample,dataset,path' > alignment.csv
    ${lines}
    mkdir -p lapa_out

    lapa \\
        --alignment alignment.csv \\
        --fasta ${fasta} \\
        --annotation ${gtf} \\
        --chrom_sizes ${sizes} \\
        --output_dir lapa_out \\
        --min_replication_rate ${params.lapa_min_replication_rate} \\
        --cluster_window ${params.lapa_cluster_window}

    lapa_analysis.py
    
    for f in ${params.lapa_clusters_file} ${params.lapa_counts_file} ${params.lapa_usage_file}; do
        if [ ! -s "lapa_out/\$f" ]; then
            echo "ERROR: LAPA no ha generado lapa_out/\$f. Ficheros generados:" >&2
            ls -la lapa_out >&2
            exit 1
        fi
    done
    """
}

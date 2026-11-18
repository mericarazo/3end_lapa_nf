process DEMULTIPLEX {
    tag "$id"
    cpus 4
    memory '16 GB'
    conda 'bioconda::cutadapt=4.9'
    publishDir "${params.outdir}/02_demultiplex", mode: params.publish_mode

    input:
    tuple val(id), path(r1), path(r2)
    path barcodes

    output:
    path '*_R1.fq.gz', emit: r1      // incluye unknown_R1.fq.gz, se filtra en el workflow
    path "${id}.demux.log", emit: log

    // El barcode esta en R2: se pasa R2 como primera lectura (el -g se aplica a ella) y
    // sus salidas se nombran *_R2 / *_R1 como en el pipeline original.
    script:
    """
    cutadapt \\
        -j ${task.cpus} \\
        -g ^file:${barcodes} \\
        -e 1 \\
        --no-indels \\
        -o {name}_R2.fq.gz \\
        -p {name}_R1.fq.gz \\
        --untrimmed-output unknown_R2.fq.gz \\
        --untrimmed-paired-output unknown_R1.fq.gz \\
        ${r2} ${r1} > ${id}.demux.log
    """
}

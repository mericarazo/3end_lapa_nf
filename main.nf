#!/usr/bin/env nextflow
/*
 * APA 3'-seq (UMI) pipeline para UN pool:
 *   cutadapt (adaptadores) -> cutadapt (demultiplex) -> umi_tools extract -> STAR (R1)
 *   -> filtro internal priming -> umi_tools dedup -> LAPA -> filtros + Fisher + proximal/distal + volcano
 */
nextflow.enable.dsl = 2

// ------------------------------------------------------------------------------------------
// Comprobacion de parametros
// ------------------------------------------------------------------------------------------
def required = ['reads', 'barcodes', 'star_index', 'genome_fasta', 'genome_sizes', 'lapa_gtf'] + (params.apa_analysis ? ['treatment', 'control'] : [])

def missing = required.findAll { !params[it] }
if (missing) {
    error "Faltan parametros obligatorios: ${missing.collect { '--' + it }.join(', ')}"
}

// ------------------------------------------------------------------------------------------
// Modulos 
// ------------------------------------------------------------------------------------------
include { ADAPTER_TRIM } from './modules/adapter_trim'
include { DEMULTIPLEX } from './modules/demultiplex'
include { UMI_EXTRACT } from './modules/umi_extract'
include { FASTQC } from './modules/fastqc'
include { STAR_ALIGN } from './modules/star_align'
include { INTERNAL_PRIMING } from './modules/internal_priming'
include { DEDUP } from './modules/dedup'
include { LAPA } from './modules/lapa'
include { APA_ANALYSIS } from './modules/apa_analysis'
include { MULTIQC } from './modules/multiqc'

// ------------------------------------------------------------------------------------------
// Workflow
// ------------------------------------------------------------------------------------------
workflow {
    barcodes = file(params.barcodes, checkIfExists: true)
    index    = file(params.star_index, checkIfExists: true)
    fasta    = file(params.genome_fasta, checkIfExists: true)
    fai      = file(params.genome_fai ?: "${params.genome_fasta}.fai", checkIfExists: true)
    sizes    = file(params.genome_sizes, checkIfExists: true)
    gtf      = file(params.lapa_gtf, checkIfExists: true)
    gtf_names = file(params.gtf_names ?: params.lapa_gtf, checkIfExists: true)

    // un unico pool: <prefijo>_1.fq.gz / <prefijo>_2.fq.gz
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true).map { id, f -> tuple(id, f[0], f[1]) }

    ADAPTER_TRIM(reads)
    DEMULTIPLEX(ADAPTER_TRIM.out.reads, barcodes)

    // una muestra por barcode (se descartan las no asignadas)
    samples = DEMULTIPLEX.out.r1
        .flatten()
        .filter { !it.name.startsWith('unknown_') }
        .map { f -> tuple(f.name.replaceFirst(/_R1\.fq\.gz$/, ''), f) }

    UMI_EXTRACT(samples)
    FASTQC(UMI_EXTRACT.out.reads)
    STAR_ALIGN(UMI_EXTRACT.out.reads, index)
    INTERNAL_PRIMING(STAR_ALIGN.out.bam, fasta, fai, sizes)
    DEDUP(INTERNAL_PRIMING.out.bam)

    lapa_in = DEDUP.out.bam
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { l ->
            def ds = l.collect { params.lapa_dataset == 'single' ? params.pool_name : it[0].replaceFirst(/^R\d+_/, '') }
            tuple(l.collect { it[0] }, ds, l.collect { it[1] }, l.collect { it[2] })
        }

    LAPA(lapa_in, fasta, fai, gtf, sizes)
    if (params.apa_analysis) {
        APA_ANALYSIS(LAPA.out.clusters, LAPA.out.counts, LAPA.out.usage, gtf_names)
    }

    MULTIQC(
        FASTQC.out.zip
            .mix(STAR_ALIGN.out.log, INTERNAL_PRIMING.out.flagstat, DEDUP.out.flagstat)
            .collect()
    )
}

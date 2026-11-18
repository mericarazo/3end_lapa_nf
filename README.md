# Pipeline APA 3'-seq (UMI) → LAPA → uso diferencial de sitios (un pool)

```
FASTQ pool (R1/R2)
  └ ADAPTER_TRIM     cutadapt -e 0.1 -m 25
    └ DEMULTIPLEX    cutadapt -g ^file:barcodes -e 1 --no-indels   (barcode en R2)
      └ UMI_EXTRACT  umi_tools extract NNNNN                       (solo R1)
        ├ FASTQC
        └ STAR_ALIGN (R1, multimap 10, mismatch 4 %, score/match 33 %, Local)
          └ INTERNAL_PRIMING  (A en ventana de 10 nt: lectura + downstream)
            └ DEDUP    umi_tools dedup
              └ LAPA   (todas las muestras juntas)
                └ APA_ANALYSIS  filtros → ∆dPUR → Fisher + BH → pares proximal/distal → volcano
MULTIQC  (FastQC + STAR + flagstat)
```

## Estructura

```
main.nf                 workflow (solo conecta los modulos)
nextflow.config         parametros y perfiles
modules/                un proceso por fichero: adapter_trim, demultiplex, umi_extract, fastqc,
                        star_align, internal_priming, dedup, lapa, apa_analysis, multiqc
bin/                    scripts que usan los modulos (internal_priming_filter.sh, apa_analysis.py)
```

## Lanzar
Ejecución en el cluster (SLURM + entornos conda):

```bash
nextflow run apa_nf -profile biogi -resume -params-file input_params.yaml
```

Hace falta el `.fai` del genoma junto al FASTA (`samtools faidx genoma.fa`).
Lanzar con `nextflow run` desde un nodo de login o con un `sbatch` que solo lance Nextflow.

Ejecución en local:

```bash
nextflow run apa_nf -profile standard -resume -params-file input_params.yaml
```

Se han creado los contenedores de docker necesarios para poder trabajar con este pipeline. Con tener docker instalado es suficiente. 


## Convenciones que debes respetar

- Los nombres del FASTA de barcodes son los nombres de muestra (`R1_NT`, `R2_TDP43`, …).
  La condición debe estar incluida en el nombre de las muestras; `--treatment` y `--control` deben coincidir con ella.
- Un pool por ejecución.

## Salidas (`--outdir`)

| Carpeta | Contenido |
|---|---|
| `01_adapter_trim`, `02_demultiplex` | JSON / logs de cutadapt (nº de lecturas por barcode) |
| `03_star` | `Log.final.out` por muestra |
| `04_internal_priming` | flagstat y lista de lecturas eliminadas por muestra |
| `05_dedup` | BAM + índice finales por muestra |
| `06_lapa` | `alignment.csv` generado y salida completa de LAPA |
| `07_apa_analysis` | `*.filter_summary.tsv`, `*.fisher_events_individual.tsv`, `*.fisher_events_paired.tsv`, `*.volcano.png/pdf` |
| `qc` | FastQC y `multiqc_report.html` |
| `pipeline_info` | timeline, report y trace de Nextflow |

## Parámetros principales

Ver `nextflow.config`. Los de filtrado y estadística son `--min_mean 5`, `--min_usage 0.0025`,
`--padj 0.05`, `--delta 0.05`, `--highlight STMN2`.

## Cosas que conviene revisar antes de fiarte de los resultados

1. **Formato del CSV de LAPA.** El pipeline escribe `sample,dataset,path`, con `dataset` = condición
   (`--lapa_dataset condition`) o un único dataset (`--lapa_dataset single`). Compáralo con la cabecera de tu
   `lapa_pool1_dedup.csv`, que no tengo, y ajusta si difiere.
2. **Nombres de los ficheros de LAPA.** El análisis espera `polyA_clusters.bed`, `lapa_counts.tsv` y
   `lapa_usage.tsv` dentro de `lapa_out/` (los que lee tu notebook). Si tu versión de LAPA los llama de otra forma,
   cámbialo con `--lapa_clusters_file / --lapa_counts_file / --lapa_usage_file`; el proceso falla con un mensaje y el
   listado de lo que sí generó.
3. **Umbral de usage.** El notebook usa `usage >= 0.0025` (0,25 %), aunque el texto de los notebooks y del póster
   dice 2,5 % o 5 %. Por defecto reproduzco el código que ejecutaste. Con `--min_usage 0.025` tendrías el filtro del
   texto.
4. **Filtro «presente en todas las muestras de una condición».** En el notebook los NaN se rellenan con 0 antes de
   aplicarlo, así que no filtraba nada (38 775 → 38 775). Por defecto se mantiene así; con
   `--require_detected_all` se exige counts > 0 en todas las réplicas de una condición.
5. **Fisher con réplicas agregadas.** Los counts de las 3 réplicas se suman por condición, así que el test no ve la
   variabilidad entre réplicas y los p-valores salen optimistas.
6. **STAR solo con R1**, y UMI extraído solo de R1 (R2 solo sirve para el barcode). Se omite `--outReadsUnmapped Fastx`
   porque no se usa.
7. **Entornos conda.** El perfil `biogi` reutiliza tus entornos (`trimming`, `STAR`, `biotools`, `QControl`, `lapa`).
   Sin ese perfil, cada proceso usa un paquete bioconda con `-profile conda`; LAPA sigue necesitando un entorno
   con `lapa` (`process.withName: 'LAPA' { conda = '/ruta/a/env' }`), y el entorno de `APA_ANALYSIS` lo crea
   Nextflow (pandas, numpy, scipy, matplotlib).

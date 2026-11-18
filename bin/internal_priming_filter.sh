#!/usr/bin/env bash
# Filtro de internal priming (version "double" del pipeline original).
# Se elimina la lectura si, en cualquiera de las dos ventanas de W nt:
#   (a) los ultimos W nt ALINEADOS de la lectura, o
#   (b) los W nt genomicos justo DOWNSTREAM de su extremo 3',
# (en orientacion del transcrito) hay > MAX_A adenosinas o una racha de > MAX_RUN A.
#
# Uso: internal_priming_filter.sh in.bam prefix W MAX_A MAX_RUN genome.fa chrom.sizes
set -euo pipefail

bam=$1; base=$2; W=$3; MAX_A=$4; MAX_RUN=$5; hg38=$6; genome_sizes=$7

# 1) ultimos W nt del alineamiento, sin soft-clip, en orientacion del transcrito
samtools view -F 4 "${bam}" | awk -v OFS='\t' -v W="${W}" '
function rc(s,   i, c, r) {
	r = ""
	for (i = length(s); i >= 1; i--) {
		c = substr(s, i, 1)
		r = r (c=="A" ? "T" : c=="T" ? "A" : c=="C" ? "G" : c=="G" ? "C" : "N")
	}
	return r
}
{
	ls=0; rs=0
	if (match($6,/^[0-9]+S/)) ls=substr($6,1,RLENGTH-1)+0
	if (match($6,/[0-9]+S$/)) rs=substr($6,RSTART,RLENGTH-1)+0
	aln=substr($10,ls+1,length($10)-ls-rs)
	if (int($2/16)%2==0) tail=substr(aln,length(aln)-W+1)
	else tail=rc(substr(aln,1,W))
	print $1, tail, $2, $6
}' > "${base}_readend.tsv"

# 2) W nt genomicos downstream del extremo 3' alineado, en orientacion del transcrito
bedtools bamtobed -i "${bam}" \
| awk 'BEGIN{OFS="\t"}
	$6=="+"{print $1,$3,$3+1,$4,$5,$6}
	$6=="-" && $2>0{print $1,$2-1,$2,$4,$5,$6}' \
| bedtools slop -i - -g "${genome_sizes}" -l 0 -r $((W-1)) -s \
| bedtools getfasta -fi "${hg38}" -bed - -s -tab -nameOnly \
> "${base}_downstream.tsv"

# 3) lecturas a eliminar
cat "${base}_readend.tsv" "${base}_downstream.tsv" \
| awk -F'\t' -v MAXA="${MAX_A}" -v MAXRUN="${MAX_RUN}" '
BEGIN{ run=sprintf("%" (MAXRUN+1) "s",""); gsub(/ /,"A",run) }
{
	id=$1; sub(/\([+-]\)$/,"",id)
	s=toupper($2)
	n=gsub(/A/,"A",s)
	if (n>MAXA || index(s,run)) print id
}' | sort -u > "${base}_remove_reads.txt"

echo "Reads to remove: $(wc -l < "${base}_remove_reads.txt")" >&2

# 4) quitar esas lecturas del BAM (samtools >= 1.12): -N selecciona, -U guarda las no seleccionadas
samtools view -b -N "${base}_remove_reads.txt" -U "${base}.ip.bam" "${bam}" > /dev/null
samtools index "${base}.ip.bam"
samtools flagstat "${base}.ip.bam" > "${base}.ip.flagstat"
rm -f "${base}_readend.tsv" "${base}_downstream.tsv"

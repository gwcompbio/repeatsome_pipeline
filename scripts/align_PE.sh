#!/bin/bash

SN="align_PE.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && ref="$1"
[[ -n "$2" ]] && fastq1="$2"
[[ -n "$3" ]] && fastq2="$3"

#--- Check that fastq files exist
[[ ! -e "$fastq1" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq1 does not exist" && exit 1
[[ ! -e "$fastq2" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq2 does not exist" && exit 1

#--- Check that reference index exists
[[ ! -e "$ref.1.bt2" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1

#--- Check file type
[[ $(head -c1 $fastq1) == '@' ]] && ftype="q" || ftype="f"

#--- Get output file names
[[ -n "$4" ]] && out="$4" || out=${fastq1%%.*}.aligned

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Ref index: $ref"
echo "[---$SN---] ($(date)) File type: $ftype"
echo "[---$SN---] ($(date)) Fastq1:    $fastq1"
echo "[---$SN---] ($(date)) Fastq2:    $fastq2"
echo "[---$SN---] ($(date)) Out:       $out"

#--- Load modules
module load bowtie2

#--- Start the timer
t1=$(date +"%s")

#--- Alignment
echo "[---$SN---] ($(date)) Starting bowtie2."
bowtie2 -p $(nproc) \
  --very-sensitive-local -k 100 --score-min L,0,1.6 \
  -${ftype} \
  -x $ref \
  -1 $fastq1 \
  -2 $fastq2 \
  -S $out.sam

echo "[---$SN---] ($(date)) bowtie2 completed successfully."

#--- Sorting and compression
echo "[---$SN---] ($(date)) Sorting and compressing alignment."
module load samtools
samtools view -uS $out.sam | samtools sort - $out
samtools index $out.bam

echo "[---$SN---] ($(date)) samtools completed successfully."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

#!/bin/bash

SN="filter_PE.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && ref="$1"
[[ -n "$2" ]] && fastq1="$2"

#--- Check that fastq files exist
[[ ! -e "$fastq1" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq1 does not exist" && exit 1

#--- Check that reference index exists
[[ ! -e "$ref.1.bt2" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1

#--- Check file type
[[ $(head -c1 $fastq1) == '@' ]] && ftype="q" || ftype="f"

#--- Get output file names
if  [[ "$ftype" == "q" ]]; then
  [[ -n "$3" ]] && out1="$3" || out1=${fastq1%.*}.filt.fq
else
  [[ -n "$3" ]] && out1="$3" || out1=${fastq1%.*}.filt.fa
fi

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Ref index: $ref"
echo "[---$SN---] ($(date)) File type: $ftype"
echo "[---$SN---] ($(date)) Fastq1:    $fastq1"
echo "[---$SN---] ($(date)) Out1:      $out1"

#--- Load modules
module load bowtie2
module load cbiC1

#--- Start the timer
t1=$(date +"%s")

#--- Align and write directly to FASTQ/A
echo "[---$SN---] ($(date)) Starting bowtie2."
bowtie2 -p $(nproc) --very-fast \
  -${ftype} \
  -x $ref \
  -U $fastq1 \
  --un $out1 > /dev/null

echo "[---$SN---] ($(date)) bowtie2 completed successfully."  

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

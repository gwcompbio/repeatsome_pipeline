#!/bin/bash

SN="qc_SE.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && fastq1="$1"
[[ -n "$2" ]] && out1="$2" || out1=${fastq1%.*}.qc.fq

#--- Check that fastq files exist
[[ ! -e "$fastq1" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq1 does not exist" && exit 1

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Fastq1:    $fastq1"
echo "[---$SN---] ($(date)) Out1:      $out1"

#--- Load modules
module load cbiC1
module load seqtk
module load fastx

#--- Start the timer
t1=$(date +"%s")

#--- Filtering
echo "[---$SN---] ($(date)) Starting cbi_fastq_filter."
cbi_fastq_filter -v $fastq1 > $out1
echo "[---$SN---] ($(date)) cbi_fastq_filter completed successfully."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

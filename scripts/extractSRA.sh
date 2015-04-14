#!/bin/bash

SN="extractSRA.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && srafile="$1"

#--- Check that SRA files exist
[[ ! -e "$srafile" ]] && echo "[---$SN---] ($(date)) FAILED: file $srafile does not exist" && exit 1

#--- Get output file names
outdir=$(dirname $srafile)

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) SRA file:  $srafile"
echo "[---$SN---] ($(date)) Outdir:    $outdir"

#--- Load modules
module load sratoolkit

#--- Start the timer
t1=$(date +"%s")

#--- Extract the files from SRA
echo "[---$SN---] ($(date)) Extracting from SRA."
fastq-dump -O $outdir --split-3 -F -R --defline-qual "+"  $srafile
echo "[---$SN---] ($(date)) Extracting complete."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

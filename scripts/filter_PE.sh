#!/bin/bash

SN="filter_PE.sh"

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
if  [[ "$ftype" == "q" ]]; then
  [[ -n "$4" ]] && out1="$4" || out1=${fastq1%.*}.filt.fq
  [[ -n "$5" ]] && out2="$5" || out2=${fastq2%.*}.filt.fq
else
  [[ -n "$4" ]] && out1="$4" || out1=${fastq1%.*}.filt.fa
  [[ -n "$5" ]] && out2="$5" || out2=${fastq2%.*}.filt.fa
fi

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Ref index: $ref"
echo "[---$SN---] ($(date)) File type: $ftype"
echo "[---$SN---] ($(date)) Fastq1:    $fastq1"
echo "[---$SN---] ($(date)) Fastq2:    $fastq2"
echo "[---$SN---] ($(date)) Out1:      $out1"
echo "[---$SN---] ($(date)) Out2:      $out2"

#--- Load modules
module load bowtie2
module load cbiC1

#--- Start the timer
t1=$(date +"%s")

#--- Align and write to SAM file
echo "[---$SN---] ($(date)) Starting bowtie2."
tmpfile=$(mktemp)
bowtie2 -p $(nproc) --very-fast \
  -${ftype} \
  -x $ref \
  -1 $fastq1 \
  -2 $fastq2 \
  --no-head > $tmpfile.sam

echo "[---$SN---] ($(date)) bowtie2 completed successfully."

#--- Extract unaligned reads from SAM file
echo "[---$SN---] ($(date)) Extracting unaligned reads."
if  [[ "$ftype" == "q" ]]; then
  cat $tmpfile.sam | grep 'YT:Z:UP' | perl -lane 'print "\@$F[0]\n$F[9]\n+\n$F[10]"' | \
    deinterleave_fastq  $out1 $out2
else
  cat $tmpfile.sam | grep 'YT:Z:UP' | perl -lane 'print ">$F[0]\n$F[9]"' | \
    deinterleave_fasta $out1 $out2
fi
rm -f $tmpfile*
echo "[---$SN---] ($(date)) reads extracted successfully."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

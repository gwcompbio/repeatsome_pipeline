#!/bin/bash

SN="repenrich_PE.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && ref="$1"
[[ -n "$2" ]] && bedfile="$2"
[[ -n "$3" ]] && setupdir="$3"
[[ -n "$4" ]] && fastq1="$4"
[[ -n "$5" ]] && fastq2="$5"

#--- Check that fastq files exist
[[ ! -e "$fastq1" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq1 does not exist" && exit 1
[[ ! -e "$fastq2" ]] && echo "[---$SN---] ($(date)) FAILED: file $fastq2 does not exist" && exit 1

#--- Check that reference index exists
#--- NOTE: this is a bowtie1 index!
[[ ! -e "$ref.1.ebwt" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1

#--- Check file type
[[ $(head -c1 $fastq1) == '@' ]] && ftype="q" || ftype="f"

#--- Get output file names
[[ -n "$6" ]] && outprefix="$6" || outprefix=$(basename $(dirname $bedfile))
[[ -n "$7" ]] && out="$7" || out=$(dirname $fastq1)/RepEnrich
mkdir -p $out

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Ref index: $ref"
echo "[---$SN---] ($(date)) BED file:  $bedfile"
echo "[---$SN---] ($(date)) Setup dir: $setupdir"
echo "[---$SN---] ($(date)) File type: $ftype"
echo "[---$SN---] ($(date)) Fastq1:    $fastq1"
echo "[---$SN---] ($(date)) Fastq2:    $fastq2"
echo "[---$SN---] ($(date)) Prefix:    $outprefix"
echo "[---$SN---] ($(date)) Out:       $out"

#--- Load modules
module load bowtie
module load bedtools
module load samtools

#--- Start the timer
t1=$(date +"%s")

#--- Alignment
echo "[---$SN---] ($(date)) Starting alignment."
bowtie --threads $(nproc) --time \
  -${ftype} \
  -m 1 \
  --sam \
  -X 600 \
  --max $out/multimap \
  --chunkmbs 512 \
  $ref \
  -1 $fastq1 \
  -2 $fastq2 \
  $out/unique.sam 2> $out/bowtie.out

echo "[---$SN---] ($(date)) Alignment completed successfully."
cat $out/bowtie.out

#--- Post-process alignment
echo "[---$SN---] ($(date)) Post-processing alignment."

samtools view -uS $out/unique.sam | samtools sort - $out/unique
samtools index $out/unique.bam

nproc=$(grep '# reads processed: ' $out/bowtie.out | awk -F' ' '{print $NF}')
nfail=$(grep '# reads that failed to align:' $out/bowtie.out | awk -F' ' '{print $(NF-1)}')
nmapped=$(( $nproc - $nfail ))

if [[ "$ftype" == "f" ]]; then
  module load cbiC1
  inlineFasta2Fastq < $out/multimap_1 > $out/multimap_1.fastq
  inlineFasta2Fastq < $out/multimap_2 > $out/multimap_2.fastq
else
  mv $out/multimap_1 $out/multimap_1.fastq
  mv $out/multimap_2 $out/multimap_2.fastq
fi

echo "[---$SN---] ($(date)) Post-processing complete."
echo "[---$SN---] ($(date)) $nmapped reads mapped"

#--- Run RepEnrich
echo "[---$SN---] ($(date)) Running RepEnrich."
module load RepEnrich

RepEnrich.py \
  $bedfile \
  $out \
  $outprefix \
  $setupdir \
  $out/multimap_1.fastq \
  --fastqfile2 $out/multimap_2.fastq \
  $out/unique.bam \
  --is_bed TRUE \
  --cpus 16 \
  --pairedend TRUE

echo "[---$SN---] ($(date)) RepEnrich completed successfully."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

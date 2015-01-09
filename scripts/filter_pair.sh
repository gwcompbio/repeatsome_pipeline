#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -p short,defq
#SBATCH -N 1

module load cbiC1
module load seqtk
module load fastx
module load bowtie2

fastq1=$(ls $prefix* | sed -n '1p')
fastq2=$(ls $prefix* | sed -n '2p')

echo "Fastq1: $fastq1"
echo "Fastq2: $fastq2"

dest=/home/bendall/fs/Projects/HERV_pilot/data/cleaned
mkdir -p $dest

tmpfile=$(mktemp)
time cbi_fastq_filter -v $fastq1 $fastq2 | deinterleave_fastq $tmpfile.1.fastq $tmpfile.2.fastq

ref=/lustre/groups/cbi/shared/References/HIV-1/NC_001802.1/Sequence/Bowtie2Index/genome
bowtie2 \
  -x $ref \
  -1 $tmpfile.1.fastq \
  -2 $tmpfile.2.fastq \
  --no-head --reorder > $tmpfile.sam

cat $tmpfile.sam | grep 'YT:Z:UP' | perl -lane 'print "\@$F[0]\n$F[9]\n+\n$F[10]"' | \
  deinterleave_fastq ${dest}/$(basename $prefix).1.fastq ${dest}/$(basename $prefix).2.fastq 

rm $tmpfile*

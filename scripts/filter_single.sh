#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -p short,defq
#SBATCH -N 1

module load cbibio
module load seqtk
module load fastx
module load bowtie2

fastq=$(ls $prefix* | sed -n '1p')

echo "Fastq: $fastq"

dest=/home/bendall/fs/Projects/HERV_pilot/data/cleaned
mkdir -p $dest

tmpfile=$(mktemp)
time cbi_fastq_filter -v $fastq > $tmpfile.fastq

ref=/lustre/groups/cbi/shared/References/HIV-1/NC_001802.1/Sequence/Bowtie2Index/genome
bowtie2 \
  -x $ref \
  -U $tmpfile.fastq \
  --no-head --reorder > $tmpfile.sam

cat $tmpfile.sam |  awk '$2 == "4" { print $0 }' | \
  perl -lane 'print "\@$F[0]\n$F[9]\n+\n$F[10]"' > ${dest}/$(basename $prefix).fastq

rm $tmpfile*

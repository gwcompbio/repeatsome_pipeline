#! /bin/bash
#SBATCH -t 4:00:00
#SBATCH -p short,defq
#SBATCH -N 1

module load bowtie2
module load cbiC1

##########################################################################################
# Required parameters
#   fastx1 - File with #1 mates
#   fastx1 - File with #2 mates
#   index  - Path to bowtie2 index of sequences to remove
##########################################################################################
[[ -z "$fastx1" ]] && echo "Reads file is not specified" && exit 1
[[ -z "$fastx2" ]] && echo "Pairs file is not specified" && exit 1
[[ -z "$index" ]] && echo "Index is not specifed" && exit 1

out1="${fastx1%.*}.filtered.${fastx1##*.}"
out2="${fastx2%.*}.filtered.${fastx2##*.}"

[[ $(head -c1 $fastx1) == '@' ]] && ftype="q" || ftype="f"

echo "Reads1:    $fastx1"
echo "Reads2:    $fastx2"
echo "Out1:      $out1"
echo "Out2:      $out2"
echo "File type: $ftype"

tmpfile=$(mktemp)

bowtie2 \
  -p $(nproc) \
  -${ftype} \
  -x $index \
  -1 $fastx1 \
  -2 $fastx2 \
  --no-head --reorder > $tmpfile.sam

if  [[ "$ftype" == "q" ]]; then
  cat $tmpfile.sam | grep 'YT:Z:UP' | perl -lane 'print "\@$F[0]\n$F[9]\n+\n$F[10]"' | \
    deinterleave_fastq  $out1 $out2
else
  cat $tmpfile.sam | grep 'YT:Z:UP' | perl -lane 'print ">$F[0]\n$F[9]"' | \
    deinterleave_fasta $out1 $out2
fi

rm $tmpfile*

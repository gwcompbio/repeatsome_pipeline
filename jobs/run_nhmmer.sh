#! /bin/bash
##########################################################################################
# Submission script for nhmmer
#   Required variables:
#     $hmmdb  - the name of the HMM database, in binary format (see hmmpress) 
#   Required variables for single-file processing:
#     $fastx  - the sequence file to be processed
#   Required variables for array processing:
#     $chunks - a file with filenames to be processed, one per line
#     (The scheduler must set $SLURM_ARRAY_TASK_ID)
##########################################################################################

echo "[---run_nhmmer.sh---] ($(date)) Starting run_nhmmer.sh."
umask 0002
module load hmmer/3.1b1

if [[ -z "$fastx" ]]; then
  #--- Check whether this is an array job
  [[ -z "$SLURM_ARRAY_TASK_ID" ]] && echo "[---run_nhmmer.sh---] ($(date)) ERROR: Not an array job." && exit 1
  #--- Check whether file list exists
  [[ -z "$chunks" ]] && echo "[---run_nhmmer.sh---] ($(date)) ERROR: no chunk file specified." && exit 1
  [[ ! -e "$chunks" ]] && echo "[---run_nhmmer.sh---] ($(date)) ERROR: chunk file \"$chunks\" not found." && exit 1
  #--- Set file to be processed
  fastx=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chunks)
fi

#--- Check whether sequence file exists
[[ ! -e "$fastx" ]] && "[---run_nhmmer.sh---] ($(date)) ERROR: sequence file \"$fastx\" does not exist." && exit 1

#--- Check whether input file is fasta or fastq
[[ $(head -c1 $fastx) == "@" ]] && seqtype="fastq" || seqtype="fasta"

#--- Count number of lines in fastx
nlines=$(wc -l < $fastx)

#--- nhmmer does not use all 16 cores, and actually works better with 8 cores
#--- Best performance is when two instances are running on a node
ncpu=$(( $(nproc) / 2 ))

tmp1=$(mktemp)
tmp2=$(mktemp)

#--- Print information
{
echo "[---run_nhmmer.sh---] ($(date)) Sequence file:     $fastx"
echo "[---run_nhmmer.sh---] ($(date)) Sequence format:   $seqtype"
echo "[---run_nhmmer.sh---] ($(date)) Lines in file:     $nlines"
echo "[---run_nhmmer.sh---] ($(date)) CPUs per instance: $ncpu"
echo "[---run_nhmmer.sh---] ($(date)) Temp output 1:     $tmp1"
echo "[---run_nhmmer.sh---] ($(date)) Temp output 2:     $tmp2"
echo "[---run_nhmmer.sh---] ($(date)) Final output:      ${fastx}.hmm_hits.out"
}

#--- Start time
t1=$(date +"%s")

#--- Run nhmmer
if [[ "$seqtype" == "fasta" ]]; then
  head -n $(( $nlines / 2 )) $fastx | nhmmscan --dfamtblout $tmp1 --notextw --cut_ga --cpu $ncpu $hmmdb - > /dev/null &
  tail -n+$(( $(( $nlines / 2 )) + 1 )) $fastx | nhmmscan --dfamtblout $tmp2 --notextw --cut_ga --cpu $ncpu $hmmdb - > /dev/null &
else
  module load cbiC1
  head -n $(( $nlines / 4 )) $fastx | inlineFastq2Fasta | nhmmscan --dfamtblout $tmp1 --notextw --cut_ga --cpu $ncpu $hmmdb - > /dev/null &
  tail -n+$(( $(( $nlines / 4 )) + 1 )) $fastx | inlineFastq2Fasta | nhmmscan --dfamtblout $tmp2 --notextw --cut_ga --cpu $ncpu $hmmdb - > /dev/null &
fi

#--- Wait for jobs to finish
wait
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---run_nhmmer.sh---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."


#--- Combine files and cleanup
cat $tmp1 $tmp2 >> ${fastx}.hmm_hits.out
rm $tmp1
rm $tmp2
echo "[---run_nhmmer.sh---] ($(date)) Complete." && exit 0

#! /bin/bash
#SBATCH -t 48:00:00
#SBATCH -p short,defq
#SBATCH -N 1

module load hmmer/3.1b1

# Set HMM database
[[ -z "$hmmdb" ]] && hmmdb=/lustre/groups/cbi/Repeatsome/data/Dfam_db/Dfam.hmm
[[ ! -e "$hmmdb.h3m" ]] && echo "[ ERROR ] HMM Database \"$hmmdb\" does not exist" && exit 1

# Set reads file
[[ -z "$reads" ]] && echo "[ ERROR ] Reads file is not specified" && exit 1
[[ ! -e "$reads" ]] && echo "[ ERROR ] Reads file \"$reads\" does not exist" && exit 1

# Set output file
[[ -z "$outfile" ]] && outfile=${reads%.*}.hmm_hits.txt

echo "HMM Database: $hmmdb"
echo "Reads:        $reads"
echo "Outfile:      $outfile"

cmd="nhmmscan --tblout ${outfile} --notextw --cut_ga --cpu $(nproc) ${hmmdb} ${reads}"
echo "Command:"
echo $cmd

time $cmd

exit 0

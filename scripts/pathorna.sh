#!/bin/bash

SN="pathorna.sh"

#--- Read command line args, if present
[[ -n "$1" ]] && samfile="$1"
[[ -n "$2" ]] && gtffile="$2"
[[ -n "$3" ]] && outdir="$3" || outdir="$PWD"
[[ -n "$4" ]] && outpre="$4" || outpre="psout"

#--- Check that files exist
[[ ! -e "$samfile" ]] && echo "[---$SN---] ($(date)) FAILED: file $samfile does not exist" && exit 1
[[ ! -e "$gtffile" ]] && echo "[---$SN---] ($(date)) FAILED: file $gtffile does not exist" && exit 1

#--- Print parameters
echo "[---$SN---] ($(date)) Starting $SN"
echo "[---$SN---] ($(date)) Samfile:    $samfile"
echo "[---$SN---] ($(date)) GTF file:   $gtffile"
echo "[---$SN---] ($(date)) Outdir:     $outdir"
echo "[---$SN---] ($(date)) Out prefix: $outpre"

#--- Load modules
module load pysam/0.8.2.1
module load pathorna

#--- Load python virtualenv
#export WORKON_HOME=$HOME/.virtualenvs
#export PROJECT_HOME=$HOME/Devel
#export VIRTUALENVWRAPPER_SCRIPT=$HOME/.local/bin/virtualenvwrapper.sh
#source $HOME/.local/bin/virtualenvwrapper_lazy.sh

#workon pathoRNA

#--- Start the timer
t1=$(date +"%s")

#--- Create directory
mkdir -p $outdir

#--- Pathoscope
echo "[---$SN---] ($(date)) Starting pathoscope_rna."
pathoscope_rna.py --verbose --outdir $outdir --maxIter 1000 --thetaPrior 200000 --exp_tag $outpre $samfile $gtffile
echo "[---$SN---] ($(date)) pathoscope_rna completed successfully."

#--- Sorting and compression
echo "[---$SN---] ($(date)) Sorting and compressing alignment."
module load samtools
sed 's/ZQ:d/ZQ:i/g' $outdir/$outpre-updated.sam | samtools view -uS - | samtools sort - $outdir/$outpre-updated
samtools index $outdir/$outpre-updated.bam

echo "[---$SN---] ($(date)) samtools completed successfully."

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

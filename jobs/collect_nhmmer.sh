#! /bin/bash
##########################################################################################
# Collecting results for nhmmer
#   Required variables:
#     $chunks  - a file with filenames to be processed, one per line
#     $outfile - filename to write results
##########################################################################################

echo "[---collect_nhmmer.sh---] ($(date)) Starting collect_nhmmer.sh."
umask 0002

dodelete=false

#--- Check whether file list exists
[[ -z "$chunks" ]] && echo "[---collect_nhmmer.sh---] ($(date)) ERROR: no chunk file specified." && exit 1
[[ ! -e "$chunks" ]] && echo "[---collect_nhmmer.sh---] ($(date)) ERROR: chunk file \"$chunks\" not found." && exit 1

#--- Check whether output file is specified
[[ -z "$outfile" ]] && echo "[---collect_nhmmer.sh---] ($(date)) ERROR: Output not specified." && exit 1

#--- Move existing outfile to backup
[[ -e "$outfile" ]] && mv ${outfile} ${outfile}.bak

#--- Print information
echo "[---collect_nhmmer.sh---] ($(date)) Chunks file: $chunks"
echo "[---collect_nhmmer.sh---] ($(date)) Output file: $outfile"


#--- This function checks whether all the output files are present
function alldone() 
{
  local _chunks=$1
  local _allfound=true;
  while read fastx; do
    [[ ! -e "${fastx}.hmm_hits.out" ]] && _allfound=false && break
  done < $_chunks
  echo $_allfound
}

#--- Stay in this loop until the output files are present
while [[ $(alldone $chunks) = false ]]; do
  echo "[---collect_nhmmer.sh---] ($(date)) Waiting for jobs to complete..."
  sleep 10
done

#--- Start time
echo "[---collect_nhmmer.sh---] ($(date)) Jobs complete. Concatenating hits files"
t1=$(date +"%s")

#--- Concatenate hits files
allfound=true
while read fastx; do
  hitfile="${fastx}.hmm_hits.out" 
  [[ ! -e "$hitfile" ]] && echo "[---collect_nhmmer.sh---] ($(date)) WARNING: \"$hitfile\" not found." && \
    allfound=false && continue
  echo "[---collect_nhmmer.sh---] ($(date)) Found output \"$hitfile\"."
  grep -v '^#' $hitfile >> $outfile
  #[[ "$first" = true ]] && cat $hitfile > $outfile || tail -n+3 $hitfile >> $outfile
  #first=false
done < $chunks

#--- Do not delete anything if all were not found
[[ "$allfound" = false ]] && echo "[---collect_nhmmer.sh---] ($(date)) WARNING: Not all output files were found. Results may be incomplete."

#--- Clean up temporary files if all were found
if [[ "$allfound" = true ]]; then
  echo "[---collect_nhmmer.sh---] ($(date)) All output files were found. Cleaning up."
  ###while read fastx; do rm -f ${fastx}*; done < $chunks

  #--- Remove temporary directory if empty
  [[ $(ls $(dirname $chunks) | wc -l) == "1" ]] && [[ -e $chunks ]] && \
    echo "[---collect_nhmmer.sh---] ($(date)) Removing temporary directory \"$(dirname $chunks)\"." && \
    rm -rf $(dirname $chunks)
fi

t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---collect_nhmmer.sh---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---collect_nhmmer.sh---] ($(date)) Complete." && exit 0


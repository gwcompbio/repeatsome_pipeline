#!/bin/bash
### From HMMER 1=fasta file in 2=results file
proc=$(nproc)
nhmmscan --tblout $2 --notextw -E 0.0005 -cut_ga -cpu $proc $1










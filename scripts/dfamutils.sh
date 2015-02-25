#! /bin/bash
###############################################################################
# Functions for converting between Dfam format table and other table format
# Example usage:
#   tbl2dfam < hits.tblout > hits.dfamtblout
#   dfam2tbl < hits.dfamtblout > hits.tblout
# 
# Command line arguments to nhmmscan for producing tables:
#   --tblout <f>       : save parseable table of per-sequence hits to file <s>
#   --dfamtblout <f>   : save table of hits to file, in Dfam format <s>
###############################################################################

function tbl2dfam() {
  grep -v '^#' | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[13]\t$F[12]\t$F[14]\t$F[4]\t$F[5]\t$F[11]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[10]"'
}

function dfam2tbl() {
  grep -v '^#' | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t-\t$F[6]\t$F[7]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[13]\t$F[8]\t$F[4]\t$F[3]\t$F[5]"'
}

#! /usr/bin/env python

import utils
from utils.scatterReads import identify_sequence_file, split_reads, submit_scatter, submit_collect

#--- Constants
# Default path to HMM database
HMM_DATABASE = "/lustre/groups/cbi/Repeatsome/data/Dfam/Dfam_db/Dfam_1.2.hmm"
# Default number of minutes for each job.
WALLTIME = 2880

""" Timing information for nhmmer
    !!! THIS INFORMATION ONLY APPLIES TO DFAM 1.2 !!!
    It takes approximately 120 minutes to process 200K reads
    nhmmer will process approximately 100K reads per hour using two instances per node (8 cpus each)
"""
# Default number of sequences per job per hour
SEQS_PER_HOUR = 100000

# Maximum walltime. In this case, max walltime to be in short queue
MAX_WALLTIME = 2880

def main(parser):
  import sys
  import os
  
  args = parser.parse_args()
  seqfile,seqtype = identify_sequence_file(args)
  outfile = args.outfile if args.outfile is not None else '%s.hmm_hits.out' % seqfile
  
  if not os.path.exists(seqfile): sys.exit('ERROR: Sequence file does not exist')

  # Set other command line args
  hmmdb = args.hmmdb
  nosubmit = args.nosubmit
  seqs_per_hour = args.seqs_per_hour 

  # Calculate number of sequences per chunk
  if args.seqs_per_chunk is None:
    minutes = args.walltime
    seqs_per_chunk = int((float(minutes) / 60) * seqs_per_hour)
  else:
    seqs_per_chunk = args.seqs_per_chunk
    minutes = int((float(seqs_per_chunk) / seqs_per_hour) * 60)

  # Print job parameters
  print >>sys.stderr, '%s%s' % ('Sequence file:'.ljust(35), seqfile)
  print >>sys.stderr, '%s%s' % ('Sequence type:'.ljust(35), seqtype)
  print >>sys.stderr, '%s%s' % ('Output file:'.ljust(35), outfile)  
  print >>sys.stderr, '%s%s' % ('Sequences (per job):'.ljust(35), seqs_per_chunk)
  print >>sys.stderr, '%s%s' % ('Walltime:'.ljust(35), minutes)
  print >>sys.stderr, '%s%s' % ('HMM Database:'.ljust(35), hmmdb)
  print >>sys.stderr, '%s%s' % ('Submit jobs?'.ljust(35), nosubmit==False)
  
  """ Process reads and submit """
  # Split the read file
  tmpdir,chunks = split_reads(seqfile, seqtype, seqs_per_chunk)
  print >>sys.stderr, '%s%s' % ('Path to sequence files:'.ljust(35), chunks[0])
  print >>sys.stderr, '\n'.join(' '*35 + _ for _ in chunks[1:])
  
  # Write the chunk file
  chunk_fofn = '%s/chunk.fofn' % tmpdir
  with open(chunk_fofn,'w') as outh:
    print >>outh, '\n'.join(chunks)
  print >>sys.stderr, '%s%s' % ('File of file names:'.ljust(35), chunk_fofn)

  # Submit the scatter jobs
  export_vars = {'chunks':chunk_fofn,'hmmdb':hmmdb,}
  scatter_jobnum = submit_scatter(utils.JOBS['run_nhmmer'], len(chunks), MAX_WALLTIME, export_vars, nosubmit)
  
  # Submit the collect job
  export_vars = {'chunks':chunk_fofn,'outfile':outfile,}  
  collect_jobnum = submit_collect(utils.JOBS['collect_nhmmer'], minutes, export_vars, scatter_jobnum, nosubmit)
  
  if nosubmit:
    import os
    print >>sys.stderr, '%s%s' % ('Detected --nosubmit, deleting temporary directory: '.ljust(35), tmpdir)    
    for f in chunks: os.unlink(f)
    os.unlink('%s/chunk.fofn' % tmpdir)
    os.rmdir(tmpdir)


if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Process reads in parallel',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--fastq', help="Fastq file of reads.")
  parser.add_argument('--fasta', help="Fasta file of reads.")

  parser.add_argument('--outfile', help="Path for final output")
  
  parser.add_argument('--seqs_per_chunk', type=int, help="Number of sequences per job.")
  parser.add_argument('--seqs_per_hour', type=int, help="Number of sequences processed per job per hour.", default=SEQS_PER_HOUR)  
  parser.add_argument('--walltime', type=int, help="Number of minutes for jobs.", default=WALLTIME)
  parser.add_argument('--hmmdb', help="path to HMM database", default=HMM_DATABASE)
  parser.add_argument('--nosubmit', action='store_true')

  main(parser)

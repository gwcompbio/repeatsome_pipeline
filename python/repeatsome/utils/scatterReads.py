#! /usr/bin/env python
import sys

VERBOSE = True

def submit_scatter(script,njobs,minutes,vardict,debug=False):
  ''' Submit array job for reads
      Parameters:
        script  - Path to script file to be submitted
        njobs   - Number of jobs to be submitted
        minutes - Number of minutes to request        
        vardict - Additional parameters to be passed to sbatch (using --export)
  '''
  import sys
  import re  
  import subprocess
  
  _args   = {'script':script,
             'njobs':njobs,
             'minutes':minutes,             
             'export':','.join('%s=%s' % t for t in vardict.iteritems())  
            }
  # Choose queue based on walltime requested
  if _args['minutes'] <= 2880:
    _args['queue'] = 'short,defq'
  else:
    _args['queue'] = 'defq'
  
  _cmd = 'sbatch -N 1 -a 1-%(njobs)d -p %(queue)s -t %(minutes)d --export %(export)s %(script)s' % _args
  print >>sys.stderr, 'Command for scatter:\n\t%s' % _cmd
  if debug:
    return "DEBUG"
  else:
    _pout = subprocess.check_output(_cmd,shell=True)
    print >>sys.stderr, _pout.strip()
    return re.match('Submitted batch job (\d+)',_pout).group(1)

def submit_collect(script,minutes,vardict,depend_jobnum,debug=False):
  ''' Submit job for collecting 
      Parameters:
        script  - Path to script file to be submitted
        njobs   - Number of jobs to be submitted
        minutes - Number of minutes to request        
        vardict - Additional parameters to be passed to sbatch (using --export)
  '''
  import sys
  import re
  import subprocess  
  
  _args   = {'script':script,
             'minutes':minutes,
             'export':','.join('%s=%s' % t for t in vardict.iteritems()),
             'depend_jobnum': depend_jobnum,
            }
  # Choose queue based on walltime requested
  if _args['minutes'] <= 2880:
    _args['queue'] = 'short,defq'
  else:
    _args['queue'] = 'defq'
              
  _cmd = 'sbatch -N 1 -p %(queue)s -t %(minutes)d --export %(export)s --depend=afterany:%(depend_jobnum)s %(script)s' % _args
  print >>sys.stderr, 'Command for collect:\n\t%s' % _cmd
  if debug:
    return "DEBUG"
  else:
    _pout = subprocess.check_output(_cmd,shell=True)
    print >>sys.stderr, _pout.strip()
    return re.match('Submitted batch job (\d+)',_pout).group(1)

def split_reads(seqfile,seqtype,chunksize):
  ''' Split read file into chunks
      Parameters:
        seqfile   - Name of sequence file
        seqtype   - Type of sequences (fastq or fasta)
        chunksize - Number of sequences per chunk
      Return value:
        _chunks  - List of paths to read files
  '''
  import sys
  import subprocess
  from glob import glob
  
  #--- Create temporary directory for split files
  _pout = subprocess.check_output("mktemp -d",shell=True)
  _tmpdir = _pout.strip('\n')
  print >>sys.stderr, '%s%s' % ('Created temporary directory:'.ljust(35), _tmpdir)

  #--- Split the file
  _nlines = int(chunksize) * 4 if seqtype == 'fastq' else int(chunksize) * 2
  retval = subprocess.call('split -l %d %s %s/seq_' % (_nlines, seqfile, _tmpdir),shell=True)
  
  #--- Get paths for chunk file
  _chunks = sorted(glob('%s/seq_*' % _tmpdir))
  return _tmpdir,_chunks

def count_reads(seqfile,seqtype):
  ''' Counts number of reads in file '''
  import subprocess
  
  _lines_per_read = 4 if seqtype == 'fastq' else 2
  _pout = subprocess.check_output("echo $(( $(wc -l < %s) / %d ))" % (seqfile,_lines_per_read),shell=True)
  return int(pout.strip())

def parse_other_args(arglist):
  if not arglist: return {}
  _argdict = {}
  if len(arglist) % 2 != 0:
    print >>sys.stderr, "[ ERROR ] Uneven number of arguments:\n%s" % arglist
    sys.exit()    
  for i in range(0,len(arglist),2):
    if not arglist[i].startswith('--'):
      print >>sys.stderr, "[ ERROR ] Flag has incorrect format: %s" % arglist[i]
      sys.exit()  
    _argdict[arglist[i].strip('--')] = arglist[i+1]
  return _argdict

def identify_sequence_file(_args):
  """ Identify sequence file and type """
  import sys
  
  if _args.fastq is None and _args.fasta is None:
    sys.exit("[ ERROR ] No sequence file was specified.")
  if _args.fastq is not None and _args.fasta is not None:
    sys.exit("[ ERROR ] Only one sequence file may be specified.")
  
  seqfile = _args.fastq if _args.fastq is not None else _args.fasta
  seqtype = 'fastq' if _args.fastq is not None else 'fasta'
  
  return seqfile,seqtype

def main(parser):
  args,oargs = parser.parse_known_args()
  export_vars = parse_other_args(oargs)
  
  seqfile,seqtype = identify_sequence_file(args)

  outfile = args.outfile if args.outfile is not None else '%s.scatter.out' % seqfile
  
  # Other variables
  seqs_per_chunk = args.seqs_per_chunk
  
  tmpdir,chunks = split_reads(seqfile, seqtype, seqs_per_chunk)
  if VERBOSE:
    print >>sys.stderr, '%s%s' % ('Path to sequence files:'.ljust(35), chunks[0])
    print >>sys.stderr, '\n'.join(' '*35 + _ for _ in chunks[1:])

  #--- Make fofn for tracking files
  chunk_fofn = args.chunk_fofn if args.chunk_fofn is not None else '%s/chunk.fofn' % tmpdir 
  with open(chunk_fofn,'w') as outh:
    print >>outh, '\n'.join(chunks)
  if VERBOSE:
    print >>sys.stderr, '%s%s' % ('File of file names:'.ljust(35), chunk_fofn)

  if args.scatter is not None:
    #--- Set the "chunks" variable
    export_vars['chunks'] = chunk_fofn
    scatter_jobnum = submit_scatter(args.scatter, len(chunks), args.minutes, export_vars)
    if args.collect is not None:
      export_vars['outfile'] = outfile
      collect_jobnum = submit_collect(args.collect, args.minutes, export_vars, scatter_jobnum)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Process reads in parallel')

  parser.add_argument('--fastq', help="Fastq file of reads.")
  parser.add_argument('--fasta', help="Fasta file of reads.")

  parser.add_argument('--outfile', help="Path for final output")
  parser.add_argument('--chunk_fofn', help="Path to save chunks.fofn")  
  parser.add_argument('--scatter', help="Path to script for processing chunks of reads.")
  parser.add_argument('--collect', help="Path to script for collecting outputs.")
  
  parser.add_argument('--minutes', type=int, default=120)
  parser.add_argument('--seqs_per_chunk', type=int, default=200000)

  #parser.add_argument('--nchunks', type=int, help="Number of chunks to split")
  #parser.add_argument('--nreads', type=int, help="Number of reads in file. Will be calculated if not specified")

  
  main(parser)

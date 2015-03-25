#! /usr/bin/env python

import sys

def identify_sequence_file(_args):
  """ Identify sequence file and type """
  if _args.fastq is None and _args.fasta is None:
    sys.exit("[ ERROR ] No sequence file was specified.")
  if _args.fastq is not None and _args.fasta is not None:
    sys.exit("[ ERROR ] Only one sequence file may be specified.")
  
  if _args.fastq is not None:
    seqfile = _args.fastq
    seqtype = 'fastq'
  elif _args.fasta is not None:
    seqfile = _args.fasta
    seqtype = 'fasta'
  else:
    sys.exit("[ ERROR ] No sequence file was specified.")

  return seqfile, seqtype  

def get_read_assignments(hitfile, selected_models=None, threshold=None):
  ''' Assign each read to one or more models according to hits file
      Parameters:
        hitfile         - File handle or file name of HMM hits (output by nhmmer)
        selected_models - List of models to allow in assignments. Default (None) is to use
                          all models
        threshold       - E-value must be less than or equal to this threshold. Default
                          (None) is to accept any e-value contained in file.
      Return Values:
        _assignments    - Dictionary with read names as keys and list of models as value.
        _models         - Set of models matched with any read
  '''
  from collections import defaultdict
  
  fh = open(hitfile,'rU') if type(hitfile) is str else hitfile
  _modelset = set(selected_models) if selected_models is not None else None
  _models = set()
  _assignments = defaultdict(set)
  evalue_idx = None
  lines = (l.strip('\n').split() for l in fh if not l.startswith('#'))
  for l in lines:
    if _modelset is not None and l[1] not in _modelset:  # Skip hit if model does not match
      continue
    if threshold is not None:
      if evalue_idx is None:
        # Find the evalue column by testing column 3. If the value is "-" then
        # this is tblout format and evalue is in column 12. Otherwise it is 
        # dfamtblout and evalue is in column 4.
        evalue_idx = 12 if l[3] == "-" else 4
      if float(l[evalue_idx]) > threshold: # Skip hit if e-value exceeds threshold
        continue
    # Model accession is in column 1
    # Read name is in column 2
    _models.add(l[1])
    _assignments[l[2]].add(l[1])
  
  return _assignments,_models

def fastq_iterator(file):
  ''' Iterates over fastq file (4 lines at a time) '''
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next(),lines.next(),lines.next()]

def fasta_iterator(file):
  ''' Iterates over fasta file (2 lines at a time) '''
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next()]

def main(parser):
  args = parser.parse_args()

  #--- Create iterator over sequence file
  seqfile, seqtype = identify_sequence_file(args)
  siter = fastq_iterator(seqfile) if seqtype == 'fastq' else fasta_iterator(seqfile)
  
  #--- Get a list of accepted models
  if args.modellist is not None:
    modellist = set([l.strip('\n').split('\t')[0] for l in open(args.modellist,'rU')])
    print >>sys.stderr, '%s%s' % ('Number of models:'.ljust(35), len(modellist))
  else:
    modellist = None

  #--- Parse hits files
  print >>sys.stderr, '%s%s' % ('Parsing hits file:'.ljust(35), args.hitfile)
  assignments,models = get_read_assignments(args.hitfile, selected_models=modellist, threshold=args.evalue)
  print >>sys.stderr, '%s%s reads after filter.' % ('Hitfile loaded:'.ljust(35), len(assignments))
  
  nreads = npass = 0
  seqprefix = args.output
  
  outh = open('%s.match_hmm.%s' % (seqprefix,seqtype),'w')
  for s1 in siter:
    nreads += 1
    readname = s1[0][1:]
    r1hit = readname in assignments
    if r1hit:
      npass += 1
      print >>outh, '\n'.join(s1)
  
  outh.close()

  print >>sys.stderr, '%s%s' % ('Reads with hits:'.ljust(35), npass)
  print >>sys.stderr, '%s%s' % ('Total reads:'.ljust(35), nreads)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Filter reads that match given model or models')

  parser.add_argument('--fastq', help="Fastq file of reads, single end")
  parser.add_argument('--fasta', help="Fasta file of reads, single end")
  parser.add_argument('--hitfile', help="HMMER output file for reads")
  
  parser.add_argument('--output', default="filteredreads", help="Output prefix.")# If --bymodel, a directory is created, otherwise writes to file")
  parser.add_argument('--evalue', type=float, help="Only accept hits <= this E-value threshold")
  
  parser.add_argument('--modellist', help="File with list of models to be included, one per line. The line may be tab delimited, and the model number is assumed to be in the first column.")
  main(parser)

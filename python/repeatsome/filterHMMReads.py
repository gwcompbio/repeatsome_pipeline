#! /usr/bin/env python

import sys
from utils.scatterReads import identify_sequence_file

def get_read_assignments(fh,threshold=None):
  ''' Assign each read to one or more models according to hits file
        _assignments - Dictionary with read names as keys and corresponding model list as value.
        _models      - Set of models matched with any read
  '''
  from collections import defaultdict
  _models = set()
  _assignments = defaultdict(list)
  lines = (l.strip('\n').split() for l in fh if not l.startswith('#'))
  for l in lines:
    read_name = l[2]
    target_accn = l[1]
    evalue = float(l[12])  
    if threshold is None or evalue <= threshold:
      _models.add(target_accn)
      _assignments[read_name].append(target_accn)
  return _assignments,_models

def open_output_files(outdir,models,seqtype):
  import os
  if not os.path.isdir(outdir): os.makedirs(outdir)
  # Create dictionary mapping model accessions to file handles
  fh_dict = {}
  for model_accn in models:
    fh_dict[model_accn] = open('%s/%s.%s' % (outdir,model_accn,seqtype),'w')
  return fh_dict

def close_output_files(fh_dict):
  # Close all open file handles  
  for k,fh in fh_dict.iteritems():
    fh.close()

def fastq_iterator(file):
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next(),lines.next(),lines.next()]

def fasta_iterator(file):
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next()]

def main(parser):
  args = parser.parse_args()

  #--- Create iterator over sequence file
  seqfile,seqtype = identify_sequence_file(args)
  seqiter = fastq_iterator(seqfile) if seqtype == 'fastq' else fasta_iterator(seqfile)
  
  #--- Get a list of accepted models
  if args.modellist is not None:
    accepted_models = set([l.strip('\n').split('\t')[0] for l in open(args.modellist,'rU')])
  else:
    accepted_models = None

  #--- Parse hits file
  assignments,models = get_read_assignments(args.hitfile, args.evalue)
  #assignments,models = get_read_assignments(open('seq2M.fastq.hmm_hits.out','rU'),None)
  
  nreads = npass = 0
  seqprefix = '.'.join(seqfile.split('.')[:-1])
  if args.bymodel:
    #--- Write seperate files for each model
    outdir = '%s.match_hmm.out' % seqprefix if args.output is None else args.output
    fh_dict = open_output_files(outdir,list(models),seqtype)
    for s in seqiter:
      nreads += 1
      if s[0][1:] in assignments:
        hasone = False
        for hit_model in assignments[s[0][1:]]:
          if accepted_models is None or hit_model in accepted_models:
            print >>fh_dict[hit_model], '\n'.join(s)
            hasone = True
        if hasone: npass += 1
    close_output_files(fh_dict)
    
  else:
    #--- Write one file with all passing reads
    outfile = '%s.match_hmm.%s' % (seqprefix,seqtype) if args.output is None else args.output
    with open(outfile,'w') as outh:
      for s in seqiter:
        nreads += 1
        if s[0][1:] in assignments:
          if accepted_models is None or any([m in accepted_models for m in assignments[s[0][1:]]]):
            npass += 1
            print >>outh, '\n'.join(s)

  print >>sys.stderr, '%s%s' % ('Total reads:'.ljust(35), nreads)
  print >>sys.stderr, '%s%s' % ('Reads with hit(s):'.ljust(35), npass)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Filter reads that match given model or models')

  parser.add_argument('--fastq', help="Fastq file of reads")
  parser.add_argument('--fasta', help="Fasta file of reads")
  parser.add_argument('--output', help="Place for output. If --bymodel, a directory is created, otherwise writes to file")
  parser.add_argument('--evalue', type=float, help="Only accept hits <= this E-value threshold")
  
  parser.add_argument('--modellist', help="File with list of models to be included, one per line. The line may be tab delimited, and the model number is assumed to be in the first column.")
  parser.add_argument('--bymodel', action='store_true' ,help="Write reads matching each model into a seperate file")
  
  parser.add_argument('hitfile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
  
  main(parser)

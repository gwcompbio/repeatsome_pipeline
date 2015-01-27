#! /usr/bin/env python

import sys

def identify_sequence_files(_args):
  """ Identify sequence file and type """
  import sys
  
  if _args.fastq1 is None and _args.fasta1 is None:
    sys.exit("[ ERROR ] No sequence file was specified.")
  if _args.fastq1 is not None and _args.fasta1 is not None:
    sys.exit("[ ERROR ] Only one sequence file may be specified.")
  
  if _args.fastq1 is not None:
    if _args.fastq2 is None:
      sys.exit("[ ERROR ] Mate pairs not specified.")
    else:
      seqfile1 = _args.fastq1
      seqfile2 = _args.fastq2
      seqtype  = 'fastq'
  elif _args.fasta1 is not None:
    if _args.fasta2 is None:
      sys.exit("[ ERROR ] Mate pairs not specified.")
    else:
      seqfile1 = _args.fasta1
      seqfile2 = _args.fasta2
      seqtype  = 'fasta'
  
  return seqfile1,seqfile2,seqtype


def get_read_assignments(hitfile, accepted_models=None, threshold=None):
  ''' Assign each read to one or more models according to hits file
        _assignments - Dictionary with read names as keys and corresponding model list as value.
        _models      - Set of models matched with any read
  '''
  from collections import defaultdict
  
  fh = open(hitfile,'rU') if type(hitfile) is str else hitfile
  _modelset = set(accepted_models) if accepted_models is not None else None
  
  _models = set()
  _assignments = defaultdict(set)
  lines = (l.strip('\n').split() for l in fh if not l.startswith('#'))
  for l in lines:
    if _modelset is not None and l[1] not in _modelset:
      continue
    evalue = float(l[12])
    if threshold is not None and evalue > threshold:
      continue
    read_name = l[2]
    target_accn = l[1]
    _models.add(target_accn)
    _assignments[read_name].add(target_accn)
  
  return _assignments,_models

"""
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
"""

def fastq_iterator(file):
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next(),lines.next(),lines.next()]

def paired_fastq_iterator(file1,file2):
  iter1 = fastq_iterator(file1)
  iter2 = fastq_iterator(file2)
  while iter1:
    yield (iter1.next(), iter2.next())

def fasta_iterator(file):
  fh = open(file,'rU') if type(file) is str else file
  lines = (l.strip() for l in fh)
  while lines:
    yield [lines.next(),lines.next()]

def paired_fasta_iterator(file1,file2):
  iter1 = fasta_iterator(file1)
  iter2 = fasta_iterator(file2)
  while iter1:
    yield (iter1.next(), iter2.next())

def main(parser):
  args = parser.parse_args()

  #--- Create iterator over sequence file
  seqfile1, seqfile2, seqtype = identify_sequence_files(args)
  piter = paired_fastq_iterator(seqfile1, seqfile2) if seqtype == 'fastq' else paired_fasta_iterator(seqfile1, seqfile2)
  
  #--- Get a list of accepted models
  if args.modellist is not None:
    modellist = set([l.strip('\n').split('\t')[0] for l in open(args.modellist,'rU')])
  else:
    modellist = None

  #--- Parse hits file
  assignments1,models1 = get_read_assignments(args.hitfile1, accepted_models=modellist, threshold=args.evalue)
  print >>sys.stderr, '%s%s reads pass filters.' % ('Read assignments 1 loaded:'.ljust(35), len(assignments1))  
  #print >>sys.stderr, 'Read assignments 1 loaded: %d' % len(assignments1)
  assignments2,models2 = get_read_assignments(args.hitfile2, accepted_models=modellist, threshold=args.evalue)  
  print >>sys.stderr, '%s%s reads pass filters.' % ('Read assignments 2 loaded:'.ljust(35), len(assignments2))    
  #print >>sys.stderr, 'Read assignments 2 loaded: %d' % len(assignments2) 
  
  nreads = npass = 0
  seqprefix = args.output
  
  outh1 = open('%s.1.match_hmm.%s' % (seqprefix,seqtype),'w')
  outh2 = open('%s.2.match_hmm.%s' % (seqprefix,seqtype),'w')
  for s1,s2 in piter:
    nreads += 1
    readname = s1[0][1:]
    r1hit = readname in assignments1
    r2hit = readname in assignments2
    if args.concordant:
      if r1hit and r2hit:
        npass += 1
        print >>outh1, '\n'.join(s1)
        print >>outh2, '\n'.join(s2)        
    else:
      if r1hit or r2hit:
        npass += 1    
        print >>outh1, '\n'.join(s1)
        print >>outh2, '\n'.join(s2)
  
  outh1.close()
  outh2.close()
  
  print >>sys.stderr, '%s%s' % ('Pairs with hits:'.ljust(35), npass)
  print >>sys.stderr, '%s%s' % ('Total pairs:'.ljust(35), nreads)  

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Filter reads that match given model or models')

  parser.add_argument('--fastq1', help="Fastq file of reads, read 1")
  parser.add_argument('--fastq2', help="Fastq file of reads, read 2")  
  parser.add_argument('--fasta1', help="Fasta file of reads, read 1")
  parser.add_argument('--fasta2', help="Fasta file of reads, read 2")
  parser.add_argument('--hitfile1', help="HMMER output file for reads 1")
  parser.add_argument('--hitfile2', help="HMMER output file for reads 2")
  
  parser.add_argument('--output', default="filteredreads", help="Output prefix. If --bymodel, a directory is created, otherwise writes to file")
  parser.add_argument('--evalue', type=float, help="Only accept hits <= this E-value threshold")
  
  parser.add_argument('--modellist', help="File with list of models to be included, one per line. The line may be tab delimited, and the model number is assumed to be in the first column.")
  parser.add_argument('--concordant', action='store_true', help="Require that both reads in a pair have a hit. Otherwise, read pair will be output if either read has a hit")  
  # parser.add_argument('--bymodel', action='store_true' ,help="Write reads matching each model into a seperate file")
    
  main(parser)

#! /usr/bin/env python
import sys

def get_read_assignments(fh,threshold=None):
  ''' Assign each read to one or more models according to hits file '''
  from collections import defaultdict
  models = set()
  assignments = defaultdict(list)
  lines = (l.strip('\n').split() for l in fh if not l.startswith('#'))
  for l in lines:
    read_name = l[2]
    target_accn = l[1]
    evalue = float(l[12])  
    if threshold is None or evalue <= threshold:
      models.add(target_accn)
      assignments[read_name].append(target_accn)
  return assignments,models

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
  
def write_reads(seqs,assignments,fh_dict,seqtype):
  from Bio import SeqIO
  for s in seqs:
    # Write read to each file
    for hit_model in assignments[s.id]:
      z = SeqIO.write(s,fh_dict[hit_model],seqtype)

def main(parser):
  args = parser.parse_args()
  if args.fastq is None and args.fasta is None:
    print >>sys.stderr, "[ ERROR ] No sequence file was specified"
    parser.print_help()
    sys.exit()
  if args.fastq is not None and args.fasta is not None:
    print >>sys.stderr, "[ ERROR ] Only one sequence file may be specified"
    parser.print_help()
    sys.exit()

  from Bio import SeqIO  
  seqfile = args.fastq if args.fastq is not None else args.fasta
  seqtype = 'fastq' if args.fastq is not None else 'fasta'
  seqs = (s for s in SeqIO.parse(seqfile,seqtype))  
  
  if args.outdir is None:
    outdir = '.'.join(seqfile.split('.')[:-1])
  else:
    outdir = args.outdir

  assignments,models = get_read_assignments(args.hitfile,args.evalue)
  fh_dict = open_output_files(outdir,list(models),seqtype)
  write_reads(seqs,assignments,fh_dict,seqtype)
  close_output_files(fh_dict)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Split reads according to HMM model')

  parser.add_argument('--fastq', help="Fastq file of reads")
  parser.add_argument('--fasta', help="Fasta file of reads")
  parser.add_argument('--outdir', help="Directory for output files")
  parser.add_argument('--evalue', type=float, help="Only accept hits <= this E-value threshold")
  parser.add_argument('hitfile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
  
  main(parser)

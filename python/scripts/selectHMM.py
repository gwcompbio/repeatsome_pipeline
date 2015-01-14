#! /usr/bin/env python
"""
Select specific HMMs from HMM file
"""
import sys

def hmmModelIter(hmm_file):
  ''' Iterate over models in an HMM file
  '''
  import gzip
  if hmm_file.endswith('gz'):
    liter = (l.strip('\n') for l in gzip.open(hmm_file,'rb') if not l.startswith('#'))
  else:
    liter = (l.strip('\n') for l in open(hmm_file,'rU') if not l.startswith('#'))
  cur = []
  while liter:
    l = liter.next()
    if l.startswith('//'):
      yield cur
      cur = []
    else:
      cur.append(l)


def hmmModelData(rec):
  ''' Extract fields from HMM record
      The Type, Class, and Superfamily information is important
  '''
  from collections import defaultdict
  recdict = defaultdict(list)
  for l in rec:
    if not l.startswith(' '):
      recdict[ l[:6].strip() ].append(l[6:].strip())
  
  retdict = {'accn':recdict['ACC'][0],
             'name':recdict['NAME'][0],
             'desc':recdict['DESC'][0],
             'length':recdict['LENG'][0],
             }
  for ctline in recdict['CT']:
    if ctline.startswith('Type'):
      retdict['type'] = ctline.split(';')[1].strip()
    elif ctline.startswith('Class'):
      retdict['class'] = ctline.split(';')[1].strip()
    elif ctline.startswith('Superfamily'):
      retdict['superfamily'] = ctline.split(';')[1].strip().strip('?')
  return retdict


def main(parser):
  args = parser.parse_args()
  
  #import re
  #types         = [re.compile('^%s$' % _) for _ in args.type.split(',')] if args.type else None
  #classes       = [re.compile('^%s$' % _) for _ in args.rep_class.split(',')]  if args.rep_class else None
  #superfamilies = [re.compile('^%s$' % _) for _ in args.superfamily.split(',')] if args.superfamily else None
  #accessions    = [re.compile('^%s$' % _) for _ in args.accessions.split(',')] if args.accessions else None
  
  from fnmatch import fnmatch, fnmatchcase
  descriptions  = args.description.split(',') if args.description else None  
  accessions    = args.accession.split(',') if args.accession else None
  names         = args.name.split(',') if args.name else None  
  superfamilies = args.superfamily.split(',') if args.superfamily else None
  classes       = args.rep_class.split(',')  if args.rep_class else None  
  types         = args.type.split(',') if args.type else None

  print >>sys.stderr, '\t'.join(['#accn','name','superfamily','class','type','length','desc'])  
  
  for record in hmmModelIter(args.hmmfile):
    hd = hmmModelData(record)
    passes = True
    #passes &= accessions is None or any([rx.match(hd['accn']) is not None for rx in accessions])
    #passes &= types is None or any([rx.match(hd['type']) is not None for rx in types])
    #passes &= classes is None or any([rx.match(hd['class']) is not None for rx in classes])
    #passes &= superfamilies is None or any([rx.match(hd['superfamily']) is not None for rx in superfamilies])

    passes &= descriptions is None or any([fnmatch(hd['desc'],pat) for pat in descriptions])    
    passes &= accessions is None or any([fnmatch(hd['accn'],pat) for pat in accessions])    
    passes &= names is None or any([fnmatch(hd['name'],pat) for pat in names])    
    passes &= superfamilies is None or any([fnmatch(hd['superfamily'],pat) for pat in superfamilies])    
    passes &= classes is None or any([fnmatch(hd['class'],pat) for pat in classes])    
    passes &= types is None or any([fnmatch(hd['type'],pat) for pat in types])

    if passes:
      print >>sys.stderr, '\t'.join(hd[f] for f in ['accn','name','superfamily','class','type','length','desc'])
      if args.hmmout is not None:
        print >>args.hmmout, '\n'.join(record)
        print >>args.hmmout, '//'

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description='Create table with data for HMM models')
  parser.add_argument('--hmmout', type=argparse.FileType('w'), help="Write HMM models to this file")  
  parser.add_argument('--description', help="Select models only with this description")  
  parser.add_argument('--accession', help="Select models only with these accessions")
  parser.add_argument('--name', help="Select models only with these names")  
  parser.add_argument('--superfamily', help="Select models only from this superfamily")
  parser.add_argument('--rep_class', help="Select models only from this class")  
  parser.add_argument('--type', help="Select models only from this type")
  parser.add_argument('hmmfile', help="name of HMM file to filter")


  
  main(parser)

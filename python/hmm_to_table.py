#! /usr/bin/env python

"""

"""

def hmmIter(hmm_file):
  ''' Iterate over records in HMM file
  '''
  import gzip
  if hmm_file.endswith('gz'):
    liter = (l.strip('\n') for l in gzip.open(hmm_file,'rb') if not l.startswith('#') and not l.startswith(' '))
  else:
    liter = (l.strip('\n') for l in open(hmm_file,'rU') if not l.startswith('#') and not l.startswith(' '))
  cur = []
  while liter:
    l = liter.next()
    if l.startswith('//'):
      yield cur
      cur = []
    else:
      cur.append(l)


def processRecord(rec):
  ''' Extract fields from HMM record
      The Type, Class, and Superfamily information is important
  '''
  from collections import defaultdict
  recdict = defaultdict(list)
  for l in rec:
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
  records = [processRecord(_) for _ in hmmIter('Dfam/Dfam.hmm.gz')]
  print "%d records processed" % len(records)

  with open('Dfam/Dfam.records.txt','w') as outh:
    for r in records:
      print >>outh, '\t'.join(r[f] for f in ['type','class','superfamily','accn','length','name','desc'])

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='Create table with data for HMM models')

  parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
  parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
  main(parser)

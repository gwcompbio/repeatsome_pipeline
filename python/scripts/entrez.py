#! /usr/bin/env python

from Bio import Entrez
Entrez.email = "bendall@gwu.edu"
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import time
import os
from glob import glob



import argparse
parser = argparse.ArgumentParser(description='Get SRA data for SRP project')
parser.add_argument('sraID')
args = parser.parse_args()

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")

def get_rundata(xmlfile):
  rundata = {}
  tree = ET.parse(xmlfile)
  root = tree.getroot()
  
  samp = root.findall('./EXPERIMENT_PACKAGE/SAMPLE')
  assert len(samp)==1
  samp = samp[0]
  rundata['sample_accession'] = samp.get('accession')
  rundata['sample_alias'] = samp.get('alias')
  rundata['sample_title'] = samp.find('TITLE').text
  for sattr in samp.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
    rundata['sample_%s' % sattr.find("TAG").text] = sattr.find("VALUE").text
  
  run = root.findall('./EXPERIMENT_PACKAGE/RUN_SET/RUN')
  assert len(run)==1
  run = run[0]
  rundata['run_accession'] = run.get('accession')
  rundata['run_alias'] = run.get('alias')
  rundata['run_bases'] = run.find('Bases').get('count')
  rundata['run_readperspot'] = run.find('Statistics').get('nreads')
  rundata['run_nspots'] = run.find('Statistics').get('nspots')
  
  return rundata

def make_url(accn):
  _url = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'
  _url += '%s/%s/%s/%s.sra' % (accn[:3],accn[:6],accn,accn)
  return _url

# Make metadata directory
if not os.path.isdir('metadata'):
  os.mkdir('metadata')

# Get the SRA run IDs list from the SRP record
sraID = args.sraID
record = Entrez.read(Entrez.esearch(db="sra",term=sraID,retmax=1000))

print >>sys.stderr, "Found %d run IDs" % len(record['IdList'])

# Fetch XML data for SRA run IDs
counter = 0
for sampid in record['IdList']:
  print >>sys.stderr, 'Fetching %s' % sampid
  pr = Entrez.efetch(db="sra",id=sampid).read()
  if pr: counter += 1
  root = ET.fromstring(pr)
  run_ids = root.findall('./EXPERIMENT_PACKAGE/RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')
  print >>sys.stderr, 'Runs: %s' % ', '.join(_.text for _ in run_ids)
  with open('metadata/%s.xml' % run_ids[0].text,'w') as outh:  
    print >>outh, prettify(root).encode('ascii','ignore')

print >>sys.stderr, "Download run data complete: %d runs" % counter
print >>sys.stderr, "Creating sample matrix file"

# Parse XML files to get run data for all runs
xmlfiles = sorted(glob('metadata/*.xml'))
allrundata = []
for f in xmlfiles:
  allrundata.append(get_rundata(f))

# Print the run metadata to sample_matrix.txt
columns = sorted(allrundata[0].keys())
with open('metadata/sample_matrix.txt','w') as outh:
  print >>outh, '\t'.join(_.replace(' ','_') for _ in columns)
  for d in allrundata:
    print >>outh, '\t'.join(d[c] if c in d else '.' for c in columns)

print >>sys.stderr, "Sample matrix file created"
print >>sys.stderr, "Creating sample URL file"

# Create file with sample URLs
lines = [l.strip('\n').split('\t') for l in open('metadata/sample_matrix.txt','rU')]
lines = lines[1:]
with open('metadata/sample_urls.txt','w') as outh:
  for l in lines:
    read_accn = l[0]
    print >>outh, make_url(read_accn)

print >>sys.stderr, "Sample URL file created"

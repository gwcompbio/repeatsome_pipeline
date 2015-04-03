#! /bin/bash

#--- Load the repeatsome module
module load repeatsome

#--- Each sample is in its own directory containing one FASTQ file
samp=SRR000001

#--- Set paths to filter database, alignment database, and annotation
filtref=/lustre/groups/cbi/Repeatsome/data/contaminants/db
alignref=/lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
annot=/lustre/groups/cbi/Repeatsome/data/annotations/knownHERV.hg19/knownHERV.hg19.gtf

#--- Run QC
rp_qc_SE ${samp}/${samp}_pass.fastq ${samp}/${samp}_qc.fastq

#--- Run filtering
rp_filter_SE $filtref ${samp}/${samp}_qc.fastq ${samp}/${samp}_qcfilt.fastq

#--- Run alignment
rp_align_SE $alignref ${samp}/${samp}_qcfilt.fastq ${samp}/aligned

#--- Run pathoscope
rp_pathorna ${samp}/aligned.sam $annot ${samp}/PS knownHERV_all

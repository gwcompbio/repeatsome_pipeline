#! /bin/bash

#--- Load the repeatsome module
module load repeatsome

#--- Each sample is in its own directory containing paired FASTQ files
samp=SRR123456

#--- Set paths to filter database, alignment database, and annotation
filtref=/lustre/groups/cbi/Repeatsome/data/contaminants/db
alignref=/lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
annot=/lustre/groups/cbi/Repeatsome/data/annotations/knownHERV.hg19/knownHERV.hg19.gtf

#--- Run QC
rp_qc_PE ${samp}/${samp}_pass_1.fastq ${samp}/${samp}_pass_1.fastq ${samp}/${samp}_qc_1.fastq ${samp}/${samp}_qc_2.fastq

#--- Run filtering
rp_filter_PE $filtref ${samp}/${samp}_qc_1.fastq ${samp}/${samp}_qc_2.fastq ${samp}/${samp}_qcfilt_1.fastq ${samp}/${samp}_qcfilt_2.fastq

#--- Run alignment
rp_align_PE $alignref ${samp}/${samp}_qcfilt_1.fastq ${samp}/${samp}_qcfilt_2.fastq ${samp}/aligned

#--- Run pathoscope
rp_pathorna ${samp}/aligned.sam $annot ${samp}/PS knownHERV_all

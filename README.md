repeatsome_pipeline
===================

### nhmmerWrapper

##### Split sequencing reads into smaller files for running nhmmscan

Running nhmmscan is embarrassingly parallel. This wrapper takes advantage of the compute
resources by splitting the file into many pieces, scattering the computation among many
nodes, and collecting the results into a single file.

    usage: nhmmerWrapper [-h] [--fastq FASTQ] [--fasta FASTA] [--outfile OUTFILE]
                         [--seqs_per_chunk SEQS_PER_CHUNK] [--walltime WALLTIME]
                         [--hmmdb HMMDB] [--nosubmit]
    
    Process reads in parallel
    
    optional arguments:
      -h, --help            show this help message and exit
      --fastq FASTQ         Fastq file of reads.
      --fasta FASTA         Fasta file of reads.
      --outfile OUTFILE     Path for final output
      --seqs_per_chunk SEQS_PER_CHUNK
                            Number of sequences per job.
      --walltime WALLTIME   Number of minutes for jobs.
      --hmmdb HMMDB         path to HMM database
      --nosubmit

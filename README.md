repeatsome_pipeline
===================

### nhmmerWrapper

##### Use nhmmer to search for reads that match HMM models

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

---

### filterHMMReads

##### Filter reads matching one or more HMM models

Input the original reads file along with an HMM hits file to extract reads that match
the HMM database. Output can be a single file with all matching reads, or one file for each
model (using --bymodel). Read matching only certain models can be selected by supplying a
file with the model names (--modellist).


    usage: filterHMMReads [-h] [--fastq FASTQ] [--fasta FASTA] [--output OUTPUT]
                          [--evalue EVALUE] [--modellist MODELLIST] [--bymodel]
                          [hitfile]
    
    Filter reads that match given model or models
    
    positional arguments:
      hitfile
    
    optional arguments:
      -h, --help            show this help message and exit
      --fastq FASTQ         Fastq file of reads
      --fasta FASTA         Fasta file of reads
      --output OUTPUT       Place for output. If --bymodel, a directory is
                            created, otherwise writes to file
      --evalue EVALUE       Only accept hits <= this E-value threshold
      --modellist MODELLIST
                            File with list of models to be included
      --bymodel             Write reads matching each model into a seperate file

---


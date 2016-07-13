# BaitSTR_type.pl
Script for use with the baitSTR pipeline.

Provide "blocks" fasta file from the extend_STR_reads script (BaitSTR), and this companion script uses lobSTR and/or BWA to 
prepare reference sequence files and align one or more sample datasets. Optionally, generate probe sequences.

Usage:
perl BaitSTR_type.pl --stem [stem] [options]

Sample usage:
perl BaitSTR_type.pl --stem my_run --full --index --target Blocks.fasta --index_prefix newIndexFiles \
  --r1 reads.r1.sample1.fastq,sample1 reads.r1.sample2.fastq,sample2
  --r2 reads.r2.sample1.fastq,sample1 reads.r2.sample2.fastq,sample2
  
Options:

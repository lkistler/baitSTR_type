# BaitSTR_type.pl
Script for use with the baitSTR pipeline.

Provide "blocks" fasta file from the extend_STR_reads script (BaitSTR), and this companion script uses lobSTR and/or BWA to 
prepare reference sequence files and align one or more sample datasets. Optionally, generate probe sequences.

Usage:
perl BaitSTR_type.pl --stem [stem] [options]

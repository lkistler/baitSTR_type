# BaitSTR_type.pl
Script for use with the baitSTR pipeline.

Provide "blocks" fasta file from the extend_STR_reads script (BaitSTR), and this companion script uses lobSTR and/or BWA to 
prepare reference sequence files and align one or more sample datasets. Optionally, you can also generate probe sequences.

Usage:
	perl BaitSTR_type.pl --stem [stem] [options]

Sample usage:
	perl BaitSTR_type.pl --stem my_run --full --index --target Blocks.fasta --index_prefix newIndexFiles \
  	--r1 reads.r1.sample1.fastq,sample1 reads.r1.sample2.fastq,sample2 \
  	--r2 reads.r2.sample1.fastq,sample1 reads.r2.sample2.fastq,sample2 \

Options:
	--index			            Create index files (requires --target)
	--align			            Perform read alignemnt (requires --r1 and --r2, and/or --SR)
	--allelotype		        Perform allelotype calling (requires --align or --bams)
	--design_probes		      Produce a fasta file of probes from qualifying input blocks
	--full			            Do alignment, allelotype, and probe design (full pipeline excluding index)
	--mem			              Use bwa-mem for alignemnt (uses lobSTR for allelotyping)
	--lobSTR		            Use lobSTR for alignment (uses lobSTR for allelotyping)
	--backtrack		          Use bwa-backtrack for alignment DEPRECATED AND NOT RECOMMENDED (skips lobSTR allelotyping and generates non-standard genotype file)
	
	--stem [str]		        Prefix for run files (will overwrite without warning)
	--index_prefix [str]	  Prefix for index files
	--target [str]		      Fasta file of extended blocks from BaitSTR output

  	Read files. Sample IDs must be given for each read file, multiple files are allowed per sample ID.
	--r1 [str]		          Forward reads in format "infile.R1.fastq,sampleID [file2,id2 file3,id3...]", requires --r2 with corresponding sample IDs
	--r2 [str]		          Reverse reads in format "infile.R2.fastq,sampleID", requires --r1 with corresponding sample IDs
	--SR [str]		          Unpaired reads in format "infile.fastq,sampleID [file2,id2 file3,id3...]"
	--bams [str]		        Compatible bam files in format "infile.bam,sampleID [file2,id2 file3,id3...]"
	
	Other options.
	--path_to_lobSTR [str]		    Full path to lobSTR binaries
	--path_to_samtools [str]	    Full path to samtools
	--path_to_bwa [str]		        Full path to BWA
	--lob_args [str]		          Additional command line arguments for running lobSTR (alignment), in quotes
	--allel_args [str]		        Additional command line arguments for running allelotype, in quotes
	--mem_args [str]		          Additional command line arguments for running bwa-mem (alignment), in quotes
	--backtrack_args [str]		    Additional command line arguments for running bwa-backtrack (alignment), in quotes
	--noise [str]			            Provide a noise model to use with allelotype (default uses PCR-free noise model)
	--samqual [n]			            Minimum read alignment quality scores ("samtools view -q [n]"), with mem and backtrack
	--gzip				                Read files are compressed using gzip (specify for lobSTR)
	--fasta				                Reads are in fasta format (default fastq; specify for lobSTR)
	--minlen [n]			            Minimum block length to include in analysis (400)
	--flank [n]			              Minimum non-repeat flank length in block (200)
	--alnflank [n]			          Minimum aligned read flank in backtrack (15)
	--mincvg [n]			            Minimum callable-read coverage in backtrack (5)
	--minallele [n]			          Minimum callable-reads per allele (2)
	--no_rmdup			              Do not perform duplicate removal following read alignment
	--help				                Help screen
	--silent			                Run silently

#Details for running BaitSTR_type.pl

Each run 

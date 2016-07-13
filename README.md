# BaitSTR_type.pl
Script for use with the baitSTR pipeline.

Provide "blocks" fasta file from the extend_STR_reads script (BaitSTR), and this companion script uses lobSTR and/or BWA to 
prepare reference sequence files and align one or more sample datasets. Optionally, you can also generate probe sequences.

Usage:

	perl BaitSTR_type.pl --stem [stem] --index_prefix [index_prefix] [options]

Sample usage:

	perl BaitSTR_type.pl --stem my_run --full --index --target Blocks.fasta --index_prefix newIndexFiles \
  	--r1 reads.r1.sample1.fastq,sample1 reads.r1.sample2.fastq,sample2 \
  	--r2 reads.r2.sample1.fastq,sample1 reads.r2.sample2.fastq,sample2 \

Options:

	--index			Create index files (requires --target)
	--align			Perform read alignemnt (requires --r1 and --r2, and/or --SR)
	--allelotype		Perform allelotype calling (requires --align or --bams)
	--design_probes		Produce a fasta file of probes from qualifying input blocks
	--full			Do alignment, allelotype, and probe design (full pipeline excluding index)
	--mem			Use bwa-mem for alignemnt (uses lobSTR for allelotyping)
	--lobSTR		Use lobSTR for alignment (uses lobSTR for allelotyping)
	--backtrack		Use bwa-backtrack for alignment DEPRECATED AND NOT RECOMMENDED (skips lobSTR allelotyping and generates non-standard genotype file)
	
	--stem [str]		Prefix for run files (will overwrite without warning)
	--index_prefix [str]	Prefix for index files
	--target [str]		Fasta file of extended blocks from BaitSTR output

  	Read files. Sample IDs must be given for each read file, multiple files are allowed per sample ID.
	--r1 [str]		Forward reads in format "infile.R1.fastq,sampleID [file2,id2 file3,id3...]", requires --r2 with corresponding sample IDs
	--r2 [str]		Reverse reads in format "infile.R2.fastq,sampleID", requires --r1 with corresponding sample IDs
	--SR [str]		Unpaired reads in format "infile.fastq,sampleID [file2,id2 file3,id3...]"
	--bams [str]		Compatible bam files in format "infile.bam,sampleID [file2,id2 file3,id3...]"
	
	Other options.
	--path_to_lobSTR [str]		Full path to lobSTR binaries
	--path_to_samtools [str]	Full path to samtools
	--path_to_bwa [str]		Full path to BWA
	--lob_args [str]		Additional command line arguments for running lobSTR (alignment), in quotes
	--allel_args [str]		Additional command line arguments for running allelotype, in quotes
	--mem_args [str]		Additional command line arguments for running bwa-mem (alignment), in quotes
	--backtrack_args [str]		Additional command line arguments for running bwa-backtrack (alignment), in quotes
	--noise [str]			Provide a noise model to use with allelotype (default uses PCR-free noise model)
	--samqual [n]			Minimum read alignment quality scores ("samtools view -q [n]"), with mem and backtrack
	--gzip				Read files are compressed using gzip (specify for lobSTR)
	--fasta				Reads are in fasta format (default fastq; specify for lobSTR)
	--minlen [n]			Minimum block length to include in analysis (400)
	--flank [n]			Minimum non-repeat flank length in block (200)
	--alnflank [n]			Minimum aligned read flank in backtrack only (15)
        --probes_per_locus [n]		Number of probes per locus [4]
        --probe_length [n]		length of probes [100]
        --probe_stagger	[n]		Stagger distance between probes [20]
	--mincvg [n]			Minimum callable-read coverage in backtrack (5)
	--minallele [n]			Minimum callable-reads per allele (2)
	--no_rmdup			Do not perform duplicate removal following read alignment
	--help				Help screen
	--silent			Run silently

#Details for running BaitSTR_type.pl

Each run requires, at a minimum:

	--stem [str]: This is the prefix used for all outfiles generated.
	--index_prefix [str]: Either a new index prefix or one referring to a previously built index

In addition, choose an alignment strategy:

	--mem uses bwa-MEM
	--lobSTR uses lobSTR alignment

And choose one or more functions:

	--full runs the entire pipeline, aligning reads, calling genotypes (vcf file output), and designing probes
	--align just aligns reads
	--allelotype just calls genotypes (using lobSTR's "allelotype" function)
	--design_probes designs probes

You must either provide an index prefix from a previous run of BaitSTR_type.pl, or create a new one. To create a new index:

	perl BaitSTR_type.pl --index --index_prefix [newPrefix] --target [blocks.fa] [...]

You can constrain the index to blocks of a minimum overall length (--minlen) or non-repeat flank (--flank).

If using --align or --full, you must provide reads, either paired (--r1 and --r2) or single (--SR).

Fastq read files are provide in a whitespace-separated list of "file,sample1 file,sample2 file,sample3...". You must provide matching --r1 and --r2 files to use paired reads, but you can provide as many read files as you like per sample. For example:

	perl BaitSTR_type.pl --index_prefix [indexPrefix] --align --mem --r1 Tatiana.R1.fastq.gz,Tatiana Oberon.R1.fastq.gz,Oberon --r2 Tatiana.R2.fastq.gz,Tatiana Oberon.R2.fastq.gz,Oberon

If you have previously built an alignment and would like to run allelotype only, provide bam files in the same way:

	perl BaitSTR_type.pl --index_prefix [indexPrefix] --allelotype --mem --bam Tatiana.bam,Tatiana Oberon.bam,Oberon
	

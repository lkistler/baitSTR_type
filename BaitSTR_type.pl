#!/usr/bin/perl

use Getopt::Long;

$help = 1 unless (@ARGV);
GetOptions (
	'help' => \$help,
	'target=s' => \$targ,
	'stem=s' => \$stem,
	'index_prefix=s' => \$index_prefix,
	'minlen=i' => \$minlen,
	'flank=i' => \$flank,
	'index' => \$index,
	'align' => \$align,
	'allelotype' => \$allel,
	'genotype' => \$genotype,
	'design_probes' => \$design_probes,
	'r1=s{1,}' => \@r1,
	'r2=s{1,}' => \@r2,
	'SR=s{1,}' => \@SR,
	'bams=s{1,}' => \@bams,
	'vcf=s' => \$vcf,
	'mincvg=i' => \$mincvg,
	'gzip' => \$gzip,
	'full_pipe' => \$full,
	'fasta' => \$fasta,
	'lob_args=s' => \$lob_args,
	'bwa_args=s' => \$bwa_args,
	'backtrack_args=s' => \$backtrack_args,
	'allel_args=s' => \$allel_args,
	'path_to_lobSTR=s' => \$path_to_lobSTR,
	'path_to_samtools=s' => \$path_to_samtools,
	'path_to_BWA=s' => \$path_to_BWA,
	'noise=s' => \$noise,
	'lobSTR' => \$lobSTR,
	'alnflank=i' => \$alnflank,
	'silent' => \$silent,
	'samqual=i' => \$samqual,
	'no_rmdup' => \$no_rmdup,
	'mem' => \$mem,
	'backtrack' => \$backtrack,
	'probes_per_locus=n' => \$probes_per_locus,
	'probe_length=n' => \$probe_length,
	'probe_stagger=i' => \$probe_stagger,
);

&psandqs;

if ($index and ($mem or $lobSTR)) {
	open IN, $targ;
	mkdir "$index_prefix.lobSTRindex";
	open MB, ">$index_prefix.lobSTRindex/lobSTR_mergedref.bed";
	open SB, ">$index_prefix.lobSTRindex/lobSTR_mergedref.targets.bed";
	open REFFA, ">$index_prefix.lobSTRindex/lobSTR_ref.fasta"; #if ($lobSTR);
	open BWAFA, ">$index_prefix.lobSTRindex/BWA_ref.fasta"; # if ($mem);
	open CHROM, ">$index_prefix.lobSTRindex/lobSTR_chromsizes.tab";
	open TAB, ">$index_prefix.lobSTRindex/lobSTR_ref_map.tab";
	open INFO, ">$index_prefix.lobSTRindex/strinfo.tab";
	print INFO "chrom\tstart\tend\tscore\tGC\tentropy\n";
	$count = 0;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /^>/) {
			$head = substr $line, 1;
			$seq = <IN>;
			chomp $seq;
			$len = length $seq;
			if ($len < $minlen) {
				next;
			}
			@d = split(/\s/, $head);
			$block = $d[0];
			$dat = $d[1];
			($motif, $cop, $beg, $fin) = split(/:/, $dat);
			if ($cop =~ /,/) {
				$cop = (split(/,/, $cop))[1];
			}
			if ($beg < $flank or ($len-$fin) < $flank) {
				next;
			}
			if ($beg < 1000) {
				$dif = 1000-$beg;
				$beg += $dif;
				$fin += $dif;
				$seq = ("N" x $dif).$seq;
			}
			if (($len-$fin) < $flank) {
				$dif = 1000-($len-$fin);
				$seq = $seq.("N" x $dif);
			}
			$beg++;
			$lo = $beg - 1000;
			$hi = $fin + 1000;
			@outblock = ($block, $lo, $hi, "$beg\_$fin\_$motif;");
			$dump = join("\t", @outblock);
			print MB "$dump\n";
			@outblock = ($block, $beg, $fin);
			$dump = join("\t", @outblock);
			print SB "$dump\n";
			$rc = reverse $motif;
			$rc =~ tr/[acgtACGT]/[tgcaTGCA]/;
			@pick = sort ($rc, $motif);
			$motif = $pick[0];
			$rc = $pick[1];
			$dump = "$beg\_$fin\_$rc\_$motif;";
			print TAB "$count\t$motif\t$dump\n";
			$seq = substr $seq, $lo, $hi-$lo+1;
			$seq = ("N" x 50).$seq.("N" x 50);
			$fahead = "$count\$$block\$$lo\$$hi";
			print REFFA "\>$fahead\n$seq\n";
			print BWAFA "\>$block\n$seq\n";
			print CHROM "$block\t$len\n";
			$gccheck = $seq;
			$gccheck =~ s/[NnGCgc]//g;
			$gc = 1-((length $gccheck)/$len);
			$flanks = (substr $seq, 0, 1000).(substr $seq, (length $seq)-1000);
			$ent = entropy($flanks);
			print INFO "$block\t$beg\t$fin\t1.0\t$gc\t$ent\n";
			$count++;
		}
	}
	close IN;
	if ($lobSTR) {
		if ($path_to_lobSTR) {
			print "Building Index. CMD: $path_to_lobSTR/lobSTRIndex index -a is $index_prefix.lobSTRindex/lobSTR_ref.fasta\n" unless ($silent);
			`$path_to_lobSTR/lobSTRIndex index -a is $index_prefix.lobSTRindex/lobSTR_ref.fasta $silent`;
		}
		else {
			print "Building Index. CMD: lobSTRIndex index -a is $index_prefix.lobSTRindex/lobSTR_ref.fasta\n" unless ($silent);
			`lobSTRIndex index -a is $index_prefix.lobSTRindex/lobSTR_ref.fasta $silent`;
		}
	}
	if ($mem) {
		print "Building Index. CMD: $bwacall index $index_prefix.lobSTRindex/BWA_ref.fasta\n" unless ($silent);
		`$bwacall index $index_prefix.lobSTRindex/BWA_ref.fasta $silent`;
	}
}
if ($index and $backtrack) {
        open IN, $targ;
        mkdir "$index_prefix.BWAindex";
        open MB, ">$index_prefix.BWAindex/$index_prefix.bed";
        open REFFA, ">$index_prefix.BWAindex/$index_prefix.fasta";
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /^>/) {
			$head = substr $line, 1;
			$seq = <IN>;
			chomp $seq;
			$len = length $seq;
			if ($len < $minlen) {
				next;
			}
			@d = split(/\s/, $head);
			$block = $d[0];
			$dat = $d[1];
			($motif, $cop, $beg, $fin) = split(/:/, $dat);
			if ($cop =~ /,/) {
				$cop = (split(/,/, $cop))[1];
			}
			if ($beg < $flank or ($len-$fin) < $flank) {
				next;
			}
			print REFFA "\>$block\n$seq\n";
			$beg++;
			print MB "$block\t$motif\t$cop\t$beg\t$fin\n";
		}
	}
	close REFFA;
	`$bwacall index -a is $index_prefix.BWAindex/$index_prefix.fasta $silent`;
}

if ($design_probes) {
	open PROBES, ">$index_prefix.probes.fasta";
	if ($backtrack) {
		open FA, "$index_prefix.BWAindex/$index_prefix.fasta";
		open BED, "$index_prefix.BWAindex/$index_prefix.bed";
		while ($line = <BED>) {
			chomp $line;
			@d = split /\s+/, $line;
			$chr = $d[0];
			$beg = $d[3];
			$fin = $d[4];
			$ping = "$beg.$fin";
			$probehold{$chr} = $ping;
		}
		close BED;
	}
	else {
		open BED, "$index_prefix.lobSTRindex/lobSTR_mergedref.targets.bed";
		open FA, "$index_prefix.lobSTRindex/BWA_ref.fasta";
		while ($line = <BED>) {
			@d = split /\s+/, $line;
			$chr = $d[0];
			$beg = $d[1];
			$fin = $d[2];
			$ping = "$beg.$fin";
			$probehold{$chr} = $ping;
		}
		close BED;
	}
	while ($line = <FA>) {
		chomp $line;
		if ($line =~ /^>/) {
			$head = substr $line, 1;
			$pull = $probehold{$head};
			($beg, $fin) = split /\./, $pull;
			$seq = <FA>;
			chomp $seq;
			$p1 = $probes_per_locus/2;
			$i = 1;
			$offset = 0;
			foreach $p (1..$p1) {
				$ping = $beg-$probe_length;
				$ping -= $offset;
				$probeseq = substr $seq, $ping, $probe_length;
				print PROBES "\>$head\tProbe$i\t$ping\n$probeseq\n";
				$offset += $probe_stagger;
				$i++;
			}
			$offset = 0;
			foreach $p (1..$p1) {
				$ping = $fin;
				$ping += $offset;
				$probeseq = substr $seq, $ping, $probe_length;
				print PROBES "\>$head\tProbe$i\t$ping\n$probeseq\n";
				$offset += $probe_stagger;
				$i++;
			}
		}
	}
	close PROBES;
}

if ($align and $mem or $lobSTR) {
	foreach $rg (keys %grps) {
		@grp = @{$grps{$rg}};
		my $cmd = "";
		$SR_wrk = "";
		$r1_wrk = "";
		$r2_wrk = "";
		undef $SRcmd;
		undef $PRcmd;
		foreach $file (@grp) {
			if ($type{$file} eq "SR") {
				$SR_wrk = $SR_wrk.",".$file;
			}
			if ($type{$file} eq "r1") {
				$r1_wrk = $r1_wrk.",".$file;
			}
			if ($type{$file} eq "r2") {
				$r2_wrk = $r2_wrk.",".$file;
			}
		}
		if ($SR_wrk) {
			$SR_wrk =~ s/^,//;
			$SRcmd = $cmd." -f $SR_wrk";
		}
		if ($r1_wrk) {
			$r1_wrk =~ s/^,//;
			$r2_wrk =~ s/^,//;
			$PRcmd = $cmd." --p1 $r1_wrk --p2 $r2_wrk";
		}

		if ($mem) {
			if ($SR_wrk) {
				$SR_wrk =~ s/,/\ /g;
				`$bwacall mem $bwa_args -aM -R "\@RG\\tID:lobSTR;sample_$rg;lib_$rg\\tLB:lib_$rg\\tSM:sample_$rg" $index_prefix.lobSTRindex/BWA_ref.fasta $SR_wrk $silent | $samcall view -F4 -q1 -Sb -o $stem.sample_$rg.SR.aligned.bam -`;
				`$samcall sort $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.SR.aligned`;
				`$samcall rmdup -S $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.SR.aligned.bam.tmp ; mv $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.SR.aligned.prermdup.bam ; mv $stem.sample_$rg.SR.aligned.bam.tmp $stem.sample_$rg.SR.aligned.bam` unless ($no_rmdup);
			}
			if ($r1_wrk) {
				$r1_wrk =~ s/,/\ /g;
				$r2_wrk =~ s/,/\ /g;
				`$bwacall mem $bwa_args -aM -R "\@RG\\tID:lobSTR;sample_$rg;lib_$rg\\tLB:lib_$rg\\tSM:sample_$rg" $index_prefix.lobSTRindex/BWA_ref.fasta $r1_wrk $r2_wrk | \
				$samcall view -F4 -q1 -Sb -o $stem.sample_$rg.PR.aligned.bam -`;
				#Split the bam file in 2 for rmdup
				`$samcall sort $stem.sample_$rg.PR.aligned.bam $stem.sample_$rg.PR.aligned`;
				unless ($no_rmdup) {
					`$samcall view -b -f2 -o $stem.sample_$rg.PR.aligned.bam.mated $stem.sample_$rg.PR.aligned.bam`;
					`$samcall view -b -F2 -o $stem.sample_$rg.PR.aligned.bam.unmated $stem.sample_$rg.PR.aligned.bam`;
					#`mv $stem.sample_$rg.PR.aligned.bam $stem.sample_$rg.PR.aligned.prermdup.bam`;
					`$samcall rmdup $stem.sample_$rg.PR.aligned.bam.mated $stem.sample_$rg.PR.aligned.bam.mated.rmdup`;
					`$samcall rmdup -S $stem.sample_$rg.PR.aligned.bam.unmated $stem.sample_$rg.PR.aligned.bam.unmated.rmdup`;
					`$samcall merge -f $stem.sample_$rg.PR.aligned.bam $stem.sample_$rg.PR.aligned.bam.mated.rmdup $stem.sample_$rg.PR.aligned.bam.unmated.rmdup`;
					`$samcall sort $stem.sample_$rg.PR.aligned.bam $stem.sample_$rg.PR.aligned`;
					`rm $stem.sample_$rg.PR.aligned.bam.mated.rmdup $stem.sample_$rg.PR.aligned.bam.unmated.rmdup`;
					`rm $stem.sample_$rg.PR.aligned.bam.mated $stem.sample_$rg.PR.aligned.bam.unmated`;
				}
			}
		}	
		if ($lobSTR) {
			if ($gzip) {
				$SRcmd = $SRcmd." --gzip" if ($SRcmd);
				$PRcmd = $PRcmd." --gzip" if ($PRcmd);
			}
			unless ($fasta) {
				$SRcmd = $SRcmd." -q" if ($SRcmd);
				$PRcmd = $PRcmd." -q" if ($PRcmd);
			}
			$SRcmd = $SRcmd." --out $stem.sample_$rg.SR" if ($SRcmd);
			$PRcmd = $PRcmd." --out $stem.sample_$rg.PR" if ($PRcmd);
			$SRcmd = $SRcmd." --rg-sample sample_$rg --rg-lib lib_$rg" if ($SRcmd);
			$PRcmd = $PRcmd." --rg-sample sample_$rg --rg-lib lib_$rg" if ($PRcmd);
			$SRcmd = $SRcmd." --index-prefix $index_prefix.lobSTRindex/lobSTR_" if ($SRcmd);
			$PRcmd = $PRcmd." --index-prefix $index_prefix.lobSTRindex/lobSTR_" if ($PRcmd);
			if ($lob_args) {
				$SRcmd = $SRcmd." $lob_args" if ($SRcmd);
				$PRcmd = $PRcmd." $lob_args" if ($PRcmd);
			}
			if ($path_to_lobSTR) {
				if ($SRcmd) {
					print "\n\nCMD: $path_to_lobSTR/lobSTR $SRcmd\n\n" unless ($silent);
					`$path_to_lobSTR/lobSTR $SRcmd $silent`;
				}
				if ($PRcmd) {
					print "\n\nCMD: $path_to_lobSTR/lobSTR $PRcmd\n\n" unless ($silent);
					`$path_to_lobSTR/lobSTR $PRcmd $silent`;
				}
			}
			else {
				if ($SRcmd) {
					print "\n\nCMD: lobSTR $SRcmd\n\n" unless ($silent);
					`lobSTR $SRcmd`;
				}
				if ($PRcmd) {
					print "\n\nCMD: lobSTR $PRcmd\n\n" unless ($silent);
					`lobSTR $PRcmd`;
				}
			}
		}
		if ($SRcmd and $PRcmd) {
			`$samcall merge -f $stem.sample_$rg.aligned.bam $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.PR.aligned.bam`;
			`rm $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.PR.aligned.bam`;
		}
		elsif ($SRcmd) {
			`mv $stem.sample_$rg.SR.aligned.bam $stem.sample_$rg.aligned.bam`;
		}
		elsif ($PRcmd) {
			`mv $stem.sample_$rg.PR.aligned.bam $stem.sample_$rg.aligned.bam`;
		}
		push @newbams, "$stem.sample_$rg.aligned.bam";
		$filekey{"$stem.sample_$rg.aligned.bam"} = $rg;
	}
	foreach $i (@newbams) {
		$pref = $i;
		$pref =~ s/\.bam$//;
		`$samcall sort $i $pref ; $samcall index $i`;
	}
}

if ($backtrack) {
	foreach $rg (keys %grps) {
		@grp = @{$grps{$rg}};
		my $cmd = "";
		undef @SR;
		undef @r1;
		undef @r2;
		undef @bam;
		foreach $file (@grp) {
			if ($type{$file} eq "SR") {
				push @SR, $file;
			}
			if ($type{$file} eq "r1") {
				push @r1, $file;
			}
			if ($type{$file} eq "r2") {
				push @r2, $file;
			}
			if ($type{$file} eq "bam") {
				push @bam, $file;
			}
		}
		undef @midbams;
		if (@SR) {
			foreach $file (@SR) {
				$filestem = $file;
				$filestem =~ s/^.+\///g;
				`$bwacall aln $bwa_args $index_prefix.lobSTRindex/BWA_ref.fasta $file $silent |\
				$bwacall samse $index_prefix.lobSTRindex/BWA_ref.fasta - $file $silent |\
				$samcall view -Sb -F4 $samqual -o $stem.$filestem.bam - $silent`;
				`$samcall sort $stem.$filestem.bam $stem.$filestem`;
				`$samcall rmdup -S $stem.$filestem.bam  $stem.$filestem.bam.tmp ; mv  $stem.$filestem.bam.tmp $stem.$filestem.bam` unless ($no_rmdup);
				push @midbams, "$stem.$filestem.bam";
			}
		}
		if (@r1) {
			foreach $i (0..(@r1-1)) {
				$f1 = $r1[$i];
				$f2 = $r2[$i];
				$f1stem = $f1;
				$f2stem = $f2;
				$f1stem =~ s/^.+\///g;
				$f2stem =~ s/^.+\///g;
				`$bwacall aln $bwa_args $index_prefix.BWAindex/$index_prefix.fasta $f1 $silent > $stem.$f1stem.sai`;
				`$bwacall aln $bwa_args $index_prefix.BWAindex/$index_prefix.fasta $f2 $silent > $stem.$f2stem.sai`;
				`$bwacall sampe $index_prefix.BWAindex/$index_prefix.fasta $stem.$f1stem.sai $stem.$f2stem.sai $f1 $f2 |\
				samtools view -Sb -F4 $samqual -o $stem.$f1stem.bam - $silent`;
				`$samcall sort $stem.$f1stem.bam $stem.$f1stem`;
				`$samcall rmdup $stem.$f1stem.bam $stem.$f1stem.bam.tmp ; mv $stem.$f1stem.bam.tmp $stem.$f1stem.bam` unless ($no_rmdup);
				`rm $stem.$f1stem.sai $stem.$f2stem.sai`;
				push @midbams, "$stem.$f1stem.bam";
			}
		}
		if (@bam) {
			foreach $file (@bam) {
				push @midbams, $file;
			}
		}
		if (@midbams == 1) {
			`mv $midbams[0] $stem.sample_$rg.aligned.bam`;
		}
		elsif (@midbams > 1) {
			$mergestring = join("\ ", @midbams);
			`$samcall merge -f $stem.sample_$rg.aligned.bam $mergestring`;
		}
		`$samcall sort $stem.sample_$rg.aligned.bam $stem.sample_$rg.aligned`;
		push @newbams, "$stem.sample_$rg.aligned.bam";
		$filekey{"$stem.sample_$rg.aligned.bam"} = $rg;
	}
	open CALLS, ">$stem.rawCalls.txt";
	print CALLS "#Contig\tMotif\trefCopies\t";
	open BED, "$index_prefix.BWAindex/$index_prefix.bed";
	while ($line = <BED>) {
		chomp $line;
		@d = split(/\s+/, $line);
		($contig, $motif, $ref, $beg, $fin) = @d;
		shift @d;
		$hold = join("\t", @d);
		$contigs{$contig} = $hold;
		push @contigMap, $contig;
	}
	if (@newbams) {
		@bams = @newbams
	}
	foreach $bam (@bams) {
		$id = $filekey{$bam};
		print CALLS "sample_$id\t";
	}
	print CALLS "\n";
	foreach $bam (@bams) {
		$. = 0;
		$id = $filekey{$bam};
		undef %gtcall;
		open BAM, "$samcall view $samqual $bam |";
		$called = "$stem.sample_$id.called.bam";
		open OUTBAM, "| $samcall view -Sb -o $called -T $index_prefix.BWAindex/$index_prefix.fasta -";
		print "Reading $bam and writing $called\n";
		$check = "#####";
		while ($line = <BAM>) {
			chomp $line;
			@d = split(/\s+/, $line);
			($block, $pos1, $cigar, $pair) = @d[2, 3, 5, 7];
			if ($. == 1) {
				$check = $block;
				undef %duphold;
				undef %duphold2;
				undef %midcall;
				undef @dupmap;
			}
			if ($block ne $check) {
				$k1 = 0;
				$k2 = 0;
				$tot = 0;
				foreach $call (keys %calls) {
					$tot++;
					if ($calls{$call}>$m1) {
						$m2 = $m1;
						$k2 = $k1;
						$m1 = $calls{$call};
						$k1 = $call;
					}
				}
				undef %calls;
				foreach $dupstat (@dupmap) {
					@arr = @{$duphold{$dupstat}};
					@arr2 = @{$duphold2{$dupstat}};
					if (@arr == 1) {
						print OUTBAM $arr[0]."\n";
						$calls{$arr2[0]}++;
					}
					else {
						$pull = index(@arr2, $k1);
						if ($pull == -1) {
							$pull = index(@arr2,$k2);
						}
						if ($pull == -1) {
							$pull = int(rand(@arr2));
						}
						$read = $arr[$pull];
						$calls{$arr2[$pull]}++;
						print OUTBAM $read."\n";
					}
							
				}
				&dumpCalls;
				$check = $block;
				undef %duphold;
				undef %duphold2;
				undef %midcall;
				undef @dupmap;
			}
			$dupstat = "$pos1.$pair";
			$hit = statread($line);
			if ($hit !~ /^$/) {
				push @dupmap, $dupstat unless (exists ($duphold{$dupstat}));
				$midcall{$line} = $hit;
				push @{$duphold{$dupstat}}, $line;
				push @{$duphold2{$dupstat}}, $hit;
			}
		}
		close BAM;
		close OUTBAM;
		undef %calls;
		$hashref = "$bam.gt";
		%{$hashref} = %gtcall;
		push @master, $hashref;
		print "Done reading $called\n";
	}
	print "Done with calling, summarizing calls\n";
	foreach $block (@contigMap) {
		($motif, $ref, $beg, $fin) = split(/\s+/, $contigs{$block});
		print CALLS "$block\t$motif\t$ref\t";
		foreach $rec (@master) {
			$call = ${$rec}{$block};
			unless ($call) {$call = "NA"};
			print CALLS $call."\t";
		}
		print CALLS "\n";
	}
}

if ($allel) {
	die "\nPlease run with \"--align\" or provide one or more bam files to run allelotype\n\n" unless (@bams or $align);
	if ($path_to_lobSTR) {
		`$path_to_lobSTR/allelotype --version 2> $stem.allvers.txt`;
		$version = `cat $stem.allvers.txt`;
		`rm $stem.allvers.txt`;
	}
	else {
		`allelotype --version 2> $stem.allvers.txt`;
		$version = `cat $stem.allvers.txt`;
		`rm $stem.allvers.txt`;
	}
	if ($version =~ /^3/) {
		$cmd = "--command classify --strinfo $index_prefix.lobSTRindex/strinfo.tab --out $stem --index-prefix $index_prefix.lobSTRindex/lobSTR_ --dont-include-flank";
	}
	if ($version =~ /^4/) {
		$cmd = "--command classify --strinfo $index_prefix.lobSTRindex/strinfo.tab --out $stem --index-prefix $index_prefix.lobSTRindex/lobSTR_ --regions $index_prefix.lobSTRindex/lobSTR_mergedref.targets.bed";
	}
	if ($mem or $bwa) {
		$cmd = $cmd." --realign --filter-clipped --min-read-end-match 10 --filter-mapq0 --max-repeats-in-ends 3";
	}
	$cmd = $cmd." --no-rmdup" if ($mem or $no_rmdup);
	
	unless ($noise) {
		&tmpnoise;
		$noise = "$stem.noisetmp";
		$modelclean = 1;
	}
	$cmd = $cmd." --noise_model $noise";
	@bams = (@bams, @newbams);
	$cmd = $cmd." --bam $bams[0]";
	shift @bams;
	foreach $i (@bams) {
		$cmd = $cmd.",$i";
	}
	if ($allel_args) {
		$cmd = "$cmd $allel_args";
	}
	if ($path_to_lobSTR) {
		print "\nCMD: $path_to_lobSTR/allelotype $cmd\n\n" unless ($silent);
		`$path_to_lobSTR/allelotype $cmd $silent`;
	}
	else {
		print "\nCMD: allelotype $cmd\n\n" unless ($silent);
		`allelotype $cmd`;
	}
	`rm $stem.noisetmp.stuttermodel $stem.noisetmp.stepmodel` if ($modelclean);
}

sub entropy {
	$flankseq = $_[0];
	foreach $i ("A", "C", "G", "T") {
		$calc{$i} = 0;
	}
	foreach $i (split("", $flankseq)) {
		$calc{$i}++ if ($i =~ /[ACGT]/);
	}
	$total = 0;
	foreach $i (values %calc) {
		$total += $i;
	}
	foreach $i (keys %calc) {
		$calc{$i} = $calc{$i}*(1/$total);
	}
	$entropy = 0;
	foreach $i (keys %calc) {
		$entropy += -1*$calc{$i}*(log($calc{$i})/log(2)) unless ($calc{$i} == 0);
	}
	return $entropy;
}

sub tmpnoise {
	open NOISE, ">$stem.noisetmp.stuttermodel";
	print NOISE "solver_type L2R_LR
nr_class 2
label -1 1
nr_feature 4
bias 0
w
-7.36733
-0.865223
0.0824449
-1.04024
4.17049
";
	close NOISE;
	open NOISE, ">$stem.noisetmp.stepmodel";
	print NOISE "0.0001
0.0166667
0.555556
0.241379
0.25
0.0001
ProbIncrease=0.221239
Period1Model 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0778761 0.233628 0.467257 0 0.132743 0.0663717 0.0221239 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Period2Model 0 0 0 0 0 0 0 0 0 0 0 0 0.0245946 0.00120957 0.184459 0.00604784 0.553378 0.00907176 0 0.0025772 0.15721 0.00171814 0.0524032 0.000343627 0.0069871 0 0 0 0 0 0 0 0 0 0 0 0
Period3Model 0 0 0 0 0 0 0 0 0 0.00406671 0.00435719 0.0116192 0.0759119 0.0542228 0.0903714 0.337386 0.120495 0.0803301 0 0.0228211 0.0342316 0.0958484 0.0256737 0.0154042 0.0215659 0.0033009 0.00123784 0.00115532 0 0 0 0 0 0 0 0 0
Period4Model 0 0 0 0 0 0 0.00150203 0.000876185 0.000468515 0.00602377 0.0697036 0.027107 0.0092239 0.0711558 0.45743 0.0889447 0.0129711 0.0333543 0 0.00947565 0.00368497 0.0252684 0.129952 0.0202147 0.00262043 0.00770084 0.0198022 0.0017113 0.000133101 0.000248916 0.000426713 0 0 0 0 0 0
Period5Model 0 0 0 0.000422567 0.00025354 0.000141982 0.000369154 0.00442985 0.0487284 0.0194914 0.00701689 0.011227 0.0785891 0.471535 0.0943069 0.0150891 0.00905347 0.0181069 0 0.00514401 0.00257201 0.00428668 0.0267917 0.133959 0.0223265 0.00318949 0.00199343 0.00553732 0.0138433 0.00125848 0.000104873 4.03359e-05 7.20284e-05 0.000120047 0 0 0
Period6Model 0.000178069 5.34296e-08 1.51409e-11 4.03824e-15 1.00939e-10 2.35486e-06 0.0510135 1.02044e-05 1.87112e-09 3.11905e-13 4.67779e-09 6.23602e-05 0.727415 7.27536e-05 6.0638e-09 4.04321e-13 2.02127e-09 6.73644e-06 0 1.91376e-06 5.74224e-10 1.14864e-13 1.72267e-09 2.06686e-05 0.206652 1.7716e-05 1.32892e-09 8.86092e-14 5.31567e-10 2.89898e-06 0.0144925 6.68995e-07 2.8676e-11 1.14723e-15 4.30139e-12 1.51789e-08 5.05878e-05
Period1Obs 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 101 22 39 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Period2Obs 0 0 0 0 0 0 0 0 0 0 0 0 4 1 13 0 55 0 25 1 9 0 2 0 10 0 0 0 0 0 0 0 0 0 0 0 0
Period3Obs 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 3 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0
Period4Obs 0 0 0 0 0 0 1 0 0 0 2 0 0 1 6 1 0 4 8 1 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Period5Obs 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 3 4 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0
Period6Obs 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
";
	close NOISE;
}

sub psandqs {
$help = 1 unless ($index or $align or $allel or $genotype or $design_probes or $full);

&helpcall if ($help);

if ($silent) {
	$silent = "2\> /dev/null";
}
else {
	$silent = "";
}

print "\nWARNING: Choosing --bwa and --lobSTR will take a while.\n\n" if ($bwa and $lobSTR);
die "\nProvide a prefix for this run with \"--stem\"\n\n" unless ($stem);

if ($full) {
	$align = 1;
	$allel = 1;
	$genotype = 1;
	$design_probes = 1;
}

die "\nProvide an index prefix with \"--index_prefix\". You can use one from a previous run, or create a new one and use \"--index\" if you also provide \"--target\".\n\n" unless ($index_prefix);
$flank = 200 unless ($flank);
$minlen = 2*$flank unless ($minlen);
$mincvg = 5 unless ($mincvg);
$loci = 5000 unless ($loci);
$probes_per_locus = 4 unless ($probes_per_locus);
$probe_length = 100 unless ($probe_length);
$probe_stagger = 20 unless ($probe_stagger);
if ($index or $align or $allel) {
	unless ($path_to_lobSTR) {
		$pathcheck = `which lobSTR`;
		chomp $pathcheck;
		unless ($pathcheck =~ /lobSTR$/) {
			die "\nlobSTR not found, please provide a full path using \"--path_to_lobSTR\"\n\n";
		}
	}
}
if ($align) {
	if ($path_to_samtools) {
		$samcall = "$path_to_samtools/samtools";
	}
	else {
		$path_to_samtools = "";
		$pathcheck = `which samtools`;
		chomp $pathcheck;
		unless ($pathcheck =~ /samtools$/) {
			die "\nsamtools not found, please provide a full path using \"--path_to_samtools\"\n\n";
		}
		$samcall = "samtools";
	}
}

if ($path_to_BWA) {
	$bwacall = "$path_to_BWA/bwa";
}
else {
	$bwacall = "bwa";
}

if ($lobSTR) {
	$count++;
}
if ($mem) {
	$count++;
}
if ($backtrack) {
	$count++;
}

if ($lobSTR and @SR) {
	unless ($gzip) {
		foreach $set (@SR) {
			$file = (split /\,/, $set)[0];
			if ($file =~ /.gz$/) {
				die "\nIt looks like $file is compressed.\n\nlobSTR requires that you specify this at the command line with \"--gzip\"\n\n";
			}
		}
	}
}
		
if ($lobSTR and @r1) {
	unless ($gzip) {
		foreach $set (@r1) {
			$file = (split /\,/, $set)[0];
			if ($file =~ /.gz$/) {
				die "\nIt looks like $file is compressed.\n\nlobSTR requires that you specify this at the command line with \"--gzip\"\n\n";
			}
		}
		foreach $set (@r2) {
			$file = (split /\,/, $set)[0];
			if ($file =~ /.gz$/) {
				die "\nIt looks like $file is compressed.\n\nlobSTR requires that you specify this at the command line with \"--gzip\"\n\n";
			}
		}
	}
}

die "\nPlease only choose one method for this run (lobSTR, mem, or backtrack).\n\n" if ($count > 1);
die "\nPlease specify --mem or --lobSTR for alignment method\n\n" if ($count == 0);

unless ($index) {
	if ($lobSTR) {
		@reqd = (
		"$index_prefix.lobSTRindex/lobSTR_chromsizes.tab", 
		"$index_prefix.lobSTRindex/lobSTR_mergedref.bed", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.amb", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.ann", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.bwt", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.pac", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.rbwt", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.rpac", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.rsa", 
		"$index_prefix.lobSTRindex/lobSTR_ref.fasta.sa", 
		"$index_prefix.lobSTRindex/lobSTR_ref_map.tab", 
		"$index_prefix.lobSTRindex/strinfo.tab");
		foreach $i (@reqd) {
			die "\nCouldn't find $i, please re-run with \"--index\"\n\n" unless (-f $i);
		}
	}
	if ($mem or $backtrack) {
		@reqd = (
		"$index_prefix.BWAindex/$index_prefix.fasta",
		"$index_prefix.BWAindex/$index_prefix.fasta.amb",
		"$index_prefix.BWAindex/$index_prefix.fasta.ann",
		"$index_prefix.BWAindex/$index_prefix.fasta.bwt",
		"$index_prefix.BWAindex/$index_prefix.fasta.pac",
		"$index_prefix.BWAindex/$index_prefix.fasta.sa",
		"$index_prefix.BWAindex/$index_prefix.bed");
		foreach $i (@reqd) {
			die "\nCouldn't find $i, please re-run with \"--index\"\n\n" unless (-f $i);
		}
		$alnflank = 15 unless ($alnflank);
	}
}

if (@r1 or @r2 or @SR) {
	print "Without choosing --align read files are useless. Ignoring.\n" unless ($align);
}
if ($design_probes) {
	die "\nFlank size isn't large enough to design your probes!\n\n" if ($probe_length+((($probes_per_locus/2)-1)*$probe_stagger) > $flank);
}
if ($genotype) {
	die "\nProvide a vcf file or use --allelotype to run --genotype\n\n" unless ($vcf or $allel);
}
if ($vcf and $allel) {
	print "\nWARNING: using --allelotype creates a new vcf file, so $vcf will be ignored. Run without --allelotype to parse an existing vcf file\n\n";
}
if ($index) {
	die "\nPlease provide a target --target fasta file for indexing. This MUST be in the format output by extend_STR_reads.\n\n" unless ($targ);
}

if ($align) {
	die "\nPlease provide read files if you would like to align reads.\n\n" unless (@r1 or @r2 or @SR);
}

die "\nYou must provide corresponding forward and reverse read files.\n\n" unless (@r2 == @r1);

foreach $i (@SR) {
	($file, $grp) = split(/,/, $i);
	die "\nCouldn't find $file\n\n" unless (-f $file);
	die "\nPlease provide a sample identifier for $file.\n\n" unless ($grp);
	@{$grps{$grp}} = () unless (exists ($grps{$grp}));
	push @{$grps{$grp}}, $file;
	$type{$file} = "SR";
	die "\n$file used more than once\n\n" if (exists ($check{$file}));
	$filekey{$file} = $grp;
	$check{$file} = 0;
}

foreach $i (@r1) {
	($file, $grp) = split(/,/, $i);
	die "\nCouldn't find $file\n\n" unless (-f $file);
	die "\nPlease provide a sample identifier for $file.\n\n" unless ($grp);
	@{$grps{$grp}} = () unless (exists ($grps{$grp}));
	push @{$grps{$grp}}, $file;
	$type{$file} = "r1";
	die "\n$file used more than once\n\n" if (exists ($check{$file}));
	$filekey{$file} = $grp;
	$check{$file} = 0;
}

foreach $i (@r2) {
	($file, $grp) = split(/,/, $i);
	die "\nCouldn't find $file\n\n" unless (-f $file);
	die "\nPlease provide a sample identifier for $file.\n\n" unless ($grp);
	@{$grps{$grp}} = () unless (exists ($grps{$grp}));
	push @{$grps{$grp}}, $file;
	$type{$file} = "r2";
	die "\n$file used more than once\n\n" if (exists ($check{$file}));
	$filekey{$file} = $grp;
	$check{$file} = 0;
}

foreach $i (@bams) {
	($file, $grp) = split(/,/, $i);
	die "\nCouldn't find $i\n\n" unless (-f $file);
	die "\nPlease provide a sample identifier for $file.\n\n" unless ($grp);
	@{$grps{$grp}} = () unless (exists ($grps{$grp}));
	push @{$grps{$grp}}, $file;
	$type{$file} = "bam";
	die "\n$file used more than once\n\n" if (exists ($check{$file}));
	$filekey{$file} = $grp;
	$check{$file} = 0;
	$i = $file;
}

if ($samqual) {
	$samqual = "-q $samqual";
}

else {
	$samqual = "";
}

if ($path_to_lobSTR) {
	`$path_to_lobSTR/allelotype --version 2> $stem.allvers.txt`;
	$version = `cat $stem.allvers.txt`;
	`rm $stem.allvers.txt`;
}

else {
	`allelotype --version 2> $stem.allvers.txt`;
	$version = `cat $stem.allvers.txt`;
	`rm $stem.allvers.txt`;
}

if ($version !~ /^4/ and $mem) {
	chomp $version;
	die "\nERROR: lobSTR version 4 is required for allelotyping a bwa alignment, you're using version $version. Please use --path_to_lobSTR to give the path to version 4.\n\n" if ($mem);
}

if ($backtrack) {
	$allel = 0;
	print "\nWARNING: Using the BWA-backtrack algorithm is not recommended, and should be considered deprecated.
Data will be passed through the entire pipeline and genotypes called without quality scores or other support.
Parameters are fixed and cannot be modified from the command line.
Please consider BWA-mem with lobSTR v.4 for best results.\n\n";
}
}

sub statread {
	$read = $_[0];
	@d = split(/\s+/, $read);
	($block, $pos1, $cigar, $pair, $seq, $qual) = @d[2, 3, 5, 6, 9, 10];
	$pos2 = $pos1 + length $seq;
	die "\nERROR: Wrong reference! $block is found in $bam, but not in $index_prefix.BWAindex/$index_prefix.bed\n\n" unless (exists ($contigs{$block}));
	($motif, $ref, $beg, $fin) = split(/\s+/, $contigs{$block});
	@cig = split(/(?<=[SMID])/, $cigar);
	if ($cig[0] =~ /[HS]/) {
		$dif = $cig[0];
		$dif =~ s/[HS]//;
		$pos1 += $dif;
		shift @cig;
	}
	if ($cig[-1] =~ /[HS]/) {
		$dif = $cig[-1];
		$dif =~ s/[HS]//;
		$pos2 -= $dif;
		pop @cig;
	}
	if (@cig == 1) {
		if ($pos1 <= $beg-$alnflank and $pos2 >= $fin+$alnflank) {
			$calls{0}++;
			push @readres, $read;
			$readcount++;
			return "0";
		}
	}
	elsif ($pos1 > $beg-$alnflank) {
		next;
	}
	else {
		while (@cig > 0) {
			$elem = $cig[0];
			($len, $type) = split(/(?=[IDM])/, $elem);
			if ($pos1 + $len < $beg - $alnflank) {
				shift @cig;
				$pos1 += $len;
				next;
			}
			elsif ($pos1 + $len < $beg) {
				undef @cig;
				next;
			}
			elsif ($pos1 + $len > $fin) {
				if ($pos1 + $len < $fin + $alnflank) {
					undef @cig;
					next;
				}
				else {
					$calls{0}++;
					push @readres, $read;
					$readcount++;
					undef @cig;
					return "0";
					next;
				}
			}
			else {
				if ($type == "M") {
					shift @cig;
					$pos1 += $len;
					$call = shift @cig;
					($len, $type) = split(/(?=[IDM])/, $call);
					if ($type =~ /D/) {
						$len = 0-$len;
					}
					$back = shift @cig;
					$back =~ s/M//;
					if ($pos1 + ($fin-$beg) + $back >= $len + $fin + $alnflank) {
						if ($len % (length $motif) == 0) {
							$calls{$len}++;
							push @readres, $read;
							$readcount++;	
							return $len;
						}
					}
				}
				last;
			}
		}
	}
}

sub pickDup {
	@panel = @_;
	foreach $read (@panel) {
		$cigar = (split(/\s+/, $read))[5];
		$pick{$cigar}++;
		$deets{$cigar} = $read;
	}
	$max = (sort {$a <=> $b} values %pick)[-1];
	foreach $read (keys %pick) {
		if ($pick{$read} == $max) {
			$return = $deets{$read};
		}
	}
	return $return;
}
	

sub dumpCalls {
	undef $sum;
	$cvg = 0;
	foreach $var (sort {$a <=> $b} keys %calls) {
		$add = "$var|".$calls{$var};
		$cvg += $calls{$var};
		if ($sum) {
			$sum = $sum.",".$add;
		}
		else {
			$sum = $add;
		}
	}
	if ($sum) {
		$gtcall{$check} = "$sum;$cvg";
	}
}

sub helpcall {
	die "

Usage: perl BaitSTR_type.pl --stem [runid] --index_prefix [index_prefix] [options]

Example:

perl BaitSTR_type.pl --index --mem --full --target [blocks.fa] --stem [run_prefix] \\
--r1 [reads.sample1.R1.fq,sampleID1 reads.sample2.R1.fq,sampleID2 ...] \\
--r2 [reads.sample1.R2.fq,sampleID1 reads.sample2.R2.fq,sampleID2 ...] \\
--SR [singleReads.sample3.fq,sampleID3 ...]
Creates an index and runs the entire pipeline (mapping through genotype calling) with defaults using BWA-mem to align reads and lobSTR defaults to call genotypes.

Requirements:
	SAMtools in PATH (or give path with --path_to_samtools)
	lobSTR in PATH to use lobSTR functionality (or give path with --path_to_lobSTR)
	BWA in PATH to use BWA functionality (or give path with --path_to_BWA)

Options:

	--index			Create index files (requires --target)
	--align			Perform read alignemnt (requires --r1 and --r2, and/or --SR)
	--allelotype		Perform allelotype calling (requires --align or --bams)
	--design_probes		Produce a fasta file of probes from qualifying input blocks
	--full			Do alignment, allelotype, genotype, and probe design (full pipeline excluding index)
	--mem			Use bwa-mem for alignemnt (uses lobSTR for allelotyping)
	--lobSTR		Use lobSTR for alignment (uses lobSTR for allelotyping)
	--backtrack		Use bwa-backtrack for alignment DEPRECATED AND NOT RECOMMENDED (skips lobSTR allelotyping and generates non-standard genotype file)
	
	--stem [str]		Prefix for run files (will overwrite without warning)
	--index_prefix [str]	Prefix for index files
	--target [str]		Fasta file of extended blocks from BaitSTR output

	Read files. Sample IDs must be given for each read file, multiple files are allowed per sample ID.
	--r1 [str]		Forward reads in format \"infile.R1.fastq,sampleID [file2,id2 file3,id3...]\", requires --r2 with corresponding sample IDs
	--r2 [str]		Reverse reads in format \"infile.R2.fastq,sampleID\", requires --r1 with corresponding sample IDs
	--SR [str]		Unpaired reads in format \"infile.fastq,sampleID [file2,id2 file3,id3...]\"
	--bams [str]		Compatible bam files in format \"infile.bam,sampleID [file2,id2 file3,id3...]\"
	
	Other options.
	--path_to_lobSTR [str]		Full path to lobSTR binaries
	--path_to_samtools [str]	Full path to samtools
	--path_to_bwa [str]		Full path to BWA
	--lob_args [str]		Additional command line arguments for running lobSTR (alignment), in quotes
	--allel_args [str]		Additional command line arguments for running allelotype, in quotes
	--mem_args [str]		Additional command line arguments for running bwa-mem (alignment), in quotes
	--backtrack_args [str]		Additional command line arguments for running bwa-backtrack (alignment), in quotes
	--noise [str]			Provide a noise model to use with allelotype (default uses PCR-free noise model)
	--samqual [n]			Minimum read alignment quality scores (\"samtools view -q [n]\"), with mem and backtrack
	--gzip				Read files are compressed using gzip (specify for lobSTR)
	--fasta				Reads are in fasta format (default fastq; specify for lobSTR)
	--minlen [n]			Minimum block length to include in analysis (400)
	--flank [n]			Minimum non-repeat flank length in block (200)
	--alnflank [n]			Minimum aligned read flank in backtrack (15)
        --probes_per_locus [n]		Number of probes per locus [4]
        --probe_length [n]		length of probes [100]
        --probe_stagger	[n]		Stagger distance between probes [20]
	--mincvg [n]			Minimum callable-read coverage in backtrack (5)
	--no_rmdup			Do not perform duplicate removal following read alignment
	--help				And here we are
	--silent			Run silently
"
}

#!/usr/bin/perl

# number of processors to be used
$jobNumbers = 240;
$runMSA = 0; # Don't run MSA
$runHMM = 0; # Don't run HMMsearch
$runEmit = 0; # Don't run HMMemit
$runOutSwiss = 0; # Don't run HMHsearch againt swissprot
$runOutEmit = 0; # Don't run HMMsearch against HMMemit
$runBlast = 0; # Don't run Blast
$runAnnotation = 0; # Don't run annotation
$runFilterScan = 0; # Don't filter the hmmscan results

# Parse the command line arguments.  
for($i = 0; $i < @ARGV; $i++) {
  if($ARGV[$i] =~ /^\-n/) { $jobNumbers = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-\-msa/) { $runMSA = 1;  next; }
  if($ARGV[$i] =~ /^\-\-hmm/) { $runHMM = 1;  next; }
  if($ARGV[$i] =~ /^\-\-emit/) { $runEmit = 1;  next; }
  if($ARGV[$i] =~ /^\-\-outswiss/) { $runOutSwiss = 1;  next; }
  if($ARGV[$i] =~ /^\-\-outemit/) { $runOutEmit = 1;  next; }
  if($ARGV[$i] =~ /^\-\-blast/) { $runBlast = 1;  next; }
  if($ARGV[$i] =~ /^\-\-annotate/) { $runAnnotation = 1;  next; }
  if($ARGV[$i] =~ /^\-\-filter/) { $runFilterScan = 1;  next; }
}


if(not $runMSA and not $runHMM and not $runEmit and not $runOutSwiss and not $runOutEmit and not $runBlast and not $runAnnotation and not $runFilterScan){
	die("Please indicate if you want to run MSA (use option --msa), run HMM build (use option --hmm), run HMM emit (use option --emit), run HMM search against SwissProt (use option --outswiss), or run HMM search against emitted sequences (use option --outemit)! Please provide one and only one of these options.");
}

# working directory
$wd = "/data1/pcl/swissprot";

# programs to be called
$muscle = "$wd/bin/muscle";

# hmm2 commands
$hmm2build = "$wd/bin/hmmer2/hmmbuild";
$hmm2search = "$wd/bin/hmmer2/hmmsearch";
$hmm2calibrate = "$wd/bin/hmmer2/hmmcalibrate";

# hmm3 commands
$hmm3build = "$wd/bin/hmmer3/hmmbuild";
$hmm3search = "$wd/bin/hmmer3/hmmsearch";
$hmm3emit = "$wd/bin/hmmer3/hmmemit";
$hmm3scan = "$wd/bin/hmmer3/hmmscan";
$blastp = "$wd/bin/ncbi-blast-2.2.23+/bin/blastp";

# filter hmmscan
$annotateGenome = "perl $wd/annotateGenome2.pl"; 

# Input and Output directories
$swissprot = "$wd/uniprot_sprot.fasta";
$seqDir = "$wd/seq"; 
$msaDir = "$wd/msa";
$hmmDir = "$wd/hmm";
$emitDir = "$wd/emit";
$outswissDir = "$wd/outswiss";
$outemitDir = "$wd/outemit";
$qsubDir = "$wd/qsub";
$blastOut = "$wd/blastout";
$blastSeq = "$wd/blastseq";
$blastDatabase = "$wd/swissProtein.fasta";
$genomeFaa = "$wd/genomes";
$genomeOutput = "$wd/scanout";
$dbhmm = "$wd/dbhmm/swiss.hmm";
$annotation = "$wd/tabanno";

# remove previous qsub scripts in the existing directories
# and create new directories
if(-d $qsubDir){ system("rm -rf $qsubDir"); }
mkdir "$qsubDir", 0755 or die("couldn't create qsub directory");



if($runMSA and (not ($runHMM or $runEmit or $runOutSwiss or $runOutEmit) )){
	if(not -d $seqDir){ die("can't find the seq directory for input"); }
	if(not -d $msaDir){ mkdir "$msaDir", 0755 or die("couldn't create hmm directory"); }

	@input = &findDirDiff($seqDir,$msaDir,"fasta","msa");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$seq = $input[$i];
			next if not $seq =~ /\.fasta$/;

			$msa = $seq;
			$msa =~ s/.fasta/.msa/;
			
			# run MSA	
			print FH "$muscle -in $seqDir/$seq -out $msaDir/$msa\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}
}
elsif($runHMM and (not ($runMSA or $runEmit or $runOutSwiss or $runOutEmit) )){
	if(not -d $msaDir){ die("can't find the msa directory for input"); }
	if(not -d $hmmDir){ mkdir "$hmmDir", 0755 or die("couldn't create msa directory"); }

	@input = &findDirDiff($msaDir,$hmmDir,"msa","hmm");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$msa = $input[$i];
			next if not $msa =~ /\.msa$/;

			$hmm = $msa;
			$hmm =~ s/.msa/.hmm/;

			# build HMM2 models
			# print FH "$hmm2build $hmmDir/$hmm $msaDir/$msa\n";
			# print FH "$hmm2calibrate $hmmDir/$hmm\n";
			
			# build HMM3 models
			print FH "$hmm3build --informat afa -o /dev/null $hmmDir/$hmm $msaDir/$msa\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}	
}
elsif($runEmit and (not ($runMSA or $runHMM or $runOutSwiss or $runOutEmit) )){
	if(not -d $hmmDir){ die("can't find the hmm directory for input"); }
	if(not -d $emitDir){ mkdir "$emitDir", 0755 or die("couldn't create emitDir directory"); }

	@input = &findDirDiff($hmmDir,$emitDir,"hmm","emit");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$hmm = $input[$i];
			next if not $hmm =~ /\.hmm$/;

			$emit = $hmm;
			$emit =~ s/.hmm/.emit/;
			
			# build HMM3 models
			print FH "$hmm3emit -N 1000 -o $emitDir/$emit $hmmDir/$hmm\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}	
}
elsif($runOutSwiss and (not ($runMSA or $runHMM or $runEmit or $runOutEmit) )){
	if(not -d $hmmDir){ die("can't find the hmm directory for input"); }
	if(not -e $swissprot ){ die("can't find the swissprot database for hmmsearch");}
	if(not -d $outswissDir){ mkdir "$outswissDir", 0755 or die("couldn't create outswissDir directory"); }

	@input = &findDirDiff($hmmDir,$outswissDir,"hmm","out");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$hmm = $input[$i];
			next if not $hmm =~ /\.hmm$/;

			$out = $hmm;
			$out =~ s/.hmm/.out/;
			
			# build HMM3 models
			print FH "$hmm3search --cpu 1 --tblout $outswissDir/$out --noali -o /dev/null $hmmDir/$hmm $swissprot\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}	
}
elsif($runOutEmit and (not ($runMSA or $runHMM or $runEmit or $runOutSwiss) )){
	if(not -d $hmmDir){ die("can't find the hmm directory for input"); }
	if(not -d $emitDir){ die("can't find the emit directory for input"); }
	if(not -d $outemitDir){ mkdir "$outemitDir", 0755 or die("couldn't create outemitDir directory"); }

	@input = &findDirDiff($hmmDir,$outemitDir,"hmm","out");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$hmm = $input[$i];
			next if not $hmm =~ /\.hmm$/;

			$out = $hmm;
			$out =~ s/.hmm/.out/;

			$emit = $hmm;
			$emit =~ s/.hmm/.emit/;
			if( -e "$emitDir/$emit" ){
				print FH "$hmm3search --cpu 1 --tblout $outemitDir/$out --noali -o /dev/null $hmmDir/$hmm $emitDir/$emit\n";
			}
			else{
				# if there is no emitted sequences, then move on to the next model
				print "Can't find the emitted database $emitDir/$emit!	Skipping model $hmm.\n" 
			}

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}	
}
elsif($runBlast and (not ($runMSA or $runHMM or $runEmit or $runOutSwiss or $runOutEmit) )){

	# Remember to format the database first:
	# makeblastdb -in swissProtein.fasta

	if(not -d $blastSeq){ die("can't find the seq directory for input"); }
	if(not -d $blastOut){ mkdir "$blastOut", 0755 or die("couldn't create blast directory"); }

	@input = &findDirDiff($blastSeq,$blastOut,"fasta","blast");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$seq = $input[$i];
			next if not $seq =~ /\.fasta$/;

			$blast = $seq;
			$blast =~ s/.fasta/.blast/;
			
			# run blast	
			print FH "$blastp -query $blastSeq/$seq -db $blastDatabase -outfmt '6 qseqid sseqid bitscore evalue length nident mismatch gaps qstart qend sstart send' -evalue 0.00001 -out $blastOut/$blast -num_descriptions 1000 -comp_based_stats 0 -use_sw_tback\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}	
}
elsif($runAnnotation and (not ($runMSA or $runHMM or $runEmit or $runOutSwiss or $runOutEmit or $runBlast) )){

	if(not -d $genomeFaa){ die("can't find the $genomeFaa directory for input"); }
	if(not -d $genomeOutput){ mkdir "$genomeOutput", 0755 or die("couldn't create $genomeOutput directory"); }
	
	@input = &findDirDiff($genomeFaa,$genomeOutput,"faa","scan");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$faa = $input[$i];
			next if not $faa =~ /\.faa$/;

			$scan = $faa;
			$scan =~ s/.faa/.scan/;
			
			# run blast	
			print FH "$hmm3scan --cpu 1 -o $genomeOutput/$scan $dbhmm $genomeFaa/$faa\n";

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}


}
elsif($runFilterScan and (not ($runMSA or $runHMM or $runEmit or $runOutSwiss or $runOutEmit or $runBlast or $runAnnotation) )){

	if(not -d $genomeOutput){ die("can't find the $genomeOutput directory for input"); }
	if(not -d $genomeFaa){ die("can't find the $genomeFaa directory for input"); }	
	if(not -d $annotation ){ mkdir "$annotation", 0755 or die("couldn't create $annotation directory"); }
	
	@input = &findDirDiff($genomeOutput,$annotation,"scan","ann");
	# calculate how many families every job should process
	$count = @input;
	$jobSize = int($count/$jobNumbers+1);
	print "Number of files to be processed: $count\n";

	$i = 0;
	while( $i < $count ) {

		$qsubTmp = "$qsubDir/qsub.$i";	
		# make a qsub script
		open FH, ">$qsubTmp";
		print FH "#PBS -m n\n";
	#	print FH "#PBS -l cput=24:00:00\n";

		for($j = 0; $j < $jobSize && $i < $count; $j++) {

			$scan = $input[$i];
			next if not $scan =~ /\.scan$/;

			$faa = $scan;
			$faa =~ s/.scan/.faa/;			
			
			$ann = $scan;
			$ann =~ s/.scan/.ann/;

			if( -e "$genomeFaa/$faa" ){
				# run annotateGenome	
				print FH "$annotateGenome -r -f $genomeFaa/$faa -o $genomeOutput/$scan -a $annotation/$ann\n";
			}
			else{
				# if there is no faa sequences, then move on to the next
				print "Can't find the faa $genomeFaa/$faa!	Skipping genome $faa.\n" 
			}
			

			$i++;
		}
		close FH;
		system "qsub $qsubTmp";
	}


}
else{
	die("Please indicate if you want to run MSA (use option --msa), run HMM build (use option --hmm), run HMM emit (use option --emit), run HMM search against SwissProt (use option --outswiss), or run HMM search against emitted sequences (use option --outemit)! Please provide one and only one of these options.");
}

# Argument: directory A; directory B; extension name for files in directory A; extension name for files in directory B
# find the files that are in directory A, but not in directory B
sub findDirDiff {
	$dirA = shift @_;
	$dirB = shift @_;
	$ExtensionA = shift @_;
	$ExtensionB = shift @_;

	opendir(DIRA, $dirA) or die "can't opendir $dirA: $!";
	opendir(DIRB, $dirB) or die "can't opendir $dirB: $!";

	%filesB = ();
	@filesDiff = ();

	while (defined($file = readdir(DIRB))) {
		next if $file =~ /^\.\.?$/;     # skip . and ..
		if( $file =~ /(\S+)\.$ExtensionB/ ){
			$baseNameB = $1;
			$filesB{$baseNameB} =  1;
		}
	}

	while (defined($file = readdir(DIRA))) {
		next if $file =~ /^\.\.?$/;     # skip . and ..
		if( $file =~ /(\S+)\.$ExtensionA/ ){
			$baseNameA = $1;
			if(not exists $filesB{$baseNameA}){
				push( @filesDiff, $file );
			}
		}
	}

	closedir(DIRA);
	closedir(DIRB);

	return @filesDiff;
}

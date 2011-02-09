
# input
$family = "familyTable.txt";
$hmmstat = "dbhmm.stat";
$emitDir = "emit";
$outswissDir = "outswiss";
$outemitDir = "outemit";


# output
$modelTable = "modelTable.txt";
$emitScoreTable = "emitScoreTable.txt";


# check if the input directories exist
if( not -d $emitDir){
	die("couldn't find $emitDir \n");
}
#if( not -d $outswissDir){
#	die("couldn't find $outswissDir \n");
#}
if( not -d $outemitDir){
	die("couldn't find $outemitDir \n");
}

#
# Parse "familyTable.txt" and "dbhmm.stat"
#
print "extracting familyTable\n";

# parsing familyTable
open FAMILY, "<$family" or die("couldn't open $family \n");
my %familyTable;
# skip the first line, which contains column headers
$row = <FAMILY>;
while ($row = <FAMILY>) {
        chomp $row;
        my @familyInfo = split(";\t", $row);
	$id = $familyInfo[0];
	$familyTable{$id} = [@familyInfo];
}

# parsing hmmstat
print "Parsing $hmmstat\n";
open HMMSTAT, "<$hmmstat" or die("couldn't open $hmmstat \n");
my %hmmstat_model_length;
my %hmmstat_eff_nseq;
while ($row = <HMMSTAT>) {
        chomp $row;
	next if($row =~ /^#/ );
	my @statInfo = split(/\s+/, $row);
	$id = $statInfo[1];
	$hmmstat_model_length{$id} = $statInfo[5];
	$hmmstat_eff_nseq{$id} = $statInfo[4];
}
close FAMILY;
close HMMSTAT;

#
# Parse emitted sequences
#
print "Parsing $outemitDir\n";
my %emitseq_score_curve;
my %emitseq_score_min;
my %emitseq_score_ave;
my %emitseq_score_max;
opendir(DIREMITOUT, $outemitDir) or die "can't opendir $outemitDir: $!";
while (defined($file = readdir(DIREMITOUT))) {
	next if $file =~ /^\.\.?$/;     # skip . and ..
	if( $file =~ /(\S+)\.out/ ){
		$id = $1;
		$sum_score = 0;
		$max_score = 0;
		$min_score = 1000000;
		$score_curve = "";
		$count = 0;
		if( not open OUT, "<$outemitDir/$file"){
			print "couldn't open emitted sequence scores from $file \n";
			next;
		}
		# collect a score every 10 rows, corresponding to one percent
		$interval = 10;
		# count how many percentages have been counted
		$curve_counter = 1;
		while ($row = <OUT>) {
			chomp $row;
			next if($row =~ /^#/ );
			@outInfo = split(/\s+/,$row);
			$score = $outInfo[5];
			$sum_score += $score;
			$count++;
			if($score > $max_score){ $max_score = $score; }
			if($score < $min_score){ $min_score = $score; }
			if( $interval == 10 ){
				$interval = 0;
				$score_curve = $score_curve . "\t" . $score;
				$curve_counter++;
			}
			$interval++;
		}
		# many models' emitted sequences may have a bit score less than the reporting cutoff
		# if there are less than 101 percentage points, pad with $min_score
		# 100%, 99%, ..., 0%
		while($curve_counter<102){
			$score_curve = $score_curve . "\t" . $min_score;
			$curve_counter++;
		}
		if($count > 0 )
		{
			$emitseq_score_ave{$id} = int ($sum_score / $count + 0.5);
			$emitseq_score_min{$id} = int ($min_score + 0.5);
			$emitseq_score_max{$id} = int ($max_score + 0.5);
			$emitseq_score_curve{$id} = $score_curve;
		}
		close OUT;

	}
}


closedir(DIREMITOUT);



# 
# Write output
#

open MODELTABLE, ">$modelTable" or die("couldn't open $modelTable \n"); 
open EMITSCORE, ">$emitScoreTable" or die("couldn't open $emitScoreTable \n"); 

print "Printing results\n";
print MODELTABLE "FamilyID;\tFamilySize;\tFamilyRecNameFull;\tFamilyRecNameEC;\tFamilyGeneName;\tModelLength;\tEffSeqN;\tMinEmitScore;\tAveEmitScore;\tMaxEmitScore;\tFamilyMember;\tFamilyGO\n";
for  $protein (sort keys %familyTable) {
	@infoArray = @{$familyTable{$protein}};

	
	print MODELTABLE "$infoArray[0];\t$infoArray[1];\t$infoArray[2];\t$infoArray[3];\t$infoArray[4];\t";
	

	if( exists $hmmstat_model_length{$protein} ) {
		$length = $hmmstat_model_length{$protein};
		print MODELTABLE "$length;\t";
	}
	else{
		print MODELTABLE "#N/A;\t";
	}

	if( exists $hmmstat_eff_nseq{$protein} ) {
		$nseq = $hmmstat_eff_nseq{$protein};
		print MODELTABLE "$nseq;\t";
	}
	else{
		print MODELTABLE "#N/A;\t";
	}

	if( exists $emitseq_score_ave{$protein} ) {
		$min_score = $emitseq_score_min{$protein};
		$ave_score = $emitseq_score_ave{$protein};
		$max_score = $emitseq_score_max{$protein};
		print MODELTABLE "$min_score;\t$ave_score;\t$max_score;\t";
	}
	else{
		print MODELTABLE "#N/A;\t#N/A;\t#N/A;\t";
	}
	

	print MODELTABLE "$infoArray[5];\t$infoArray[6]\n";
	
	if( exists $emitseq_score_curve{$protein} ) {
		print EMITSCORE "$protein$emitseq_score_curve{$protein}\n";
	}
}

close MODELTABLE;
close EMITSCORE;









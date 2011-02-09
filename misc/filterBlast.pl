
# usage: -t swissTable.txt -f swissProtein.fasta -o familyTable.txt -d seq
$table = "swissTable.txt";
$blastDir = "blastout";
$blastFiltered = "swissBlast.txt";
$blastFilteredIndex = "swissBlast.index";

$minDomainSize = 100;
# disable filtering based on coverage
# $minDomainSize = 1000000;

# Parse the command line arguments.  
for($i = 0; $i < @ARGV-1; $i++) {
  if($ARGV[$i] =~ /^\-t/) { $table = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-b/) { $blastDir = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-o/) { $blastFiltered = $ARGV[$i+1]; $i++; next; }
  
}

open TABLE, "<$table" or die("couldn't open $table \n");

#
# Parse "swissTable.txt" and "swissProtein.fasta"
# Populate %swissTable
#
print "extracting swissTable and swissProtein\n";
my %swissTable;
# skip the first line, which contains column headers
$row = <TABLE>;
while ($row = <TABLE>) {
        chomp $row;
#        ($id, $OC0, $OC1, $RecNameFull, $RecNameEC, $GeneName, $FastaPosition, $proteinLength) = split(";\t", $_);
        my @swissInfo = split(";\t", $row);
	$id = $swissInfo[0];
	$length = $swissInfo[7];
	$swissTable{$id} = $length;
}

if(not -d $blastDir){ die("can't find the blastDir directory for input"); }
opendir(DIR, $blastDir) or die "can't opendir $blastDir: $!";
open BLAST, ">$blastFiltered" or die("couldn't open $blastFiltered \n");
# open BLASTINDEX, ">$blastFilteredIndex" or die("couldn't open $blastFilteredIndex \n");


%blastMatrix;
while (defined($file = readdir(DIR))) {
	next if $file =~ /^\.\.?$/;     # skip . and ..
	if( $file =~ /(\S+)\.blast/ ){
		open FILE, "<$blastDir/$file" or die("couldn't open $file \n");
		my %blastTable;
		print "working on $file ...\n";
		while ($row = <FILE>) {
			chomp $row;
			($qseqid, $sseqid, $bitscore, $evalue, $length, $nident, $mismatch, $gaps, $qstart, $qend, $sstart, $send) = split("\t", $row);			

			# skip self vs self alignments
			next if( $qseqid eq $sseqid );

			if ( exists $swissTable{$qseqid} ){
				$qseqLength = $swissTable{$qseqid};
			}
			else {
				print "can't find $qseqid in swissTable.txt\n";
			}
			if ( exists $swissTable{$sseqid} ){
				$sseqLength = $swissTable{$sseqid};
			}
			else {
				print "can't find $sseqid in swissTable.txt\n";
			}

			$maxLength = $qseqLength;
			if($qseqLength < $sseqLength){
				$maxLength = $sseqLength;
			}

			$PercentIdentity = $nident / $maxLength;
			$AlignmentCoverage = $length / $maxLength;


			# determine if there is an extra domain in either N or C termini of the two proteins
			$noExtraDomain = 0;
			if( $qstart < $minDomainSize and $sstart < $minDomainSize and ($qseqLength-$qend) < $minDomainSize and ($sseqLength - $send) < $minDomainSize ){
				$noExtraDomain = 1;
			}

			if($evalue > 0 ){
				$logEvalue = -log($evalue) / log(10);
				$logEvalue = int($logEvalue)
			}
			else{
				$logEvalue = 1000;
			}
			


		# filter based on PercentIdentity; cutoff = 0.6
			# if($PercentIdentity > 0.6){
		# filter based on alignment coverage and e evalue
			# if( $AlignmentCoverage > 0.8 and $logEvalue > 5 ){

		# filter based on any extra domain and e evalue
			if( $noExtraDomain and $logEvalue > 5 ){

			#	@seqPair = sort($qseqid, $sseqid);
			#	$seqPairKey = $seqPair[0]."*".$seqPair[1];

				# a pair of sequences may have multiple alignments
				if(exists $blastTable{$qseqid}{$sseqid}){
					if($blastTable{$qseqid}{$sseqid} < $logEvalue ){
						$blastTable{$qseqid}{$sseqid} = $logEvalue;
					}
				}
				else{
					$blastTable{$qseqid}{$sseqid} = $logEvalue;
				}
			}
		}
		for  $query (sort keys %blastTable) {
			#	$Position = tell(BLAST);
			#	print BLASTINDEX "$query#$Position\n";
			#	print BLAST "$query#";
			$hits = "";
			for $subject (sort keys %{$blastTable{$query}}){
				$hits .= "$subject=$blastTable{$query}{$subject};";
			}
			$blastMatrix{$query} = $hits;
		}
	}

}

# use the larger score between A->B and B->A
print "Re-calculating scores ...\n";
&resetBlastScore;

# print the results
print "Printing ...\n";
for  $query (sort keys %blastMatrix) {
	@hitList = split(/;/, $blastMatrix{$query});
	$newHits = "";
	foreach $hit ( @hitList ){
		($subject, $score) = split(/=/, $hit);
		# only show the match at the query that is less than the subject
		if( $query lt $subject ){
			$newHits .= "$subject=$score;";
		}
	}
	if( $newHits ne "" ){
		print BLAST "$query#$newHits\n";
	}
}

sub resetBlastScore{
	
	my $seq0 = shift;
	my $seq1 = shift;
	my $similarityScore0 = 0;
	my $similarityScore1 = 0;

	# score of seq0 -> seq1
	$line = $blastMatrix{$seq0};
	if($line =~ /$seq1=(\d+);/ ){
		$similarityScore0 = $1;
	}
	else{
		# seq0 didn't hit seq1
		# add seq1
		$blastMatrix{$seq0} .= "$seq1=0;";
	}

	# score of seq1 -> seq0
	$line = $blastMatrix{$seq1};	
	if($line =~ /$seq0=(\d+);/ ){
		$similarityScore1 = $1;
	}
	else{
		# seq1 didn't hit seq0
		# add seq0
		$blastMatrix{$seq1} .= "$seq0=0;";
	}

	# re-set the smaller one.
	if( $similarityScore0 < $similarityScore1 ){
		$blastMatrix{$seq0} =~ s/$seq1=(\d+);/$seq1=$similarityScore1;/;
	}
	else{
		$blastMatrix{$seq1} =~ s/$seq0=(\d+);/$seq0=$similarityScore0;/;
	}

}

sub getBlastScore{
	
	my $seq0 = shift;
	my $seq1 = shift;
	my $similarityScore0 = 0;
	my $similarityScore1 = 0;

	# score of seq0 -> seq1
	$line = $blastMatrix{$seq0};
	if($line =~ /$seq1=(\d+);/ ){
		$similarityScore0 = $1;
	}

	# score of seq1 -> seq0
	$line = $blastMatrix{$seq1};	
	if($line =~ /$seq0=(\d+);/ ){
		$similarityScore1 = $1;
	}

	# return the bigger one.
	if( $similarityScore0 > $similarityScore1 ){
		return $similarityScore0;
	}
	else{
		return $similarityScore1;
	}

}

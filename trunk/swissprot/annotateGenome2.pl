#!/usr/bin/env perl

$hmmscanBin = "/data1/pcl/swissprot/bin/hmmer3/hmmscan";
$fasta = "/data1/pcl/swissprot/rpal/protein_gene_revised.faa";
$dbhmm = "/data1/pcl/swissprot/dbhmm/swiss.hmm";
$hmmscanOutput = "/data1/pcl/swissprot/rpal/hmmscan_output.txt";
$hmmscanTable = "/data1/pcl/swissprot/rpal/hmmscan_table.txt";
$modelFile = "/data1/pcl/swissprot/modelTable.txt";
$annotation = "/data1/pcl/swissprot/rpal/swiss_annotation.txt";
$emitScoreTable = "/data1/pcl/swissprot/emitScoreTable.txt";

$runHmmscan = 1;

# Parse the command line arguments.  
# usage: -f fasta_genome -d hmm_database -m model_table -a annotation_output
for($i = 0; $i < @ARGV; $i++) {
  if($ARGV[$i] =~ /^\-f/) { $fasta = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-d/) { $dbhmm = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-m/) { $modelFile = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-a/) { $annotation = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-o/) { $hmmscanOutput = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-t/) { $hmmscanTable = $ARGV[$i+1]; $i++; next; } 
  if($ARGV[$i] =~ /^\-e/) { $emitScoreTable = $ARGV[$i+1]; $i++; next; } 
  if($ARGV[$i] =~ /^\-r/) { $runHmmscan = 0;next; } # don't run Hmmscan; use the previous results
}

open MODEL, "<$modelFile" or die("couldn't open $modelFile \n");
my %modelTable;
my %modelLength;
# skip the first line, which contains column headers
$row = <MODEL>;
while ($row = <MODEL>) {
        chomp $row;
        my @modelInfo = split(";\t", $row);
	$id = $modelInfo[0];
	# FamilyRecNameFull; FamilyRecNameEC; FamilyGeneName; FamilyGO; MinEmitScore"
	$modelTable{$id} = [($modelInfo[2],$modelInfo[3],$modelInfo[4],$modelInfo[11],$modelInfo[7])];
	$modelLengthTable{$id} = $modelInfo[5];
}
close MODEL;

open EMITSCORE,  "<$emitScoreTable" or die("couldn't open $emitScoreTable \n");
my %emitScoreTable;
while ($row = <EMITSCORE>) {
        chomp $row;
        my @rowInfo = split("\t", $row);
	$id = $rowInfo[0];
	shift(@rowInfo);
	$emitScoreTable{$id} = [(@rowInfo)];
}
close EMITSCORE;


# do hmmscan

if($runHmmscan){
#	print "$hmmscanBin --tbl $hmmscanTable -o $hmmscanOutput $dbhmm $fasta\n";
#	system("$hmmscanBin --tbl $hmmscanTable -o $hmmscanOutput $dbhmm $fasta");
	print "$hmmscanBin -o $hmmscanOutput $dbhmm $fasta\n";
	system("$hmmscanBin -o $hmmscanOutput $dbhmm $fasta");
}


print "parsing hmmscan output\n";
my %HitModel;
my %HitEvalue;
my %HitScore;
my %HitHmmCoverage;
my %HitAliCoverage;
my %QueryLength;

# 
# parse the output from the Tabular output
#

#open HMMTABLE, "<$hmmscanTable" or die("couldn't open $hmmscanTable \n");
#while ($row = <HMMTABLE>) {
#        chomp $row;
#	next if($row =~ /^#/ );
#	my @rowInfo = split(/\s+/, $row);
#	$model = $rowInfo[0];
#	$query = $rowInfo[2];
#	$Evalue = $rowInfo[4];
#	$score = $rowInfo[5];
#	if( exists $HitModel{$query} ){
#		if( $HitScore{$query} < $score ){
#			$HitModel{$query} = $model;
#			$HitEvalue{$query} = $Evalue;
#			$HitScore{$query} = $score;
#		}
#	}
#	else {
#		$HitModel{$query} = $model;
#		$HitEvalue{$query} = $Evalue;
#		$HitScore{$query} = $score;
#	}
#}
#close HMMTABLE;

# parse the full output
open HMMOUT, "<$hmmscanOutput" or die("couldn't open $hmmscanOutput \n");
while ($row = <HMMOUT>) {
        chomp $row;
	next if($row =~ /^#/ );

	if($row =~ /^Query:\s+(\S+)\s+\[L=(\d+)\]/ ){
		$query = $1;
		$queryLength = $2;
		if(exists $HitModel{$query}){
			print "redundant scan results for $query\n";
			next;
		}
		$QueryLength{$query} = $queryLength;
		$HitScore{$query} = 0;

		my %modelScore = ();
		my %modelEvalue = ();
		my %modelHmmCoverage = ();
		my %modelAliCoverage = ();

		# get the full-sequence e-value and bit score
		while($row = <HMMOUT>){
			# move to this row
			if($row =~ /E-value  score  bias    E-value  score  bias    exp  N  Model    Description/){
				last;
			}
		}
		# skip a line
		$row = <HMMOUT>;
		# read the list up to "inclusion threshold" or to the end
		while($row = <HMMOUT>){
			chomp $row;
			# parse this row
#    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
#    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
#   4.1e-147  490.4  16.1   6.1e-131  437.0   6.1    2.0  2  FIXG    
			my @rowInfo = split(/\s+/, $row);
			if(scalar @rowInfo >= 9){
				$modelScore{$rowInfo[9]}= $rowInfo[2];
				$modelEvalue{$rowInfo[9]}= $rowInfo[1];

			}
			else{
				last;
			}
		}
		
		# loop thru all models' alignments
		while($row = <HMMOUT>){
			chomp $row;
			next if($row =~ /^#/ );
			if($row =~ /^>>\s+(\S+)/ ){
				$model = $1;
				$bestDomainScore = 0;
				$bestDomainEvalue = 1;
				@hmmfrom = ();
				@hmmto = ();
				@alifrom = ();
				@alito = ();
				# skill two rows
				$row = <HMMOUT>;
				$row = <HMMOUT>;
				# loop thru all domains
				while($row = <HMMOUT>){
					# parse the row
##       score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
# ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#   1 !   33.9   0.1   6.7e-12     5e-09      31     221 ..      65     270 ..      51     285 .. 0.80	

					if($row =~ /\s+\d+\s+[!?]\s+(\S+)\s+/ ){
						my @rowInfo = split(/\s+/, $row);
						push(@hmmfrom, $rowInfo[7]);
						push(@hmmto, $rowInfo[8]);
						push(@alifrom, $rowInfo[10]);
						push(@alito, $rowInfo[11]);

						# save the best scoring domain info 
						if($bestDomainScore < $rowInfo[3] ){
							$bestDomainScore = $rowInfo[3];
							$bestDomainEvalue = $rowInfo[6];
							$bestDomainHmmCoverage = $rowInfo[8] - $rowInfo[7] + 1;
							$bestDomainAliCoverage = $rowInfo[11] - $rowInfo[10] + 1;
						}
					}
					else{
						last;
					}

				}

				if(exists($modelScore{$model})){
					if(&isAlignmentOverlapped(@hmmfrom, @hmmto) or &isAlignmentOverlapped(@alifrom, @alito)){
						# use the best domain score/evalue/coverage
						$modelScore{$model} = $bestDomainScore;
						$modelEvalue{$model} = $bestDomainEvalue;
						$modelHmmCoverage{$model} = $bestDomainHmmCoverage;
						$modelAliCoverage{$model} = $bestDomainAliCoverage;
					}
					else{
						# use the full-sequence score/evalue/coverage
						$modelHmmCoverage{$model} = &computeCoverage(\@hmmfrom, \@hmmto);
						$modelAliCoverage{$model} = &computeCoverage(\@alifrom, \@alito);
					}
				}
			}
			elsif($row =~ /^\/\//){
				last;
			}
		}
		# complete current query sequence;
		for $model (keys %modelScore) {
			if( $modelScore{$model} > $HitScore{$query} ){

				$HitModel{$query} = $model;
				$HitScore{$query} = $modelScore{$model};
				$HitEvalue{$query} = $modelEvalue{$model};
				$HitHmmCoverage{$query} = $modelHmmCoverage{$model};
				$HitAliCoverage{$query} = $modelAliCoverage{$model};

			}

		}

	}
}
close HMMOUT;



# calculate the bit score percentile
# the scores should correspond to 101 percentages points: 100%, 99%, ..., 0%.
my %HitEmitScorePercentile;
for  $protein (sort keys %HitScore) {
	$score = $HitScore{$protein};
	$model = $HitModel{$protein};
	if( exists $emitScoreTable{$model} ){
		@scorePercentile = @{$emitScoreTable{$model}};
		$index = 0;
		# test up to 1%
		$PercentileNumber = 100;
		while( $index < $PercentileNumber ){
			if( $scorePercentile[$index] < $score ){
				last;
			}
			$index++;
		}
		$HitEmitScorePercentile{$protein} = ($PercentileNumber-$index)/$PercentileNumber;
		# if the score is less than 1% and greater than minimum score, assign 0.1 to the percentile
		if($index == $PercentileNumber && $score > $scorePercentile[100]){
			$HitEmitScorePercentile{$protein} = 0.1
		}
	}
}


# get the list of proteins in the FASTA file
open FASTA, "<$fasta" or die("couldn't open $fasta \n");
@ProteinList = ();
%Protein2GenbankDescription = ();
%Protein2Species = ();
while ($row = <FASTA>){
	chomp $row;
	if($row =~ /^>/ ){
		# example:
		# >gi|170016304|ref|YP_001727225.1| Cu2+-Cu+-Ag+-P-type ATPase [Leuconostoc citreum KM20]

		$row =~ s/^>//;
		my @rowInfo = split(/\s+/, $row);
		$proteinName = $rowInfo[0];
		push(@ProteinList, $proteinName);

		$desc = "NA";
		$species = "NA";
		$row =~ s/^\Q$proteinName\E\s+//;
		if( $row =~ /(.+)\s+\[(.+)\]/ ){
			$desc = $1;
			$species = $2;
		}
		$Protein2GenbankDescription{$proteinName} = $desc;
		$Protein2Species{$proteinName} = $species;
	}
}

close FASTA;

# write results
open ANNOT, ">$annotation" or die("couldn't open $annotation \n");
print ANNOT "Species;\tProtein;\tHit;\tEvalue;\tScore;\tDescription;\tEC;\tGeneName;\tGO;\tMinEmitScore;\tEmitScorePercentile;\tPercentModelCoverage;\tPercentSeqCoverage;\tGenbankAnnotation\n";

foreach $protein (@ProteinList) {

	$model = "NA";
	$Evalue = "NA";
	$Score = "NA";
	$modelInfo = "NA;\tNA;\tNA;\tNA;\tNA";	

	$ScorePercentile = "NA";
	$PctModelCoverage = "NA";
	$PctSeqCoverage = "NA";
	if( exists $HitModel{$protein} ){
		$model = $HitModel{$protein};
		$Evalue = $HitEvalue{$protein};
		$Score = $HitScore{$protein};
		$hmmCoverage = $HitHmmCoverage{$protein};
		$aliCoverage = $HitAliCoverage{$protein};
		$queryLength = $QueryLength{$protein};
		$modelLength = 0;
 
		if( exists $HitEmitScorePercentile{$protein} ){
			$ScorePercentile = $HitEmitScorePercentile{$protein};
		}
		if(exists $modelTable{$model} ){
			@array = @{$modelTable{$model}};
			$modelInfo = join(";\t", @array);
			$modelLength = $modelLengthTable{$model};
		}

		if($queryLength != 0){
			$PctSeqCoverage = $aliCoverage / $queryLength;
			$PctSeqCoverage = sprintf("%.2f", $PctSeqCoverage);
			
		}

		if($modelLength != 0){
			$PctModelCoverage = $hmmCoverage / $modelLength;
			$PctModelCoverage = sprintf("%.2f", $PctModelCoverage);
		}

	}
	$species = $Protein2Species{$protein};
	$GenbankDesc = $Protein2GenbankDescription{$protein};
	print ANNOT "$species;\t$protein;\t$model;\t$Evalue;\t$Score;\t$modelInfo;\t$ScorePercentile;\t$PctModelCoverage;\t$PctSeqCoverage;\t$GenbankDesc\n";
}


close ANNOT;

sub isAlignmentOverlapped(){
	my @from = shift;
	my @to = shift;

	if( scalar @from != scalar @to ){
		print "Error: invalid alignment\n";
		return 0; # false
	}

	if( scalar @from == 1 ){
		return 0; # false
	}

	$maxOverlapAllowed = 10;

	for ($i = 0; $i < (scalar @from) - 1; $i++) {
		for( $j = $i + 1; $j < (scalar @from); $j++ ){
			if( $to[$i] - $from[$j] > $maxOverlapAllowed or $to[$j] - $from[$i] > $maxOverlapAllowed ){
				return 1; # true
			}
		}
	}
	return 0; # false
}

sub computeCoverage(){
	my @from = @{ $_[0] };
	my @to   = @{ $_[1] };

	if( scalar @from != scalar @to ){
		print "Error: invalid alignment\n";
		return 0;
	}

	my %coveredSequence = ();
	for ($i = 0; $i < (scalar @from) ; $i++) {
		for( $j = $from[$i]; $j <= $to[$i]; $j++ ){
			$coveredSequence{$j} = 0;
		}
	}

	return scalar keys %coveredSequence;
}


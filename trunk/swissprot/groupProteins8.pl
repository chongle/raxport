
# usage: -t swissTable.txt -f swissProtein.fasta -o familyTable.txt -d seq
$table = "swissTable.txt";
$fasta = "swissProtein.fasta";
$blastFile = "swissBlast.txt";
$family = "familyTable.txt";
$profamily = "profamilyTable.txt";
$familyRecNameFullFile = "familyNames.txt";
$seqDir = "seq";
$pass = 0;

$BlastCutoff1 = 5;
$BlastCutoff2 = 10;

# Parse the command line arguments.  
for($i = 0; $i < @ARGV-1; $i++) {
  if($ARGV[$i] =~ /^\-t/) { $table = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-f/) { $blastFile = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-o/) { $family = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-d/) { $seqDir = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-x/) { $pass = $ARGV[$i+1]; $i++; next; }
}

if( $pass != 1 and $pass != 2 ){
	die("Please use either option -x 1 for the first pass and option -x 2 for the second pass\n");
}

#
# Parse "swissTable.txt" and populate %swissTable
# Parse "swissBlast.txt" and populate %blastMatrix 
#
%swissTable;
%blastMatrix;
&populateInput;

# This is the information for a family: $familyInfo{$protein} = [ $FamilyRecNameFull, $FamilyRecNameEC, $FamilyGeneName ]; 
%familyInfo;
# This is the list of proteins in a family: protein -> an array of IDs;
%familyTable;

# group proteins and populate %familyInfo and %familyTable
if($pass == 1){
	&groupProteins1;
}
elsif($pass == 2){
	&groupProteins2;
}


#
# Print familyTable.txt and proFamilyTable.txt based on %familyInfo
# For $pass ==2, Print a FASTA file for every family in the seqDir directory
#
&printResults;



#
# filter proteins before grouping
#
sub filterProteins {

# The parameter array @swissInfo include $protein, $OC0, $OC1, $RecNameFull, $RecNameEC, $GeneName; $FastaPosition; $proteinLength

	my $myProtein = shift;
	my $OC0 = shift;
	my $OC1 = shift;
	my $myRecNameFull = shift;
	my $myRecNameEC = shift;
	my $myGeneName = shift;
	my $myFastaPosition = shift;
	my $myProteinLength = shift;

	if($myProteinLength < 30 ){
		return 0;
	}

	# remove viruses, animals and plants
#	if($OC0 eq "Viruses" or $OC1 eq "Metazoa" or $OC1 eq "Viridiplantae" ){
#		return 0;
#	}
	# remove "uncharacterized" proteins
	if($myRecNameFull =~ /[Uu]ncharacterized/ ){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}
	# remove "UPF" proteins
	if($myRecNameFull =~ /UPF\d+\s/ ){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	# remove "X kDa proteins/antigens"
	if($myRecNameFull =~ /[0-9]*\.?[0-9]+\skDa/ && $myRecNameFull =~ /protein|antigen/ ){
		if( $myRecNameEC eq "NA" && $myGeneName eq "NA" ){
	#		print "Discarding $myRecNameFull\n";
			return 0;
		}
	}

	# remove "HLA class I histocompatibility antigen "
	if($myRecNameFull =~ /histocompatibility antigen,/ ){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	# remove "domain/repeat/motif-containing protein", because these proteins only share a domain 
	if($myRecNameFull =~ /-containing protein/){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	# remove zinc finger
	if($myRecNameFull =~ /[zZ]inc finger/ or $myRecNameFull =~ /F-box/  or $myRecNameFull =~ /[aA]nkyrin repeat/ or $myRecNameFull =~ /Homeobox/ ){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	# remove starting with "Transmembrane protein" 
	if($myRecNameFull =~ /^Transmembrane protein/ ){
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	# matches ~300 proteins
	if($myRecNameFull =~ /^Maf-like protein/ ){ 
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	if($myRecNameFull =~ /^Ig / ){ 
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	if($myRecNameFull =~ /^Gene\s+.*\s+protein$/ ){ 
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	if($myRecNameFull =~ /^Protein FAM\S+$/ ){ 
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	if($myRecNameFull =~ /^Olfactory receptor/ ){ 
	#	print "Discarding $myRecNameFull\n";
		return 0;
	}

	

	return 1;

}

# sub filterFamilies {
#	($arraySize, $FamilyRecNameFull, $FamilyRecNameEC, $FamilyGeneName, $idString) = @_;
#	if($arraySize == 1 && $FamilyRecNameFull =~ /family/ ) {return 0;}
#	if($arraySize == 1 && $idString =~ /_ARATH/) {return 0;}	
#	if($arraySize == 1 && $FamilyRecNameFull =~ /^Metalloprotease /) {return 0;}	
#	return 1;	
#}

sub findMostCommonItem {
	my @List = @_;
	# count how many times an item occurs
	my %ItemCount;
	foreach $item (@List) {
		if ( exists $ItemCount{$item} ){
			$ItemCount{$item}++;
		}
		else{
			$ItemCount{$item}=1;
		}
	}
	# find MostCommonItem that is not "NA"
	$MostCommonItem = "NA";
	$MaxCount = 0;

	while (($item, $count) = each(%ItemCount)){
		if( $count > $MaxCount && $item ne "NA" ){
			$MaxCount = $count;
			$MostCommonItem = $item;
		}
	}
	return $MostCommonItem;
}

sub isCompatible{
	($name0, $name1) = @_;

	# change both to lower case
	$name0 = lc($name0);
	$name1 = lc($name1);

	# Remove uninformative words, including protein, proteins, homolog, isoform, putative, hypothetical, probable, uncharacterized, unknown.
#	$name0 =~ s{(protein|proteins|complex|homolog|isoform|like|related)}{}g;
#	$name1 =~ s{(protein|proteins|complex|homolog|isoform|like|related)}{}g;

	# replace every non-alphanumeric character with a space
#	$name0 =~ s{[^a-z0-9]}{ }g;
#	$name1 =~ s{[^a-z0-9]}{ }g;

	# replace symbols, including -, \, /, with a space
#	$name0 =~ s{(/|-|\\)}{ }g;
#	$name1 =~ s{(/|-|\\)}{ }g;
	
	# Remove all white spaces
#	$name0 =~ s{\s+}{}g;
#	$name1 =~ s{\s+}{}g;

	if( $name0 eq $name1) {
		return 1;
	}
	else {
		return 0;
	}

}

sub formatRecNameFull{

	($name) = @_;

	$name =~ s{, (chloroplastic|chromoplast|mitochondrial|peroxisomal|cytosolic|cytoplasmic|organellar chromatophore|apicoplast|cyanelle|plastid|plasmid|chromosomal|amyloplastic|glycosomal|endoplasmic reticulum)/(chloroplastic|chromoplast|mitochondrial|peroxisomal|cytosolic|cytoplasmic|organellar chromatophore|apicoplast|cyanelle|plastid|plasmid|chromosomal|amyloplastic|glycosomal|endoplasmic reticulum)}{}g;
	$name =~ s{, (chloroplastic|chromoplast|mitochondrial|peroxisomal|cytosolic|cytoplasmic|organellar chromatophore|apicoplast|cyanelle|plastid|plasmid|chromosomal|amyloplastic|glycosomal|endoplasmic reticulum)(|\s+isoform|\s+isozyme)}{}g;

	if ( $name =~ s{^(Putative|Probable|Hypothetical)\s+}{} ) {
		$name =~ s/^([a-z])/uc($1)/e;
	}

	# the previous change will change "Probable tRNA" to "TRNA", which is corrected here.
	$name =~ s/^TRNA/tRNA/;

	# remove "from XXXX";
	# example: change "DNA-invertase from lambdoid prophage e14" to "DNA-invertase" 
	$name =~ s/\s+from.+$//;

	# remove "XXX homolog"
	$name =~ s/\s+homolog(\s+\d+|)$//;

	return $name;

}
# split family based on just Blast similarity
sub splitFamilty{
	my $currentProtein = shift;
	print "$currentProtein\n";
	my @idList = @{$familyTable{$currentProtein}};
	my $arraySize = @idList;

#	my %IdName = ();
#	foreach $id (@idList) {
#		@swissInfo = @{$swissTable{$id}};
#		$IdName{$id} = $swissInfo[3]; # RecNameFull
#	}

	# build adjacency matrix
	my %adjacency;
	for ($i = 0; $i < $arraySize; $i++){
		my @neighbors = ();
		$adjacency{$idList[$i]} = [@neighbors]
	}

	for ($i = 0; $i < $arraySize - 1; $i++){
		for( $j = $i + 1; $j < $arraySize; $j++ ){
			
		#	$nameCompatibility = &isCompatible($IdName{$idList[$i]}, $IdName{$idList[$j]});

			# if the names are compatible
		#	if( $nameCompatibility == 1 ){
				# get the similarity score
				my $similarityScore = getBlastScore($idList[$i],$idList[$j]);
				if( $similarityScore > $BlastCutoff1 ){
					push(@{$adjacency{$idList[$i]}}, $idList[$j]);
					push(@{$adjacency{$idList[$j]}}, $idList[$i]);
				}
		#	}

		}	
	}
	my @componentList = &breathFirstSearch(\%adjacency);
	return @componentList;
}


sub mergeRecNameFull{
	# a hash of ID => Name
	my (%IdName) = @_;
	# a hash of Name => an array of IDs
	my %MergedNameID;
	# a hash of Name => count of IDs
	my %MergedNameCount;
	for  $Id (sort keys %IdName) {
		$Name = $IdName{$Id};
		if ( exists $MergedNameID{$Name} ){
			# Appending a new element to one in the arrays
			push(@{$MergedNameID{$Name}}, $Id);	
			$MergedNameCount{$Name}++;
		}
		else{
			# Adding an array to the hash
			@newArray = ($Id);
			$MergedNameID{$Name} = [@newArray];
			$MergedNameCount{$Name} = 1;
		}
	}

	my $size = scalar(keys %MergedNameID);

	if($size == 1)
	{
		return %MergedNameID;
	}

	# sort the keys of the hash by the number of proteins they contains
	my @SortedName = sort { $MergedNameCount{$b} <=> $MergedNameCount{$a} } keys %MergedNameCount;

	my $i;
	my $j;
	for ($i = 0; $i < $size; $i++){
		for($j = $i+1; $j < $size; $j++) {
			# $name0 has more proteins than $name1
			my $name0 = $SortedName[$i];
			my $name1 = $SortedName[$j];
			# if these two names are not deleted in previous iterations
			if( exists $MergedNameID{$name0} && exists $MergedNameID{$name1} ){
				# test if these two names are similar enough
				if( isCompatible($name0,$name1) ){
					# if yes, merge the protein list of name1 to the protein list of name0
					push(@{$MergedNameID{$name0}}, @{$MergedNameID{$name1}});
					# and delete name1 from the hash
					delete($MergedNameID{$name1});	
				}
			}
		}
	}

	return %MergedNameID
}

sub breathFirstSearch{
    	
	my $params = shift;

    	my %adjacency = %$params;

	my @componentList;
	# Clear visited hash
	foreach(keys %adjacency) { $visited{$_} = 0; }

	foreach(keys %adjacency) {
	  next if($visited{$_} == 1);

	  # Add next element to the queue
	  push @queue, $_;
	  $component = $_;

	  # Mark element as visited
	  $visited{$_} = 1;

	  # While the queue isn't empty, add all unvisited
	  # nodes to queue and mark them as visited.
	  while(@queue > 0) {
	    $cur = shift @queue;
	    @inf = @{$adjacency{$cur}};
	    for($i = 0; $i < @inf; $i++) {
	      next if($inf[$i] eq $cur);
	      next if($visited{$inf[$i]} == 1);
	      $visited{$inf[$i]} = 1;
	      $component = $component . "#" . $inf[$i];
	      push @queue, $inf[$i];
	    }
	  }
	  push @componentList, $component;
#	  print "$component\n";
	}

	return @componentList;
}


# retrieve the blast score between two proteins
#
sub getBlastScore{
	my $seq0 = shift;
	my $seq1 = shift;
	my $similarityScore = 0;
	if($seq0 lt $seq1){
		# score of seq0 -> seq1	
		$line = $blastMatrix{$seq0};
		if($line =~ /$seq1=(\d+);/ ){
			$similarityScore = $1;
		}		
	}
	else{
		# score of seq1 -> seq0
		$line = $blastMatrix{$seq1};	
		if($line =~ /$seq0=(\d+);/ ){
		       $similarityScore = $1;
		}
	}
	return $similarityScore;
}

# populate input swissTable and blastMatrix
sub populateInput{

	# Parse "swissTable.txt" and populate %swissTable
	print "extracting swissTable and swissProtein\n";
	open TABLE, "<$table" or die("couldn't open $table \n");
	# skip the first line, which contains column headers
	$row = <TABLE>;
	while ($row = <TABLE>) {
		chomp $row;
	#       $id, $OC0, $OC1, $RecNameFull, $RecNameEC, $GeneName, $FastaPosition, $proteinLength, $GO
		my @swissInfo = split(";\t", $row);
		$id = $swissInfo[0];
		($protein, $species) = split("_", $id);
		$swissInfo[0] = $protein;
		# format the RecNameFull
		$swissInfo[3] = formatRecNameFull( $swissInfo[3] );
		if(filterProteins( @swissInfo )){
			# a hash of arrays
			# The key of the hash is the id
			# The array include $protein, $OC0, $OC1, $RecNameFull, $RecNameEC, $GeneName; $FastaPosition; $proteinLength; $GO
			$swissTable{$id} = [@swissInfo];
		}
	}
	close TABLE;

	# 
	# parse "swissBlast.txt" and populate %blastMatrix
	#
	print "parse swissBlast.txt and populate %blastMatrix\n";
	open BLAST, "<$blastFile" or die("couldn't open $blastFile \n");
	while ($row = <BLAST>) {
		chomp $row;
		($seq, $hits) = split("#", $row);
		$blastMatrix{$seq} = $hits;
	}
	close BLAST;


}

# print results:
sub printResults{
	
	print "Printing results\n";
	open FAMILY, ">$family" or die("couldn't open $family \n");
	open PROFAMILY, ">$profamily" or die("couldn't open $profamily \n"); 

	if($pass == 2){
		open FASTA, "<$fasta" or die("couldn't open $fasta \n");

		# if $seqDir doesn't exist, creat it.
		# if it does, delete all files under it
		if(-d $seqDir){
			system("rm -rf $seqDir");
		}
		mkdir "$seqDir", 0755 or die("couldn't create seq directory");
	}


	print FAMILY "FamilyID;\tFamilySize;\tFamilyRecNameFull;\tFamilyRecNameEC;\tFamilyGeneName;\tFamilyMember;\tFamilyGO\n";
	for  $protein (sort keys %familyTable) {
		@idArray = @{$familyTable{$protein}};
		
		$arraySize = @idArray;
		@infoArray = @{$familyInfo{$protein}};
		$FamilyRecNameFull = $infoArray[0];
		$FamilyRecNameEC = $infoArray[1];
		$FamilyGeneName = $infoArray[2];
		$FamilyGOterm = $infoArray[3];
		$idString = join(",", @idArray);

		# Print family info		
		print FAMILY "$protein;\t$arraySize;\t$FamilyRecNameFull;\t$FamilyRecNameEC;\t$FamilyGeneName;\t($idString);\t$FamilyGOterm\n";
		print PROFAMILY "$protein;\t$arraySize;\t$FamilyRecNameFull;\t$FamilyRecNameEC;\t$FamilyGeneName;\t$FamilyGOterm\n";


		foreach $id (@idArray) {
			@swissInfo = @{$swissTable{$id}};
			print PROFAMILY "$id;\t##;\t$swissInfo[3];\t$swissInfo[4];\t$swissInfo[5];\t$swissInfo[8]\n";
		}

		if($pass == 2 ){
			# Creat family FASTA
			open PROTEINFASTA, ">./$seqDir/$protein.fasta" or die("couldn't open $protein.fasta \n");
			foreach $id (@idArray) {
				@swissInfo = @{$swissTable{$id}};			
				$position = $swissInfo[6];
				seek(FASTA, $position, 0);
				$line = <FASTA>;
				chomp $line;
				if($line eq ">$id"){
					$sequence = <FASTA>;
					chomp $sequence;
				}
				else{
					print "Couldn't find the sequence for protein $id .\n";
				}
				print PROTEINFASTA ">$id\n";
				print PROTEINFASTA "$sequence\n";
			}
			close PROTEINFASTA;
		}
		
	}

	if($pass == 2){
		close FASTA;
	}

	close FAMILY;
	close PROFAMILY;

}

sub groupProteins1{
	#
	# Group proteins by their name
	# Populate %familyTable
	#
	print "Making familyTable\n";
	for  $id (keys %swissTable) {
		@swissInfo = @{$swissTable{$id}};
		$protein = $swissInfo[0];
		if ( exists $familyTable{$protein} ){
			# Appending a new element to one in the arrays
			push(@{$familyTable{$protein}}, $id);	
		}
		else{
			# Adding an array to the hash
			@newArray = ($id);
			$familyTable{$protein} = [@newArray];
		}
	}

	# 
	# Check the consistency of RecNameFull in a family
	# Remove inconsistant member to a separate family
	#
	print "Checking familyTable\n";
	for  $protein (sort keys %familyTable) {
		@idArray = @{$familyTable{$protein}};
		my %IdName = ();
		foreach $id (@idArray) {
			@swissInfo = @{$swissTable{$id}};
			$IdName{$id} = $swissInfo[3]; # RecNameFull
		}

		my %MergedNameID = mergeRecNameFull(%IdName);
		my $size = scalar(keys %MergedNameID);
		if($size > 1 ){
			delete($familyTable{$protein});
			my $i = 0;
			for $name (keys %MergedNameID) {
				$familyTable{$protein."[".$i."]"} = $MergedNameID{$name};
				$i++;
			}
		}
	}

	# 
	# merge families
	#

	print "Merging families\n";
	&mergeFamilies;


	#
	# Find the most common names in this family
	# Populate %familyInfo
	#
	print "Extracting familyInfo\n";
#	my %familyInfo;
	for  $protein (sort keys %familyTable) {
		@idArray = @{$familyTable{$protein}};
		$arraySize = @idArray;
		@RecNameFull = ();
		@RecNameEC = ();
		@GeneName = ();
		@GOarray = ();
		foreach $id (@idArray) {
			@swissInfo = @{$swissTable{$id}};
			push(@RecNameFull, $swissInfo[3]);
			push(@RecNameEC, $swissInfo[4]);
			push(@GeneName, $swissInfo[5]);
			push(@GOarray, $swissInfo[8]);
		}

		# find family info
		$FamilyRecNameFull = &findMostCommonItem(@RecNameFull);
		$FamilyRecNameEC = &findMostCommonItem(@RecNameEC);
		$FamilyGeneName = &findMostCommonItem(@GeneName);
		$FamilyGOterm = &findConsensusGO(@GOarray);


		# modify full rec name
		foreach $name (@RecNameFull){
			if( not isCompatible( $FamilyRecNameFull, $name ) ){
				# this is a merged protein
				$compatibity2 = isCompatible2( $FamilyRecNameFull, $name);
				if( $compatibity2 ne "false" and $compatibity2 ne $FamilyRecNameFull ){
					$FamilyRecNameFull = $compatibity2;
					last;
				}
			}
		}

		$familyInfo{$protein} = [ $FamilyRecNameFull, $FamilyRecNameEC, $FamilyGeneName, $FamilyGOterm ];
	}

}

sub groupProteins2{

	# read in the list of curated protein names
	# and populate the keys of $familytable
	print "Reading $familyRecNameFullFile\n";
	open NAME, "<$familyRecNameFullFile" or die("couldn't open $familyRecNameFullFile \n");
	while ($row = <NAME>) {
		chomp $row;
		# make sure the protein names are all unique
		if( not exists $familyTable{$row} ){
			@newArray = ();	
			# use the protein name as the key for $familyTable	
			$familyTable{$row} = [@newArray];
		}	
	}
	close NAME;
	

	print "building familtyTable\n";
	for  $id (sort keys %swissTable) {
		@swissInfo = @{$swissTable{$id}};
		$currentFullName = $swissInfo[3];

		if(exists $familyTable{$currentFullName} ){
			push(@{$familyTable{$currentFullName}}, $id);
		}
		else{
			$matchScore = 0;
			$matchedKey = "";

			for $fullNameKey (keys %familyTable){
				if($currentFullName =~ /^\Q$fullNameKey\E/i){
					if($matchScore < length($fullNameKey)){
						$matchScore = length($fullNameKey);
						$matchedKey = $fullNameKey;
					}
					
				}
			}
			if ( $matchScore > 0 ){
				# Appending a new element to one in the arrays
				push(@{$familyTable{$matchedKey}}, $id);
			}
			else{
				print "can't find a family for $id +++ $currentFullName\n";
			}
		}
	}



	# 
	# Check the consistency of RecNameFull in a family
	# Remove inconsistant member to a separate family
	#
	print "splitting families\n";
	for  $protein (sort keys %familyTable) {
		my @componentList = &splitFamilty($protein); 
		my $size = scalar(@componentList);
		if($size > 1 ){
			delete($familyTable{$protein});
			my $i = 0;
			foreach $newFamily (@componentList) {
				$familyTable{$protein."=".$i} = [split("#", $newFamily)]; 
				$i++;
			}
		}
	}

	#
	# Find the most common names in this family
	# Populate %familyInfo
	#
	print "Extracting familyInfo\n";
	for  $protein (sort keys %familyTable) {
		@idArray = @{$familyTable{$protein}};
		$arraySize = @idArray;
		@RecNameFull = ();
		@RecNameEC = ();
		@GeneName = ();
		@GOarray = ();
		foreach $id (@idArray) {
			@swissInfo = @{$swissTable{$id}};
			push(@RecNameFull, $swissInfo[3]);
			push(@RecNameEC, $swissInfo[4]);
			push(@GeneName, $swissInfo[5]);
			push(@GOarray, $swissInfo[8]);
		}

		# find family info
		if($arraySize == 1){
			$FamilyRecNameFull = $RecNameFull[0];
		}
		else{
			$commonFamilyRecNameFull = &findMostCommonItem(@RecNameFull);
			$SameName = "true";
			foreach $name (@RecNameFull){
				if(not $name eq $commonFamilyRecNameFull){
					$SameName = "false";
					last;
				}
			}
			if($SameName eq "true"){
				$FamilyRecNameFull = $commonFamilyRecNameFull;
			}
			elsif( $protein =~ /^(.+)=/){
				$FamilyRecNameFull = $1;
			}
			else{
				$FamilyRecNameFull = $protein;
			}
		}
		$FamilyRecNameEC = &findMostCommonItem(@RecNameEC);
		$FamilyGeneName = &findMostCommonItem(@GeneName);
		$FamilyGOterm = &findConsensusGO(@GOarray);

		$familyInfo{$protein} = [ $FamilyRecNameFull, $FamilyRecNameEC, $FamilyGeneName, $FamilyGOterm ];

	}


	# replace the key to the ID
	my @oldKey = sort keys %familyTable;
	for  $protein (@oldKey) {
		@idArray = @{$familyTable{$protein}};
		@ProteinID = ();
		foreach $id (@idArray) {
			@swissInfo = @{$swissTable{$id}};
			push(@ProteinID, $swissInfo[0]);
		}
	       
		$FamilyKey = &findMostCommonItem(@ProteinID);
	       
	       
		if( exists $familyTable{$FamilyKey} ){
			$i = 1;
			while ( exists $familyTable{"$FamilyKey$i"} ){
				$i++;
			}
			$FamilyKey .= $i;
		}

		# The delete operator returns the value being deleted.
		$familyInfo{$FamilyKey} = delete $familyInfo{$protein};			
		$familyTable{$FamilyKey} = delete $familyTable{$protein};			
	}

}

sub mergeFamilies{
	my @familyNames = sort keys %familyTable;
	my $arraySize = @familyNames;

	# build adjacency matrix
	my %adjacency;
	for ($i = 0; $i < $arraySize; $i++){
		my @neighbors = ();
		$adjacency{$familyNames[$i]} = [@neighbors]
	}

	my @RecNameFullArray = ();
	my @family;
	for ($i = 0; $i < $arraySize; $i++){
		@family = @{$familyTable{$familyNames[$i]}};
		@RecNameFull = ();
		foreach $id (@family) {
			@swissInfo = @{$swissTable{$id}};
			push(@RecNameFull, $swissInfo[3]);
		}
		$FamilyRecNameFull = &findMostCommonItem(@RecNameFull);
		push(@RecNameFullArray, $FamilyRecNameFull);
	
	}

	
	print "start building adj matrix\n";

	my @family0;
	my @family1;
	for ($i = 0; $i < $arraySize - 1; $i++){
		for( $j = $i + 1; $j < $arraySize; $j++ ){

			@family0 = @{$familyTable{$familyNames[$i]}};
			@family1 = @{$familyTable{$familyNames[$j]}};

			# don't merge two large families
	#		if( scalar @family0 >= 3 and scalar @family0 >= 3 ){
	#			next;
	#		}

			$RecNameFull0 = $RecNameFullArray[$i];
			$RecNameFull1 = $RecNameFullArray[$j];

		#	print "$familyNames[$i] === $familyNames[$j]\n";
			if(&isCompatible2($RecNameFull0, $RecNameFull1) ne "false"){

# commented out the blast score filtering				
	#			$similarityScore = 0;
	#			$similarEnough = 0;
	#			foreach $id0 (@family0) {
	#				foreach $id1 (@family1) {
	#					$similarityScore = getBlastScore($id0,$id1);
	#					if( $similarityScore > $BlastCutoff1 ){
	#						$similarEnough = 1;
	#						last;
	#					}
	#				}
	#				if($similarEnough){
	#					last;
	#				}
	#			}
	#			if( $similarEnough ){
					push(@{$adjacency{$familyNames[$i]}}, $familyNames[$j]);
					push(@{$adjacency{$familyNames[$j]}}, $familyNames[$i]);
	#			}

			}
		}	
	}

	print "bfs\n";

	my @componentList = breathFirstSearch(\%adjacency);
	print "merging\n";
	foreach $component (@componentList) {
		print "$component aa\n";
		@mergedFamilies = split("#", $component); 
		@mergedFamiliesIdList = ();

		my $size = scalar(@mergedFamilies);
		if($size > 1 ){
			
			foreach $family (@mergedFamilies) {
				push(@mergedFamiliesIdList, @{$familyTable{$family}});
				delete($familyTable{$family});
			}
			$familyTable{$component} = [@mergedFamiliesIdList];
		}

	}
}


sub isCompatible2{
	my $name0 = shift;
	my $name1 = shift;

	if($name0 eq $name1){
		return $name0;
	}

	my $LastWord = "";
	# \Q  tells Perl where to start escaping special characters, and the \E  tells it where to stop
	if( $name0 =~ /^\Q$name1\E(\s+|-)(\S+)$/ ){
		$LastWord = $2;
		if($LastWord =~ /^\d+$/ or $LastWord =~ /^[IVX]+$/){
			return $name1;
		}
	}
	elsif( $name1 =~ /^\Q$name0\E(\s+|-)(\S+)$/ ){
		$LastWord = $2;
		if($LastWord =~ /^\d+$/ or $LastWord =~ /^[IVX]+$/){
			return $name0;
		}
	}

	$name0String = $name0;
	$name1String = $name1;
	$name0String =~ s/-/ /g;
	$name1String =~ s/-/ /g;
	@name0Words = split(/\s+/, $name0String);
	@name1Words = split(/\s+/, $name1String);

	my $LastWord0 = pop(@name0Words);
	my $LastWord1 = pop(@name1Words);

	if($LastWord0 =~ /^\d+$/ or $LastWord0 =~ /^[IVX]+$/){
		if($LastWord1 =~ /^\d+$/ or $LastWord1 =~ /^[IVX]+$/){
			if($name0Words[-1] ne "subunit" and $name1Words[-1] ne "subunit" ){
				
				$name0String = join(" ",  @name0Words);
				$name1String = join(" ",  @name1Words);
				if($name0String eq $name1String){
					# remove last word
					$name0 =~ s/(\s+|-)+$LastWord0$//;		
					return $name0;
				}
			}
		}
	}


	# false
	return "false";
}

sub findConsensusGO{
	my @FamilyGoTerms = @_;
	my $proteinCount = scalar @FamilyGoTerms;
	my %termFreq = ();
	foreach $proteinGoTerms (@FamilyGoTerms){
		next if($proteinGoTerms eq "NA");
		@termList = split(':', $proteinGoTerms);
		# remove "GO" in the "GO:0000000".
		shift(@termList);
		foreach $term (@termList){
			if( exists($termFreq{$term})){
				$termFreq{$term}++;
			}
			else{
				$termFreq{$term}=1;
			}
		}
	}

	my @ConsensusGO;
	for  $term (sort keys %termFreq) {
		if( $termFreq{$term} * 2 >= $proteinCount ){
			push(@ConsensusGO,$term);
		}
	}
	if( scalar(@ConsensusGO) > 0 ){
		unshift(@ConsensusGO, "GO");
		return join(':',@ConsensusGO);
	}
	else{
		return "NA";
	}
}

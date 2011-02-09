#!/usr/bin/perl

# move genomes from inputDir to outputDir
$inputDir = "E:\\Research\\swissprot\\HAMAP_Genomes\\hamap_genomes";
$output = "E:\\Research\\swissprot\\HAMAP_Genomes\\hamap.txt";
$characterFile = "E:\\Research\\swissprot\\HAMAP_Genomes\\genomeCharacters.txt";

if(not -d $inputDir){ die("can't find the $inputDir directory for input"); }

opendir(INDIR, $inputDir) or die "can't opendir $inputDir";
#open TABLE, ">$output" or die("couldn't open $output \n");

%speciesTable;
%infoTable;
%flagellaTable;
while (defined($file = readdir(INDIR))) {
	next if $file =~ /^\.\.?$/;
	if($file =~ /(\S+)\.html$/){
		$organismCode = $1;
		open WEBPAGE, "<$inputDir/$file" or die("couldn't open $file \n");
		$species = "NA";
		$info = "NA";
		$flagellaPresence = "NA";
		while( $row = <WEBPAGE> ){
			chomp $row;
			if($row =~ /<div class="section">Species:&#160;&#160;(.+)<\/div>/){
			#	print "$row\n";
				$species = $1;
			}
			elsif( $row =~ /\Q<td valign="top"><strong>Description:<\/strong><\/td>\E/ ){
				$row = <WEBPAGE>;
				$row =~ /<td>(.+)<\/td>/;
				$info = $1;
			}
			elsif($row =~ /Presence of flagella:/){
				$row = <WEBPAGE>;
				$row = <WEBPAGE>;
				$row =~ /^\s+(.+)<br>/;
				$flagellaPresence = $1
			}
		}
#		print TABLE "$organismCode;\t$species;\t$flagellaPresence;\t$info\n";
		$speciesTable{$organismCode} = $species;
		$infoTable{$organismCode} = $info;
		$flagellaTable{$organismCode} = $flagellaPresence;
	}
}

closedir(INDIR);
#close( TABLE );

open CHAR, ">$characterFile" or die("couldn't open $characterFile \n");
print CHAR "Code\tSpecies\tFlagella\tGram\tDesc\n";
for  $organismCode (sort keys %speciesTable) {
	print CHAR "$organismCode\t";

	# print species
	if(exists $speciesTable{$organismCode}){
		$species = $speciesTable{$organismCode};
		print CHAR "$species\t";
	}
	else{
		print CHAR "NA\t";

	}	

	# print flagella
	if(exists $flagellaTable{$organismCode}){
		$flagellaPresence = $flagellaTable{$organismCode};
		print CHAR "$flagellaPresence\t";
	}
	else{
		print CHAR "NA\t";

	}

	# print characters
	$info = $infoTable{$organismCode};

	# get Gram positive/negative
	if($info =~ /Gram-negative/ ){
		print CHAR "negative\t"

	}
	elsif($info =~ /Gram-positive/){
		print CHAR "positive\t"
	}
	else{
		print CHAR "NA\t";
	}

	# print info
	print CHAR "$info\n";
}

close CHAR;

#!/usr/bin/env perl

use SWISS::Entry;
use SWISS::OCs;
use SWISS::GNs;

# usage: -i uniprot_sprot.dat -t swissTable.txt -f swissProtein.fasta
$input = "uniprot_sprot.dat";
$table = "swissTable.txt";
$fasta = "swissProtein.fasta";
$blastseq = "blastseq";

# Parse the command line arguments.  
for($i = 0; $i < @ARGV-1; $i++) {
  if($ARGV[$i] =~ /^\-i/) { $input = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-t/) { $table = $ARGV[$i+1]; $i++; next; }
  if($ARGV[$i] =~ /^\-f/) { $fasta = $ARGV[$i+1]; $i++; next; }
}

open INPUT, "<$input" or die("couldn't open $index \n");
open TABLE, ">$table" or die("couldn't open $table \n");
open FASTA, ">$fasta" or die("couldn't open $fasta \n");

# make a directory to store sequence files for blast
if(-d $blastseq){ system("rm -rf $blastseq"); }
mkdir "$blastseq", 0755 or die("couldn't create blastseq directory");

# set separator
local $/ = "\n//\n";

print TABLE "ID;\tOC0;\tOC1;\tRecNameFull;\tRecNameEC;\tGeneName;\tFastaPosition;\tProteinLength;\tGO\n";
$fasta_size = 100;
$counter = 0;
$seq_file_id = 1;

open BLAST_SEQ, ">$blastseq/Query_$seq_file_id.fasta" or die("couldn't open blast seq file $seq_file_id \n");

# Read an entire record at a time
while (<INPUT>){
  # Read the entry
  $entry = SWISS::Entry->fromText($_);
  
  next if $entry->isFragment;     # skip entries that are fragments
  

  $id = "NA";
  @OCs = ("NA", "NA");
  $RecNameFull = "NA";
  $RecNameEC = "NA";
  $GeneName = "NA";
  $sequence = "NA";
  $GoTerms = "";

  # get the primary ID of each entry.
  $id = $entry->ID;
  $sequence = $entry->SQ;

  # get OC
  $i = 0;
  foreach my $oc ($entry->OCs->elements) {
	  $OCs[$i] = $oc;
	  $i++;
  }


  # get GO terms from DRs
  my @terms = ();
  foreach my $DR ($entry->DRs->list ){
	  my @DRarray = @$DR;
	  foreach my $xref (@DRarray){
		  my @itemArray = @$xref;
		  if($itemArray[0] eq "GO" ){
			  $itemArray[1] =~ /GO:(\d+)/;
			  push(@terms, $1);
		  }
	  }

  }

  if(scalar(@terms) > 0 ){
	  $GoTerms = join(':', @terms);
	  $GoTerms = "GO:$GoTerms";
  }
  else{
	  $GoTerms = "NA";
  }

  # get RecName
  foreach my $de ($entry->DEs->elements) {
	  if( $de->category eq "RecName" && $de->type eq  "Full")
	  {
		  $RecNameFull = $de->text;
	  }
	  if( $de->category eq "RecName" && $de->type eq  "EC")
	  {
		  $RecNameEC = $de->text;
	  }	  
  }
  

  # get the first gene name
  foreach my $geneGroup ($entry->GNs->elements) {
	foreach my $gn ($geneGroup->Names->elements) {
		$GeneName = $gn->text;
		last;
	}
	last;
  }

  if(not -d $blastseq){ die("can't find the blastseq directory for input"); }

  if( filterProteins($id, $OCs[0], $OCs[1], $RecNameFull, $RecNameEC, $GeneName) ){
  	$FastaPosition = tell(FASTA);
	$proteinLength = length($sequence);
  	print TABLE "$id;\t$OCs[0];\t$OCs[1];\t$RecNameFull;\t$RecNameEC;\t$GeneName;\t$FastaPosition;\t$proteinLength;\t$GoTerms\n";
  	print FASTA ">$id\n";
  	print FASTA "$sequence\n";

	if( $counter < $fasta_size ){
		print BLAST_SEQ ">$id\n"; 
  		print BLAST_SEQ "$sequence\n";
		$counter++;
	}
	else{
		close BLAST_SEQ;
		$seq_file_id++;
		open BLAST_SEQ, ">$blastseq/Query_$seq_file_id.fasta" or die("couldn't open blast seq file $seq_file_id \n");
		print BLAST_SEQ ">$id\n"; 
  		print BLAST_SEQ "$sequence\n";
		$counter = 1;
	}
   }

}

close TABLE;
close FASTA;
close INPUT;
close BLAST_SEQ;

sub filterProteins {

	my @mySwissEntry = @_;
#       @swissEntry = $protein, $OC0, $OC1, $RecNameFull, $RecNameEC, $GeneName
	my $OC0 = $mySwissEntry[1];
	my $OC1 = $mySwissEntry[2];
	my $myRecNameFull = $mySwissEntry[3];
	my $myRecNameEC = $mySwissEntry[4];
	my $myGeneName = $mySwissEntry[5];

	# remove viruses, animals and plants
#	if($OC0 eq "Viruses" or $OC1 eq "Metazoa" or $OC1 eq "Viridiplantae" ){
#		return 0;
#	}

	if($OC0 eq "Viruses" or $OC0 eq "Eukaryota" ){
		return 0;
	}

#	if($OC0 eq "Viruses"){
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

	return 1;

}




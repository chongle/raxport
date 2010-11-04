#!/usr/bin/perl

@directories = ();
@options = ();
$Raxport = "Raxport.exe";

# Parse the command line arguments.  
for($i = 0; $i < @ARGV; $i++) {
	if($ARGV[$i] =~ /\-\-/){
		push(@options, $ARGV[$i]);
	}
	elsif($ARGV[$i] =~ /\\/){
		push(@directories, $ARGV[$i]);
	}
	else{
		
	}
}

$optionText = join(' ', @options);

foreach $dir (@directories){
	opendir(DIR, $dir) or die "can't opendir $dir: $!";
	while (defined($file = readdir(DIR))) {
		if( $file =~ /(\S+)\.RAW/ or $file =~ /(\S+)\.raw/ ){
			print "\n\n$Raxport $optionText -f $dir\\$file \n\n";
			system("$Raxport $optionText -f $dir\\$file");
		}
	}

}

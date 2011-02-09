#!/usr/bin/perl

# move genomes from inputDir to outputDir
$inputDir = "/data1/pcl/genomes";
$outputDir = "/data1/pcl/swissprot/genomes";

if(not -d $inputDir){ die("can't find the $inputDir directory for input"); }
if(-d $outputDir){ system("rm -rf $outputDir"); }
mkdir "$outputDir", 0755 or die("couldn't create $outputDir directory");

opendir(INDIR, $inputDir) or die "can't opendir $inputDir";
while ( ($subdir = readdir(INDIR))) {
	next if $subdir =~ /^\.\.?$/;     # skip . and ..
	if( -d "$inputDir/$subdir" ){
		print "$subdir\n";
		@fileList = ();
		opendir(SUB, "$inputDir/$subdir") or die "can't open $subdir";
		while (defined($file = readdir(SUB))) {
			next if $file =~ /^\.\.?$/;
			if($file =~ /(\S+)\.faa/){
				push(@fileList, "$inputDir/$subdir/$file");
			}
		}
		if( scalar @fileList > 0 ){
			$files = join(' ', @fileList);
			system("cat $files > $inputDir/$subdir/$subdir.faa");
			system("mv $inputDir/$subdir/$subdir.faa $outputDir");

		}
		else{
			print "can't find any faa file under $subdir \n";
		}
		closedir(SUB)
	}
}
closedir(INDIR);


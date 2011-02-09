$dbhmmDir = "/data1/pcl/swissprot/dbhmm";
if(-d $dbhmmDir){
	system("rm -rf $dbhmmDir");
}
mkdir "$dbhmmDir", 0755 or die("couldn't create qsub directory");
system("find hmm -name \"*.hmm\" | xargs cat > $dbhmmDir/swiss.hmm");
system("/data1/pcl/swissprot/bin/hmmer3/hmmpress $dbhmmDir/swiss.hmm");
system("/data1/pcl/swissprot/bin/hmmer3/hmmstat $dbhmmDir/swiss.hmm > dbhmm.stat");

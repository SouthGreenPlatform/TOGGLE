#!/usr/bin/perl -w

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################

#Will test if samTools module work correctly works correctly
use strict;
use warnings;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);

########################################
#use of samtools modules ok
########################################
use_ok('toolbox') or exit;
use_ok('samTools') or exit;

can_ok( 'samTools','samToolsFaidx');
can_ok( 'samTools','samToolsIndex');
#can_ok( 'samTools','samToolsSort');
#can_ok( 'samTools','mergeHeader');
#can_ok( 'samTools','samToolsMerge');
can_ok( 'samTools','samToolsView');
#can_ok( 'samTools','mergeHeader');
#can_ok( 'samTools','samToolsIdxstats');
#can_ok( 'samTools','samToolsDepth');
##can_ok( 'samTools','samToolsFlagstat');

use toolbox;
use samTools;


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"samTools\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf samTools_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot remove the previous log files with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/samtoolsTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $creatingDirCom \n$!\n");

my $originalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $originalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastaRef in the test directory with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested

my $originalBamFile="../DATA/expectedData/RC3.PICARDTOOLSSORT.bam";
my $bamFile="$testingDir/RC3.PICARDTOOLSSORT.bam";
my $bamFileCopyCom="cp $originalBamFile $bamFile";
system($bamFileCopyCom) and die("ERROR: $0 :Cannot copy the file $originalBamFile in the test directory with the command $bamFileCopyCom\n$!\n");

toolbox::readFileConf("software.config.txt");


################################################################################################
###Samtools faidx
################################################################################################
is(samTools::samToolsFaidx($fastaRef),1,'Ok for samToolsFaidx running');

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
my $observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools faidx output structure');


################################################################################################
##Samtools index
################################################################################################
is(samTools::samToolsIndex($bamFile),1,'Ok for samToolsIndex running');

###Checking the correct structure for the output file using md5sum
$expectedMD5sum = "8e1c53324486a6455790a26dfdc5a464";
$observedMD5sum=`md5sum $bamFile.bai`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools index output structure');


################################################################################################
##Samtools view
################################################################################################
my $bamFileOut="$testingDir/RC3.SAMTOOLSVIEW.bam";

my %optionsRef = ("-h" => '', "-b" => '', "-F" => "0*02");
my $optionRef = \%optionsRef; 
is(samTools::samToolsView($bamFile, $bamFileOut, $optionRef),1,'Ok for samToolsView running');


###Checking the correct structure for the output file using md5sum
$expectedMD5sum = "af10b544a1fa63c5544115f56445a59c";
$observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools view output structure');

exit;

################################################################################################
##Samtools sort
################################################################################################
#is(samTools::samToolsSort($bamFile),1,'Ok for samToolsSort running');
#
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "118c12f23985225eee198927007c2e73";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/RC3.PICARDTOOLSSORT.bam`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools sort output structure');



################################################################################################
###Samtools merge
################################################################################################
#my @bamFiles=('../DATA/expectedData/RC3.SAMTOOLSVIEW.bam','../DATA/expectedData/RC3.PICARDTOOLSSORT.bam');
#my $headerExtractCommand="samtools view -H ../DATA/expectedData/RC3.SAMTOOLSVIEW.bam > ../DATA-TEST/samtoolsTestDir/headerFile.sam";  #Extracting header for the following test
#TODO: {
#        system($headerExtractCommand) and die ("\nCannot launch the header extract command: $!\n Aborting tests\n");
#        $optionsHachees=$configInfos->{'samtools merge'};
#        is(samTools::samToolsMerge(\@bamFiles,"$testingDir/out.bam",'../DATA-TEST/samtoolsTestDir/headerFile.sam',$optionsHachees),1,'Ok for samToolsMerge running');
#        
#    }  
###Verifying if the output files are existing for sort
#my $expectedOutputMerge="$testingDir/out.bam";
#is(toolbox::existsFile($expectedOutputMerge),1,'Ok for samToolsMerge produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "d326235da0035dbe76e8214cadb46f8f";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/out.bam`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools merge output structure');
#
#
#################################################################################################
###mergeHeader
#################################################################################################
###Running
#TODO: {
#        is(samTools::mergeHeader(\@bamFiles,"$testingDir/tested_header.txt"),1,'Ok for mergeHeader running');   
#}
###Verifying if the output files are existing for sort
#my $expectedOutputMergeHeader="$testingDir/tested_header.txt";
#is(toolbox::existsFile($expectedOutputMergeHeader),1,'Ok for mergeHeader produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "3911b0e41d2336eba54c973e6a97c66a";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/tested_header.txt`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the mergeHeader output structure');

################################################################################################
##Samtools idxstats
################################################################################################
#is(samTools::samToolsIdxstats($bamFile,"$testingDir/samIdx.txt"),1,'Ok for samtools Idx Stats running');
###Verifying if the output files are existing for sort
#my $expectedOutputIdxstats="$testingDir/samIdx.txt";
#is(toolbox::existsFile($expectedOutputIdxstats),1,'Ok for samtools Idx Stats produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "689274eaea49281bff09c0924d200df4";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/samIdx.txt`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools Idx Stats output structure');

################################################################################################
##Samtools Depth
################################################################################################
#is(samTools::samToolsDepth(\@bamFiles,"$testingDir/depth.txt"),1,'Ok for samtools Depth running');
###Verifying if the output files are existing for sort
#my $expectedOutputDepth="$testingDir/depth.txt";
#is(toolbox::existsFile($expectedOutputDepth),1,'Ok for samtools Depth produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "f6147d7552230e774787bfc073627f9f";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/depth.txt`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools Depth output structure');

################################################################################################
##Samtools Flagstat
################################################################################################
#is(samTools::samToolsFlagstat($bamFile,"$testingDir/flagstats.txt"),1,'Ok for samtools Flagstats running');
###Verifying if the output files are existing for sort
#my $expectedOutputFlag="$testingDir/flagstats.txt";
#is(toolbox::existsFile($expectedOutputFlag),1,'Ok for samtools Flagstats produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = "cb54a20b65967bbd0c2a15dbcfb4a122";
#$observedMD5sum=`md5sum ../DATA-TEST/samtoolsTestDir/flagstats.txt`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools Flagstats output structure');

exit;
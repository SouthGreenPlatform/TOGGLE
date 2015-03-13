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
use_ok('tophat') or exit;

can_ok( 'tophat','indexRef');
can_ok( 'tophat','tophatRun');

use toolbox;
use tophat;


#######################################
#Creating the IndividuSoft.txt file for samTools_test.t
#######################################
my $creatingCommand="echo \"TEST\ntophatTEST\" > individuSoft.txt";
system($creatingCommand) and die ("\nCannot create the individuSoft.txt file needed for logging the test: $!\nAborting\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf TEST_tophatTEST_log.*";
system($cleaningCommand) and die ("\nCannot clean the previous log files for this test: $!\nAborting\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/tophatTestDir";
toolbox::makeDir($testingDir,1);                                    #Allows to have a working directory for the tests

my $OriginalFastaRef="../DATA/RC1/Reference.fasta";
my $fastaRef="$testingDir/Reference.fa";
my $refCopyCom="cp $OriginalFastaRef $fastaRef";
system($refCopyCom) and die ("\nCannot copy the Reference for test:$!\nAborting\n");     #Now we have a ref to be tested

my $originalFastq1="../DATA/RC1/0_PAIRING_FILES/RC1_1.fastq";
my $fastq1="$testingDir/RC1_1.fastq";
my $fastq1CopyCom="cp $originalFastq1 $fastq1";
system($fastq1CopyCom) and die ("\nCannot copy the FASTQ File 1 for test:$!\nAborting\n");     #Now we have a fastq1 to be tested

my $originalFastq2="../DATA/RC1/0_PAIRING_FILES/RC1_2.fastq";
my $fastq2="$testingDir/RC1_2.fastq";
my $fastq2CopyCom="cp $originalFastq1 $fastq1";
system($fastq2CopyCom) and die ("\nCannot copy the FASTQ File 2 for test:$!\nAborting\n");     #Now we have a fastq2 to be tested

########################################
#reading config
########################################
toolbox::readFileConf("software.config.txt");

########################################
#Tests
########################################

#################################################################################################
###tophat indexRef
##Running
is(tophat::indexRef($fastaRef),1,'Ok for tophat indexRef running');
###Verifying if the output files are existing for sort
my $expectedOutputFilesNumber="7
";
is(`ls $testingDir/ | grep -c Reference`,$expectedOutputFilesNumber,'Ok for tophat indexRef produced files');
####Checking the correct structure for the output file using md5sum
my $expectedMD5sum = 'a14c09818383990a1b6d0ad3b8259007  ../DATA-TEST/tophatTestDir/Reference.fa
c5b3489bd75784efff7cabdd9ae7fcd7  ../DATA-TEST/tophatTestDir/Reference.fa.1.bt2
a53657f51ff9d1604cbdbdb22b321b15  ../DATA-TEST/tophatTestDir/Reference.fa.2.bt2
64de2859464956db84c5d97442c8e310  ../DATA-TEST/tophatTestDir/Reference.fa.3.bt2
f6796769163bff9ea37faa188b99090f  ../DATA-TEST/tophatTestDir/Reference.fa.4.bt2
138cb44fa8d50632e7b5e9c96ba610c5  ../DATA-TEST/tophatTestDir/Reference.fa.rev.1.bt2
b3763aee9a0909174e73b08f05b5c52c  ../DATA-TEST/tophatTestDir/Reference.fa.rev.2.bt2
';
my $observedMD5sum=`md5sum ../DATA-TEST/tophatTestDir/Reference*`;#md5sum values observed for the current files produced
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the tophat indexRef output structure');


#################################################################################################
###tophat running
##Running
is(tophat::tophatRun($fastaRef,$fastq1,$fastq2),1,'Ok for tophat tophatRun running');
###Verifying if the output files are existing for sort
#$expectedOutputFilesNumber="7
#";
#is(`ls $testingDir/ | grep -c Reference`,$expectedOutputFilesNumber,'Ok for tophat tophat Run produced files');
####Checking the correct structure for the output file using md5sum
#$expectedMD5sum = 'c5b3489bd75784efff7cabdd9ae7fcd7  ../DATA-TEST/tophatTestDir/Reference.1.bt2
#';
#$observedMD5sum=`md5sum ../DATA-TEST/tophatTestDir/`;#md5sum values observed for the current files produced
#is($observedMD5sum,$expectedMD5sum,'Ok for the content of the tophat Run output structure');


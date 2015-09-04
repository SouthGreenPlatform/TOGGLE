#!/usr/bin/perl -w

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################

#Will test if picardsTools module works correctly

use strict;
use warnings;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use FindBin qw($Bin);
use lib qw(../Modules/);

########################################
#use of samtools modules ok
########################################
use_ok('toolbox') or exit;
use_ok('picardTools') or exit;
can_ok( 'picardTools','picardToolsMarkDuplicates');
can_ok( 'picardTools','picardToolsCreateSequenceDictionary');
can_ok( 'picardTools','picardToolsSortSam');

use toolbox;
use picardTools;

toolbox::readFileConf("software.config.txt");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"picardTools\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCommand\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf picardTools_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCommand \n$!\n");

#########################################
#Remove the files and directory created by the previous test
#########################################
$cleaningCommand="rm -Rf ../DATA-TEST/picardtoolsTestDir";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous test dir with the command $cleaningCommand \n$!\n");

########################################
#Creation of test directory
########################################
my $testingDir="../DATA-TEST/picardtoolsTestDir";
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

##########################################
#picardToolsCreateSequenceDictionary test
##########################################
### Input files
my $originalRefFile = "../DATA/expectedData/Reference.fasta";     # Ref fasta file 
my $refFile = "$testingDir/Reference.fasta";         # Ref fasta file for test
my $refFileCopyCom = "cp $originalRefFile $refFile"; # command to copy the original Ref fasta file into the test directory
system ($refFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalRefFile in the test directory with the command $refFileCopyCom\n$!\n");    # RUN the copy command

my $refFileDict = "$testingDir/Reference.dict";      # output of this module

### TEST OF FUNCTION
is(picardTools::picardToolsCreateSequenceDictionary($refFile,$refFileDict),1,'Test for picardTools::picardToolsCreateSequenceDictionary');     # test if picardTools::picardToolsCreateSequenceDictionary works

### TEST OF STRUCTURE
my $nbOfLineExpected= "952";
my $nbOfLineObserved= `wc -l $refFileDict`;
my @nameless=split /\s/, $nbOfLineObserved;

is_deeply($nameless[0],$nbOfLineExpected,'Test for the lines number of the output file of picardTools::picardToolsCreateSequenceDictionary');

my $firstM5Expected= "bedc1338f03b37384785c231069eae0e";
my $firstM5Observed= `head -n2 $refFileDict`;
chomp $firstM5Observed;
my @m5Observed=split /M5:/, $firstM5Observed;

is_deeply($m5Observed[1],$firstM5Expected,'Test for the MD5 value in the first line of the output file of picardTools::picardToolsCreateSequenceDictionary');

my $lastM5Expected= "872df605abe21dcfe7cfcc7f4d491ea1";
my $lastM5Observed= `tail -n 1 $refFileDict`;
chomp $lastM5Observed;
@m5Observed=split /M5:/, $lastM5Observed;

is_deeply($m5Observed[1],$lastM5Expected,'Test for the MD5 value in the last line of the output file of picardTools::picardToolsCreateSequenceDictionary');




##########################################
#picardToolsSortSam test
##########################################
###### SINGLE ######
## Input files test for single analysis
my $originalSamFile = "../DATA/expectedData/RC3.BWASAMSE.sam";        # original SAM file
my $samFile = "$testingDir/RC3.BWASAMSE.sam";                                 # SAM file of test
my $samFileCopyCom = "cp $originalSamFile $samFile";                            # command to copy the original Ref fasta file into the test directory
system ($samFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalSamFile in the test directory with the command $samFileCopyCom\n$!\n");    # RUN the copy command

my $bamFileOut = "$testingDir/RC3Single.PICARDTOOLSSORT.bam";

my %optionsRef = ("SORT_ORDER" => "coordinate","VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
my $optionRef = \%optionsRef;                           # Ref of the hash


#### TEST OF FUNCTION
is(picardTools::picardToolsSortSam($samFile,$bamFileOut,$optionRef),1,'Test for picardTools::picardToolsSortSam single');  # test if picardTools::picardToolsSortSam works

#### TEST OF STRUCTURE
my $md5sumExpected = "5b4ca218eb2b99bd8d127c8e26756bf0"";
my $md5sumObserved = `md5sum $bamFileOut`;
@nameless = split (" ", $md5sumObserved);           # to separate the structure and the name of file
$md5sumObserved = $nameless[0];                        # just to have the md5sum result

is_deeply ($md5sumObserved,$md5sumExpected, 'Test for the structure of the output file of picardTools::picardToolsSortSam for single');    # test if the structure of the output file is ok



####### PAIR ######
#### Input files test for pair analysis
$originalSamFile = "../DATA/expectedData/RC3.BWASAMPE.sam";        # original SAM file
$samFile = "$testingDir/RC3.BWASAMPE.sam";                                 # SAM file of test
$samFileCopyCom = "cp $originalSamFile $samFile";                            # command to copy the original Ref fasta file into the test directory
system ($samFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalSamFile in the test directory with the command $samFileCopyCom\n$!\n");    # RUN the copy command

$bamFileOut = "$testingDir/RC3.PICARDTOOLSSORT.bam";

%optionsRef = ("SORT_ORDER" => "coordinate","VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionRef = \%optionsRef;                           # Ref of the hash


#### TEST OF FUNCTION
is(picardTools::picardToolsSortSam($samFile,$bamFileOut,$optionRef),1,'Test for picardTools::picardToolsSortSam pair');  # test if picardTools::picardToolsSortSam works

#### TEST OF STRUCTURE
$md5sumExpected = "b8a74a771fed5eb41fe185b77c70e3f8";
$md5sumObserved = `md5sum $bamFileOut`;
@nameless = split (" ", $md5sumObserved);           # to separate the structure and the name of file
$md5sumObserved = $nameless[0];                        # just to have the md5sum result

is_deeply ($md5sumObserved,$md5sumExpected, 'Test for the structure of the output file of picardTools::picardToolsSortSam for pair');    # test if the structure of the output file is ok



###########################################
##picardToolsMarkDuplicates test
###########################################
my $originalBamFile = "../DATA/expectedData/RC3.GATKINDELREALIGNER.bam";        # original BAM file
my $bamFile = "$testingDir/RC3.GATKINDELREALIGNER.bam";                                 # BAM file of test
my $bamFileCopyCom = "cp $originalBamFile $bamFile";                            # command to copy the original Ref fasta file into the test directory
system ($bamFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalBamFile in the test directory with the command $bamFileCopyCom\n$!\n");    # RUN the copy command

$bamFileOut = "$testingDir/RC3.PICARDTOOLSMARKDUPLICATES.bam";
my $duplicatesFileOut = "$testingDir/RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates";

%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionRef = \%optionsRef;                           # Ref of the hash


#### TEST OF FUNCTION
is(picardTools::picardToolsMarkDuplicates($bamFile, $bamFileOut, $duplicatesFileOut, $optionRef),1,'Test for picardTools::picardToolsMarkDuplicates');  # test if picardTools::picardToolsMarkDuplicates works

#### TEST OF STRUCTURE
my $expectedBam = "ed71faf044357f8bd36eeb30c500100f";
my $observedBam = `md5sum $bamFileOut`;
@nameless = split (" ", $observedBam);           # to separate the structure and the name of file
$observedBam = $nameless[0];                        # just to have the md5sum result

is_deeply ($observedBam, $expectedBam, 'Test for BAM file of picardTools::picardToolsMarkDuplicates');      # test if the structure of BAM file is ok

my $observedDup = `less $duplicatesFileOut`;
like($observedDup, qr/## METRICS CLASS/, 'Test for duplicates file of picardTools::picardToolsMarkDuplicates');     # test if the structure of duplicates file is ok


exit;

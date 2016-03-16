#!/usr/bin/perl -w

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
# Version 1 written by Cecile Monat, Ayite Kougbeadjo, Christine Tranchant, Cedric Farcy, Mawusse Agbessi, Maryline Summo, and Francois Sabot
# Version 2 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Enrique Ortega-Abboud, Julie Orjuela-Bouniol, Sebastien Ravel, Souhila Amanzougarene, and Francois Sabot
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

#Will test if picardsTools module works correctly
use strict;
use warnings;

use Test::More 'no_plan'; 
use Test::Deep;
use Data::Dumper;
use FindBin qw($Bin);
use lib qw(../Modules/);

########################################
#use of picardTools modules ok
########################################
use_ok('toolbox') or exit;
use_ok('picardTools') or exit;
can_ok( 'picardTools','picardToolsMarkDuplicates');
can_ok( 'picardTools','picardToolsCreateSequenceDictionary');
can_ok( 'picardTools','picardToolsSortSam');

use toolbox;
use picardTools;

my $configInfos = toolbox::readFileConf("software.config.txt");

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/picardtoolsTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $expectedData="../../DATA/expectedData/";

########################################
#Creation of test directory
########################################
my $makeDirCmd = "mkdir $testingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");
chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"picardtools\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf picardtools_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

##########################################
#picardToolsCreateSequenceDictionary test
##########################################

my $refFile = "Reference.fasta";  
my $originalRefFile = $expectedData."/".$refFile;    

my $lnCmd = "ln -s $originalRefFile ."; # command to copy the original Ref fasta file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalRefFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command

my $refFileDict = "Reference.dict";      # output of this module

### TEST OF FUNCTION
is(picardTools::picardToolsCreateSequenceDictionary($refFile,$refFileDict),1,'picardTools::picardToolsCreateSequenceDictionary... Running');     # test if picardTools::picardToolsCreateSequenceDictionary works

### TEST OF STRUCTURE
my $nbOfLineExpected= "952";
my $nbOfLineObserved= `wc -l $refFileDict`;
my @nameless=split /\s/, $nbOfLineObserved;

is_deeply($nameless[0],$nbOfLineExpected,'picardTools::picardToolsCreateSequenceDictionary... Test for the lines number of the output file');

my $firstM5Expected= "bedc1338f03b37384785c231069eae0e";
my $firstM5Observed= `head -n2 $refFileDict`;
chomp $firstM5Observed;
my @m5Observed=split /M5:/, $firstM5Observed;

is_deeply($m5Observed[1],$firstM5Expected,'picardTools::picardToolsCreateSequenceDictionary... Test for the MD5 value in the first line of the output file');

my $lastM5Expected= "872df605abe21dcfe7cfcc7f4d491ea1";
my $lastM5Observed= `tail -n 1 $refFileDict`;
chomp $lastM5Observed;
@m5Observed=split /M5:/, $lastM5Observed;

is_deeply($m5Observed[1],$lastM5Expected,'picardTools::picardToolsCreateSequenceDictionary... Test for the MD5 value in the last line of the output file');


##########################################
#picardToolsSortSam test
##########################################
###### SINGLE ######
## Input files test for single analysis
my $samFile = "RC3.BWASAMSE.sam";              # SAM file of test
my $originalSamFile = $expectedData."/".$samFile;        # original SAM file
$lnCmd = "ln -s $originalSamFile .";                  # command to copy the original Ref fasta file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalSamFile in the test directory with the command $lnCmd\n$!\n");

my $bamFileOut = "RC3Single.PICARDTOOLSSORT.bam";

my %optionsRef = ("SORT_ORDER" => "coordinate","VALIDATION_STRINGENCY" => "SILENT");   
my $optionRef = \%optionsRef;                           # Ref of the hash

#### TEST OF FUNCTION
is(picardTools::picardToolsSortSam($samFile,$bamFileOut,$optionRef),1,'picardTools::picardToolsSortSam... Running single');  

#### TEST OF STRUCTURE
my $md5sumExpected = "22e0135ae3488cf16fdb095283ac91c4";
my $md5sumObserved = `md5sum $bamFileOut`;
@nameless = split (" ", $md5sumObserved);           # to separate the structure and the name of file
$md5sumObserved = $nameless[0];                        # just to have the md5sum result

is_deeply ($md5sumObserved,$md5sumExpected, 'picardTools::picardToolsSortSam... Test for the structure of the output file for single');    # test if the structure of the output file is ok



####### PAIR ######
#### Input files test for pair analysis
$samFile = "RC3.BWASAMPE.sam";            # SAM file of test
$originalSamFile = $expectedData."/".$samFile;        # original SAM file
$lnCmd = "ln -s $originalSamFile .";                            # command to copy the original Ref fasta file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalSamFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command

$bamFileOut = "RC3.PICARDTOOLSSORT.bam";

%optionsRef = ("SORT_ORDER" => "coordinate","VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionRef = \%optionsRef;                           # Ref of the hash


#### TEST OF FUNCTION
is(picardTools::picardToolsSortSam($samFile,$bamFileOut,$optionRef),1,'picardTools::picardToolsSortSam... Running pair');  # test if picardTools::picardToolsSortSam works

#### TEST OF STRUCTURE
$md5sumExpected = "7e5a7dc36c0f0b599cc158c599c9913d";
$md5sumObserved = `md5sum $bamFileOut`;
@nameless = split (" ", $md5sumObserved);           # to separate the structure and the name of file
$md5sumObserved = $nameless[0];                        # just to have the md5sum result

is_deeply ($md5sumObserved,$md5sumExpected, 'picardTools::picardToolsSortSam... Test for the structure of the output file for pair');    # test if the structure of the output file is ok




###########################################
##picardToolsMarkDuplicates test
###########################################
my $bamFile = "RC3.GATKINDELREALIGNER.bam";                         # BAM file of test
my $originalBamFile = $expectedData."/".$bamFile;        # original BAM file

$lnCmd = "ln -s $originalBamFile .";                            # command to copy the original Ref fasta file into the test directory
system ($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalBamFile in the test directory with the command $lnCmd\n$!\n");    # RUN the copy command

$bamFileOut = "RC3.PICARDTOOLSMARKDUPLICATES.bam";
my $duplicatesFileOut = "RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates";

%optionsRef = ("VALIDATION_STRINGENCY" => "SILENT");        # Hash containing informations
$optionRef = \%optionsRef;                           # Ref of the hash


#### TEST OF FUNCTION
is(picardTools::picardToolsMarkDuplicates($bamFile, $bamFileOut, $duplicatesFileOut, $optionRef),1,'picardTools::picardToolsMarkDuplicates... Running');  # test if picardTools::picardToolsMarkDuplicates works

#### TEST OF STRUCTURE
my $expectedBam = "63aa7627a3658ad513351fa73f5d8f93";
my $observedBam = `md5sum $bamFileOut`;
@nameless = split (" ", $observedBam);           # to separate the structure and the name of file
$observedBam = $nameless[0];                        # just to have the md5sum result

is_deeply ($observedBam, $expectedBam, 'picardTools::picardToolsMarkDuplicates... Test for BAM file');      # test if the structure of BAM file is ok

my $observedDup = `less $duplicatesFileOut`;
like($observedDup, qr/## METRICS CLASS/, 'picardTools::picardToolsMarkDuplicates... Test for duplicates file');     # test if the structure of duplicates file is ok


exit;

# AF
#PicardToolsValidateSamFile
#PicardToolsCleanSam
#PicardToolsSamFormatConverter
#PicardToolsAddOrReplaceGroup
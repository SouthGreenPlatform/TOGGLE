#!/usr/bin/perl

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

#Will test if samTools module work correctly works correctly
use strict;
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
#can_ok( 'samTools','samToolsIdxstats');
#can_ok( 'samTools','samToolsDepth');
can_ok( 'samTools','samToolsFlagstat');

use toolbox;
use samTools;


my $configInfos = toolbox::readFileConf("software.config.txt");

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/samToolsTestDir";
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
my $creatingCmd="echo \"samTools\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf samTools_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

########################################
#initialisation and setting configs
########################################
my $fastaRef="Reference.fasta";
my $originalFastaRef=$expectedData."/".$fastaRef;
my $lnCmd="ln -s $originalFastaRef .";
system($lnCmd) and die ("ERROR: $0 : Cannot copy the file $originalFastaRef in the test directory with the command $lnCmd\n$!\n");     #Now we have a ref to be tested

my $bamFile="RC3.PICARDTOOLSSORT.bam";
my $originalBamFile=$expectedData."/".$bamFile;
$lnCmd="ln -s $originalBamFile .";
system($lnCmd) and die("ERROR: $0 :Cannot copy the file $originalBamFile in the test directory with the command $lnCmd\n$!\n");



################################################################################################
###Samtools faidx
################################################################################################
is(samTools::samToolsFaidx($fastaRef),1,'samTools::samToolsFaidx... running');

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
my $observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsFaidx... Ok for the content of the samtools faidx output structure');


################################################################################################
##Samtools index
################################################################################################
is(samTools::samToolsIndex($bamFile),1,'samTools::samToolsIndex... running');

###Checking the correct structure for the output file using md5sum
$expectedMD5sum = "29bed7c8c70c24cd84a439d3fc473474";
$observedMD5sum=`md5sum $bamFile.bai`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsIndex... Ok for the content of the samtools index output structure');


################################################################################################
##Samtools view
################################################################################################
my $bamFileOut="RC3.SAMTOOLSVIEW.bam";

my %optionsRef = ("-h" => '', "-b" => '', "-F" => "0*02");
my $optionRef = \%optionsRef; 
is(samTools::samToolsView($bamFile, $bamFileOut, $optionRef),1,'samTools::samToolsView... Ok for samToolsView running');


###Checking the correct structure for the output file using md5sum
$expectedMD5sum = "c5db29f185507f5433f0c08163a2dc57";
$observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsView...Ok for the content of the samtools view output structure');


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
is(samTools::samToolsFlagstat($bamFile,"RC3.SAMTOOLSFLAGSTAT.txt"),1,'samTools::samToolsFlagStat... running');
my $expectedOutputFlag="RC3.SAMTOOLSFLAGSTAT.txt";

####Checking the correct structure for the output file using md5sum
$expectedMD5sum = "c08ff58a41733e3e1ab782ca22653397";
$observedMD5sum=`md5sum RC3.SAMTOOLSFLAGSTAT.txt`;	# structure of the test file
@withoutName = split (" ", $observedMD5sum);    				# to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];     						# just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsFlagStat... Ok for the content of the samtools Flagstats output structure');

exit;

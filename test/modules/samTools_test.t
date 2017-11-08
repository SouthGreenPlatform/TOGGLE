#!/usr/bin/perl

###################################################################################################################################
#
# Copyright 2014-2017 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

#Will test if samTools module work correctly works correctly
use strict;
use warnings;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);

########################################
#use of samtools modules ok
########################################
use_ok('localConfig') or exit;
use_ok('samTools') or exit;

can_ok( 'samTools','samToolsFaidx');
can_ok( 'samTools','samToolsIndex');
can_ok( 'samTools','samToolsSort');
can_ok( 'samTools','mergeHeader');
can_ok( 'samTools','samToolsMerge');
can_ok( 'samTools','samToolsView');
can_ok( 'samTools','samToolsIdxstats');
can_ok( 'samTools','samToolsDepth');
can_ok( 'samTools','samToolsFlagstat');
can_ok( 'samTools','samToolsMpileUp');

use localConfig;
use samTools;

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/samToolsTestDir";
my $cleaningCmd="rm -Rf $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

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

my $bankData="$toggle/data/Bank/";
my $bamData="$toggle/data/testData/samBam/";

# COPY reference file
my $fastaRef="referenceIrigin.fasta";
my $originalFastaRef=$bankData."/referenceIrigin.fasta";
my $copyCmd= "cp $originalFastaRef $fastaRef";           # command to copy the original fasta file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalFastaRef in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command



################################################################################################
###Samtools faidx
################################################################################################

# execution test
is(samTools::samToolsFaidx($fastaRef),1,'samTools::samToolsFaidx');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('individuSoft.txt','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o');
#
is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsFaidx - output list');

# expected content test

my $expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
my $observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsFaidx - output structure');



################################################################################################
##Samtools index
################################################################################################

#Input files
my $originalBam="$bamData/oneBam/RC3-SAMTOOLSVIEW.bam";
my $bamFile="RC3.bam";
$copyCmd= "cp $originalBam $bamFile";           # command to copy the original BAM file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot link the file $originalBam in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

#execution test
is(samTools::samToolsIndex($bamFile),1,'samTools::samToolsIndex');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','RC3.bam','RC3.bam.bai','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsIndex - output list');

# expected output structure
$expectedMD5sum = "ab22a5d60dfc20bc5d4e54608e62095c";
$observedMD5sum=`md5sum $bamFile.bai`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsIndex - output structure');

################################################################################################
##Samtools view
################################################################################################

#Output file
my $bamFileOut="RC3-SAMTOOLSVIEW.bam";
my %optionsRef = ("-h" => '', "-b" => '', "-F" => "0*02");
my $optionsHachees = \%optionsRef; 

#execution test
is(samTools::samToolsView($bamFile, $bamFileOut, $optionsHachees),1,'samTools::samToolsView');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsView - output list');

# expected output structure
$expectedMD5sum = "c5db29f185507f5433f0c08163a2dc57";
$observedMD5sum=`md5sum $bamFileOut`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsView - output structure');


################################################################################################
##Samtools sort
################################################################################################

#Output file
$bamFileOut = "RC3-SAMTOOOLSSORT.bam";

#execution test
is(samTools::samToolsSort($bamFile, $bamFileOut),1,'samTools::samToolsSort');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsSort - output list');

# expected output structure
$expectedMD5sum = "c5db29f185507f5433f0c08163a2dc57";
$observedMD5sum=`md5sum RC3-SAMTOOOLSSORT.bam`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsSort - output structure');


################################################################################################
###Samtools merge
################################################################################################

#Input files
my @bamFiles=("$originalBam","$originalBam");
my $headerExtractCommand="samtools view -H $originalBam > headerFile.sam";  #Extracting header for the following test
system($headerExtractCommand) and die ("\nCannot launch the header extract command: $!\n Aborting tests\n");


#Outputfile
$bamFileOut = "MergeBam.bam";

#execution test
is(samTools::samToolsMerge(\@bamFiles,$bamFileOut,'headerFile.sam'),1,'samTools::samToolsMerge');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsMerge - output list');

# expected output structure
my $expectedLineNumber = "3996";
my $observedLineNumber=`samtools view $bamFileOut | wc -l`;# structure of the test file
chomp $observedLineNumber;
is($observedLineNumber,$expectedLineNumber,'samTools::samToolsMerge - output structure');


#################################################################################################
###mergeHeader
#################################################################################################

#Output file

my $header = "testedHeader.txt";

#execution test

is(samTools::mergeHeader(\@bamFiles,$header),1,'samTools::mergeHeader');   

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samTools_TEST_log.e','samTools_TEST_log.o','testedHeader.txt');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::mergeHeader - output list');

# expected output structure
$expectedMD5sum = "d4e9ac5401144ced034f7f09ba0622a4";
$observedMD5sum=`md5sum $header`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::mergeHeader - output structure');


################################################################################################
##Samtools idxstats
################################################################################################

#output file
my $statsFile = "samIdxStat.txt";

#execution test
is(samTools::samToolsIdxstats($bamFile,$statsFile),1,'samTools::samToolsIdxstats');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samIdxStat.txt','samTools_TEST_log.e','samTools_TEST_log.o','testedHeader.txt');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsIdxstats - output list');

# expected output structure
$expectedMD5sum = "7e6925990265c3ad6ee54dad91cf9682";
$observedMD5sum=`md5sum $statsFile`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsIdxstats - output structure');



################################################################################################
##Samtools Depth
################################################################################################

#outputfile
my $depthFile = "depth.txt";

#execution test
is(samTools::samToolsDepth(\@bamFiles,$depthFile),1,'samTools::samToolsDepth');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('depth.txt','headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samIdxStat.txt','samTools_TEST_log.e','samTools_TEST_log.o','testedHeader.txt');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsIdxstats - output list');

# expected output structure
$expectedMD5sum = "9ae50587f19328ce6c4246f308e941a9";
$observedMD5sum=`md5sum $depthFile`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsDepth - output structure');


################################################################################################
##Samtools Flagstat
################################################################################################

#output file
my $statsBamFile = "RC3-SAMTOOLSFLAGSTAT.txt";

#execution test
is(samTools::samToolsFlagstat($bamFile,$statsBamFile),1,'samTools::samToolsFlagStat');


# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('depth.txt','headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSFLAGSTAT.txt','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samIdxStat.txt','samTools_TEST_log.e','samTools_TEST_log.o','testedHeader.txt');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsIdxstats - output list');

# expected output structure
$expectedMD5sum = "0cbd99e7e43e2a515b53f61c696d4949";
$observedMD5sum=`md5sum $statsBamFile`;	# structure of the test file
@withoutName = split (" ", $observedMD5sum);    				# to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];     						# just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsFlagStat - output structure');


################################################################################################
##Samtools MpileUp
################################################################################################

#output file
my $mpileupFile = "RC3-SAMTOOLSMPILEUP.mpileup";

#execution test
is(samTools::samToolsMpileUp(\@bamFiles,$mpileupFile),1,'samTools::samToolsMpileUp');


# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('depth.txt','headerFile.sam','individuSoft.txt','MergeBam.bam','RC3.bam','RC3.bam.bai','RC3-SAMTOOLSFLAGSTAT.txt','RC3-SAMTOOLSMPILEUP.mpileup','RC3-SAMTOOLSVIEW.bam','RC3-SAMTOOOLSSORT.bam','referenceIrigin.fasta','referenceIrigin.fasta.fai','samIdxStat.txt','samTools_TEST_log.e','samTools_TEST_log.o','testedHeader.txt');

is_deeply(\@observedOutput,\@expectedOutput,'samTools::samToolsMpileUp - output list');

# expected output structure
$expectedMD5sum = "0247045dd18fa40cf0a75b612e1c484d";
$observedMD5sum=`md5sum $mpileupFile`;	# structure of the test file
@withoutName = split (" ", $observedMD5sum);    				# to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];     						# just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'samTools::samToolsMpileUp - output structure');

exit;

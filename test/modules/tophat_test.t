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

use strict;
use warnings;

#Will test if tophat works correctly
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../../modules/);

########################################
#Test of the use of gatk modules
########################################
use_ok('localConfig') or exit;
use_ok('tophat') or exit;

can_ok( 'tophat','bowtieBuild');
can_ok( 'tophat','bowtie2Build');
can_ok( 'tophat','tophat2');

use localConfig;
use tophat;

my $expectedData="$toggle/data/expectedData/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/tophatTestDir";
my $creatingDirCom="rm -rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"tophat\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -rf tophat_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


##########################################
### Test for tophat::bowtieBuild
##########################################

# input file
my $fastaRefIni=$expectedData."/referenceRNASeq.fa";
my $fastaRef="referenceRNASeq.fa";

#copy fasta reference into test directory where the index will be created
my $copyCommand="cp $fastaRefIni .";
system ($copyCommand) and die "ERROR: $0: Cannot copy the refence file with the command $copyCommand \n$!\n";

# execution test
is(tophat::bowtieBuild($fastaRef),$fastaRef,'tophat::bowtieBuild');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('individuSoft.txt','referenceRNASeq.fa','referenceRNASeq.fa.1.ebwt','referenceRNASeq.fa.2.ebwt','referenceRNASeq.fa.3.ebwt','referenceRNASeq.fa.4.ebwt','referenceRNASeq.fa.rev.1.ebwt','referenceRNASeq.fa.rev.2.ebwt','tophat_TEST_log.e','tophat_TEST_log.o');
is_deeply(\@observedOutput,\@expectedOutput,'tophat::bowtieBuild - output list');

# expected output content
my $expectedMD5sum="167cd0622bda91392673aaf207255d2b";
my $observedMD5sum=`md5sum $fastaRef.1.ebwt`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content 1.ebwt');

$expectedMD5sum="dd30c97b610f5dc53cf6f02123fcf807";
$observedMD5sum=`md5sum $fastaRef.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content 2.ebwt');

$expectedMD5sum="8aa5c56a0ba0b0ab7e9e7f3fb7ee4a76";
$observedMD5sum=`md5sum $fastaRef.3.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content 3.ebwt');

$expectedMD5sum="a4ebbf39ff457e410253b571ee79088d";
$observedMD5sum=`md5sum $fastaRef.4.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content 4.ebwt');

$expectedMD5sum="4bd4f23dc8b98a5dc4b56d7f4d89a9b5";
$observedMD5sum=`md5sum $fastaRef.rev.1.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content rev.1.ebwt');

$expectedMD5sum="619322d189d42f4eaede8aaaedf9890e";
$observedMD5sum=`md5sum $fastaRef.rev.2.ebwt`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtieBuild - output content rev.2.ebwt');





###############################################################################################
##tophat::bowtie2Build
###############################################################################################

my %optionsHachees = ();                # Hash containing informations
my $optionHachees = \%optionsHachees;   # Ref of the hash

# execution test
is(tophat::bowtie2Build($fastaRef,$optionHachees),$fastaRef, 'tophat::bowtie2Build');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('individuSoft.txt','referenceRNASeq.fa','referenceRNASeq.fa.1.bt2','referenceRNASeq.fa.1.ebwt','referenceRNASeq.fa.2.bt2','referenceRNASeq.fa.2.ebwt','referenceRNASeq.fa.3.bt2','referenceRNASeq.fa.3.ebwt','referenceRNASeq.fa.4.bt2','referenceRNASeq.fa.4.ebwt','referenceRNASeq.fa.rev.1.bt2','referenceRNASeq.fa.rev.1.ebwt','referenceRNASeq.fa.rev.2.bt2','referenceRNASeq.fa.rev.2.ebwt','tophat_TEST_log.e','tophat_TEST_log.o');
##print Dumper(\@observedOutput);
is_deeply(\@observedOutput,\@expectedOutput,'tophat::bowtie2Build - output list');

# expected output content
###Checking the correct structure for the output file using md5sum
$expectedMD5sum="4e0329b55cd2a67490ef96b7b2567e5d";
$observedMD5sum=`md5sum $fastaRef.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content 1.bt2');

$expectedMD5sum="4ad45a523ecaef6310bb6f7f608eb311";
$observedMD5sum=`md5sum $fastaRef.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content 2.bt2');

$expectedMD5sum="8aa5c56a0ba0b0ab7e9e7f3fb7ee4a76";
$observedMD5sum=`md5sum $fastaRef.3.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content 3.bt2');

$expectedMD5sum="a4ebbf39ff457e410253b571ee79088d";
$observedMD5sum=`md5sum $fastaRef.4.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content 4.bt2');

$expectedMD5sum="beeb0f44030b631af6f182f4a85f045d";
$observedMD5sum=`md5sum $fastaRef.rev.1.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content rev.1.bt2');

$expectedMD5sum="df753d8d0522ee01e2479f99a04525dd";
$observedMD5sum=`md5sum $fastaRef.rev.2.bt2`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::bowtie2Build - output content rev.2.bt2');





#################################################################################################
####tophat::tophat2
#################################################################################################

#input file
my $gffRef=$expectedData."/referenceRNASeq.gff3";

my $fastqFile1=$expectedData."/RNASeq_1.fastq";
my $fastqFile2=$expectedData."/RNASeq_2.fastq";

#tophat option 
%optionsHachees = ( "-i" => "30",
                    "-I" => "20000",
                    "-a" => "8",
                    "-m" => "1",
                    "--no-coverage-search" => '',
                    "-g" => "10",
                    "--bowtie-n" => '',
                    "--library-type" => 'fr-secondstrand',
                    "--microexon-search" => '');
$optionHachees = \%optionsHachees;      # Ref of the hash

my $tmpDir=`pwd`;
chomp $tmpDir;
$tmpDir.='/tophatOut';

# execution test
is(tophat::tophat2($tmpDir, $fastaRef, $fastqFile1, $fastqFile2, $gffRef, $optionHachees),1,'tophat::tophat2');

# expected output test
$observedOutput = `ls $tmpDir`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RNASeq.accepted_hits.bam','RNASeq.align_summary.txt','RNASeq.deletions.bed','RNASeq.insertions.bed','RNASeq.junctions.bed','RNASeq.logs','RNASeq.prep_reads.info','RNASeq.unmapped.bam');
is_deeply(\@observedOutput,\@expectedOutput,'tophat::tophat2 - output list');


# expected output test
$expectedMD5sum="974ed083118b23e373ec072cc91aaf54";
$observedMD5sum=`md5sum $tmpDir/RNASeq.accepted_hits.bam`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'tophat::tophat2- output content bam');


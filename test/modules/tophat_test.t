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
use_ok('bowtie') or exit;

can_ok( 'tophat','tophat2');

use localConfig;
use tophat;
use bowtie; #I know but cannot do something else

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

my $bankData="$toggle/data/Bank/";
my $testData="$toggle/data/testData/rnaseq/";

##########################################
### Test for bowtie::bowtieBuild
##########################################

# input file

my $fastaRefIni=$bankData."/referenceRnaseq.fa";
my $fastaRef="referenceRnaseq.fa";

#copy fasta reference into test directory where the index will be created
my $copyCommand="cp $fastaRefIni ./$fastaRef";
system ($copyCommand) and die "ERROR: $0: Cannot copy the $fastaRefIni file with the command $copyCommand \n$!\n";

# execution test
is(bowtie::bowtieBuild($fastaRef),$fastaRef,'bowtie::bowtieBuild');



###############################################################################################
##bowtie::bowtie2Build
###############################################################################################

my %optionsHachees = ();                # Hash containing informations
my $optionHachees = \%optionsHachees;   # Ref of the hash

# execution test
is(bowtie::bowtie2Build($fastaRef,$optionHachees),$fastaRef, 'bowtie::bowtie2Build');



#################################################################################################
####tophat::tophat2
#################################################################################################

# input file
my $gffRef="$bankData/referenceRnaseqGFF.gff3";

my $fastqFile1=$testData."/pairedOneIndividu/RNASeq_1.fastq";
my $fastqFile2=$testData."/pairedOneIndividu/RNASeq_2.fastq";

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
my $observedOutput = `ls $tmpDir`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('RNASeq.accepted_hits.bam','RNASeq.align_summary.txt','RNASeq.deletions.bed','RNASeq.insertions.bed','RNASeq.junctions.bed','RNASeq.logs','RNASeq.prep_reads.info','RNASeq.unmapped.bam');
is_deeply(\@observedOutput,\@expectedOutput,'tophat::tophat2 - output list');


# expected output test
#$expectedMD5sum="2d3a8ce4123320066b059133be256876";
#$observedMD5sum=`md5sum $tmpDir/RNASeq.accepted_hits.bam`;# structure of the test file
#@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
#$observedMD5sum = $withoutName[0];       # just to have the md5sum result
#is($observedMD5sum,$expectedMD5sum,'tophat::tophat2- output content bam');



my $expectedEndLine="HWI-D00393:103:C6KCUANXX:1:1101:1904:31492	147	EG5_Chr1	68172906	50	118M	=	68172845	-179	AATGGATCAGGCTATAATAATGGAAGCTTGAGTTATGAAAATGGTGAGAGTAACTTTGGTTTGGGCACCAGCAGTTTTGGAAGAAGTGGTCAAACTGGAGCTGCTAACACTTCTCTTA	FFFFFFFFFFFFFFFFFFFFFFFB/FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:118	YT:Z:UU	XS:A:+	NH:i:1
";
my $observedEndLine=`samtools view $tmpDir/RNASeq.accepted_hits.bam | tail -1`  or die ("ERROR: $0 : Cannot execute: samtools view $tmpDir/RNASeq.accepted_hits.bam | tail -1  \n$!\n");
is($observedEndLine,$expectedEndLine,'tophat::tophat2- output content bam - output endFile');

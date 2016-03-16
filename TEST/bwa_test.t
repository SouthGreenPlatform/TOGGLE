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

use strict;
use warnings;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../Modules/);
use Data::Dumper;

########################################
#use of bwa module ok
########################################
use_ok('toolbox') or exit;
use_ok('bwa') or exit;
can_ok( 'bwa','bwaIndex');
can_ok('bwa','bwaAln');
can_ok('bwa','bwaSampe');
can_ok('bwa','bwaSamse');
can_ok('bwa','bwaMem');

use toolbox;
use bwa;

my $configInfos = toolbox::readFileConf("software.config.txt");

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/bwaTestDir";
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
my $creatingCmd="echo \"bwa\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf bwa_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

########################################
#Picking up data for tests
########################################
my $fastaRef="Reference.fasta";
my $originalFastaRef=$expectedData."/Reference.fasta";
my $lnCmd="ln -s $originalFastaRef .";
system($lnCmd) and die ("ERROR: $0 : Cannot copy the Reference for test with the Command $lnCmd \n$!\n");

my $fastqFile1="RC3_1.REPAIRING.fastq";
my $originalFastqFile1=$expectedData."/".$fastqFile1;
$lnCmd="ln -s $originalFastqFile1 .";
system($lnCmd) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile1 for test with the command $lnCmd \n$!\n");    #The sequences are copied for testing

my $fastqFile2="RC3_2.REPAIRING.fastq";
my $originalFastqFile2=$expectedData."/".$fastqFile2;
$lnCmd="ln -s $originalFastqFile2 .";
system($lnCmd) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile2 for test with the command $lnCmd \n$!\n");    #The sequences are copied for testing

my $forwardSaiFileIn="RC3_1.BWAALN.sai";
my $reverseSaiFileIn="RC3_2.BWAALN.sai";
my $sampeFileOut="RC3.BWASAMPE.sam";

#######################################################################################################
####Test for bwa index running
#######################################################################################################
##Test for running
my $optionsHachees=$configInfos->{'BWA index'};
is(bwa::bwaIndex($fastaRef,$optionsHachees),1,'bwa::bwaIndex... running');
###Verify if output are correct for bwa index
my @expectedOutput=('./bwa_TEST_log.e','./bwa_TEST_log.o','./individuSoft.txt','./RC3_1.REPAIRING.fastq','./RC3_2.REPAIRING.fastq','./Reference.fasta','./Reference.fasta.amb','./Reference.fasta.ann','./Reference.fasta.bwt','./Reference.fasta.pac','./Reference.fasta.sa');
my @observedOutput=toolbox::readDir(".");

is_deeply(@observedOutput,\@expectedOutput,'bwa::bwaIndex');

###Test for correct file value of bwa index using a md5sum file control -  work through the different bwa versions
my $expectedMD5sum='b86728bb71903f8641530e61e9687b59  Reference.fasta.amb
a51b6a5152f51b13833a40fe609474ea  Reference.fasta.ann
e4fdc0af9540ee8365e7e324fc5c0cc3  Reference.fasta.bwt
03454e7242900c436d9d7126f492e4d5  Reference.fasta.pac
9243bf066de0cc18aa0d3813f174cae8  Reference.fasta.sa
'; #Expected values for the files produced by the bwa index
my $observedMD5sum=`md5sum Reference.fasta.*`;#md5sum values observed for the current files produced
is($observedMD5sum,$expectedMD5sum,'bwa::bwaIndex... Test for the content of the bwa index output');

#######################################################################################################
###Test for bwa Aln running
#######################################################################################################

$optionsHachees=$configInfos->{'BWA aln'};
is (bwa::bwaAln($fastaRef,$fastqFile1,$forwardSaiFileIn,$optionsHachees),'1',"bwa::bwaAln... Test for bwa Aln running for forward");
is (bwa::bwaAln($fastaRef,$fastqFile2,$reverseSaiFileIn,$optionsHachees),'1',"bwa::bwaAln... Test for bwa Aln running for reverse");

###Verify if output are correct for bwa Aln
@expectedOutput=('./bwa_TEST_log.e','./bwa_TEST_log.o','./individuSoft.txt','./RC3_1.BWAALN.sai','./RC3_1.REPAIRING.fastq','./RC3_2.BWAALN.sai','./RC3_2.REPAIRING.fastq','./Reference.fasta','./Reference.fasta.amb','./Reference.fasta.ann','./Reference.fasta.bwt','./Reference.fasta.pac','./Reference.fasta.sa');
@observedOutput=toolbox::readDir(".");

is_deeply(@observedOutput,\@expectedOutput,'bwa::bwaAln... Test for bwa aln output files');
exit;

###Test for correct file value of bwa aln using a md5sum - BE CAREFUL, the sum changes based on the version of BWA!!
TODO: {
    local $TODO = "The file structure depends of the version of BWA in use. Here, tested for bwa 0.7.9a";
    
    $expectedMD5sum ='1be54c4d8d37870cbd78d93cee30b26f  RC3_1.BWAALN.sai
53ae4174a5bd3cdf2958a59614cf0eb1  RC3_2.BWAALN.sai
';
    $observedMD5sum=`md5sum *.sai`;#md5sum values observed for the current files produced
    is($observedMD5sum,$expectedMD5sum,'bwa::bwaAln... est for the content of the bwa aln output');
}

########################################################################################################
####Test for bwa sampe
#######################################################################################################
is(bwa::bwaSampe($sampeFileOut,$fastaRef,$forwardSaiFileIn,$reverseSaiFileIn,$fastqFile1,$fastqFile2,"RC3"),'1',"bwa::bwaSampe... Test for bwa sampe running");
####Verify if output are correct for sampe
@expectedOutput=('./bwa_TEST_log.e','./bwa_TEST_log.o','./individuSoft.txt','./RC3_1.BWAALN.sai','./RC3_1.REPAIRING.fastq','./RC3_2.BWAALN.sai','./RC3_2.REPAIRING.fastq','./RC3.BWASAMPE.sam','./Reference.fasta','./Reference.fasta.amb','./Reference.fasta.ann','./Reference.fasta.bwt','./Reference.fasta.pac','./Reference.fasta.sa');
@observedOutput=toolbox::readDir($testingDir);
is_deeply(@observedOutput,\@expectedOutput,'bwa::bwaSampe... Test for Output file ok for bwa sampe');


###Test for correct file value of bwa sampe
#GREP command result
my $grepResult=`grep -c "XT:A:U" RC3.BWASAMPE.sam`;
chomp $grepResult;
is($grepResult,1704,'Test for the result of bwa sampe');


####################################################################################################################
#####Test for bwa Samse
########################################################################################################
my $fastqFile="RC3.REPAIRING.fastq";
my $originalFastqFile=$expectedData."/".$fastqFile;

$lnCmd="ln -s $originalFastqFile .";
system($lnCmd) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile for test with the command $lnCmd \n$!\n");    #The sequences are copied for testing

my $singleSaiFileIn="RC3.BWAALN.sai";
my $samseFileOut="RC3.BWASAMSE.sam";

is (bwa::bwaAln($fastaRef,$fastqFile,$singleSaiFileIn,$optionsHachees),'1',"bwa::bwaAln... Test for bwa Aln running for single");
is (bwa::bwaSamse($samseFileOut,$fastaRef,$singleSaiFileIn,$fastqFile,"RC3"),'1',"bwa::bwaSamse... Test for bwa samse running");

####Verify if output are correct for samse
@expectedOutput=('./bwa_TEST_log.e','./bwa_TEST_log.o','./individuSoft.txt','./RC3_1.BWAALN.sai','./RC3_1.REPAIRING.fastq','./RC3_2.BWAALN.sai','./RC3_2.REPAIRING.fastq','./RC3.BWAALN.sai','./RC3.BWASAMPE.sam','./RC3.BWASAMSE.sam','./RC3.REPAIRING.fastq','./Reference.fasta','./Reference.fasta.amb','./Reference.fasta.ann','./Reference.fasta.bwt','./Reference.fasta.pac','./Reference.fasta.sa');
@observedOutput=toolbox::readDir(".");

is_deeply(\@expectedOutput,@observedOutput,'Test for bwa samse output files');

####Test for correct file value of bwa samse - 
$grepResult= `grep -c "XT:A:U" RC3.BWASAMSE.sam`;
chomp $grepResult;
is($grepResult,1,'bwa::bwaSamse... Test for the content result of bwa samse');


########################################################################################################
#####Test for bwa mem single
########################################################################################################

##Running test
is (bwa::bwaMem($samseFileOut,$fastaRef,$fastqFile1,"","truc"),'1',"bwa::bwaMem... Test for bwa mem running single");

###Verify if output are correct for mem single
@expectedOutput=('./bwa_TEST_log.e','./bwa_TEST_log.o','./individuSoft.txt','./RC3_1.REPAIRING.fastq','./RC3_2.REPAIRING.fastq','./RC3.BWASAMPE.sam','./Reference.fasta','./Reference.fasta.amb','./Reference.fasta.ann','./Reference.fasta.bwt','./Reference.fasta.pac','./Reference.fasta.sa');
my @outPut=toolbox::readDir(".");
is_deeply(@outPut,\@expectedOutput,'bwa::bwaMem... Test for bwa mem single output files');

##Output value test
$grepResult= `grep -c 'LOC' RC3.BWASAMPE.sam`;
chomp $grepResult;
is($grepResult,717,'bwa::bwaMem... Test for the result of bwa mem single');

exit;

### ADD BWA MEM PAIRED TEST


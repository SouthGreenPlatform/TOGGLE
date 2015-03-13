#!/usr/bin/perl -w

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot
#
###################################################################################################################################

use strict;

#Will test if bwa works correctly
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../Modules/);
use Data::Dumper;

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"bwa\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf bwa_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

#########################################
#Remove the files and directory created by the previous test
#########################################
$cleaningCommand="rm -Rf ../DATA-TEST/bwaTestDir";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous test dir with the command $cleaningCommand \n$!\n");

########################################
#Creation of test directory
########################################
my $testingDir="../DATA-TEST/bwaTestDir";
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

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

toolbox::readFileConf("software.config.txt");

########################################
#Picking up data for tests
########################################
my $originalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $originalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference for test with the command $refCopyCom \n$!\n");
  #Now we have a ref to be tested
my $originalFastqFile1="../DATA/expectedData/RC3_1.REPAIRING.fastq";
my $originalFastqFile2="../DATA/expectedData/RC3_2.REPAIRING.fastq";
my $fastqFile1="$testingDir/RC3_1.REPAIRING.fastq";
my $fastqFile2="$testingDir/RC3_2.REPAIRING.fastq";
my $seqCopyCom1="cp $originalFastqFile1 $fastqFile1";
my $seqCopyCom2="cp $originalFastqFile2 $fastqFile2";
system($seqCopyCom1) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile1 for test with the command $seqCopyCom1 \n$!\n");    #The sequences are copied for testing
system($seqCopyCom2) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile2 for test with the command $seqCopyCom2 \n$!\n");    #The sequences are copied for testing
my $forwardSaiFileIn="$testingDir/RC3_1.BWAALN.sai";
my $reverseSaiFileIn="$testingDir/RC3_2.BWAALN.sai";
my $sampeFileOut="$testingDir/RC3.BWASAMPE.sam";

#######################################################################################################
####Test for bwa index running
#######################################################################################################
##Test for running
my $optionsHachees=$configInfos->{'BWA index'};
is(bwa::bwaIndex($fastaRef,$optionsHachees),1,'Test for bwaIndex running');
###Verify if output are correct for bwa index
my @expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
my @outPut=toolbox::readDir($testingDir);

is_deeply(@outPut,\@expectedOutput,'Test for the output files produced by bwa index');


###Test for correct file value of bwa index using a md5sum file control -  work through the different bwa versions
my $expectedMD5sum='b86728bb71903f8641530e61e9687b59  ../DATA-TEST/bwaTestDir/Reference.fasta.amb
a51b6a5152f51b13833a40fe609474ea  ../DATA-TEST/bwaTestDir/Reference.fasta.ann
e4fdc0af9540ee8365e7e324fc5c0cc3  ../DATA-TEST/bwaTestDir/Reference.fasta.bwt
03454e7242900c436d9d7126f492e4d5  ../DATA-TEST/bwaTestDir/Reference.fasta.pac
9243bf066de0cc18aa0d3813f174cae8  ../DATA-TEST/bwaTestDir/Reference.fasta.sa
'; #Expected values for the files produced by the bwa index
my $observedMD5sum=`md5sum ../DATA-TEST/bwaTestDir/Reference.fasta.*`;#md5sum values observed for the current files produced
is($observedMD5sum,$expectedMD5sum,'Test for the content of the bwa index output');


#######################################################################################################
###Test for bwa Aln running
#######################################################################################################

$optionsHachees=$configInfos->{'BWA aln'};
is (bwa::bwaAln($fastaRef,$fastqFile1,$forwardSaiFileIn,$optionsHachees),'1',"Test for bwa Aln running for forward");
is (bwa::bwaAln($fastaRef,$fastqFile2,$reverseSaiFileIn,$optionsHachees),'1',"Test for bwa Aln running for reverse");

###Verify if output are correct for bwa Aln
@expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
@outPut=toolbox::readDir($testingDir);

is_deeply(@outPut,\@expectedOutput,'Test for bwa aln output files');

###Test for correct file value of bwa aln using a md5sum - BE CAREFUL, the sum changes based on the version of BWA!!
TODO: {
    local $TODO = "The file structure depends of the version of BWA in use. Here, tested for bwa 0.7.9a";
    
    $expectedMD5sum ='1be54c4d8d37870cbd78d93cee30b26f  ../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai
53ae4174a5bd3cdf2958a59614cf0eb1  ../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai
';
    $observedMD5sum=`md5sum ../DATA-TEST/bwaTestDir/*.sai`;#md5sum values observed for the current files produced
    is($observedMD5sum,$expectedMD5sum,'Test for the content of the bwa aln output');
    
    
}

########################################################################################################
####Test for bwa sampe
#######################################################################################################
is(bwa::bwaSampe($sampeFileOut,$fastaRef,$forwardSaiFileIn,$reverseSaiFileIn,$fastqFile1,$fastqFile2,"RC3"),'1',"Test for bwa sampe running");
####Verify if output are correct for sampe
@expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
@outPut=toolbox::readDir($testingDir);
is_deeply(@outPut,\@expectedOutput,'Test for Output file ok for bwa sampe');


###Test for correct file value of bwa sampe
#GREP command result
my $grepResult=`grep -c "XT:A:U" ../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam`;
chomp $grepResult;
is($grepResult,1704,'Test for the result of bwa sampe');

######Remove the file created
#unlink ('../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam');

####################################################################################################################
#####Test for bwa Samse
########################################################################################################
my $originalFastqFile3="../DATA/expectedData/RC3.REPAIRING.fastq";
my $fastqFile3="$testingDir/RC3.REPAIRING.fastq";
my $seqCopyCom3="cp $originalFastqFile3 $fastqFile3";
system($seqCopyCom3) and die ("ERROR: $0 : Cannot copy the Fastq file $fastqFile3 for test with the command $seqCopyCom3 \n$!\n");    #The sequences are copied for testing

my $singleSaiFileIn="$testingDir/RC3.BWAALN.sai";
my $samseFileOut="$testingDir/RC3.BWASAMSE.sam";

is (bwa::bwaAln($fastaRef,$fastqFile3,$singleSaiFileIn,$optionsHachees),'1',"Test for bwa Aln running for single");
is (bwa::bwaSamse($samseFileOut,$fastaRef,$singleSaiFileIn,$fastqFile3,"RC3"),'1',"Test for bwa samse running");

####Verify if output are correct for samse
@expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/RC3.BWASAMSE.sam','../DATA-TEST/bwaTestDir/RC3.REPAIRING.fastq','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
@outPut=toolbox::readDir($testingDir);

is_deeply(\@expectedOutput,@outPut,'Test for bwa samse output files');

####Test for correct file value of bwa samse - 
$grepResult= `grep -c "XT:A:U" ../DATA-TEST/bwaTestDir/RC3.BWASAMSE.sam`;
chomp $grepResult;
is($grepResult,1,'Test for the content result of bwa samse');

######Remove the file created
#unlink('../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai');

exit;

########################################################################################################
#####Test for bwa mem single
########################################################################################################

##Running test
is (bwa::bwaMem($samseFileOut,$fastaRef,$fastqFile1,"","truc"),'1',"Test for bwa mem running single");

###Verify if output are correct for mem single
@expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
@outPut=toolbox::readDir($testingDir);
is_deeply(@outPut,\@expectedOutput,'Test for bwa mem single output files');

##Output value test

$grepResult= `grep -c 'LOC' ../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam`;
chomp $grepResult;
is($grepResult,717,'Test for the result of bwa mem single');

####Remove the file created
unlink('../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/RC3_1.BWAALN.sai','../DATA-TEST/bwaTestDir/RC3_2.BWAALN.sai');

###################################################################################################################
####Test for bwa mem paired
#######################################################################################################

## Running test
#is (bwa::bwaMem($sampeFileOut,$fastaRef,$fastqFile1,$fastqFile2,"truc"),'1',"Test for bwa mem running paired");

#### Verify if output are correct for mem paired
#@expectedOutput=('../DATA-TEST/bwaTestDir/RC3_1.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3_2.REPAIRING.fastq','../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam','../DATA-TEST/bwaTestDir/Reference.fasta','../DATA-TEST/bwaTestDir/Reference.fasta.amb','../DATA-TEST/bwaTestDir/Reference.fasta.ann','../DATA-TEST/bwaTestDir/Reference.fasta.bwt','../DATA-TEST/bwaTestDir/Reference.fasta.pac','../DATA-TEST/bwaTestDir/Reference.fasta.sa');
#@outPut=toolbox::readDir($testingDir);

#print Dumper(\@outPut);

#is_deeply(@outPut,\@expectedOutput,'Test for bwa mem paired output files');

##Output value test
#$grepResult= `grep -c 'LOC' ../DATA-TEST/bwaTestDir/RC3.BWASAMPE.sam`;
#homp $grepResult;
#is($grepResult,1436,'Test for the result of bwa mem paired');

exit;




#/usr/bin/perl

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################


#Will test if the modue fastqc work correctly

use strict;
use warnings;

use Test::More  'no_plan';
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);

########################################
#use of fastqc module ok
########################################
use_ok('toolbox') or exit;
use_ok('fastqc') or exit;
can_ok( 'fastqc','execution');
can_ok('fastqc','parse');

use toolbox;
use fastqc;
toolbox::readFileConf("software.config.txt");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"fastqc\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCommand\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf fastqc_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCommand \n$!\n");

#########################################
#Remove the files and directory created by the previous test
#########################################
$cleaningCommand="rm -Rf ../DATA-TEST/fastqcTestDir";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous test dir with the command $cleaningCommand \n$!\n");

########################################
#Creation of test directory
########################################
my $testingDir="../DATA-TEST/fastqcTestDir";
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

########################################
#Input files
########################################
my $originalFastqcFile = "../DATA/expectedData/RC3_2.fastq";     # fastq file 
my $fastqcFile = "$testingDir/RC3_2.fastq";                             # fastq file for test
my $fastqcFileCopyCom = "cp $originalFastqcFile $fastqcFile";           # command to copy the original fastq file into the test directory
system ($fastqcFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $fastqcFileCopyCom\n$!\n");    # RUN the copy command

##########################################
#Fastqc exec test
##########################################
my $executionTest = "$testingDir/execution";
$makeDirCom = "mkdir $executionTest";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

is(fastqc::execution($fastqcFile,$executionTest),1,'Test for fastqc::execution');     # test if fastqc::execution works

my @expectedOutput = ('../DATA-TEST/fastqcTestDir/execution/RC3_2_fastqc.zip',
                      '',
                      '../DATA-TEST/fastqcTestDir/execution/RC3_2_fastqc:',
                      'fastqc_data.txt',
                      'fastqc_report.html',
                      'Icons',
                      'Images',
                      'summary.txt');

my @observedOutput = toolbox::readDir($executionTest);
##DEBUG print "ICI :\n"; print Dumper(@observedOutput);
is_deeply(@observedOutput,\@expectedOutput,'Test for output file of fastqc::execution');        # test if the observed output of fastqc::execution is ok

#########################################
#Fastqc  parse test
#########################################   
my $expectedOutput= {
            'Filename' => 'RC3_2.fastq',
            'File type' => 'Conventional base calls',
            'Encoding' => 'Sanger / Illumina 1.9',
            'Total Sequences' => '1000',    
            'Filtered Sequences' => '0',       
            'Sequence length' => '38-101',      
            '%GC' => '42'};

my $observedoutput=fastqc::parse($executionTest);      # test if fastqc::parse works and will give the observed output in the same time
##DEBUG print Dumper($observedoutput);
is_deeply($observedoutput,$expectedOutput,'Test for fastqc::parse and it\'s output');      # test if the observed output of fastqc::parse is ok

exit;
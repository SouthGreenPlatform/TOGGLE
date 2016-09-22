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

#Will test if toolbox.pm works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Data::Dumper;
use lib qw(../Modules/);

use localConfig;


########################################
#Test of the use of toolbox modules
########################################

use_ok('toolbox');
can_ok('toolbox','exportLog');
can_ok('toolbox','checkFile');
can_ok('toolbox','readFile');   #remove test for return 0
can_ok('toolbox','writeFile');  #remove test for return 0
can_ok('toolbox','sizeFile');
can_ok('toolbox','existsFile');
can_ok('toolbox','existsDir');
can_ok('toolbox','makeDir');
can_ok('toolbox','readDir');
#can_ok('toolbox','readDir2');
can_ok('toolbox','readFileConf');
can_ok('toolbox','extractPath');
can_ok('toolbox','extractOptions');
can_ok('toolbox','extractName');
can_ok('toolbox','run');
can_ok('toolbox','checkNumberLines');
can_ok('toolbox','checkFormatFastq');
can_ok('toolbox','addInfoHeader');
can_ok('toolbox','checkSamOrBamFormat');
can_ok('toolbox','changeDirectoryArbo');
can_ok('toolbox','extractHashSoft');
can_ok('toolbox','checkInitialDirContent');
can_ok('toolbox','checkVcfFormat');
can_ok('toolbox','checkFormatFasta');
can_ok('toolbox','relativeToAbsolutePath');



use toolbox;

my $expectedData="../../DATA/expectedData/";
my $configFile="../../SNPdiscoveryPaired.config.txt";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/toolboxTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

my $makeDirCmd = "mkdir pairingDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"toolbox\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf toolbox_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");




#######################################
#charge needed test files into variables
#######################################

#Fastq file
my $fastqFile=$expectedData."RC3_1.fastq";

#Sam file
my $samFile=$expectedData."RC3.BWASAMPE.sam";

#Bam file
my $bamFile=$expectedData."RC3.PICARDTOOLSSORT.bam";

#VCF file
my $vcfFile=$expectedData."GATKHAPLOTYPECALLER.vcf";

#VCF file non readable
my $chmodVcfFile=$expectedData."test-nonreadrigth.vcf";

#File empty
my $emptyFile=$expectedData."empty-file.vcf";
my $createFileCmd="touch $emptyFile";
system($createFileCmd) and die ("ERROR: $0 : Cannot create the empty file with the command $createFileCmd\n$!\n"); +

#Fasta files
my $reference=$expectedData."correctReference.fasta";
my $wrongFasta=$expectedData."wrongReference.fasta";


########################################
# toolbox::exportLog tests
########################################
# Test if individuSoft.txt has been created
my $file="individuSoft.txt";
my $got=(-e $file)?1:0;
my $expectedBool=1;
is($got,$expectedBool,"toolbox::exportLog - exist $file?");

# Test if the log is written in toolbox_TEST_log.o
# Test if toolbox_TEST_log.o has been created
my $expected="INFO:  toolbox_test.t\n";
toolbox::exportLog($expected,1);

my $file_log="toolbox_TEST_log.o";
$got=(-e $file_log)?1:0;
is($got,$expectedBool,"toolbox::exportLog - exist $file_log?");

$got=`head -n1 $file_log`;
is($got,$expected,"toolbox::exportLog - Log in $file_log");


# Test if the log is written in toolbox_TEST_log.e
# Test if toolbox_TEST_log.e has been created
$expected="WARNING:  toolbox_test.t\n";
toolbox::exportLog($expected,2);
my $file_error="toolbox_TEST_log.e";
$got=(-e $file_error)?1:0;
is($got,$expectedBool,"toolbox::exportLog - exist $file_error?");

$got=`head -n1 $file_error`;
is($got,$expected,"toolbox::exportLog - Log in $file_error");


########################################
# TEST FILES
########################################
########################################
#toolbox::checkFile
########################################
is (toolbox::checkFile($configFile),'1','toolbox::checkFile - OK');

########################################
#toolbox::existsFile
########################################
is (toolbox::existsFile($configFile),'1','toolbox::existsFile - return 1');
is (toolbox::existsFile('beurk.txt',0),'0','toolbox::existsFile - return 0');

########################################
#toolbox::readFile
########################################
is (toolbox::readFile($configFile),'1','toolbox::readFile - return 1');

#my $chmodCmd="chmod -r $chmodVcfFile";
#system($chmodCmd) and die ("\nCannot change the right of the vcf file for test:$!\nAborting\n");
#is (toolbox::readFile($chmodVcfFile),'0','toolbox::readFile - return 0');

########################################
#toolbox::writeFile test TODO to verify
########################################
is (toolbox::writeFile($configFile),'1','toolbox::writeFile - return 1');
#
#my $chmodCmd="chmod -w $chmodVcfFile";
#system($chmodCmd) and die ("\nCannot change the right of the vcf file for test:$!\nAborting\n");
#is (toolbox::writeFile($chmodVcfFile),'0','toolbox::writeFile - return 0');

########################################
#toolbox::sizeFile 
########################################
ok (toolbox::sizeFile($configFile) == 1,'toolbox::sizeFile - return 1');
ok (toolbox::sizeFile($emptyFile) == 0,'toolbox::sizeFile - return 0');


########################################
#Directory test
########################################
########################################
#toolbox::existsDir 
########################################
is (toolbox::existsDir('.'),'1','toolbox::existsDir - return 1');
is (toolbox::existsDir('beurk',0),'0','toolbox::existsDir - return 0');

########################################
#toolbox::makeDir 
########################################
is (toolbox::makeDir('test_dir'),'1','toolbox::makedir - created directory OK');
is (toolbox::existsDir('test_dir'),'1','toolbox::existsDir - return 1');
system ("rm -Rf test_dir") and die ("\nCannot remove the test_dir for test:$!\nAborting\n");;

########################################
#toolbox::readDir test with a directory name as argument
########################################
my $listCom = `ls $expectedData/*`;
chomp $listCom;
my @listExpected = split /\n/, $listCom;
my @listObserved = toolbox::readDir($expectedData);
is_deeply(\@listExpected,@listObserved,'toolbox::readDir - just directory');

########################################
#toolbox::readDir test with a directory name and a format as arguments
########################################
$listCom = `ls $expectedData/*fastq`;
chomp $listCom;
@listExpected = split /\n/, $listCom;
@listObserved = toolbox::readDir($expectedData,'fastq');
is_deeply(\@listExpected,@listObserved,'toolbox::readDir - just fastq files');

########################################
#toolbox::readDir2 test with a directory name
########################################
#$listCom = `ls ../../DATA/expectedData/*`;
#chomp $listCom;
#@listExpected = split /\n/, $listCom;
#@listObserved = toolbox::readDir2('../../DATA/expectedData');
#is_deeply(\@listExpected,@listObserved,'toolbox::readDir2... just a directory');
#
########################################
#toolbox::readDir2 test with a directory name and a part of filename as arguments
########################################
#$listCom = `ls ../../DATA/expectedData/RC*`;
#chomp $listCom;
#@listExpected = split /\n/, $listCom;
#@listObserved = toolbox::readDir2('../../DATA/expectedData','RC');
#is_deeply(\@listExpected,@listObserved,'toolbox::readDir2... just on filename');


########################################
#Path test
########################################
########################################
#toolbox::extractPath
########################################
my @expectedList=("toto","/home/username/");
my @testList=toolbox::extractPath('/home/username/toto');
is_deeply (\@expectedList,\@testList,'toolbox::extractPath - OK');

########################################
#toolbox::extractName Extract name from path test
########################################
is (toolbox::extractName($samFile),'RC3.BWASAMPE.sam','toolbox::extractName - OK');


########################################
##Sequence count test
########################################
########################################
#toolbox::checkNumberLines test with a fastq file
########################################
my $count = (toolbox::checkNumberLines($fastqFile))/4;
is($count,'1000',"toolbox::checkNumberLines - OK");


########################################
#File Format test
########################################
########################################
#toolbox::checkFormatFastq
########################################
is(toolbox::checkFormatFastq($fastqFile),'1', 'toolbox::checkFormatFastq - OK');

########################################
#toolbox::checkSamOrBamFormat
########################################
is (toolbox::checkSamOrBamFormat($samFile),'1', 'toolbox::checkSamOrBamFormat - sam format');
is (toolbox::checkSamOrBamFormat($bamFile),'2', 'toolbox::checkSamOrBamFormat - bam format');

########################################
#toolbox::checkFormatFasta
########################################
is (toolbox::checkFormatFasta($reference),'1','toolbox::checkFormatFasta - Format Ok');
is (toolbox::checkFormatFasta($wrongFasta),'0','toolbox::checkFormatFasta - Format not Ok, warnings send');


########################################
#Config file test
########################################

#Reading of config file infos
my $configInfos = toolbox::readFileConf($configFile);
#checking if $configInfos exists
is (ref($configInfos),'HASH','toolbox::readFileConf - the reference returned is a HASH');

#checking how many software configs
my @listOfSoftwares=keys (%$configInfos);#Soft are BWA and samtoolsView
##DEBUG foreach my $key(@listOfSoftwares){print "$key\n";}
my $numberOfSoft= scalar (@listOfSoftwares); #expecting 17
my $command='grep "^\\\$" '.$configFile.' -c';
##DEBUG print "DEBUG: $0: Number of softwares returned by grep command: $command\n";
my $numberOfSoftExpected=`$command`;
chomp $numberOfSoftExpected;
is ($numberOfSoft, $numberOfSoftExpected, 'toolbox::readFileConf - the number of software to configure');

#checking for info retrieval, directly, ie data structure
##DEBUG print "DEBUG: $0: ".$configInfos->{"samToolsView"}{"-F"}." \n";
is ($configInfos->{"samToolsView"}{-f},'0x02','toolbox::readFileConf - samToolsView infos retrieval');
isnt ($configInfos->{"samTooView"}{-F},'0x02','toolbox::readFileConf - samToolsView infos retrieval');

#checking for info extract
my $optionLine=toolbox::extractOptions($configInfos->{"bwaAln"}," ");
is ($optionLine =~ m/-n 5/,'1','toolbox::extractOptions - is'); #Test as an Test form because of randomness of hash sorting, to be sure of controlling the data

$optionLine=toolbox::extractOptions($configInfos->{"bwaln"}," ");
isnt ($optionLine =~ m/-n 5/ && $optionLine =~ m/-e -1/,'1','toolbox::extractOptions - option test'); 



########################################
#toolbox::addInfoHeader test
########################################

#is (toolbox::addInfoHeader($bamFile, 'addInfoHeaderTest'),'1','Test for addInfoHeader');

#my $observedMD5sum = `md5sum $bamFile`; #md5sum values observed for the generated bam

#Check if the two bam are the same
#is_deeply ($observedMD5sum, "1546666d4335961d81254e91951cac6c  ../DATA-TEST/toolboxTestDir/RC3.PICARDTOOLSSORT.bam\n", "Test for addInfoHeader result");


########################################
#toolbox::changeDirectoryArbo test
########################################
my $directory = "./TEST/";
my $newDirectory = toolbox::changeDirectoryArbo($directory,'0');
is_deeply($newDirectory, "./TEST/0_PAIRING_FILES", "toolbox::changeDirectoryArbo");

is(toolbox::changeDirectoryArbo($directory,'8'),undef,"toolbox::changeDirectoryArbo");


########################################
#toolbox::extractHashSoft test
########################################
# Get option for bwa index
my $hashConfig ={
                    "-n" => "5"
                };

my $testHashConfig = toolbox::readFileConf($configFile);
my $softInfos = toolbox::extractHashSoft($testHashConfig, "bwaAln");
cmp_deeply($hashConfig, $softInfos, 'toolbox::extractHashSoft - bwaAln');


########################################
#toolbox::run command test
########################################

#testing rendering ie return 1
my $testCom="date +%D >> log.txt"; # print the date in the log, format MM/DD/YYYY
my $returnValue=toolbox::run($testCom);
ok ($returnValue== 1, 'toolbox::run - return value');

#testing correct behaviour
my $date=`date +%D`; # The previous test will print the date in the log, format MM/DD/YYYY
chomp $date;
my $endOfLog=`tail -n 1 log.txt `; #The last line of log is always "Cmd Done", so pick up the two last and keep the n-1 line
chomp $endOfLog;
ok($date eq $endOfLog,'toolbox::run - command behaviour');

########################################
#toolbox::checkVcfFormat test TODO add test negatif
########################################
is (toolbox::checkVcfFormat($vcfFile),'1','toolbox::checkVcfFormat - vcf file');
#isnt (toolbox::checkVcfFormat($samFile),'1','Test for checkVcfFormat - sam file');

########################################
#toolbox::checkInitialDirContent test TODO add test negatif
########################################
$makeDirCmd = "mkdir initialDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");

my $copyCmd= "cp $fastqFile initialDir/";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFile in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

is (toolbox::checkInitialDirContent('initialDir'),'0','toolbox::checkInitialDirContent - OK');


########################################
#toolbox::relativeToAbsolutePath
########################################
#is (toolbox::relativeToAbsolutePath('./'),'0','toolbox::relativeToAbsolutePath - OK');

$optionLine=toolbox::relativeToAbsolutePath('./');
print $optionLine;
is ($optionLine =~ m/INFOS/,'1','toolbox::relativeToAbsolutePath - OK');



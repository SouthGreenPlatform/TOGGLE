#/usr/bin/perl


###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################

#Will test if toolbox works correctly
#Modified version by Marilyne version3
#test du 13/08/14

use strict;
use warnings;
use Test::More 'no_plan'; #tests => 19; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Data::Dumper;
use Test::Deep;
use lib qw(../Modules/);

use localConfig;
my $configFile='software.config.txt';


########################################
#use of toolbox module and functions ok
########################################

use_ok('toolbox');
can_ok('toolbox','exportLog');
can_ok('toolbox','checkFile');
can_ok('toolbox','readFile');
can_ok('toolbox','writeFile');
can_ok('toolbox','sizeFile');
can_ok('toolbox','existsFile');
can_ok('toolbox','existsDir');
can_ok('toolbox','makeDir');
can_ok('toolbox','readDir');
can_ok('toolbox','readDir2');
can_ok('toolbox','readFileConf');
can_ok('toolbox','extractPath');
can_ok('toolbox','extractOptions');
can_ok('toolbox','extractName');
can_ok('toolbox','run');
can_ok('toolbox','checkFormatFastq');
can_ok('toolbox','addInfoHeader');
can_ok('toolbox','checkSamOrBamFormat');
can_ok('toolbox','changeDirectoryArbo');
can_ok('toolbox','extractHashSoft');
can_ok('toolbox','checkInitialDirContent');
can_ok('toolbox','checkVcfFormat');
can_ok('toolbox','transferDirectoryFromMasterToNode'); # TEST A IMPLEMENTER
can_ok('toolbox','transferDirectoryFromNodeToMaster'); # TEST A IMPLEMENTER

use toolbox;


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



#########################################
#Remove the files and directory created by the previous test
#########################################
my $testingDir="../DATA-TEST/toolboxTestDir";
$cleaningCommand="rm -Rf ../DATA-TEST/$testingDir"; 
system ("rm -Rf $testingDir") and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCommand \n$!\n");


########################################
#Creation of test directory
########################################
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");


#######################################
#Copy needed test files into testing directory
#######################################

#Fastq file
my $originalFastqFile="../DATA/expectedData/RC3_1.fastq";
my $fastqFile="$testingDir/RC3_1.fastq";
my $copyCommand="cp $originalFastqFile $fastqFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the fastq file with the command $copyCommand\n$!\n"); 

#Sam file
my $originalSamFile="../DATA/expectedData/RC3.BWASAMPE.sam";
my $samFile="$testingDir/RC3.BWASAMPE.sam";
$copyCommand="cp $originalSamFile $samFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the sam file with the command $copyCommand\n$!\n"); 

#Bam file
my $originalBamFile="../DATA/expectedData/RC3.PICARDTOOLSSORT.bam";
my $bamFile="$testingDir/RC3.PICARDTOOLSSORT.bam";
$copyCommand="cp $originalBamFile $bamFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the bam file with the command $copyCommand\n$!\n"); 

#VCF file
my $OriginalVcfFile="../DATA/expectedData/GATKHAPLOTYPECALLER.vcf";
my $vcfFile="$testingDir/GATKHAPLOTYPECALLER.vcf";
$copyCommand="cp $OriginalVcfFile $vcfFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the vcf file with the command $copyCommand\n$!\n"); 

#VCF file non readable
my $chmodVcfFile="$testingDir/test-nonreadrigth.vcf";
$copyCommand="cp $OriginalVcfFile $chmodVcfFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the vcf file with the command $copyCommand\n$!\n"); 

#File empty
my $emptyFile="$testingDir/empty-file.vcf";
my $createFileCommand="touch $emptyFile";
system($createFileCommand) and die ("ERROR: $0 : Cannot create the empty file with the command $createFileCommand\n$!\n"); +

#Fasta files
my $originalReference = "../DATA/expectedData/correctReference.fasta";
my $reference = "$testingDir/correctReference.fasta";
$copyCommand=" cp $originalReference $reference";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the fasta file with the command $copyCommand\n$!\n");
my $originalWrongFasta = "../DATA/expectedData/wrongReference.fasta";
my $wrongFasta = "$testingDir/wrongReference.fasta";
$copyCommand=" cp $originalWrongFasta $wrongFasta";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the fasta file with the command $copyCommand\n$!\n");


########################################
# exportLog tests
########################################
# Test if individuSoft.txt has been created
my $file="individuSoft.txt";
my $got=(-e $file)?1:0;
my $expectedBool=1;
is($got,$expectedBool,"Test for exportLog... exist $file?");

# Test if the log is written in toolbox_TEST_log.o
# Test if toolbox_TEST_log.o has been created
my $expected="INFO:  toolbox_test.t\n";
toolbox::exportLog($expected,1);

my $file_log="toolbox_TEST_log.o";
$got=(-e $file_log)?1:0;
is($got,$expectedBool,"Test for exportLog... exist $file_log?");

$got=`head -n1 $file_log`;
is($got,$expected,"Test for exportLog... Log in $file_log");


# Test if the log is written in toolbox_TEST_log.e
# Test if toolbox_TEST_log.e has been created
$expected="WARNING:  toolbox_test.t\n";
toolbox::exportLog($expected,2);
my $file_error="toolbox_TEST_log.e";
$got=(-e $file_error)?1:0;
is($got,$expectedBool,"Test for exportLog... exist $file_error?");

$got=`head -n1 $file_error`;
is($got,$expected,"Test for exportLog... Log in $file_error");


########################################
#File infos tests
########################################

#checkFile
is (toolbox::checkFile($configFile),'1','Test for checkFile');
#is (toolbox::checkFile('beurk.txt'),'0','Test for checkFile'); Gestion des erreurs

#existsFile
is (toolbox::existsFile($configFile),'1','Test for existsFile... return 1');
is (toolbox::existsFile('beurk.txt',0),'0','Test for existsFile... return 0');

#readFile test
is (toolbox::readFile($configFile),'1','Test for readFile... return 1');

my $chmodCommand="chmod -r $chmodVcfFile";
system($chmodCommand) and die ("\nCannot change the right of the vcf file for test:$!\nAborting\n");
is (toolbox::readFile($chmodVcfFile),'0','Test for readFile... return 0');

#writeFile test TODO to verify
is (toolbox::writeFile($configFile),'1','Test for writeFile... return 1');

$chmodCommand="chmod -w $chmodVcfFile";
system($chmodCommand) and die ("\nCannot change the right of the vcf file for test:$!\nAborting\n");
is (toolbox::writeFile($chmodVcfFile),'0','Test for writeFile... return 0');

#sizeFile test
ok (toolbox::sizeFile($configFile) == 1,'Test for sizeFile... return 1');
ok (toolbox::sizeFile($emptyFile) == 0,'Test for sizeFile... return 0');

########################################
#Directory test
########################################

#existsDir
is (toolbox::existsDir($testingDir),'1','Test for existsDir... return 1');
is (toolbox::existsDir('beurk',0),'0','Test for existsDir... return 0');

#makeDir
is (toolbox::makeDir($testingDir.'/test_dir'),'1','Test for makedir');
is (toolbox::existsDir($testingDir.'/test_dir'),'1','Test for existsDir... return 1');
system ("rm -Rf ".$testingDir."/test_dir") and die ("\nCannot remove the test_dir for test:$!\nAborting\n");;

#readDir test with a directory name as argument
my $listCom = `ls ../DATA/expectedData/*`;
chomp $listCom;
my @listExpected = split /\n/, $listCom;
my @listObserved = toolbox::readDir('../DATA/expectedData');
is_deeply(\@listExpected,@listObserved,'Test for readDir... just directory');

#readDir test with a directory name and a format as arguments
$listCom = `ls ../DATA/expectedData/*fastq`;
chomp $listCom;
@listExpected = split /\n/, $listCom;
@listObserved = toolbox::readDir('../DATA/expectedData','fastq');
is_deeply(\@listExpected,@listObserved,'Test for readDir... just fastq files');

#readDir2 test with a directory name
$listCom = `ls ../DATA/expectedData/*`;
chomp $listCom;
@listExpected = split /\n/, $listCom;
@listObserved = toolbox::readDir2('../DATA/expectedData');
is_deeply(\@listExpected,@listObserved,'Test for readDir2... just a directory');

#readDir2 test with a directory name and a part of filename as arguments
$listCom = `ls ../DATA/expectedData/RC*`;
chomp $listCom;
@listExpected = split /\n/, $listCom;
@listObserved = toolbox::readDir2('../DATA/expectedData','RC');
is_deeply(\@listExpected,@listObserved,'Test for readDir2... just on filename');

########################################
#Path test
########################################

#extractPath test
my @expectedList=("toto","/home/username/");
my @testList=toolbox::extractPath('/home/username/toto');
is_deeply (\@expectedList,\@testList,'Test for extractPath');

# Extract name from path test
is (toolbox::extractName($samFile),'RC3.BWASAMPE','Test for extractName');

########################################
#File Format test
########################################

#checkFormatFastq
is(toolbox::checkFormatFastq($fastqFile),'1', 'Test for checkFormatFastq');

#checkSamOrBamFormat
is (toolbox::checkSamOrBamFormat($samFile),'1', 'Test for checkSamOrBamFormat... sam format');
is (toolbox::checkSamOrBamFormat($bamFile),'2', 'Test for checkSamOrBamFormat... bam format');

#dnaFastaFormatValidator
is (toolbox::checkFormatFasta($reference),'1','Test for checkFormatFasta... Format Ok');
is (toolbox::checkFormatFasta($wrongFasta),'0','Test for checkFormatFasta... Format not Ok, warnings send');

########################################
#Config file test
########################################

#Reading of config file infos
toolbox::readFileConf($configFile);

#checking if $configInfos exists
is (ref($configInfos),'HASH','Test for readFileConf... the reference returned is a HASH');

#checking how many software configs
my @listOfSoftwares=keys (%$configInfos);#Soft are BWA and samtoolsView
##DEBUG foreach my $key(@listOfSoftwares){print "$key\n";}
my $numberOfSoft= scalar (@listOfSoftwares); #expecting 17
my $command='grep "^\\\$" '.$configFile.' -c';
##DEBUG print "DEBUG: $0: Number of softwares returned by grep command: $command\n";
my $numberOfSoftExpected=`$command`;
chomp $numberOfSoftExpected;
ok ($numberOfSoft == $numberOfSoftExpected, 'Test for readFileConf... the number of software to configure');

#checking for info retrieval, directly, ie data structure
is ($configInfos->{"samtools view pair"}{-F},'0x02','Test for readFileConf... samtools view infos retrieval');
isnt ($configInfos->{"samToo view pair"}{-F},'0x02','Test for readFileConf... samToo view infos retrieval');

#checking for info extract
my $optionLine=toolbox::extractOptions($configInfos->{"BWA aln"}," ");
is ($optionLine =~ m/-n 5/ && $optionLine =~ m/-e -1/,'1','Test for extractOptions... is'); #Test as an Test form because of randomness of hash sorting, to be sure of controlling the data

$optionLine=toolbox::extractOptions($configInfos->{"BWA ln"}," ");
isnt ($optionLine =~ m/-n 5/ && $optionLine =~ m/-e -1/,'1','Test for extractOptions... option test'); 

$optionLine=toolbox::extractOptions($configInfos->{"BwA aln"}," ");
isnt ($optionLine =~ m/-n 5/ && $optionLine =~ m/-e -1/,'1','Test for extractOptions... case test');


########################################
#addInfoHeader test
########################################

#is (toolbox::addInfoHeader($bamFile, 'addInfoHeaderTest'),'1','Test for addInfoHeader');

#my $observedMD5sum = `md5sum $bamFile`; #md5sum values observed for the generated bam

#Check if the two bam are the same
#is_deeply ($observedMD5sum, "1546666d4335961d81254e91951cac6c  ../DATA-TEST/toolboxTestDir/RC3.PICARDTOOLSSORT.bam\n", "Test for addInfoHeader result");


########################################
#changeDirectoryArbo test
########################################
my $directory = "./TEST/";
my $newDirectory = toolbox::changeDirectoryArbo($directory,'0');
is_deeply($newDirectory, "./TEST/0_PAIRING_FILES", "Test for changeDirectoryArbo");

is(toolbox::changeDirectoryArbo($directory,'8'),undef,"Test for changeDirectoryArbo");

########################################
#extractHashSoft test
########################################

# Get option for bwa index
my $hashConfig ={
                    "-a" => "is"
                };

my $testHashConfig = toolbox::readFileConf($configFile);
my $softInfos = toolbox::extractHashSoft($testHashConfig, "BWA index");
cmp_deeply($hashConfig, $softInfos, 'Test for extractHashSoft... BWA');

# Get no option for picardTools createSequenceDictionary
$hashConfig =   {
                    " " => " "
                };

$testHashConfig = toolbox::readFileConf($configFile);
$softInfos = toolbox::extractHashSoft($testHashConfig, "picardTools createSequenceDictionary");
cmp_deeply($hashConfig, $softInfos, 'Test for extractHashSoftPicard');
    
# Get option for a tool that doen't exist
$testHashConfig = toolbox::readFileConf($configFile);
$softInfos = toolbox::extractHashSoft($testHashConfig, "picarTools sortsam pair");
cmp_deeply(undef, $softInfos, 'Test for extractHashSoft... picarTools sortsam pair');

########################################
#Run command test
########################################

#testing rendering ie return 1
my $testCom="date +%D >> log.txt"; # print the date in the log, format MM/DD/YYYY
my $returnValue=toolbox::run($testCom);
ok ($returnValue== 1, 'Test for toolbox::run return value');

#testing correct behaviour
my $date=`date +%D`; # The previous test will print the date in the log, format MM/DD/YYYY
chomp $date;
my $endOfLog=`tail -n 1 log.txt `; #The last line of log is always "Command Done", so pick up the two last and keep the n-1 line
chomp $endOfLog;
ok($date eq $endOfLog,'Test for toolbox::run command behaviour');

########################################
#checkVcfFormat test TODO add test negatif
########################################
is (toolbox::checkVcfFormat($vcfFile),'1','Test for checkVcfFormat... vcf file');
#isnt (toolbox::checkVcfFormat($samFile),'1','Test for checkVcfFormat... sam file');





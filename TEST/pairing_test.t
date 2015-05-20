#/usr/bin/perl -w

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################

#Will test if pairing.pm works correctly

use strict;
use warnings;
use Test::More  'no_plan';
use Data::Dumper;

use lib qw(../Modules/);

my $configFile='software.config.txt';

########################################
#use of pairing module ok
########################################

use_ok('pairing');

use pairing;

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"pairing\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCommand\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf pairing_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
#######################################
my $testingDir="../DATA-TEST/pairingTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom: \n$!\n");


########################################
#Input files
########################################
my $originalFastqcFile = "../DATA/expectedData/RC*_1.fastq";     # fastq file 
my $fastqcFileCopyCom = "cp $originalFastqcFile $testingDir";           # command to copy the original fastq file into the test directory
system ($fastqcFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $fastqcFileCopyCom\n$!\n");    # RUN the copy command

$originalFastqcFile = "../DATA/expectedData/RC*_2.fastq";     # fastq file 
$fastqcFileCopyCom = "cp $originalFastqcFile $testingDir";           # command to copy the original fastq file into the test directory
system ($fastqcFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $fastqcFileCopyCom\n$!\n");    # RUN the copy command


########################################
#extractName
########################################
my $expectName1=("RC3_1");
my $expectRG1=("RC3");

my ($obsName1, $obsRG1)=pairing::extractName('$testingDir/RC3_1.fastq');
is_deeply($obsName1,$expectName1,'Test for pairing::extractName... individu RC3_1');
is_deeply($obsRG1,$expectRG1,'Test for pairing::extractName... RG RC3');

my $expectName2=("RC3_2");
my $expectRG2=("RC3");

my ($obsName2, $obsRG2)=pairing::extractName('$testingDir/RC3_2.fastq');
is_deeply($obsName2,$expectName2,'Test for pairing::extractName... individu RC3_2');
is_deeply($obsRG2,$expectRG2,'Test for pairing::extractName... RG RC3');



########################################
#pairRecognition
########################################
my $expectedOutput={
          '@H3:C39R6ACXX:3:1101:1215:1877' => {
                                                         'ReadGroup' => 'RC1',
                                                         'forward' => $testingDir.'/RC1_1.fastq',
                                                         'reverse' => $testingDir.'/RC1_2.fastq'
                                                       },
          '@H2:C381HACXX:5:1101:1359:1908' => {
                                                 'ReadGroup' => 'RC3',
                                                 'forward' => $testingDir.'/RC3_1.fastq',
                                                 'reverse' => $testingDir.'/RC3_2.fastq'
                                               },
          '@H3:C39R6ACXX:3:1101:1192:1848' => {
                                                   'ReadGroup' => 'RC2',
                                                   'forward' => $testingDir.'/RC2_1.fastq',
                                                } 
        };

my $observedoutput=pairing::pairRecognition($testingDir);
##DEBUG print "pairRecognition Expected :\n"; print Dumper ($expectedOutput);print "pairRecognition Observed:\n"; print Dumper ($observedoutput);
is_deeply($expectedOutput,$observedoutput,'Test for pairing::pairRecognition');



#########################################
##createDirPerCouple
#########################################
my $checkValue3=pairing::createDirPerCouple(pairing::pairRecognition($testingDir),$testingDir);
is ($checkValue3,1,'Test for pairing::createDirPerCouple... running');

# Filetree expected
my $expectedFileTree = 
        [
            '../DATA-TEST/pairingTestDir/RC1:',
            'RC1_1.fastq',
            'RC1_2.fastq',
            '',
            '../DATA-TEST/pairingTestDir/RC2:',
            'RC2_1.fastq',
            '',
            '../DATA-TEST/pairingTestDir/RC3:',
            'RC3_1.fastq',
            'RC3_2.fastq'
        ];

my $observedFileTree=toolbox::readDir($testingDir);
##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($expectedFileTree,$observedFileTree,'Test for pairing::pairRecognition... Filetree created');


########################################
#Input file for repairing
########################################
#$originalFastqcFile = "../DATA/expectedData/RC3Single.fastq";     # fastq file 
#$fastqcFileCopyCom = "cp $originalFastqcFile $testingDir";           # command to copy the original fastq file into the test directory
#system ($fastqcFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $fastqcFileCopyCom\n$!\n");    # RUN the copy command

########################################
#repairing 
########################################
$originalFastqcFile = "../DATA/expectedData/*.CUTADAPT.fastq";     # fastq file 
$fastqcFileCopyCom = "cp $originalFastqcFile $testingDir";           # command to copy the original fastq file into the test directory
system ($fastqcFileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalFastqcFile in the test directory with the command $fastqcFileCopyCom\n$!\n");    # RUN the copy command


#Check if running
my $checkValue=pairing::repairing( $testingDir.'/RC3_1.CUTADAPT.fastq',$testingDir.'/RC3_2.CUTADAPT.fastq',$testingDir);
is ($checkValue,'1','Test for pairing::repairing... running');

#Check if working
my $numberOfLinesObserved=`wc -l $testingDir/RC3_Single/RC3Single.fastq`;
chomp $numberOfLinesObserved;
is ($numberOfLinesObserved,'4 '.$testingDir.'/RC3_Single/RC3Single.fastq','Test for pairing::repairing... single file');

#Check if the files created are the same as planned
my $diffForward=`diff -q $testingDir/RC3_1.REPAIRING.fastq ../DATA/expectedData/RC3_1.REPAIRING.fastq`;
is ($diffForward,'','Test for pairing::repairing... forward file');

#Check if the files created are the same as planned
my $diffReverse=`diff -q $testingDir/RC3_2.REPAIRING.fastq ../DATA/expectedData/RC3_2.REPAIRING.fastq`;
is ($diffReverse,'','Test for pairing::repairing... reverse file');

exit;
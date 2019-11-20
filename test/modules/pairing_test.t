###################################################################################################################################
#
# Copyright 2014-2018 IRD-CIRAD-INRA-ADNid
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






######################################################################################################################################
######################################################################################################################################
## COMMON MODULE TEST HEADER
######################################################################################################################################
######################################################################################################################################

use strict;
use warnings;
use Data::Dumper;

use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;

# Load localConfig if primary test is successful
use_ok('localConfig') or exit;
use localConfig;


########################################
# Extract automatically tool name and sub name list
########################################
my ($toolName,$tmp) = split /_/ , $0;
my $subFile=$toggle."/modules/".$toolName.".pm";
my @sub = `grep "^sub" $subFile`or die ("ERROR: $0 : Cannot extract automatically sub name list by grep command \n$!\n");


########################################
#Automatically module test with use_ok and can_ok
########################################

use_ok($toolName) or exit;
eval "use $toolName";

foreach my $subName (@sub)
{
    chomp ($subName);
    $subName =~ s/sub //;
    can_ok($toolName,$subName);
}

#########################################
#Preparing test directory
#########################################
my $testDir="$toggle/dataTest/$toolName"."TestModule";
my $cmd="rm -Rf $testDir ; mkdir -p $testDir";
system($cmd) and die ("ERROR: $0 : Cannot execute the test directory $testDir ($toolName) with the following cmd $cmd\n$!\n");
chdir $testDir or die ("ERROR: $0 : Cannot go into the test directory $testDir ($toolName) with the chdir cmd \n$!\n");


#########################################
#Creating log file
#########################################
my $logFile=$toolName."_log.o";
my $errorFile=$toolName."_log.e";
system("touch $testDir/$logFile $testDir/$errorFile") and die "\nERROR: $0 : cannot create the log files $logFile and $errorFile: $!\nExiting...\n";

######################################################################################################################################
######################################################################################################################################








######################################################################################################################################
######################################################################################################################################
# SPECIFIC PART OF MODULE TEST
######################################################################################################################################
######################################################################################################################################
########################################
##### pairing::extractName
########################################
my $fastqPairedData="$toggle/data/testData/fastq/pairedTwoIndividusIrigin/";
my $fastqSingleData="$toggle/data/testData/fastq/singleOneIndividuIrigin/";
my $fastqDir="/fastqDir";
mkdir $fastqDir;
my $mkdirCmd= "mkdir  $testDir$fastqDir";           # command to copy the original fastq file into the test directory
system ($mkdirCmd) and die ("ERROR: $0 : Cannot copy the repertory $testDir$fastqDir in the test directory with the command $mkdirCmd\n$!\n");    # RUN the copy command


my $expectName1=("irigin3_1");
my $expectRG1=("irigin3");

my ($obsName1, $obsRG1)=pairing::extractName('irigin3_1.fastq');
is_deeply($obsName1,$expectName1,'pairing::extractName - individu irigin3_1');
is_deeply($obsRG1,$expectRG1,'pairing::extractName - RG irigin3');

my $expectName2=("irigin3_2");
my $expectRG2=("irigin3");

my ($obsName2, $obsRG2)=pairing::extractName('irigin3_2.fastq');
is_deeply($obsName2,$expectName2,'pairing::extractName - individu irigin3_2');
is_deeply($obsRG2,$expectRG2,'pairing::extractName - RG irigin3');


########################################
##### pairing::pairRecognition
########################################

# input file
my $checkFastq = 1;

# copy paired and single files into $testDir
my $fastqFilePaired = $fastqPairedData."irigin*_*.fastq";     # fastq file
my $copyCmd= "cp $fastqFilePaired $testDir$fastqDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFilePaired in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

my $fastqFileSingle = $fastqSingleData."irigin*_*.fastq";     # fastq file
$copyCmd= "cp $fastqFileSingle  $testDir$fastqDir";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFileSingle in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command


my $expectedOutput={
          '@H3:C39R6ACXX:3:1101:1215:1877' => {
                                                         'ReadGroup' => 'irigin1',
                                                         'forward' => "$testDir$fastqDir/irigin1_1.fastq",
                                                         'reverse' => "$testDir$fastqDir/irigin1_2.fastq"
                                                       },
          '@H2:C381HACXX:5:1101:1359:1908' => {
                                                 'ReadGroup' => 'irigin3',
                                                 'forward' => "$testDir$fastqDir/irigin3_1.fastq",
                                                 'reverse' => "$testDir$fastqDir/irigin3_2.fastq"
                                               },
          '@H3:C39R6ACXX:3:1101:1192:1848' => {
                                                   'ReadGroup' => 'irigin2',
                                                   'forward' => "$testDir$fastqDir/irigin2_1.fastq",
                                                }
        };

my $listFiles_ref=toolbox::readDir("$testDir$fastqDir");
my $observedOutput=pairing::pairRecognition($listFiles_ref,$checkFastq);
##DEBUG print "pairRecognition Expected :\n"; print Dumper ($expectedOutput);print "pairRecognition Observed:\n"; print Dumper ($observedOutput);
is_deeply($observedOutput,$expectedOutput,'pairing::pairRecognition - output list');


########################################
##### pairing::createDirPerCouple
########################################

my $checkValue3=pairing::createDirPerCouple($observedOutput,"$testDir$fastqDir");
is ($checkValue3,1,'pairing::createDirPerCouple - running');

# Filetree expected
my $expectedFileTree = "$testDir$fastqDir:
irigin1
irigin2
irigin3

$testDir$fastqDir/irigin1:
irigin1_1.fastq
irigin1_2.fastq

$testDir$fastqDir/irigin2:
irigin2_1.fastq

$testDir$fastqDir/irigin3:
irigin3_1.fastq
irigin3_2.fastq
";

my $observedFileTree = `ls -R $testDir$fastqDir`;

##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($observedFileTree,$expectedFileTree,'pairing::pairRecognition - Filetree created');

########################################
##### pairing::repairing
########################################

# input file
my $rmDirCmd= "rm -r $testDir$fastqDir";           # command to copy the original fastq file into the test directory
system ($rmDirCmd) and die ("ERROR: $0 : Cannot removed $testDir$fastqDir in the test directory with the command $rmDirCmd\n$!\n");    # RUN the rm command

 $fastqPairedData="$toggle/data/testData/fastq/pairingRepairing/";
# copy paired and single files into $testDir
$fastqFilePaired = $fastqPairedData."irigin*_*.fastq";     # fastq file
$copyCmd= "cp $fastqFilePaired ./";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFilePaired in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

$fastqFileSingle = $fastqSingleData."irigin*_*.fastq";     # fastq file
$copyCmd= "cp $fastqFileSingle ./";           # command to copy the original fastq file into the test directory
system ($copyCmd) and die ("ERROR: $0 : Cannot copy the file $fastqFileSingle in the test directory with the command $copyCmd\n$!\n");    # RUN the copy command

#Check if running

my $checkValue=pairing::repairing('irigin3_1.CUTADAPT.fastq','irigin3_2.CUTADAPT.fastq',".",$checkFastq);
is ($checkValue,'1','pairing::repairing - running');

#Check if files created
$expectedFileTree = ".:
irigin2_1.fastq
irigin3_1.CUTADAPT.fastq
irigin3_1.REPAIRING.fastq
irigin3_2.CUTADAPT.fastq
irigin3_2.REPAIRING.fastq
irigin3_Single
pairing_log.e
pairing_log.o

./irigin3_Single:
irigin3Single.fastq
";

$observedFileTree = `ls -R`;

##DEBUG print "Expected: \n"; print Dumper ($expectedFileTree);print "Observed: \n"; print Dumper ($observedFileTree);
is_deeply($observedFileTree,$expectedFileTree,'pairing::pairRecognition - Filetree created');

#Check if working
my $numberOfLinesObserved=`wc -l irigin3_Single/irigin3Single.fastq`;
chomp $numberOfLinesObserved;
is ($numberOfLinesObserved,'8 '.'irigin3_Single/irigin3Single.fastq','pairing::repairing - single file');

#Check if the files created are the same as planned
my $diffForward=`diff -q irigin3_1.REPAIRING.fastq $toggle/data/testData/toolbox/irigin3_1.REPAIRING.fastq`;
is ($diffForward,'','pairing::repairing - forward file');

#Check if the files created are the same as planned
my $diffReverse=`diff -q irigin3_2.REPAIRING.fastq $toggle/data/testData/toolbox/irigin3_2.REPAIRING.fastq`;
is ($diffReverse,'','pairing::repairing - reverse file');

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

my $expectedData="$toggle/data/expectedData/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="$toggle/dataTest/snpEffTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf snpEff_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");




########################################
#Picking up data for tests
########################################

my $database="Osa1"; #need a currently formatted db for snpEff, in the test DATA-TEST/data folder
my $nonAnnotatedVcf="$toggle/data/expectedData/nonAnnotated.vcf"; #name of the original vcf file
my $annotatedVcf="$toggle/data/expectedData/annotated.vcf"; #expected file
my $outputVcf=$testingDir."/output.vcf";
my $originalGff = "$toggle/data/expectedData/genes.gff";
my $gff = $testingDir."/genes.gff";
my $copyCommand=" cp $originalGff $gff";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the gff file with the command $copyCommand\n$!\n");
my $originalReference = "$toggle/data/expectedData/sequences.fa";
my $reference=$testingDir."/sequences.fa";
$copyCommand=" cp $originalReference $reference";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the fasta file with the command $copyCommand\n$!\n");
my $name = "testBuild";
my $originalConfigFile = "$toggle/data/expectedData/snpEff.config";
my $configFile = $testingDir."/snpEff.config";
$copyCommand=" cp $originalConfigFile $configFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the configFile file with the command $copyCommand\n$!\n");

#######################################################################################################
####Test for snpeff snpeffAnnotation running
#######################################################################################################
##Test for running
#my $optionsHachees=$configInfos->{'snpEff annotator'};
my $optionsHachees="";

is(snpeff::snpeffAnnotation($nonAnnotatedVcf,$database,$outputVcf,$optionsHachees),'1','Test for snpeffAnnotation running');

###Verify if output are correct for snpeffAnnotation
my @expectedOutput=("$toggle/dataTest/snpeffTestDir/genes.gff","$toggle/dataTest/snpeffTestDir/output.vcf","$toggle/dataTest/snpeffTestDir/sequences.fa");
my @outPut=toolbox::readDir($testingDir);
is_deeply(@outPut,\@expectedOutput,'Test for the output files produced by snpeffAnnotation');


###Test for correct file value of snpeffAnnotation using a md5sum file control -  work through the different snpeff versions
my $diffResult=`diff $outputVcf $annotatedVcf`; #Differences between the tested output and the expected one ?
chomp $diffResult;
is($diffResult,"",'Test for the content of the snpeffAnnotation output');

#######################################################################################################
####Test for snpeff dbcreator running
#######################################################################################################
##Test for running
#$optionsHachees=$configInfos->{'snpEff build'};
#is(snpeff::dbCreator($gff,$reference,$name,$optionsHachees),'1','Test for dbCreator running');

###Verify if output are correct for snpeffAnnotation
#@expectedOutput=('$toggle/dataTest/snpeffTestDir/output.vcf');
#@outPut=toolbox::readDir($testingDir);
#is_deeply(@outPut,\@expectedOutput,'Test for the output files produced by snpeffAnnotation');


###Test for correct file value of snpeffAnnotation using a md5sum file control -  work through the different snpeff versions
#my $diffResult=`diff $outputVcf $annotatedVcf`; #Differences between the tested output and the expected one ?
#chomp $diffResult;
#is($diffResult,"",'Test for the content of the snpeffAnnotation output');


exit;

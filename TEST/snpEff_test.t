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

#Will test if snpeff works correctly
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../Modules/);
use Data::Dumper;

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"snpeff\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf snpeff_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

#########################################
#Remove the files and directory created by the previous test
#########################################
$cleaningCommand="rm -Rf ../DATA-TEST/snpeffTestDir";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous test dir with the command $cleaningCommand \n$!\n");

########################################
#Creation of test directory
########################################
my $testingDir="../DATA-TEST/snpeffTestDir";
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

########################################
#use of snpeff module ok
########################################
use_ok('toolbox') or exit;
use_ok('snpeff') or exit;
can_ok( 'snpeff','snpeffAnnotation');
can_ok( 'snpeff', 'dbCreator');

use toolbox;
use snpeff;

toolbox::readFileConf("software.config.txt");

########################################
#Picking up data for tests
########################################
my $database="Osa1"; #need a currently formatted db for snpEff, in the test DATA-TEST/data folder
my $nonAnnotatedVcf="../DATA/expectedData/nonAnnotated.vcf"; #name of the original vcf file
my $annotatedVcf="../DATA/expectedData/annotated.vcf"; #expected file
my $outputVcf=$testingDir."/output.vcf";
my $originalGff = "../DATA/expectedData/genes.gff";
my $gff = $testingDir."/genes.gff";
my $copyCommand=" cp $originalGff $gff";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the gff file with the command $copyCommand\n$!\n");
my $originalReference = "../DATA/expectedData/sequences.fa";
my $reference=$testingDir."/sequences.fa";
$copyCommand=" cp $originalReference $reference";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the fasta file with the command $copyCommand\n$!\n");
my $name = "testBuild";
my $originalConfigFile = "../DATA/expectedData/snpEff.config";
my $configFile = $testingDir."/snpEff.config";
$copyCommand=" cp $originalConfigFile $configFile";
system($copyCommand) and die ("ERROR: $0 : Cannot copy the configFile file with the command $copyCommand\n$!\n");

#######################################################################################################
####Test for snpeff snpeffAnnotation running
#######################################################################################################
##Test for running
my $optionsHachees=$configInfos->{'snpEff annotator'};
#
#is(snpeff::snpeffAnnotation($nonAnnotatedVcf,$database,$outputVcf,$optionsHachees),'1','Test for snpeffAnnotation running');
#
####Verify if output are correct for snpeffAnnotation
#my @expectedOutput=('../DATA-TEST/snpeffTestDir/genes.gff','../DATA-TEST/snpeffTestDir/output.vcf','../DATA-TEST/snpeffTestDir/sequences.fa');
#my @outPut=toolbox::readDir($testingDir);
#is_deeply(@outPut,\@expectedOutput,'Test for the output files produced by snpeffAnnotation');
#
#
####Test for correct file value of snpeffAnnotation using a md5sum file control -  work through the different snpeff versions
#my $diffResult=`diff $outputVcf $annotatedVcf`; #Differences between the tested output and the expected one ?
#chomp $diffResult;
#is($diffResult,"",'Test for the content of the snpeffAnnotation output');

#######################################################################################################
####Test for snpeff dbcreator running
#######################################################################################################
##Test for running
$optionsHachees=$configInfos->{'snpeff build'};
is(snpeff::dbCreator($gff,$reference,$name,$optionsHachees),'1','Test for dbCreator running');

###Verify if output are correct for snpeffAnnotation
#@expectedOutput=('../DATA-TEST/snpeffTestDir/output.vcf');
#@outPut=toolbox::readDir($testingDir);
#is_deeply(@outPut,\@expectedOutput,'Test for the output files produced by snpeffAnnotation');


###Test for correct file value of snpeffAnnotation using a md5sum file control -  work through the different snpeff versions
#my $diffResult=`diff $outputVcf $annotatedVcf`; #Differences between the tested output and the expected one ?
#chomp $diffResult;
#is($diffResult,"",'Test for the content of the snpeffAnnotation output');


exit;

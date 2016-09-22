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
use_ok('fastqc');
can_ok('fastqc','execution');

use fastqc;

my $expectedData="../../DATA/expectedData/";

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/fastqcTestDir";
my $cleaningCmd="rm -Rf $testingDir && mkdir $testingDir"; 
system ($cleaningCmd) and die ("ERROR: $0 : Cannot clean or create the test directory with the command $cleaningCmd \n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCmd="echo \"fastqc\nTEST\" > individuSoft.txt";
system($creatingCmd) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCmd\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
$cleaningCmd="rm -Rf fastqc_TEST_log.*";
system($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCmd \n$!\n");

##########################################
#Fastqc::execution test
##########################################

# input file
my $fastqFile = "RC3_2.fastq";
my $fastqPath = $expectedData.$fastqFile;      
my $fastqcDir = "fastqcOut";
my $makeDirCmd = "mkdir $fastqcDir";
system ($makeDirCmd) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCmd\n$!\n");

# execution test
is(fastqc::execution($fastqPath,$fastqcDir),1,'fastqc::execution');     

# expected output test
my $expectedOutput = 'fastqcOut/RC3_2_fastqc.zip';
my $observedOutput = `ls $fastqcDir`;
chomp($observedOutput);
$observedOutput = $fastqcDir."/".$observedOutput;
is($observedOutput,$expectedOutput,'fastqc::execution - output list'); 

# expected content test
my $observedContent=`unzip -l $observedOutput | tail -n1`;
my $validContent = ( $observedContent =~ m/19 files/);
is($validContent,1,'fastqc::execution - output content');
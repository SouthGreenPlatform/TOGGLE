#!/usr/bin/perl -w

###################################################################################################################################
#
# Copyright 2014 IRD-CIRAD
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
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################


use strict;

#Will test if tophat works correctly
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use lib qw(../Modules/);
use Data::Dumper;

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"tophat\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");


#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf tophat_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");

########################################
#initialisation and setting configs
########################################
my $testingDir="../DATA-TEST/tophatTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

my $OriginalFastaRef="../DATA/expectedData/Reference.fasta";
my $fastaRef="$testingDir/Reference.fasta";
my $refCopyCom="cp $OriginalFastaRef $fastaRef";
system($refCopyCom) and die ("ERROR: $0 : Cannot copy the Reference $OriginalFastaRef with the command $refCopyCom\n$!\n");     #Now we have a ref to be tested


########################################
#use of module ok
########################################
use_ok('toolbox') or exit;
use_ok('tophat') or exit;
can_ok( 'tophat','bowtieBuild');

use toolbox;
use tophat;

################################################################################################
###tophat::bowtieBuild
################################################################################################
my $expectedIndexPrefix="Reference";
my $observedIndexPrefix=is(tophat::bowtieBuild($fastaRef),1, 'OK for bowtieBuild RUNNING');
is($expectedIndexPrefix,$observedIndexPrefix,'OK for prefix index');

exit;

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
my $observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools faidx output structure');

################################################################################################
###tophat::bowtieBuild
################################################################################################
my $expectedIndexPrefix="Reference";
my $observedIndexPrefix=is(tophat::bowtieBuild($fastaRef),1, 'OK for bowtieBuild RUNNING');
is($expectedIndexPrefix,$observedIndexPrefix,'OK for prefix index');

###Checking the correct structure for the output file using md5sum
my $expectedMD5sum="4b9a4431e72c9db7e5c1f2153eba9fe7";
my $observedMD5sum=`md5sum $fastaRef.fai`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'Ok for the content of the samtools faidx output structure');

exit;
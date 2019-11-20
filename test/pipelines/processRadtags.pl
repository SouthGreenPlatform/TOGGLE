#!/usr/bin/env perl

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

use strict;
use warnings;
use Test::More 'no_plan';
use Test::Deep;
use fileConfigurator;
use localConfig;
use lib qw(./);
use sedToggle;

#####################
## PIPELINE TEST
#####################

#####################
## TOGGLE radseq processRadtags Single
#####################

my $datafastqRadseq = "$toggle/data/testData/radseq/single/";
my $keyfileSingle = "$toggle/data/testData/radseq/keyfileTestSingle";

print "\n\n#################################################\n";
print "#### TEST processRadtags  (single) / no SGE mode\n";
print "#################################################\n";

# Copy file config
my $fileconfRadseq="$toggle/exampleConfigs/radseqSingle.config.txt";
my $fileconfRadseqNoSGE="radseqSingle.config.txt";

my $cmd="cp $fileconfRadseq $fileconfRadseqNoSGE";
system($cmd) and die ("#### ERROR COPY CONFIG FILE: $cmd\n");     # Copy into TEST

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/processRadtags-Single-noSGE";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $runCmd = "toggleGenerator.pl -c ".$fileconfRadseqNoSGE." -d ".$datafastqRadseq." -k ".$keyfileSingle." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for processRadtags no SGE mode";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;

my @expectedOutput = ('33-16_fastqc.html','33-16_fastqc.zip','A239_fastqc.html','A239_fastqc.zip','A272_fastqc.html','A272_fastqc.zip','A554_fastqc.html','A554_fastqc.zip','A619_fastqc.html','A619_fastqc.zip','A632_fastqc.html','A632_fastqc.zip','A654_fastqc.html','A654_fastqc.zip','A659_fastqc.html','A659_fastqc.zip','A680_fastqc.html','A680_fastqc.zip','A682_fastqc.html','A682_fastqc.zip','B103_fastqc.html','B103_fastqc.zip','B104_fastqc.html','B104_fastqc.zip','B109_fastqc.html','B109_fastqc.zip','B10_fastqc.html','B10_fastqc.zip','B73Htrhm_fastqc.html','B73Htrhm_fastqc.zip','B76_fastqc.html','B76_fastqc.zip','B77_fastqc.html','B77_fastqc.zip','B84_fastqc.html','B84_fastqc.zip','B97_fastqc.html','B97_fastqc.zip','C103_fastqc.html','C103_fastqc.zip','CH701-30_fastqc.html','CH701-30_fastqc.zip','CI28A_fastqc.html','CI28A_fastqc.zip','CI3A_fastqc.html','CI3A_fastqc.zip','CI66_fastqc.html','CI66_fastqc.zip','CI-7_fastqc.html','CI-7_fastqc.zip','CM105_fastqc.html','CM105_fastqc.zip','CM174_fastqc.html','CM174_fastqc.zip','CML14_fastqc.html','CML14_fastqc.zip','CML157Q_fastqc.html','CML157Q_fastqc.zip','CML220_fastqc.html','CML220_fastqc.zip','CML228_fastqc.html','CML228_fastqc.zip','CML238_fastqc.html','CML238_fastqc.zip','CML258_fastqc.html','CML258_fastqc.zip','CML277_fastqc.html','CML277_fastqc.zip','CML281_fastqc.html','CML281_fastqc.zip','CML311_fastqc.html','CML311_fastqc.zip','CML323_fastqc.html','CML323_fastqc.zip','CML332_fastqc.html','CML332_fastqc.zip','CML333_fastqc.html','CML333_fastqc.zip','CML45_fastqc.html','CML45_fastqc.zip','CML52_fastqc.html','CML52_fastqc.zip','CML91_fastqc.html','CML91_fastqc.zip','CML92_fastqc.html','CML92_fastqc.zip','CO106_fastqc.html','CO106_fastqc.zip','CO125_fastqc.html','CO125_fastqc.zip','DE-2_fastqc.html','DE-2_fastqc.zip','DE-3_fastqc.html','DE-3_fastqc.zip','EMPTY_fastqc.html','EMPTY_fastqc.zip','EP1_fastqc.html','EP1_fastqc.zip','F6_fastqc.html','F6_fastqc.zip','F7_fastqc.html','F7_fastqc.zip','H91_fastqc.html','H91_fastqc.zip','H95_fastqc.html','H95_fastqc.zip','I137TN_fastqc.html','I137TN_fastqc.zip','I29_fastqc.html','I29_fastqc.zip','IDS28_fastqc.html','IDS28_fastqc.zip','K55_fastqc.html','K55_fastqc.zip','Ki11_fastqc.html','Ki11_fastqc.zip','Ki14_fastqc.html','Ki14_fastqc.zip','L317_fastqc.html','L317_fastqc.zip','M14_fastqc.html','M14_fastqc.zip','NC230_fastqc.html','NC230_fastqc.zip','NC250_fastqc.html','NC250_fastqc.zip','NC262_fastqc.html','NC262_fastqc.zip','NC290A_fastqc.html','NC290A_fastqc.zip','NC298_fastqc.html','NC298_fastqc.zip','NC300_fastqc.html','NC300_fastqc.zip','NC318_fastqc.html','NC318_fastqc.zip','NC320_fastqc.html','NC320_fastqc.zip','NC328_fastqc.html','NC328_fastqc.zip','NC338_fastqc.html','NC338_fastqc.zip','NC356_fastqc.html','NC356_fastqc.zip','NC364_fastqc.html','NC364_fastqc.zip','NC368_fastqc.html','NC368_fastqc.zip','Oh40B_fastqc.html','Oh40B_fastqc.zip','Oh43E_fastqc.html','Oh43E_fastqc.zip','Oh603_fastqc.html','Oh603_fastqc.zip','Os420_fastqc.html','Os420_fastqc.zip','Pa875_fastqc.html','Pa875_fastqc.zip','R109B_fastqc.html','R109B_fastqc.zip','R168_fastqc.html','R168_fastqc.zip','SC213R_fastqc.html','SC213R_fastqc.zip','Sg18_fastqc.html','Sg18_fastqc.zip','tripsacum_fastqc.html','tripsacum_fastqc.zip','Tzi25_fastqc.html','Tzi25_fastqc.zip','Va22_fastqc.html','Va22_fastqc.zip','Va35_fastqc.html','Va35_fastqc.zip','Va59_fastqc.html','Va59_fastqc.zip','Va99_fastqc.html','Va99_fastqc.zip','W117Ht_fastqc.html','W117Ht_fastqc.zip','W153R_fastqc.html','W153R_fastqc.zip','W182B_fastqc.html','W182B_fastqc.zip','W22_fastqc.html','W22_fastqc.zip','WD_fastqc.html','WD_fastqc.zip');


# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - processRadtags (no SGE) - list ');

# expected output content
my $observedContent=`unzip -l $testingDir/finalResults/33-16_fastqc.zip | tail -n1`;
my $validContent = ( $observedContent =~ m/20 files/);
is($validContent,1,'toggleGenerator - processRadtags (no SGE) - output content');

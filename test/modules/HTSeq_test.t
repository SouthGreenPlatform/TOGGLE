#!/usr/bin/perl

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

########################################
#HTSeq::htseqCount bam file
########################################
my $bankData="$toggle/data/Bank/";
my $testData="$toggle/data/testData/";

# input file
my $gffRef="$bankData/referenceRnaseqGFF.gff3";
my $bamIn="$testData/samBam/oneBam/RC3-SAMTOOLSVIEW.bam";
my $bam="accepted_hits.SAMTOOLSSORT.bam";

#copy bam reference into test directory where the index will be created
my $copyCommand="cp $bamIn ./$bam";
system ($copyCommand) and die "ERROR: $0: Cannot copy the bam file with the command $copyCommand \n$!\n";


#htseq option
my %optionsHachees = (
                      "-r" => "name",
                      "-s" => "no",
                      "-t" => "mRNA",
                      "-m" => "union",
                      "-i" => "ID",
                      );       

my $optionHachees = \%optionsHachees;                           # Ref of the hash

#outputfile
my $htseqcountFile="accepted_hits.HTSEQCOUNT.txt";
is(HTSeq::htseqCount($bam, $htseqcountFile,$gffRef, $optionHachees),1,'HTSeq::htseqCount (bam file)');

# expected output test
my $observedOutput = `ls`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('accepted_hits.HTSEQCOUNT.txt','accepted_hits.SAMTOOLSSORT.bam','accepted_hits.SAMTOOLSSORT.sam','htseq_log.e','htseq_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'HTSeq::htseqCount - output list');

# expected content test
my $expectedMD5sum="a97aaf22fa76469ba1ec630429600a5e";
my $observedMD5sum=`md5sum $htseqcountFile`;# structure of the test file
my @withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'HTSeq::htseqCount - output content');



########################################
#HTSeq::htseqCount sam file
########################################

#input file
my $samIni="$testData/samBam/oneSam/RC3-SAMTOOLSVIEW.sam";
my $sam="accepted_hits.SAMTOOLSSORT.sam";

#copy fasta reference into test directory where the index will be created
my $rmCommand="rm $bam $sam";
system ($rmCommand) and die "ERROR: $0: Cannot remove bam and sam files generated precedently with the command $rmCommand \n$!\n";

#copy fasta reference into test directory where the index will be created
$copyCommand="cp $samIni ./$sam";
system ($copyCommand) and die "ERROR: $0: Cannot copy the refence file with the command $copyCommand \n$!\n";

#htseq option
%optionsHachees = (
                      "-r" => "name",
                      "-s" => "no",
                      "-t" => "mRNA",
                      "-m" => "union",
                      "-i" => "ID",
                      );       

$optionHachees = \%optionsHachees;                           # Ref of the hash

#outputfile
$htseqcountFile="accepted_hits.SAM-HTSEQCOUNT.txt";
is(HTSeq::htseqCount($sam, $htseqcountFile,$gffRef, $optionHachees),1,'HTSeq::htseqCount (sam file)');

# expected output test
$observedOutput = `ls`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('accepted_hits.HTSEQCOUNT.txt','accepted_hits.SAM-HTSEQCOUNT.txt','accepted_hits.SAMTOOLSSORT.sam','htseq_log.e','htseq_log.o');

is_deeply(\@observedOutput,\@expectedOutput,'HTSeq::htseqCount - output list');

# expected content test
$expectedMD5sum="a97aaf22fa76469ba1ec630429600a5e";
$observedMD5sum=`md5sum $htseqcountFile`;# structure of the test file
@withoutName = split (" ", $observedMD5sum);     # to separate the structure and the name of the test file
$observedMD5sum = $withoutName[0];       # just to have the md5sum result
is($observedMD5sum,$expectedMD5sum,'HTSeq::htseqCount - output content');

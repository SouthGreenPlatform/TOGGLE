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

#####################
## TOGGLE samtools sortsam
#####################

#Input data
my $dataOneBam = "$toggle/data/testData/samBam/oneBamUnsorted/";
my $dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

print "\n\n#################################################\n";
print "#### TEST SAMtools sort / no SGE mode\n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/samToolsSort-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("samToolsSort");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");


my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools sort";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('unsorted.SAMTOOLSSORT.bam');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) sorting list ');

# expected output content
$observedOutput=`samtools view $testingDir/finalResults/unsorted.SAMTOOLSSORT.bam | head -n 2 | cut -f4`; # We pick up only the position field
chomp $observedOutput;
my @position = split /\n/, $observedOutput;
$observedOutput= 0;
$observedOutput = 1 if ($position[0] < $position[1]); # the first read is placed before the second one
is($observedOutput,"1", 'toggleGenerator - One Bam (no SGE) sorting content ');

#####################
## TOGGLE samtools MpileUp
#####################

#Input data
$dataOneBam = "$toggle/data/testData/samBam/oneBam/";

print "\n\n#################################################\n";
print "#### TEST SAMtools MpileUp\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/samToolsMpileUp";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("samToolsMpileUp");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");


$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools MpileUp";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RC3-SAMTOOLSVIEW.SAMTOOLSMPILEUP.mpileup');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) MpileUp list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/RC3-SAMTOOLSVIEW.SAMTOOLSMPILEUP.mpileup`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,"155888 $toggle/dataTest/samToolsMpileUp/finalResults/RC3-SAMTOOLSVIEW.SAMTOOLSMPILEUP.mpileup", 'toggleGenerator - One Bam (no SGE) MpileUp content ');

#####################
## TOGGLE samtools depth
#####################

#Input data
my $dataMultipleBam = "$toggle/data/testData/samBam/twoBamsIrigin/";
$dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

print "\n\n#################################################\n";
print "#### TEST SAMtools depth\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/samToolsDepth-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("samToolsDepth");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");


$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataMultipleBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools depth";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('irigin1-PICARDTOOLSMARKDUPLICATES.SAMTOOLSDEPTH.depth', 'irigin3-PICARDTOOLSMARKDUPLICATES.SAMTOOLSDEPTH.depth');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - samtools depth list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/irigin1-PICARDTOOLSMARKDUPLICATES.SAMTOOLSDEPTH.depth`; # We pick up only the position field
chomp $observedOutput;
@position = split /\s/, $observedOutput;
$observedOutput= $position[0];
is($observedOutput,"6172", 'toggleGenerator - samtools depth content ');

#####################
## TOGGLE samtools merge
#####################

#Input data
$dataMultipleBam = "$toggle/data/testData/samBam/twoBamsIrigin/";
$dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

print "\n\n#################################################\n";
print "#### TEST SAMtools merge\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/samToolsMerge-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("1000","samToolsMerge");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");


$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataMultipleBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools merge";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('multipleAnalysis.SAMTOOLSMERGE.bam');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - samtools merge list ');

# expected output content
$observedOutput=`samtools view -h $testingDir/finalResults/multipleAnalysis.SAMTOOLSMERGE.bam | wc -l`; # We pick up only the position field
chomp $observedOutput;
@position = split /\s/, $observedOutput;
$observedOutput= $position[0];
is($observedOutput,"4964", 'toggleGenerator - samtools merge content ');

#####################
## TOGGLE samtools idxstats
#####################

#Input data
$dataMultipleBam = "$toggle/data/testData/samBam/twoBamsIrigin/";
$dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

print "\n\n#################################################\n";
print "#### TEST SAMtools idxstats\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/samToolsIdxstats-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("samToolsIdxStats");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");


$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataMultipleBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools idxstats";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('irigin1-PICARDTOOLSMARKDUPLICATES.idxstats','irigin3-PICARDTOOLSMARKDUPLICATES.idxstats');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - samtools idxstats list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/irigin1-PICARDTOOLSMARKDUPLICATES.idxstats`; 
chomp $observedOutput;
@position = split /\s/, $observedOutput;
$observedOutput= $position[0];
is($observedOutput,"952", 'toggleGenerator - samtools idxstats content ');


#####################
## TOGGLE samtools faidx
#####################

#Input data
my $dataFasta = "$toggle/data/testData/fasta/TGICL/";
$dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";

print "\n\n#################################################\n";
print "#### TEST SAMtools faidx\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/samToolsFaidx-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("samToolsFaidx");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");
##Adding command to faidx fileConf
#my $addCom = "echo \"\n\$samtoolsfaidx\n\TRINITY_DN75358_c0_g1_i1:1-100\n\" >> blockTestConfig.txt";
#system($addCom);


$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataFasta." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for samtools faidx";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('contig.faidx.fasta');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - samtools faidx list ');

# expected output content
$observedOutput=`cat $testingDir/finalResults/contig.faidx.fasta`; 
chomp $observedOutput;
my $expectedOutput = ">TRINITY_DN75358_c0_g1_i1:1-100
ATCACACGAACACATACGCTCTTCCTTTCTCTGAAATAAGAAGGAGGAGAAGAAGAAGGT
CATGGGGTTCACAAAACTCGCCCTAGTTCTGGCAATAGCA";
is($observedOutput,$expectedOutput, 'toggleGenerator - samtools faidx content ');

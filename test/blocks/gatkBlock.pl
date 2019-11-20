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
## PATH for datas test
#####################

# references files
my $dataRefIrigin = "$toggle/data/Bank/referenceIrigin.fasta";
my $dataRefArcad = "$toggle/data/Bank/referenceArcad.fasta";
my $dataOneBam = "$toggle/data/testData/samBam/oneBam/";



print "\n\n#################################################\n";
print "#### TEST gatk BaseRecalibrator\n";
print "#################################################\n";

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/gatkBaseRecalibrator-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
my @listSoft = ("gatkBaseRecalibrator");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

my $runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for gatkBaseRecalibrator";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
my $observedOutput = `ls $testingDir/finalResults`;
my @observedOutput = split /\n/,$observedOutput;
my @expectedOutput = ('RC3-SAMTOOLSVIEW.GATKBASERECALIBRATOR.tableReport');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) gatkBaseRecalibrator file list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKBASERECALIBRATOR.tableReport`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,"6273 $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKBASERECALIBRATOR.tableReport", 'toggleGenerator - One Bam (no SGE) gatkBaseRecalibrator content');


print "\n\n#################################################\n";
print "#### TEST gatk UnifiedGenotyper\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/gatkUnifiedGenotyper-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("gatkUnifiedGenotyper");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for gatkUnifiedGenotyper";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RC3-SAMTOOLSVIEW.GATKUNIFIEDGENOTYPER.vcf','RC3-SAMTOOLSVIEW.GATKUNIFIEDGENOTYPER.vcf.idx');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) gatkUnifiedGenotyper file list ');

# expected output content
$observedOutput=`grep -v "#" $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKUNIFIEDGENOTYPER.vcf`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,"2233572	145	.	A	G	54.74	.	AC=2;AF=1.00;AN=2;DP=2;Dels=0.00;ExcessHet=3.0103;FS=0.000;HaplotypeScore=0.0000;MLEAC=2;MLEAF=1.00;MQ=49.84;MQ0=0;QD=27.37;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:82,6,0", 'toggleGenerator - One Bam (no SGE) gatkUnifiedGenotyper content');

print "\n\n#################################################\n";
print "#### TEST gatkIndelRealigner\n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/gatkIndelRealigner-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("samToolsIndex","gatkRealignerTargetCreator","gatkIndelRealigner");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for gatkIndelRealigner";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RC3-SAMTOOLSVIEW.GATKINDELREALIGNER.bai','RC3-SAMTOOLSVIEW.GATKINDELREALIGNER.bam');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) gatkIndelRealigner file list ');

# expected output content
$observedOutput=`samtools view $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKINDELREALIGNER.bam | wc -l`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,"1998", 'toggleGenerator - One Bam (no SGE) gatkIndelRealigner content');



print "\n\n#################################################\n";
print "#### TEST gatkprintReads \n";
print "#################################################\n";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/gatkPrindReads-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("gatkBaseRecalibrator","gatkPrintReads");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneBam." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for gatkPrindReads";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('RC3-SAMTOOLSVIEW.GATKPRINTREADS.bai','RC3-SAMTOOLSVIEW.GATKPRINTREADS.bam');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Bam (no SGE) gatkPrindReads file list ');

# expected output content
$observedOutput=`wc -l $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKPRINTREADS.bam`; # We pick up only the position field
chomp $observedOutput;
is($observedOutput,"473 $testingDir/finalResults/RC3-SAMTOOLSVIEW.GATKPRINTREADS.bam", 'toggleGenerator - One Bam (no SGE) gatkPrindReads content');








#####################
## TOGGLE VCF singleVCF gatkSelectVariants
#####################

my $dataOneVcf = "$toggle/data/testData/vcf/singleVCF";

print "\n\n#################################################\n";
print "#### TEST one VCF / no SGE mode\n";
print "#################################################\n";

# Copy file config
my $fileVcf="../exampleConfigs/vcf.config.txt";

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/gatkSelectVariants-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

#Creating config file for this test
@listSoft = ("gatkSelectVariants");
fileConfigurator::createFileConf(\@listSoft,"blockTestConfig.txt");

$runCmd = "toggleGenerator.pl -c blockTestConfig.txt -d ".$dataOneVcf." -r ".$dataRefIrigin." -o ".$testingDir;
print "\n### Toggle running : $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for One Vcf no SGE mode";

# check final results
print "\n### TEST Ouput list & content : $runCmd\n";
$observedOutput = `ls $testingDir/finalResults`;
@observedOutput = split /\n/,$observedOutput;
@expectedOutput = ('GATKVARIANTFILTRATION.GATKSELECTVARIANT.vcf','GATKVARIANTFILTRATION.GATKSELECTVARIANT.vcf.idx');

# expected output test
is_deeply(\@observedOutput,\@expectedOutput,'toggleGenerator - One Vcf (no SGE) list ');

# expected output content
$observedOutput=`tail -n 1 $testingDir/finalResults/GATKVARIANTFILTRATION.GATKSELECTVARIANT.vcf`;
chomp $observedOutput;
my $expectedOutput="2290182	1013	.	A	G	42.74	PASS	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=29.00;MQ0=0;QD=21.37;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:70,6,0";
is($observedOutput,$expectedOutput, 'toggleGenerator - One Vcf (no SGE) content ');

#rm idx file on $dataOneVcf
$cleaningCmd="rm -f $dataOneVcf/GATKVARIANTFILTRATION.vcf.idx";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove GATKVARIANTFILTRATION.vcf.idx with the command $cleaningCmd \n$!\n");

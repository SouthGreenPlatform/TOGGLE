#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2017 IRD-CIRAD-INRA-ADNid
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

# fasta files
my $dataFasta = "$toggle/data/testData/fasta/TGICL";
my $dataFastqpairedOneIndividuPacaya = "$toggle/data/testData/fastq/assembly/pairedOneIndivuPacaya";
my $reference = "$toggle/data/testData/fasta/TGICL/contig_tgicl.fasta";

print "\n\n#################################################\n";
print "#### TEST TGICL Assembly\n";
print "#################################################\n";

#Creating config file for this test
my @listSoft = ("tgicl");
fileConfigurator::createFileConf(\@listSoft,"$toggle/test/blocks/blockTestConfig.txt");

# Remove files and directory created by previous test
my $testingDir="$toggle/dataTest/tgiclPacaya-noSGE-Blocks";
my $cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

my $runCmd = "toggleGenerator.pl -c $toggle/test/blocks/blockTestConfig.txt -d ".$dataFasta." -o ".$testingDir;

print "\n### $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for TGICL-Pacaya";


# check final results
#########################################################
###  Test for tgiclRun 
#########################################################

SKIP:
{
    # get hostname 
    my $host = `hostname`;
    chomp $host;
    
    #DEBUG diag("\n\n----".$host."---\n\n");

    #if tests wasn't running on master, the 3 tests following will be skipped
    skip "No tgicl test on node", 3 if ($host !~/^master0.alineos.net$/ );
    
    # expected content test
    my $cmd = 'grep -c "^>" '.$testingDir.'/finalResults/all_contigs.fasta';
    my $expectedAnswer="74"; # tested on master 19-10-2016 --- OK
    my $observedAnswer=`$cmd`;
    chomp($observedAnswer);

    is($observedAnswer,$expectedAnswer,'tgicl::tgiclRun- output content');
}



#####################
## TOGGLE assembly pairedOneIndividuPacaya
#####################
print "\n\n#################################################\n";
print "#### TEST Trinity assembly pairedOneIndividuPacaya (one individu) / no SGE mode\n";
print "#################################################\n";

## rm readcount in initialDir if exist
if (-e "$toggle/data/testData/fastq/assembly/pairedOneIndivuPacaya/g02L5Mapped_R1.fq.readcount")
{ 
     `rm $toggle/data/testData/fastq/assembly/pairedOneIndivuPacaya/*.readcount`;
}
if (-e "$toggle/data/testData/fastq/assembly/pairedOneIndivuPacaya/g02L5Mapped_R1.readcount")
{ 
     `rm $toggle/data/testData/fastq/assembly/pairedOneIndivuPacaya/*.readcount`;
}

#Creating config file for this test
@listSoft = ("trinity");
fileConfigurator::createFileConf(\@listSoft,"$toggle/test/blocks/blockTestConfig.txt");
#Output file
my $readGroup = 'g02L5Mapped'; ## à remplacer par le readGroup

# Remove files and directory created by previous test
$testingDir="$toggle/dataTest/pairedOneIndividuPacaya-noSGE-Blocks";
$cleaningCmd="rm -Rf $testingDir";
system ($cleaningCmd) and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCmd \n$!\n");

$runCmd = "toggleGenerator.pl -c $toggle/test/blocks/blockTestConfig.txt -d ".$dataFastqpairedOneIndividuPacaya." -o ".$testingDir;

print "\n### $runCmd\n";
system("$runCmd") and die "#### ERROR : Can't run TOGGLE for pairedOneIndividuPacaya";

# expected output test
my $observedOutput = `ls $testingDir"/finalResults"`;
my @observedOutput = split /\n/,$observedOutput;

my @expectedOutput = ($readGroup.'_both.fa',$readGroup.'_both.fa.ok',$readGroup.'_both.fa.read_count',$readGroup.'_inchworm.K25.L25.DS.fa',$readGroup.'_inchworm.K25.L25.DS.fa.finished',$readGroup.'_inchworm.kmer_count',$readGroup.'_jellyfish.kmers.fa',$readGroup.'_jellyfish.kmers.fa.histo',$readGroup.'_left.fa.ok',$readGroup.'_partitioned_reads.files.list',$readGroup.'_partitioned_reads.files.list.ok',$readGroup.'_recursive_trinity.cmds',$readGroup.'_recursive_trinity.cmds.completed',$readGroup.'_recursive_trinity.cmds.ok',$readGroup.'_right.fa.ok');

is_deeply(\@observedOutput,\@expectedOutput,'trinity::trinityRun - output list - One paired bank');
 


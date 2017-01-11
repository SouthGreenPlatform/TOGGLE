#!/usr/bin/perl

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

#Will test if tgicl works correctly
use strict;
use warnings;
use Test::More 'no_plan'; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;
use Test::Warn;

use Data::Dumper;
use lib qw(../Modules/);

########################################
#Test of the use of gatk modules
########################################
use_ok('toolbox') or exit;
use_ok('tgicl') or exit;

can_ok('tgicl','tgiclRun');

use toolbox;
use tgicl;

#########################################
#Remove files and directory created by previous test
#########################################
my $testingDir="../DATA-TEST/tgiclTestDir";
my $creatingDirCom="rm -Rf $testingDir ; mkdir -p $testingDir";                                    #Allows to have a working directory for the tests
system($creatingDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingDirCom\n$!\n");

chdir $testingDir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $testingDir\"\n$!\n");

my $outdir="./outputDir";
my $creatingSubDirCom="mkdir -p $outdir";                                    #Allows to have a working directory for the tests
system($creatingSubDirCom) and die ("ERROR: $0 : Cannot execute the command $creatingSubDirCom\n$!\n");
chdir $outdir or die ("ERROR: $0 : Cannot go into the new directory with the command \"chdir $outdir\"\n$!\n");

#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"tgicl\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf tgicl_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");


##########################################
#tgiclCreateSequenceDictionary test
##########################################


## Added by CD to resolve the bug with tgicl run (chdir into this sub)
my $path=`pwd`;  
chomp $path;
##########################################################################

#Input file
my $reference = $path."/../../../DATA/testData/fasta/TGICL/contig_tgicl.fasta";  
#my $reference = "/home/adiall/TOGGLE-DEV/DATA/testData/fasta/TGICL/contig_tgicl.fasta";  ### TODO: remplacer par chemin relatif mais bug dans module
#Output file
#my $readGroup = 'g02L5'; ## à remplacer par le readGroup

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
    
    # running test
    is(tgicl::tgiclRun($outdir,$reference),1,'tgicl::tgiclRun');

    # expected output test
    my $observedOutput = `ls`;
    my @observedOutput = split /\n/,$observedOutput;
    #
    my @expectedOutput = ('all_contigs.fasta','contig_tgicl.fasta.cidx','contig_tgicl.fasta_clusters','contig_tgicl.fasta.nhr','contig_tgicl.fasta.nin','contig_tgicl.fasta.nsq','contig_tgicl.fasta.singletons','err_tgicl_contig_tgicl.fasta.log','formatdb.log','hitsort_001.Z','individuSoft.txt','masked.lst','singletons.fasta','tgicl_contig_tgicl.fasta.log','tgicl_TEST_log.e','tgicl_TEST_log.o');
    #
    is_deeply(\@observedOutput,\@expectedOutput,'tgicl::tgiclRun - output list');
    
    # expected content test
    my $cmd = 'grep -c "^>" all_contigs.fasta';
    my $expectedAnswer="74"; # tested on master 19-10-2016 --- OK
    my $observedAnswer=`$cmd`;
    chomp($observedAnswer);

    is($observedAnswer,$expectedAnswer,'tgicl::tgiclRun- output content');
}
 
1;

#exit;
#__END__

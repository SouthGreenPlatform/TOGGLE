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

##########################################
#tgiclCreateSequenceDictionary test
##########################################


#Input file
my $reference = "$toggle/data/testData/fasta/TGICL/contig_tgicl.fasta";  

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
    is(tgicl::tgiclRun($testDir,$reference),1,'tgicl::tgiclRun');

    # expected output test
    my $observedOutput = `ls`;
    my @observedOutput = split /\n/,$observedOutput;
    #
    my @expectedOutput = ('all_contigs.fasta','contig_tgicl.fasta.cidx','contig_tgicl.fasta_clusters','contig_tgicl.fasta.nhr','contig_tgicl.fasta.nin','contig_tgicl.fasta.nsq','contig_tgicl.fasta.singletons','err_tgicl_contig_tgicl.fasta.log','formatdb.log','hitsort_001.Z','masked.lst','singletons.fasta','tgicl_contig_tgicl.fasta.log','tgicl_log.e','tgicl_log.o');
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

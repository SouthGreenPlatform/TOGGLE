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

#Will test if localConfig is Ok

use strict;
use warnings;
use Test::More 'no_plan'; 
use Data::Dumper;
use Test::Deep;
use lib qw(../../modules/);


########################################
#use of localConfig module ok
########################################

use_ok('localConfig');

use localConfig;

########################################
#Capture of location of localConfig module ok
########################################

# We need to verify if the import pass weel, thus we parse the localConfig.pm to check what is expected to be recovered
#Is the comment enoughly explicit ?
open LOCALCONFIG, "../../modules/localConfig.pm";

my %dictLocation;
while (<LOCALCONFIG>)
{
	chomp $_;
	if(/our \$/)
	{
		chop $_;
		my @line = split ("our ", $_);
		$line[1] =~ s/\s=\s/=/g; # remove space $java = "/....";
		$line[1] =~ s/\\//g; # remove \ $java = "/....";
		$line[1] =~ s/;\s*$//g; # remove ;\s in the end line
		@line = split("=", $line[1]);
		$line[1]=~s/"//g if defined $line[1];
		$dictLocation{$line[0]} = $line[1];
	}
}
close LOCALCONFIG;


# Verify 0 = 0;
is($java,"/usr/java/jre1.8.0_111//bin/java -jar","java infos");

is($bwa,$dictLocation{"\$bwa"},"bwa infos");

is($picard,"/usr/java/jre1.8.0_111//bin/java -jar /usr/local/picard-tools-2.5.0//picard.jar","picard infos");

is($samtools,$dictLocation{"\$samtools"},"samtools infos");

is($GATK,'/usr/java/jre1.8.0_111//bin/java -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar',"GATK infos");

is($fastqc,$dictLocation{"\$fastqc"},"fastqc infos");

is($cutadapt,$dictLocation{"\$cutadapt"},"cutadapt infos");

is($fastxTrimmer,$dictLocation{"\$fastxTrimmer"},"fastx-Trimmer infos");

is($tophat2,$dictLocation{"\$tophat2"},"Tophat2 infos");

is($bowtie2Build,$dictLocation{"\$bowtie2Build"},"bowtie2Build infos");

is($bowtieBuild,$dictLocation{"\$bowtieBuild"},"bowtieBuild infos");

is($htseqcount,$dictLocation{"\$htseqcount"},"HTseq-count infos");

is($cufflinks,$dictLocation{"\$cufflinks"},"cufflinks infos");

is($cuffdiff,$dictLocation{"\$cuffdiff"},"cuffdiff infos");

is($cuffmerge,$dictLocation{"\$cuffmerge"},"cuffmerge infos");

is($snpEff,'/usr/java/jre1.8.0_111//bin/java -jar /usr/local/snpEff-4.3//snpEff.jar',"snpEff infos");



######################################
# Testing the correct location of bwa
######################################


`$bwa 2> /tmp/out.txt`; #We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
my $line;
while (<OUT>) {
    $line=<OUT>;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"Program: bwa (alignment via Burrows-Wheeler transformation)\n","Test for bwa location");

######################################
# Testing the correct location of samtools
######################################

`$samtools 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=<OUT>;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"Program: samtools (Tools for alignments in the SAM format)\n","Test for samtools location");

######################################
# Testing the correct location of cutadapt
######################################

`$cutadapt 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    next; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"cutadapt: error: At least one parameter needed: name of a FASTA or FASTQ file.\n","Test for cutadapt location");

######################################
# Testing the correct location of java
######################################

#Be careful, we redirect only the STDIN, some java send it in STDOUT
`$java 2>/tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    next unless $line =~ /java \[-options\] class/;
    $line =~ s/Syntaxe /Usage/; # For French locale usage
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"Usage: java [-options] class [args...]\n","Test for java location");

######################################
# Testing the correct location of GATK
######################################

`$GATK 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"##### ERROR ------------------------------------------------------------------------------------------\n","Test for GATK location");

######################################
# Testing the correct location of FASTQC
######################################

`$fastqc -h > /tmp/out.txt`;#We works with the STDOUT output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=<OUT>;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"            FastQC - A high throughput sequence QC analysis tool\n","Test for FASTQC location");

######################################
# Testing the correct location of fastx_trimmer
######################################

`$fastxTrimmer -h > /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    chomp $line;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"usage: fastx_trimmer [-h] [-f N] [-l N] [-t N] [-m MINLEN] [-z] [-v] [-i INFILE] [-o OUTFILE]","Test for fastx_trimmer location");

######################################
# Testing the correct location of bowtie-build
######################################

`$bowtieBuild 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    $line.=<OUT>;
    chomp $line;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"No input sequence or sequence file specified!
Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>","Test for bowtie-build location");


######################################
# Testing the correct location of bowtie2-build
######################################

`$bowtie2Build 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    $line.=<OUT>;
    chomp $line;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"No input sequence or sequence file specified!
Bowtie 2 version 2.2.9 by Ben Langmead (langmea\@cs.jhu.edu, www.cs.jhu.edu/~langmea)","Test for bowtie2-build location");

######################################
# Testing the correct location of tophat2
######################################

`$tophat2 2> /tmp/out.txt`;#We works with the STDERR output
open(OUT,"<", "/tmp/out.txt");
while (<OUT>) {
    $line=$_;
    $line.=<OUT>;
    chomp $line;
    last; 
}
close OUT;
unlink("/tmp/out.txt");

is($line,"tophat: 
TopHat maps short sequences from spliced transcripts to whole genomes.","Test for tophat2 location");




######################################
# Testing the correct location of cufflinks
######################################
#
#`$cufflinks/cufflinks 2> /tmp/out.txt`;#We works with the STDERR output
#open(OUT,"<", "/tmp/out.txt");
#while (<OUT>) {
#    <OUT>;
#    <OUT>;
#    $line=<OUT>;
#    last; 
#}
#close OUT;
#unlink("/tmp/out.txt");
#
#is($line,"Usage:   cufflinks [options] <hits.sam>\n","Test for cufflinks location");

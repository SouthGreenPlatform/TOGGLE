#!/usr/bin/perl -w
use strict;

#Will test if local config works correctly
use warnings;

use Test::More tests=>9; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Test::Deep;

use Data::Dumper;

use lib qw(../Modules/);

# Test we can use ../Modules
use_ok('localConfig') or exit;

use localConfig;

# We need to verify if the import pass weel, thus we parse the localConfig.pm to check what is expected to be recovered
#Is the comment enoughly explicit ?
open LOCALCONFIG, "../Modules/localConfig.pm";

my %dictLocation;
while (<LOCALCONFIG>){
	chomp $_;
	if(/our \$/){
		chop $_;
		my @line = split ("our ", $_);
		@line = split(" = ", $line[1]);
		$line[1]=~s/"//g;
		$dictLocation{$line[0]} = $line[1];
		
		
	}
	
}
close LOCALCONFIG;


# Verify 0 = 0;
is($java,$dictLocation{"\$java"},"Ok for java location bro");
is($bwa,$dictLocation{"\$bwa"},"Ok for bwa location");
is($picard,$dictLocation{"\$picard"},"Ok for picard location");
is($samtools,$dictLocation{"\$samtools"},"Ok for samtools location");
is($GATK,$dictLocation{"\$GATK"},"Ok for GATK location");
is($fastqc,$dictLocation{"\$fastqc"},"Ok for fastqc location");
#is($cuffdir,$dictLocation{"\$cuffdir"},"Ok for cufflink directory location");
is($cutadapt,$dictLocation{"\$cutadapt"},"Ok for cutadapt location");



__END__
my $currentJava=`which java`;
chomp $currentJava;
print $currentJava."\n$java\n";
ok ($java =~ m/$currentJava/, 'Ok for java location');

my $currentBwa=`which bwa`;
chomp $currentBwa;
print $currentBwa;
ok ($java =~ m/$currentBwa/, 'Ok for bwa location');

my $currentSamtools=`which samtools`;
chomp $currentSamtools;
ok ($java =~ m/$currentSamtools/, 'Ok for samtools location');

my $currentFastqc=`which fastqc`;
chomp $currentFastqc;
ok ($java =~ m/$currentFastqc/, 'Ok for fastqc location');

my $currentPicardTools=`locate picard-tools | head -n 1`;
chomp $currentPicardTools;
ok ($java =~ m/$currentPicardTools/, 'Ok for PicardTools location');

my $currentGatk=`locate GenomeAnalysisTK | grep "/usr/local" | head -n 1`;
chomp $currentGatk;
ok ($java =~ m/$currentGatk/, 'Ok for Gatk location');

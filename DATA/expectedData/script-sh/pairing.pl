#!/usr/bin/perl
 
use Data::Dumper;

my ($forward,$reverse, $forwardOut, $reverseOut, $singleOut)=@ARGV;

open IN, $forward;
open MATEF, ">$forwardOut";
open MATER, ">$reverseOut";
open SINGLE, ">$singleOut";

my %hashforward;
while (<IN>)
	{
	my $line = $_;
	chomp $line;
	$line =~ s/\/\d$//;
	my $next = <IN>;
	$next .= <IN>;
	$next .= <IN>;
	$hashforward{$line}=$next;
	}
close IN;
#print Dumper(\%hashforward);
#exit;

open IN2,"$reverse";

while (<IN2>)
	{
	my $line = $_;
	chomp $line;
	$line =~ s/\/\d$//;
	my $next = <IN2>;
	$next .= <IN2>;
	$next .= <IN2>; 
	#print $line;
	if (exists $hashforward{$line})
		{
		my $out = $line."/1\n".$hashforward{$line};
		print MATEF $out;
		my $out2 = $line."/2\n".$next;
		print MATER $out2;
		delete $hashforward{$line};
		}
	else
		{
		my $out2 = $line."/2\n".$next;
		print SINGLE $out2;
		}
	}

foreach my $remain (keys %hashforward)
	{
	my $out = $remain."/1\n".$hashforward{$remain};
	print SINGLE $out;
	}

close IN2;
close SINGLE;
close MATEF;
close MATER;

exit;
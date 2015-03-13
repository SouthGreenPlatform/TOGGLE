#/usr/bin/perl

use strict;
use warnings;
use Test::More 'no_plan'; #tests => 19; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .


use lib qw(../../Modules/);
use toolbox;

my $configFile='software.config.txt';

my $testCom="cat Log"; # print the date in the log, format MM/DD/YYYY
my $returnValue=toolbox::run($testCom);
ok ($returnValue== 0, 'Ok for toolbox::run return value');
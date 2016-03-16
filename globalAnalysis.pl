#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
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
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;

##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my ($nomprog)=$0=~/([^\/]+)$/;
unless ($#ARGV>=0)                    # if no argument given
{

  print <<"Mesg";

  perldoc $nomprog display the help

Mesg

  exit;
}

my %param = @ARGV;                     # get the parameters 
if (not defined($param{'-d'}) or not defined($param{'-c'}) or not defined($param{'-r'}) or not defined ($param{'-o'}))
{
  print <<"Mesg";

  ERROR: Parameters -d and -c and -r and -o are required.
  perldoc $nomprog display the help

Mesg
  exit;
}

if (defined $param{'-g'})
{
  system("toggleGenerator.pl -d $param{'-d'} -c $param{'-c'} -r $param{'-r'} -o $param{'-o'} -g $param{'-g'}  ") and die "Can't execute toggleGenerator.pl -d $param{'-d'} -c $param{'-c'} -r $param{'-r'} -o $param{'-o'}";

}
else
{
  system("toggleGenerator.pl -d $param{'-d'} -c $param{'-c'} -r $param{'-r'} -o $param{'-o'} ") and die "Can't execute toggleGenerator.pl -d $param{'-d'} -c $param{'-c'} -r $param{'-r'} -o $param{'-o'}";  
}

exit;

=head1 Name

globalAnalysis.pl 

=head1 Usage


globalAnalysis.pl-d DIR -c FILE -r FILE -o DIR -g FILE 

=head1 Required Arguments

      -d DIR    The directory containing initial files
      -c FILE   The configuration file
      -r FILE   The reference sequence (fasta)
      -o DIR    The directory containing output files
      
=head1 For RNAseq analysis

      -g FILE   The gff file containing reference annotations

=head1  Authors

Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

Copyright 2014-2015 IRD-CIRAD-INRA-ADNid

=cut
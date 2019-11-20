package sniplay;

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
use Data::Dumper;

use localConfig;
use toolbox;


################################################################################################
# sub sniplay::ped2fasta => to run ped2fasta
################################################################################################
# arguments :
# 	- pedIn : the PED (pedigree) file to analyze
#	- fastaOut : the fasta output file (mandatory)
################################################################################################
# return boolean :
#	- 1 if ped2fasta has runned correctly
#	- else 0
################################################################################################
sub ped2fasta
{
	toolbox::exportLog("ERROR: sniplay::ped2fasta : should get at least two arguments\n",0) if (@_ < 2);
	my ($pedIn,$fastaOut)=@_; 
	
	if (toolbox::sizeFile($pedIn)==1)
	{

		##DEBUG toolbox::exportLog("INFOS: fastqc::execution : running\n", 1);
	
		my %IUPAC =
		(
        '00'=> "?",
        'AA'=> "A",
        'CC'=> "C",
        'GG'=> "G",
        'TT'=> "T",
        'AG'=> "R",
        'GA'=> "R",
        'CT'=> "Y",
        'TC'=> "Y",
        'TG'=> "K",
        'GT'=> "K",
        'CG'=> "S",
        'GC'=> "S",
        'AT'=> "W",
        'TA'=> "W",
        'AC'=> "M",
        'CA'=> "M",
		);
	
		open(my $O,">$fastaOut");
		my $letter = "";
		#print $fastaOut;exit;
		open(my $P,$pedIn) or die "File does not exist";
		while(<$P>)
		{
			my $line = $_;
			$line =~s/\r//g;
			$line =~s/\n//g;
			my @infos = split("\t",$_);
			my $ind = $infos[0];
			print $O ">$ind\n";
			for (my $i = 6; $i <= $#infos; $i= $i+2)
			{
				my $code = $infos[$i].$infos[$i+1];
				if ($IUPAC{$code}){
					$letter = $IUPAC{$code};
				}
				print $O $letter;
			}
			print $O "\n";
		}
		close($P);
		close($O);

		return 1;
	}	
	else
	{
            toolbox::exportLog("ERROR: sniplay::ped2fasta : ABORTED\n",0);
        }
}
################################################################################################
# END sub sniplay::ped2fasta
################################################################################################






1;

=head1 NOM

package I<sniplay> 

=head1 SYNOPSIS

	use sniplay;

	sniplay::ped2fasta($pedFile,$fastaOut);                                                           

=head1 DESCRIPTION

This module is a set of functions related to sniplay software,  L<http://sniplay.southgreen.fr/>


=head2 Functions


=head3 ped2fasta()

This function allows to convert Ped to fasta alignement format.
2 arguments are required : the ped filename and the fasta file created.

Return 1 if ped2fasta has ran correctly else 0.

Example : 
C<sniplay::ped2fasta($pedFile,$fastaOut);> 	

=head1 AUTHORS

 Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut

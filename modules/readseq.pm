package readseq;

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

my %hash = (
1=>".ig",
2=>".gb",
3=>".nbrf",
4=>".embl",
5=>".gcg",
6=>".strider",
7=>".fitch",
8=>".fasta",
9=>".zuker",
10=>".olsen",
11=>".phylip2",
12=>".phylip",
13=>".seq",
14=>".pir",
15=>".msf",
16=>".asn",
17=>".nexus",
18=>".pretty",
19=>".xml",
20=>".blast",
21=>".scf",
22=>".aln",
23=>".fff",
24=>".gff",
25=>".ace"

);

################################################################################################
# sub readseq::readseq => to run readseq
################################################################################################
# arguments :
# 	- input : the input file to analyze
#	- output : the output file (mandatory)
################################################################################################
# return boolean :
#	- 1 if readseq has runned correctly
#	- else 0
################################################################################################
sub readseq
{
	my($fileIn,$fileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($fileIn)==1)
	{ ##Check if entry file exist and is not empty
	
          my $options="";
		  my $fmt = "";
          if ($optionsHachees)
          {
               $options=toolbox::extractOptions($optionsHachees);
          }
		  if  ($options !~ m/-f /)
		  {
				toolbox::exportLog("ERROR: readseq::readseq: no format number provided. ***$optionsHachees***\n",0);
		  }
		  my $extension;
		  if ($options =~ m/-f (\d+)/){
			$fmt = $1;
			$extension = $hash{$fmt};
		  }
		  
		  
		  
          my $command=$readseqjar." $options -o $fileOut$extension $fileIn";
          
          #Execute command
          if(toolbox::run($command)==1)
          {
               return 1;#Command Ok
          }
          else
          {
               toolbox::exportLog("ERROR: readseq::readseq : Uncorrectly done\n",0);
               return 0;#Command not Ok
          }
     }
     else
     {
        toolbox::exportLog("ERROR: readseq::readseq : The file $fileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
}



1;

=head1 NOM

package I<readseq> 

=head1 SYNOPSIS

	use readseq;

	readseq::readseq($input,$output);                                                           

=head1 DESCRIPTION

This module is a set of functions related to readseq software,  L<https://www.cog-genomics.org/readseq2>


=head2 Functions


=head3 readseq()

This function allows to convert between alignement format.
2 arguments are required : the filename and the file created.

Return 1 if readseq has ran correctly else 0.

Example : 
C<readseq::readseq($vcfFile,$pedOut);> 	

=head1 AUTHORS

 Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut

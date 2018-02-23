package plink;

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
# sub plink::vcf2ped => to run vcf2ped
################################################################################################
# arguments :
# 	- VCFIn : the VCF file to analyze
#	- pedOut : the ped output file (mandatory)
################################################################################################
# return boolean :
#	- 1 if vcf2ped has runned correctly
#	- else 0
################################################################################################
sub vcf2ped
{
	my($vcfFileIn,$pedFileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($vcfFileIn)==1)
	{ ##Check if entry file exist and is not empty
	
          #Check if the format is correct
          if (checkFormat::checkFormatVcf($vcfFileIn)==0)
          {#The file is not a VCF file
               toolbox::exportLog("ERROR: plink::vcf2ped : The file $vcfFileIn is not a VCF file\n",0);
               return 0;
          }
          

		  my $options="";
          
          if ($optionsHachees)
          {
               $options=toolbox::extractOptions($optionsHachees);
          }
          
          
          my $command=$plink." --recode --allow-extra-chr ".$options." --out ".$pedFileOut." --vcf ".$vcfFileIn;
          
          #Execute command
          if(toolbox::run($command)==1)
          {
               return 1;#Command Ok
          }
          else
          {
               toolbox::exportLog("ERROR: plink::vcf2ped : Uncorrectly done\n",0);
               return 0;#Command not Ok
          }
     }
     else
     {
        toolbox::exportLog("ERROR: plink::vcf2ped : The file $vcfFileIn is uncorrect\n",0);
        return 0;#File not Ok
     }
}
################################################################################################
# END sub plink::vcf2ped
################################################################################################






1;

=head1 NOM

package I<plink> 

=head1 SYNOPSIS

	use plink;

	plink::vcf2ped($vcf,$pedoutput);                                                           

=head1 DESCRIPTION

This module is a set of functions related to plink software,  L<https://www.cog-genomics.org/plink2>


=head2 Functions


=head3 vcf2ped()

This function allows to convert VCF to PED format.
2 arguments are required : the VCF filename and the ped file created.

Return 1 if vcf2ped has ran correctly else 0.

Example : 
C<plink::vcf2ped($vcfFile,$pedOut);> 	

=head1 AUTHORS

 Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut

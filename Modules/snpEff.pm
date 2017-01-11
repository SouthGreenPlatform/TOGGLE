package snpEff;

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

use strict;
use warnings;
use localConfig;
use toolbox;
use Data::Dumper;

sub snpEffAnnotation
{ #Will annotate a VCF file based on a already prepared db for SNPeff
    
    my ($vcf,$outputVcf, $optionsHachees)=@_;
    # the vcf file $vcf will be annotated using the infos from the already formatted database provided as option, and the resulting vcf will be outputted in $outputName
    #The options may be the gatk compatibility (-o gatk) and the config file positions
    
    my $options="";

    if ($optionsHachees)
	{
        $options=toolbox::extractOptions($optionsHachees);		##Get given options
    }
    my $command=$snpEff." ann ".$options." ".$vcf." > ".$outputVcf;
    
    #Execute command
    if(toolbox::run($command)==1)		##The command should be executed correctly (ie return) before exporting the log
    {
        return 1;
    }
    else
    {
        toolbox::exportLog("ERROR: snpEff::snpEffAnnotation : ABORTED\n",0);		# snpEffAnnotation has not been correctly done
        return 0;
    }
}


1;

=head1 NAME

    Package I<snpEff> 

=head1 SYNOPSIS

	use snpEff;
    
	snpEff::snpEffAnnotation($vcf,$outputVcf, $optionsHachees)
	
=head1 DESCRIPTION

    Package SNPeff (Cingolani et al, 2012, http://http://snpeff.sourceforge.net/ ) is a software package for refining SNP/Indel information based on annotation.

=head2 FUNCTIONS

=head3 snpEff::snpEffAnnotation

This module is the core function of SNPeff
It will annotate a given VCF using an annotation database already formatted for SNPeff.
The resulting VCF will have the annotation information (and effects of SNP) within it.


=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut

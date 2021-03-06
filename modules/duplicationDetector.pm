package duplicationDetector;

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
use localConfig;
use toolbox;
use Data::Dumper;
use checkFormat;

##duplicationDetector running
sub execution
{
    my($vcfFileIn,$fileOut,$optionsHachees)=@_;
    
    if (checkFormat::checkFormatVcf($vcfFileIn) != 1) #No mate file, sequences are single ends
    {
        toolbox::exportLog("ERROR: duplicationDetector::execution: THE file $vcfFileIn is not a correctly formatted VCF.", 0);
    }
    
    #Options
    my $options = "";
    if ($optionsHachees)
    {
        $options=toolbox::extractOptions($optionsHachees);		##Get given options
    }
	
    #Command
    #Basic command line
    my $command = $duplicationDetector." ".$options." -i ".$vcfFileIn." -o ".$fileOut." ";
    
    #Execute command
    if(toolbox::run($command)==1)		## if the command has been executed correctly, export the log
    {
        return 1;
    }
    else
    {
        toolbox::exportLog("ERROR: duplicationDetector::execution : ABORTED\n",0);
        return 0;
    }

}



1;
=head1 NAME

    Package I<duplicationDetector>

=head1 SYNOPSIS

	use duplicationDetector
	
	duplicationDetector::execution
	

=head1 DESCRIPTION

    Package DuplicationDetector ( http:// ) is a software package for detecting duplications through the analysis of a VCF file

=head2 FUNCTIONS

=head3 duplicationDetector::execution

This module will identify the duplication in a given VCF file

=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD, ADNid and South Green development platform
Written by Christine Tranchant, Julie Orjuela, Sebastien Ravel and Francois Sabot

=head1 SEE ALSO

L<http://toggle.southgreen.fr/>

=cut

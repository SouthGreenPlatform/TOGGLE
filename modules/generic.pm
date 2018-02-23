package generic;

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


sub generic {
	my($fileIn,$fileOut,$optionsHachees)=@_;
	if (toolbox::sizeFile($fileIn)==1)
	{ ##Check if entry file exist and is not empty
		
		my $options="";
		
		if ($optionsHachees)
		{
			$options=toolbox::extractOptions($optionsHachees);
		}
		else
		{
				toolbox::exportLog("ERROR: generic::generic : No options and commands provided, we cannot execute any command\n",0);
				return 0;
		}
		
		#The generic command system will transform the FILEIN text by the correct FILENAME
		$options =~ s/FILEIN/$fileIn/i;
		
		#The generic command system will transform the FILEOUT text by the correct FILENAME
		$options =~ s/FILEOUT/$fileOut/i;
		my $command=$options;
		
		
		#Execute command
		if(toolbox::run($command)==1)
		{
			return 1;#Command Ok
		}
		else
		{
			toolbox::exportLog("ERROR: generic::generic : Uncorrectly done\n",0);
			return 0;#Command not Ok
		}
	}
	else
	{
	   toolbox::exportLog("ERROR: generic::generic : The file $fileIn is uncorrect\n",0);
	   return 0;#File not Ok
	}
}

1;
=head1 NAME

    Package I<generic>

=head1 SYNOPSIS

	use generic;

	generice::generic($fileIn, $fileOut, $optionsHachees)

=head1 DESCRIPTION

    This package allows the user to launch ANY command in the TOGGLe pipeline. The FILEIN text in the command will be replaced by the $fileIn value, as well as the FILEOUT text by the $fileOut value
	
	/!\ NO FORMAT CONTROL ASSOCIATED /!\

=head2 FUNCTIONS

=head3 generic::generic

	launch any command. requires $fileIn, $fileOut, $optionsHachees


=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Christine Tranchant, Cecile Monat, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://toggle.southgreen.fr/>

=cut
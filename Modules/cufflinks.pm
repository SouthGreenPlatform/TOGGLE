package cufflinks;

###################################################################################################################################
#
# Copyright 2014 IRD-CIRAD
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
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################

use strict;
use warnings;

use localConfig;
use toolbox;
use Data::Dumper;



sub execution
{

    my ($cuffliksDirOut, $bamFileIn,$annotGffFile,$optionsHachees)=@_;

    if ((toolbox::sizeFile($bamFileIn)==1) and (toolbox::sizeFile($annotGffFile)==1))             ##Check if the bamfile and annotGffFile exist and are not empty
    {   


        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.configRnaSeq
        ## DEBUGG toolbox::exportLog("DEBUG: cufflinks::execution : execution option equals to $options",1);

        my $command=$cufflinks.$options." "." -o ".$cuffliksDirOut." -g ".$annotGffFile." ".$bamFileIn;				# command line to execute cutadapt
 
        ##DEBUG
        toolbox::exportLog("INFOS: cufflinks::execution : $command\n",1);

        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: kufflinks : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: cufflinks : ABBORTED\n",0);
            return 0;
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: cufflinks::execution : Problem with the files\n",0);
        return 0;
    }
    
    
}


1;

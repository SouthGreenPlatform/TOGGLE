package eaUtils;

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
use checkFormat;
use Data::Dumper;
use Switch;

################################################################################################
# sub fastq-stats:stats on fastq file
################################################################################################
# arguments : fastq file to analyze
# Returns boolean (1 if the execution is correctly done else 0)
################################################################################################
sub fastqStats
{
    my($fastqFileIn,$fileOut,$optionsHachees)=@_;
    if (toolbox::sizeFile($fastqFileIn)==1 and checkFormat::checkFormatFastq($fastqFileIn)==1)             ##Check if the fastqfile exist and is not empty
    {
        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.config

        my $command=$fastqStats." ".$options." ".$fastqFileIn." > ".$fileOut; ## Command initialization
        
        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: eaUtils::fastqStats : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: eaUtils:fastqStats : ABORTED\n",0);
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: eaUtils::fastqStats : Problem with the file $fastqFileIn (does not exist or invalid fastq format)\n",0);
    }
    
    
}


1;
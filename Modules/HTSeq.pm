package HTSeq;

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

use localConfig;
use toolbox;
use Data::Dumper;



sub htseqCount
{
    my ($bamFileIn,$htseqFileOut,$annotGffFile,$optionsHachees)=@_;
    if (toolbox::sizeFile($bamFileIn)==1 and toolbox::sizeFile($annotGffFile)==1)            ##Check if the bamfileIn exist and is not empty
    {
        
        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.config
        ## DEBUG toolbox::exportLog("DEBUG: HTSeq::htseqcount : htseqcount option equals to options $options",1);
        
        my $command;
        my ($bamFileOut,$readGroup) = pairing::extractName($bamFileIn);
        
        # Verify if $bamFileIn is a bam or sam
        if ($bamFileIn =~ m/.sam$/) # if the file type is sam
            {
                $command="$htseqcount $options $bamFileIn $annotGffFile > $htseqFileOut"; ## Command initialization
            }
        if ($bamFileIn =~ m/.bam$/) # if the file type is bam, samtools view is used. HTseq-count does not count reads if bam is used.
            {
                #$command=samTools::samToolsView($bamFileIn,$bamFileOut, XX?); # pbm avec $optionsHachees
                $command="samtools view -h  $bamFileIn"; ## Command initialization to samtools view
                if (toolbox::run($command)==1)
                {
                    toolbox::exportLog("INFOS: samTools::samToolsView : correctly done\n",1);
                    $command="$htseqcount $options $bamFileIn $annotGffFile > $htseqFileOut" ;
                }
                else
                {
                    toolbox::exportLog("ERROR: samTools::samToolsView : ABORTED\n",0);
                }
            }
        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: HTSeq::htseqCount : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: HTSeq::htseqCount : ABBORTED\n",0);
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: HTSeq::htseqCount : Problem with the file $bamFileIn (don't exist ou invalid bam format)\n",0);
    }
}


1;



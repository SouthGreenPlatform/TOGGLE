package cufflinks;

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



sub execution
{

    my ($cuffliksDirOut, $bamFileIn,$annotGffFile,$optionsHachees)=@_;

    if ((toolbox::sizeFile($bamFileIn)==1) and (toolbox::sizeFile($annotGffFile)==1))             ##Check if the bamfile and annotGffFile exist and are not empty
    {   


        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.configRnaSeq
        ## DEBUGG toolbox::exportLog("DEBUG: cufflinks::execution : execution option equals to $options",1);

        my $command=$cufflinks.$options." "." -o ".$cuffliksDirOut." -g ".$annotGffFile." ".$bamFileIn;	# command line to execute cutadapt
 
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

sub cuffmerge
{
 
    my ($cuffmergeDirOut, $RefFastaFileIn,$gffFile,$assemblyGTFlist,$optionsHachees)=@_;

    if ((toolbox::sizeFile($RefFastaFileIn)==1) and (toolbox::sizeFile($gffFile)==1) and (toolbox::sizeFile($assemblyGTFlist)==1))             ##Check if the bamfile and annotGffFile exist and are not empty
    {   


        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.configRnaSeq
        ## DEBUGG toolbox::exportLog("DEBUG: cufflinks::cuffmerge : cuffmerge option equals to $options",1);

        my $command=$cuffmerge.$options." -o ".$cuffmergeDirOut." -p 8 "." -s ".$RefFastaFileIn." -g ".$gffFile." ".$assemblyGTFlist;				# command line to execute cutadapt
 
        ##DEBUG
        toolbox::exportLog("INFOS: cufflinks::cuffmerge : $command\n",1);

        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: cuffmerge : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: cuffmerge : ABBORTED\n",0);
            return 0;
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: cufflinks::cuffmerge : Problem with the files\n",0);
        return 0;
    }
    
    
}


sub cuffdiff

{
 
    my ($cuffdiffDirOut,$assemblyGTF,$optionsHachees,@bamFileIn)=@_;
    
    
    foreach my $bam (@bamFileIn) {
        if (toolbox::checkSamOrBamFormat($bam)==0)
        
        {
            toolbox::exportLog("ERROR : cuffdiff : $bam is neither sam nor bam\n",1);
            return 0;
        }
        
    }

    if ((toolbox::sizeFile($assemblyGTF)==1))             ##Check if the annotGffFile exist and are not empty
    {   


        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.configRnaSeq
        ## DEBUGG toolbox::exportLog("DEBUG: cufflinks::cuffmerge : cuffmerge option equals to $options",1);

        my $command=$cuffdiff.$options." -o ".$cuffdiffDirOut." -p 8 ".$assemblyGTF." ".@bamFileIn;				# command line to execute cutadapt
 
        ##DEBUG
        toolbox::exportLog("INFOS: cufflinks::cuffdiff : $command\n",1);

        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: cuffdiff : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: cuffdiff : ABBORTED\n",0);
            return 0;
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: cufflinks::cuffdiff : Problem with the files\n",0);
        return 0;
    }
    
    
}


1;

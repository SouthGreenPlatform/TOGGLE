package fastxToolkit;

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


################################################################################################
# sub fastxTrimmer : to remove the n first pb of sequences in  fastq file
################################################################################################
# arguments : fastq file to analyze
# Returns boolean (1 if the execution is correctly done else 0)
################################################################################################
sub fastxTrimmer
{
    my($fastqFileIn,$fastqFileOut,$optionsHachees)=@_;
    if (toolbox::sizeFile($fastqFileIn)==1 and toolbox::checkFormatFastq($fastqFileIn)==1)             ##Check if the fastqfile exist and is not empty
    {
        my $options=toolbox::extractOptions($optionsHachees, " ");  ##Get given options by software.config
        ## DEBUG toolbox::exportLog("DEBUG: fastxToolkit::fastxTrimmer : fastxTrimmer option equals to $options",1);
        my $command=$fastxTrimmer." ".$options." -i ".$fastqFileIn." -o ".$fastqFileOut; ## Command initialization
        
        # Command is executed with the run function (package toolbox)
        if (toolbox::run($command)==1)
        {
            toolbox::exportLog("INFOS: fastxToolkit::fastxTrimmer : correctly done\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: fastxToolkit:fastxTrimmer : ABBORTED\n",0);
        }
        
    }
    else
    {
        toolbox::exportLog("ERROR: fastxToolkit::fastxTrimmer : Problem with the file $fastqFileIn (don't exist ou invalid fastq format)\n",0);
    }
    
    
}
################################################################################################
# END sub fastxTrimmer 
################################################################################################

1;

=head1 NAME

package I<fastxToolkit> 

=head1 SYNOPSIS

	use fastxToolkit;

	fastxToolkit::fastxTrimmer($fastqFilename);
	
=head1 DESCRIPTION

The FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/index.html) is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.


=head2 Functions

=head3 fastxToolkit::fastxTrimmer()

This function execute the fasttx_trimmer program that remove the n first pb of sequences in  a fastq file.
One argument is required : the fastq filename
- lines to print in a log or error file
- boolean value to control the type of the log and of the error. 0 : print logs into error log file and die.
1 : print logs into log file. 2 : print logs into error log file and the program continue.

The parameter returned is a boolean : 1 if the analysys has been correctly done else 0.

Example : 
<fastxToolkit::fastxTrimmer($fastqFilename);>


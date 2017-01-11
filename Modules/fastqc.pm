package fastqc;

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
use Data::Dumper;

use localConfig;
use toolbox;


################################################################################################
# sub fastqc::execution => to run fastqc
################################################################################################
# arguments :
# 	- fileIn : the fastq file to analyze
#	- dirOut : the directory that contains all the files generated by fastqc
################################################################################################
# return boolean :
#	- 1 if fastqc has runned correctly
#	- else 0
################################################################################################
sub execution
{
	toolbox::exportLog("ERROR: fastqc::execution : should get at least two arguments\n",0) if (@_ < 2);
	my ($filein,$dirOut,$optionsHachees)=@_; 
	
	##DEBUG toolbox::exportLog("INFOS: fastqc::execution : running\n", 1);
	
	my $options="";
        if ($optionsHachees)
	{
            $options=toolbox::extractOptions($optionsHachees);		##Get given options
        }

	
    	my $cmd_line=$fastqc." ".$options." --noextract -o $dirOut $filein"; ## Force option noextract (Do not uncompress the output file after creating it)
	
	 if(toolbox::run($cmd_line)==1)		## if the command has been excuted correctly, export the log
	{
            return 1;
        }
	else
	{
            toolbox::exportLog("ERROR: fastqc::execution : ABORTED\n",0);
        }
}
################################################################################################
# END sub fastqc::execution
################################################################################################






1;

=head1 NOM

package I<fastqc> 

=head1 SYNOPSIS

	use fastqc;

	fastqc::execution($fastqFile,$fastqcDir);                                                           

=head1 DESCRIPTION

This module is a set of functions related to fastqc software,  L<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>


=head2 Functions


=head3 execution()

This function execute the fastqc software to analyze a fastq file and a directory (that contains the files generated by fastqc) is created.
2 arguments are required : the fastq filename and the directory name created.

Return 1 if fastqc has ran correctly else 0.

Example : 
C<fastqc::execution($fqFile,$fastqcDir);> 	

=head1 AUTHORS

 Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut

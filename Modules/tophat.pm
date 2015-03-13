package tophat;



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
use lib qw(.);
use localConfig;
use toolbox;
use Data::Dumper;
###################################
##############################################
##TOPHAT
##
####Create Index
sub indexRef
    {
    my ($fastaRef) = @_;
    my $prefixRef = $fastaRef;
    #if ($fastaRef  =~ /^(.*)\.fa/)
    #    { $prefixRef=$1; }
        
    my $command= "$bowtieIndex -f $fastaRef $prefixRef ";
      
    toolbox::run($command);    
    }
    

####Run Tophat
sub tophatRun
    {
    my ($fastaRef, $fastqPair1, $fastqPair2,$annotGff)=@_;
    #Picking up library name
    my $library;
    if ($fastqPair1  =~ /^(.*)\.fastq/)
        { $library=$1; }
    #$fastaRef =~ s/\.fa$//; # to avoid the multiple .fa.fa for file... I hate tophat
    
    #Generating options
    my $options = "-o $library";
    if (defined $annotGff)
        {# A GFF of already known transcript has been submitted as an option
        #Checking if annotGFF is a true GFF
        
        #Adding to options
        $options.=" -G $annotGff";
        }
    
    
    #Checking if fastqPair1 and 2 are truely fastq
    
    #Checking if fastaRef has been already indexed
    
    
    #Running command
    my $command= "$tophat $options $fastaRef $fastqPair1 $fastqPair2";

    toolbox::run($command); 
    }
1;
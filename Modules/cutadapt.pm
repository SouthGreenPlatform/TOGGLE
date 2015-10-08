package cutadapt;



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
use Data::Dumper;

use localConfig;
use toolbox;


################################################################################################
# sub cutadapt::createConfFile => create the configuration file specific to cutadapt
################################################################################################
# ex : cutadapt $(<cutadapt.conf) MID54.fastq > MID54.cleaned.fastq
################################################################################################
# arguments :
# 	- fileAdaptator : file of adaptator (one sequence by line)
#	- fileConf : the name of the cutadapt configuration file created by this sub
#	- optionref : the reference of option hash that contains the cutadapt options.
# Theses options will be written in the cutadapt configuration file 
################################################################################################
# return boolean :
#	- 1 if the file has been correctly created
#	- else 0
################################################################################################
sub createConfFile
{
    
    my ($fileAdaptator, $fileConf, $optionref)=@_;							# recovery of arguments
    
    open(CONF, ">", $fileConf) or toolbox::exportLog("ERROR: cutadapt::createConfFile : Can't open the configuration file $fileConf $!\n",0); 	# opening the configuration file to fill
    open(ADAPTATOR, "<", $fileAdaptator) or toolbox::exportLog("ERROR: cutadapt::createConfFile : Cannot open the adaptator file $fileAdaptator $!\n",0);	# opening the adaptators file
    while (my $seq=<ADAPTATOR>)
    {
	next if ($seq=~m/^$/);										# next if empty line
	chomp $seq;											# remove "\n" at the end of the line
        print CONF "-b $seq\n";										# print in the configuration file the "-b" parameter and the adaptators sequence cooresponding
        
        $seq = reverse $seq;										# do the reverse adaptators sequence
        $seq =~ tr/AaCcGgTt/TtGgCcAa/;									# do the complement adaptators sequence
        print CONF "-b $seq\n";										# print in the configuration file the "-b" parameter and the reverse adaptators sequence cooresponding
    }
    
    ##DEBUG : print Data::Dumper::Dumper(\%$optionref);
    my %optionsRef = %$optionref;
    foreach my $parameter (keys %optionsRef)
    {
	print CONF "$parameter $optionsRef{$parameter}\n";							# print in the configuration file the parameter and the options cooresponding to
    }

    close CONF;
    close ADAPTATOR;
    
    if (toolbox::sizeFile($fileConf))
    {									# Check if the configuration file is not empty
	toolbox::exportLog("INFOS: cutadapt::createConfFile : the configuration file $fileConf has been correctly created\n",1);
	return 1;
    }
    else
    {
	toolbox::exportLog("ERROR: cutadapt::createConfFile : A probleme has occured during the creation of the configuratin file $fileConf\n",0);
    }
    
}
################################################################################################
# END sub cutadpt::createConfFile 
################################################################################################



################################################################################################
# sub cutadapt::exectution => to run cutdapt
################################################################################################
# arguments :
# 	- fileIn : the fastq file with the sequences to clean
#	- fileConf : the name of the cutadapt configuration file created by the sub
# cutadapt::createConfFile
#	- fileOut : the fastq file generated by cutadapt with the sequences cleaned
################################################################################################
# No parameters returned, parameters will be returned and errors managed by toolbox::run
################################################################################################
sub execution
{

    toolbox::exportLog("ERROR: cutadapt::execution : should get exactly three arguments\n",0) if (@_ !=3 );	# Check if the number of arguments is good

    my ($fileIn,$fileConf,$fileOut)=@_;	 # recovery of arguments

    ## Lancement de cutadapt        
    open(ADAPT, "<", $fileConf) or toolbox::exportLog("ERROR: cutadapt::execution : Cannot open the file $fileConf $!\n",0);
    my $adaptors=" ";
    while (<ADAPT>)
    {
	chomp($_);
	$adaptors=$adaptors.$_." ";
    }
    close ADAPT;
    
    my $cmd_line=$cutadapt." ".$adaptors." ".$fileIn." -o ".$fileOut;				# command line to execute cutadapt

    toolbox::run($cmd_line);									# tool to execute the command line
}
################################################################################################
# END sub cutadapt::exectution => to run cutdapt
################################################################################################

1;


=head1 NAME

package I<cutadapt> 

=head1 SYNOPSIS

	use cutadapt;

	cutadapt::createConfFile($adaptatorFile,$cutadaptFileConf,$option_prog{'cutadapt'});

	cutadapt::execution($fasqFilename,$cutadaptFileConf,$fastqCleanedFile);
	
=head1 DESCRIPTION

This module is a set of functions related to cutadapt software, L<https://code.google.com/p/cutadapt>


=head2 Functions

=head3 cutadapt::createConfFile()

This function create the configuration file used by cutadapt. It is required 3 arguments :
- a file with the adaptators sequences
- the name of the configuration file generated by this fonction and used by cutadapt
- a hash containing the options used by cutadapt and written in the file created.

This fonction parse the adaptators sequences file, generate the reverse sequences and
groups all the sequences in the configuration file with the options given by the hash

This function returns 1 if the file have been created correctly else 0.

Example : 
C<cutadapt::createConfFile($adaptatorFile,$cutadaptConf,$option_prog{'cutadapt'});>


=head3 cutadapt::execution()

This function execute the cutadapt software and generate a cleaned sequences file (without adaptators and minimal size).

It is required 3 arguments in input :
- the sequence file (fastq format),
- the configuration file generated by the fonction cutadapt::createConfFile
- the name of file (fastq format) generated with the sequences cleaned 

Return 0,1,2 with the sub toolbox::run

Example: 
C<cutadapt::exec($fqFile,$cutadaptConf,$fastqCutadapt, $logFile) ;> 	


=head1 AUTHORS

 Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>

=cut

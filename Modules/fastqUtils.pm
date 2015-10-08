package fastqUtils;



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
use Data::Translate;#To convert ASCII to decimal, required!
use Data::Dumper;



#########################################
#checkEncodeByASCIIcontrol
#Check the FASTQ format of a given file. 
#########################################
sub checkEncodeByASCIIcontrol
{ 
    my ($fileName)=@_;		# recovery of informations
    if (toolbox::checkFormatFastq($fileName)== 1)		# check if the file you give is a FASTQ file
    {
	toolbox::exportLog("INFOS: fastqUtils::checkEncodeByASCIIcontrol : The file $fileName is a fastq file\n",1);
    }
    else
    {
	toolbox::exportLog("ERROR: fastqUtils::checkEncodeByASCIIcontrol : The file $fileName is not a fastq file\n",0);
    }
    my $translator = new Data::Translate;		# ???
    return 0 unless toolbox::readFile($fileName);		#Check if file readable. Stop if not

    open(IN, "<", $fileName) or toolbox::exportLog("ERROR: fastqUtils::checkEncodeByASCIIcontrol : Cannot open the file $fileName $!\n",0);		# open the file in, and if it can't, return and error message
    <IN>;		#skip the line of the name of the sequence
    <IN>;		#skip the line of the sequence
    <IN>;		#skip the line of informations 

    my $qualityASCII=<IN>; #picking up the quality line
    chomp $qualityASCII;
    
    my @listASCII=split //,$qualityASCII;		# split each symbol
    ##DEBUG print "@listASCII","\n";
    my $s;		#Requested for Data::Translate to function
    my $phred33Control = 1;		# variable of control, will be change if PHRED is not 33
    while (@listASCII)		# for each quality symbol of the sequence ... 
    {
	my (@ASCII) = shift @listASCII;
	my $Decimal;
	($s,$Decimal)=$translator->a2d(@ASCII);		# ...translate from ascii to decimal value
	$phred33Control = 0 if $Decimal > 74;		#Phred33 maximal value is 34 + 40 = 73, if this digit is reached then the file is not in PHRED 33, so changing the variable of control
    }
    return $phred33Control;		#Return 1 if PHRED 33, O if not
}
#########################################
#changeEncode
#will change a FASTQ PHRED format in another
#########################################
#TODO: create a FASTQ line/sequence parser ?
sub changeEncode
{
    my ($fileIn,$fileOut,$formatInit,$formatOut)=@_;#The format must be numerical (eg 33 or 64)
    if (toolbox::checkFormatFastq($fileIn)== 1)		# check if the file you give is a FASTQ file
    {
	toolbox::exportLog("INFOS: fastqUtils::changeEncode : The file $fileIn is a fastq file\n",1);
    }
    else
    {
	toolbox::exportLog("ERROR: fastqUtils::changeEncode : The file $fileIn is not a fastq file\n",0);
    }
    
    open (IN,"<",$fileIn) or toolbox::exportLog("ERROR: fastqUtils::changeEncode : Cannot open the file $fileIn $!\n",0); #Can read file
    open(OUT,">", $fileOut) or toolbox::exportLog("ERROR: fastqUtils::changeEncode : Cannot create the file $fileOut\n$!\n",0); #Create the output and verify if any error
    while (my $line = <IN>)	#Pick up the Sequence Name Line from the infile
    {
	$line .= <IN>; # Add the IUPAC line
	$line .= <IN>; # Add the '+' line
	my $qualityLine = <IN>;# recovery of quality line
	chomp $qualityLine;
	
	my $finalQuality;
	if ($formatInit == 64 and $formatOut == 33) 	#From Illumina old PHRED 64 to Sanger PHRED 33
	{
	    $finalQuality=convertLinePHRED64ToPHRED33($qualityLine);# execution of the conversion
	}
	elsif ($formatInit == 33 and $formatOut == 64)	#From Sanger PHRED 33 to Illumina old PHRED 64 
	{
	    $finalQuality=convertLinePHRED33ToPHRED64($qualityLine);# execution of the conversion
	}  
	$line.=$finalQuality."\n";
	print OUT $line; #Outputting in the outfile
    }
    close IN;
    close OUT;

    toolbox::exportLog("INFOS: fastqUtils::changeEncode : The file $fileIn has been re-encoded from a PHRED $formatInit scale to a PHRED $formatOut scale in $fileOut\n",1);
    return 1;
}

#########################################
#convertLinePHRED64ToPHRED33
#From a PHRED 64 quality line, will convert in PHRED33
#########################################
sub convertLinePHRED64ToPHRED33
{
    my ($initialQuality)=@_;		# recovery of informations
    my @listOfQuality = split //, $initialQuality;		# split each symbol
    foreach (@listOfQuality)		# for each quality symbol 
    {
	tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;		# convert the scale from PHRED64 to PHRED33
    }
    my $finalQuality = join("",@listOfQuality);		# rejoin each symbol to make a quality line
    return $finalQuality;		# return the new quality line
}
#########################################
#convertLinePHRED33ToPHRED64
#From a PHRED 33 quality line, will convert in Phred 64
#########################################
sub convertLinePHRED33ToPHRED64
{
    my ($initialQuality)=@_;		# recovery of informations
    my @listOfQuality = split //, $initialQuality;		# split each symbol
    foreach (@listOfQuality)		# for each quality symbol
    {
	$_  = chr (ord($_)+31);		#convert the scale from PHRED33 to PHRED64
    }
    my $finalQuality = join("",@listOfQuality);		# rejoin each symbol to make a quality line
    return $finalQuality;		# return the new quality line
}
1;

=head1 NAME

    Package I<fastqUtils> 

=head1 SYNOPSIS

        use fastqUtils;
    
        fastqUtils::checkEncodeByASCIIcontrol ($fileName);
    
        fastqUtils::changeEncode ($fileIn,$fileOut,$formatInit,$formatOut);
	
	fastqUtils::convertLinePHRED64ToPHRED33 ($initialQuality);
	
	fastqUtils::convertLinePHRED33ToPHRED64 ($initialQuality);

=head1 DESCRIPTION

Package fastqUtils is a set of modules which deals with issues relative to FASTQ format files

=head2 FUNCTIONS
 

=head3 fastqUtils::checkEncodeByASCIIcontrol

This module check the FASTQ format of a given file
It takes only one argument, the file you want to know the encoding



=head3 fastqUtils::changeEncode

This module change a FASTQ PHRED format in another
it takes four arguments: the file you want to change encoding, the name of the outpur file, the PHRED format of your initial file, the PHRED format you want in the output file



=head3 fastqUtils::convertLinePHRED64ToPHRED33

This module change a PHRED 64 quality line in PHRED 33
It takes only one argument, the initial quality of your file



=head3 fastqUtils::convertLinePHRED33ToPHRED64

This module change a PHRED 33 quality line in PHRED 64
It takes only one argument, the initial quality of your file



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut

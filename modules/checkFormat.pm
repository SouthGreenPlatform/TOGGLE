package checkFormat;

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
use Exporter;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use lib qw(.);
use localConfig;
use toolbox;



################################################################################################
# sub checkFormatFastq => check if a file is really a FASTQ file
################################################################################################
# arguments : filename to analyze
# Returns boolean (1 if the format is fastq else 0)
################################################################################################

sub checkFormatFastq
{

    my $notOk = 0;                  # counter of error(s)
    my ($fileToTest) = @_;          # recovery of file to test

    #Checking the beginning and end structure
    my ($beginLines, $endLines);
    if ($fileToTest =~ m/gz$/)
    { # The file is in gzipped format
        #using zcat command for head and tail
        $beginLines = `zcat $fileToTest | head -n 4`;
        $endLines = `zcat $fileToTest | tail -n 4`;
    }
    else
    {
        $beginLines = `head -n 4 $fileToTest`;
        $endLines = `tail -n 4 $fileToTest`;
    }
    chomp $beginLines;
    chomp $endLines;

    if ($beginLines !~ m/^@/ and $endLines !~ m/^@/) # The file is not containing a 4 lines sequence in FASTQ format
    {
        toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : Number of lines is not a multiple of 4 in file $fileToTest, thus not a FASTQ file.\n",0);
    }


    open (my $inputHandle, $fileToTest) or toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : Cannot open the file $fileToTest\n$!\n",0); # open the file to test

    my  @linesF1=();
    my $comp=0;
    my $countlines=0;
    my $stop=0;

    #If $fileToTest is in gzip format
    if($fileToTest =~ m/\.gz$/)
    {
        $inputHandle = new IO::Uncompress::Gunzip $inputHandle or toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : Cannot open the gz file $fileToTest: $GunzipError\n",0);
    }

    while ((my $line = <$inputHandle>))                                           # scanning file and stocking in an array the four lines of a read.
    {
        chomp $line;
        $countlines++;

        if ($comp<3)
        {
            $comp++;
            push (@linesF1,$line);
        }
        else                                                            # Completing block, treatment starts here.
        {
            $stop++;
            push (@linesF1,$line);

            my $i=0;
            while ( ($i<=$#linesF1) and ($notOk <=1))                 # treatment of a block containing four lines of file and stop if 20 errors found.
            {

                my $idLine=$linesF1[$i];
                my $fastaLine=$linesF1[$i+1];
                my $plusLine=$linesF1[$i+2];
                my $qualityLine=$linesF1[$i+3];
                my $nbIDLine=$countlines-3;
                my $nbLineFasta=$countlines-2;
                my $nbPlusLine=$countlines-1;
                my $nbQualityLine=$countlines;

                if (($idLine=~m/^$/) and ($plusLine=~m/^$/))            # if the ID line and the "plus" line are not empty ...
                {
                    toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : The file $fileToTest is not a FASTQ file => The ID infos line $nbIDLine is not present.\n",0);
                    $notOk++;                                           # one error occured, so count it
                }

                elsif ( (($idLine=~m/^\@.*/) or ($idLine=~m/^\>.*/) ) and ($plusLine=~m/^\+/) )   # if ID ligne is not empty and it starts by "@" or ">" and the
                # plus line has a "+", the block of four lines ($i to $i+3) is traited.
                {
                    if ( length($fastaLine) == length($qualityLine) )   # comparing the fasta line and the quality line lengths.
                    {
                        my @fasta = split //, $fastaLine;
                        foreach my $nucleotide (@fasta)
                        {
                            if ($nucleotide!~m/A|T|G|C|a|g|t|c|N|n/)    # checking nucleotides in the fasta line.
                            {
                                toolbox::exportLog ("ERROR: checkFormat::checkFormatFastq : Not basic IUPAC letter, only ATGCNatgcn characters are allowed: unauthorized characters are in the line $nbLineFasta of $fileToTest.\n",0);
                                $notOk++;
                            }
                        }
                    }
                    else 												# error if fasta line length and quality line length are differents.
                    {
                        toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : Fasta line $nbLineFasta has not the same length than quality line $nbQualityLine in file $fileToTest.\n",0);
                        $notOk++;
                    }
                }

                else													#error if the ID line do not start with @ or >.
                {
                    toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : ID line has to start with @ or > in line $nbIDLine of file $fileToTest.\n",0);
                    $notOk++;
                }
                $i=$i+4; 												# jumping to next read.
            }

            last if ($stop==20000);                                    # stoping treatment if 50000 reads were analysed.

            undef @linesF1;
            $comp=0;
        }
        next;
    }

    if (($notOk == 0))                    						# if any error occured in the file, the format is right.
    {
        toolbox::exportLog("INFOS: checkFormat::checkFormatFastq : The file $fileToTest is a FASTQ file.\n",1);
        return 1;
    }
    else                                						# if one or some error(s) occured on the file, the fastq format is not right.
    {
        toolbox::exportLog("ERROR: checkFormat::checkFormatFastq : Invalid FASTQ requirements in file $fileToTest.\n",0);
    }

    close $inputHandle;

}
################################################################################################
# END sub checkFormatFastq
################################################################################################

################################################################################################
# sub checkFormatSamOrBam => verifying the SAM/BAM format based on samtools view system
# samtools view gave an error in case of non SAM or BAM format
################################################################################################
# arguments : filename to analyze
# Returns boolean (1 if the fileformat is sam, 2 bam and 0 neither bam or sam)
################################################################################################
sub checkFormatSamOrBam
{

    my ($samFile)=@_;

    toolbox::existsFile($samFile); # Will check if the submitted file exists

    #Is the file sam of bam through the binary mode? Requested for the -S option in samtools view
    my ($inputOption,$binary);
    if (-B $samFile) #The file is a binary BAM file
    {
        $inputOption = ""; #no specific option in samtools view requested
        $binary = 1; # the file is binary
    }
    else #the file is a not binary SAM file
    {
        $inputOption = " -S ";#-S mean input is SAM
        $binary = 0; # the file is not binary
    }
    my $checkFormatCommand="$samtools view $inputOption $samFile -H > /dev/null";
    # The samtools view will ask only for the header to be outputted (-H option), and the STDOUT is redirected to nowher (>/dev/null);
    my $formatValidation=toolbox::run($checkFormatCommand,"noprint");

    if ($formatValidation == 1)                    # if no error occured in extracting header, ok
    {
        ##DEBUG toolbox::exportLog("INFOS: toolbox::checkSamOrBamFormat : The file $samFile is a SAM/BAM file\n",1);
        return 1 if $binary == 0;# Return 1 if the file is a SAM
        return 2 if $binary == 1;# Return 2 if the file is a BAM
    }
    else                                # if one or some error(s) occured in extracting header, not ok
    {
        toolbox::exportLog("ERROR: checkFormat::checkFormatSamOrBam : The file $samFile is not a SAM/BAM file\n",0);
        return 0;
    }
}
################################################################################################
# END checkFormatSamOrBam
################################################################################################


################################################################################################
# sub checkFormatVcf=> checks if the fileformat is VCF
################################################################################################
# arguments : file name to analyze
# Returns 1 if the file format is validated else 0 (error are managed by toolbox::exportLog)
################################################################################################
sub checkFormatVcf
{
    #Inspired from The vcftools vcf-validator
    my ($file)=@_;

    #Check if we can read the file
    my $rightToRead = toolbox::readFile($file);
    if ($rightToRead == 0)
    {
        toolbox::exportLog("ERROR: toolbox::checkFormatVcf : The file $file is not readable\n",0);
        return 0;
    }

    #Parsing the file
    my @header;#List to gather the header
    my @listOfFields;
    open(my $inputHandle, "<",$file) or toolbox::exportLog("ERROR: checkFormat::checkFormatVcf : Cannot open the file $file\n$!\n",0);

    # if the input file is a gz file
    if($file =~ m/\.gz$/)
    {
        $inputHandle = new IO::Uncompress::Gunzip $inputHandle or toolbox::exportLog("ERROR: checkFormat::checkFormatVcf : Cannot open the gz file $file: $GunzipError\n",0);
    }
    while (my $line=<$inputHandle>)
    {
        chomp $line;

        ##DEBUG print $line."\n";
        if ($line =~ m/^#/)
        {
           push (@header, $line);
        }
        else
        {
            @listOfFields = split /\t/, $line;
        }
    }

    #Check if the first line of the header is including the version
    my $versionLine=defined $header[0]? $header[0]:undef;
    toolbox::exportLog("ERROR: checkFormat::checkFormatVcf : There is no header of the file $file\n",0) unless (defined $versionLine); #Thrown an error if the $version cannot be obtained (misformatted line)

    my @version=split /VCFv/,$versionLine;
    exportLog("ERROR: checkFormat::checkfFormatVcf : Cannot evaluate the VCF version of the file $file file\n",0) if (scalar(@version)==0); #Thrown an error if the $version cannot be obtained (misformatted line)
    ##DEBUG print "DEBUG: $0: vcf version $versionLine : ".scalar(@version)." : $version[1] \n";
    eval ($version[1] == $version[1]) or exportLog("ERROR: checkFormat::checkFormatVcf : Cannot obtain the VCF version of $file\n",0); #Verifying if the value obtained is numerical.

    # Check the first line format as recommanded
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Individual1
    #Chr1    228069  .       A       .       31.23   LowQual AN=2;DP=1;MQ=37.00;MQ0=0        GT:DP   0/0:1
    if (scalar @listOfFields < 10) #Less than the 10 minimum fields
    {
        toolbox::exportLog("ERROR: checkFormat::checkFormatVcf : The file $file is misformatted (less than 10 colums) ".Dumper(\@listOfFields)."\n",0);
    }

    #Check if the second field (the position) is numerical
    eval ($listOfFields[1] == $listOfFields [1]) or toolbox::exportLog("ERROR: checkFormat::checkFormatVcf : Cannot confirm that $file is a VCF.\nAborting.\n",0); #Verifying if numerical. Die if not

    close $inputHandle;

    return 1; #Return correct if all check are Ok
}
################################################################################################
# END sub checkFormatVcf
################################################################################################

################################################################################################
# sub checkFormatFasta => check if a given file is a fasta or not
################################################################################################
# arguments :
# 	- the fasta file
# returns :
#	- a boolean 1 for Ok, 0 for not.
#	- send a warning to the log if not a fasta file
################################################################################################

sub checkFormatFasta{
    my ($file)=@_;

    #Check if we can read the file
    my $rightToRead = toolbox::readFile($file);
    if ($rightToRead == 0)
    {
        toolbox::exportLog("ERROR: checkFormat::checkFormatFasta : The file $file is not readable\n",0);
        return 0;
    }

    #Opening file
    open(FILE, "<", $file) or toolbox::exportLog("ERROR: checkFormat::checkFormatFasta : Cannot open the file $file\n$!\n",0);

    #Reading file
    my $initiator = 0; #We need to scan the very first line, that begins with '>', but this is the only one beginning with that we need
    my %errors; # Will score all the error types
    my $lineNumber=0;#Recording line position
    while (my $line = <FILE>)
    {
         chomp $line;
         $lineNumber++;#Counter for scoring lines with errors
         #Checking first line, that must be '>NAME'
         if ($initiator == 0)
         {
             #very first line
             $initiator = 1; #no need for any other lines beginning with ##
             next if $line =~ m/^>/; #Correct first line, let's go to the following ones
             #First line is not correct
             $errors{$lineNumber} = "Not correctly formatted, must be a header name such as '>NAME'.";
             #Will stop immediatly as the file is not correct at all
             last;
         }
 #Other lines
         next if $line =~ m/^$/;# Representing empty lines
         if ($line =~ m/^>/) #A new sequence beginning
         {
             next;
         }
         else # The sequence is read
         {
             my $modifiedLine = $line;
             #Normal IUPAC Code
             $modifiedLine =~ s/[A|T|G|C|a|g|t|c|N|n]//g;
             #Degenerated IUPAC
             $modifiedLine =~ s/[R|r|Y|y|S|s|W|w|K|k|M|m|B|b|D|d|H|h|V|v]//g;
             #DEBUG : print "\n**",$modifiedLine,"**\n";
             #DEBUG : exit;
             next if length($modifiedLine eq "");
             $errors{$lineNumber} = "Not IUPAC letter, only ATGCNatgcn and RYSWKMBDHVryswkmbdhv characters are allowed: unauthorized characters are $modifiedLine.";
         }
 #Check if there are too much errors...
         if (scalar (keys %errors) > 19)
         {
           #Will stop after 20 errors
           my $tooMuchErrors = "WARNING : checkFormat::checkFormatFasta : Too much errors in the Fasta file $file, only the first 20 errors are shown..."; #Will be printed before all errors
           toolbox::exportLog($tooMuchErrors,2);
           last;
         }
         next;

     }
    close FILE;

    #Checking if there are errors
    my $numberOfErrors = scalar(keys %errors);

    #Sending the warning message, if any
    if ($numberOfErrors)#There are errors, we will send a warning to the system
	{
        #Formatting errors for the warning
        my $warningLog;
        $warningLog.="WARNING : checkFormat::checkFormatFasta : There are $numberOfErrors errors detected for file $file\n";
        foreach my $individualError (sort {$a <=> $b} keys %errors)#Sorting error per line number
        {
            $warningLog.="WARNING : checkFormat::checkFormatFasta : Line $individualError does not respect FASTA standard : $errors{$individualError}\n";
        }
        #Sending the warning
        toolbox::exportLog($warningLog,2);
        #returning for failure
        return 0;
    }

    #returning ok
    toolbox::exportLog("INFOS: checkFormat::checkFormatFasta : The file $file is a FASTA file\n",1);
    return 1;
}


################################################################################################
# END sub fastaFormatValidator
################################################################################################

1;

=head1 NAME

package I<checkFormat>

=head1 SYNOPSIS

	use checkFormat;

	checkFormat::::checkFormatFastq($fastqFile);


=head1 DESCRIPTION

This module aims to provide a set of functions related to check file format


=head2 Functions


=head3 checkFormat::checkFormatFastq()

This function checks if a file is really a FASTQ file.
Returns  boolean (1 if the format is fastq else 0)

Example :
C<checkFormat::checkFormatFastq("/data/projects/oryza/RC1.fastq"); >

=head3 checkFormat::checkFormatSamOrBam()

This function verifies the SAM/BAM format based on samtools view system. One argument is required : filename to analyze.
Returns  boolean (1 if the fileformat is sam, 2 bam and 0 neither bam or sam)

Example :
C<checkFormat::checkFormatSamOrBam($samFile); >

=head3 checkFormat::checkFormatVcf()

This function checks if the format of the given as argument is VCF.
Returns 1  if the file format is validated else 0 (an error is maneged by toolbox::exportLog)

Example :
C<checkFormat::checkFormatVcf($file); >

=head3 checkFormat::checkFormatFasta()

This function checks if a given file is a true DNA fasta file
The only required argument is the filemane. The file can be plain text
Returns a 1 for success, and a 0 for failure. Will send a warning to the log in case of failure.
Will return a maximum of 20 errors.
Will stop immediatly if the first line is misformatted

Example :
C<( checkFormat::checkFormatFasta($fastaFile);); >


=head1 AUTHORS

Cecile Monat, Ayite Kougbeadjo, Julie Orjuela-Bouniol, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

L<http://www.southgreen.fr/>


=cut
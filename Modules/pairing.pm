package pairing;

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

#For gz files
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#########################################
#pairRecognition
#from a large set of files in a folder, will recognize forward, reverse and associate them
#A single file will be alone in its subhash, but can be named forward or reverse (generally forward)
    #$VAR1 = \{
    #            '@HWUSI-EAS454_0001:1:1:15:303#0' => {
    #                                                   'ReadGroup' => 'single',
    #                                                   'forward' => '../DATA-TEST/Files_for_pairing_test/single.fastq'
    #                                                 },
    #            '@CJP75M1:362:C20PVACXX:7:1101:1496:2086' => {
    #                                                           'ReadGroup' => 'first_forward',
    #                                                           'forward' => '../DATA-TEST/Files_for_pairing_test/first_forward.fastq',
    #                                                           'reverse' => '../DATA-TEST/Files_for_pairing_test/first_reverse.fastq'
    #                                                         },
    #            '@HWUSI-EAS454_0001:1:1:15:911#0' => {
    #                                                   'ReadGroup' => 'second_forward_single',
    #                                                   'forward' => '../DATA-TEST/Files_for_pairing_test/second_forward_single.fastq'
    #                                                 },
    #            '@HWUSI-EAS454_0001:1:1:15:301#0' => {
    #                                                   'ReadGroup' => 'second_forward_forwardRepaired',
    #                                                   'forward' => '../DATA-TEST/Files_for_pairing_test/second_forward_forwardRepaired.fastq',
    #                                                   'reverse' => '../DATA-TEST/Files_for_pairing_test/second_reverse_reverseRepaired.fastq'
    #                                                 }
    #          };
#########################################
sub pairRecognition
{
    
    my ($folder)=@_;		# recovery of informations
    my %pairs;
    use Data::Dumper;
    
    #Reading the files in the folder
    my $listFiles_ref=toolbox::readDir($folder);
    my @listFiles=@{$listFiles_ref};
    
    foreach my $currentFile (@listFiles)		# for each file
    {  
	my $checkFastq=toolbox::checkFormatFastq($currentFile);		# check the file is a fastq one
	if ($checkFastq == 0)		# If the file is not a Fastq File, we cannot consider it
	{
	    next; # go to the next sequence
	}
	
	# If the file is a fastq one:
	my $firstLineComplete=`head -n 1 $currentFile`;		#Fetching the first line to obtain the ID sequence
	#if the file is a gz Compressed file
	$firstLineComplete = `zcat $currentFile | head -n 1` if ($currentFile =~ m/gz$/);
	chomp $firstLineComplete;
	my $namingConvention;		#Infos for the type of modification
	
	#We need 3 infos: the coding convention ($namingConvention), the name of the pair ($sequenceName) and the strand of the current mate of the pair ($typeOfStrand)
	
	#Picking up the name of the pair and its coding convention 
	my $sequenceName=$firstLineComplete; 
	$sequenceName =~ s/\/.$// and $namingConvention = 1; # old agreement /1 /2
	$sequenceName =~ s/\s\d:\w:\d:\w*$// and $namingConvention = 2; # new agreement 1:N:[A-Z]1 to 10
	
	#Picking up the strand of the current mate of the pair 
	my $typeOfStrand=$firstLineComplete;
	$typeOfStrand =~ s/.+\/(.)$/$1/ if $namingConvention == 1;
	$typeOfStrand =~ s/.+\s(\d:\w:\d:\w*)$/$1/ if $namingConvention == 2;
	
	#Converting the end of Line in forward and reverse
	my $nameOfStrand = "unknown"; # May correspond to the single sequence
	$nameOfStrand = "forward" if $typeOfStrand =~ m/^1/;
	$nameOfStrand = "reverse" if $typeOfStrand =~ m/^2/;
	
	#Completing the hash
	$pairs{$sequenceName}{$nameOfStrand}=$currentFile;
	
	
	if  (not defined($pairs{$sequenceName}{"ReadGroup"})) #If the current file is the first one of the pair, we do not have yet infos on the readGroup to be associated with
	{
	    my ($tag, $readGroup) = extractName($currentFile); #picking up the info of readGroup based on the name of the file
	    $pairs{$sequenceName}{"ReadGroup"}=$readGroup; #Adding the readGroup to the name
	}
    }
    
    #Exporting log   
    return (\%pairs);
    
}
#########################################
#createDirPerCouple
#From a hash (reference to) of paired sequences, will constrct and organize couples in separated folders, named with their RG
#########################################
sub createDirPerCouple
{
    my ($hashOfPairs,$targetDirectory)=@_;		# recovery of informations
    my @listOfSequences = keys %{$hashOfPairs}; #Pick up the names of first sequences (key values of the reference)
    foreach my $couple (@listOfSequences)		# for each pair
    {
	if ($couple=~ /^@/)		# Allows the directory to be created only for a FASTQ file pair		 
	{
	    #Extract infos
	    my $forwardFile=$hashOfPairs->{$couple}{"forward"};		#Picking up the forwardFile name
	    $forwardFile=$hashOfPairs->{$couple}{"unknown"} unless $forwardFile; #forwardFile name may not exist and the file is single
	    
	    my $reverseFile;
	    $reverseFile=$hashOfPairs->{$couple}{"reverse"} if exists $hashOfPairs->{$couple}{"reverse"}; #Reverse may not exist if the sequence is single
	    
	    my $ReadGroup=$hashOfPairs->{$couple}{"ReadGroup"};		# We need the ReadGroup to create the directory
	    $ReadGroup=$targetDirectory."/".$ReadGroup;		# We add the full path to the ReadGroup that will create the full directory name
	  
	    #Creating subfolder based on ReadGroup name
	    my $makeDirCheck=toolbox::makeDir($ReadGroup,0); #Will not erase if created
	    if ($makeDirCheck != 1)		#If not able not create the directory...
	    {
		toolbox::exportLog("ERROR: pairing::createDirPerCouple : Cannot create the folder $ReadGroup\n",0);		# ... return an error message
		return 0;
	    }
	    #Everything is Ok, we copy the files in the correct place and remove the files from their original place
	    my $forwardCopyCommand="mv $forwardFile ".$ReadGroup."/."; # removing will be effective ONLY if copy is ok
	    toolbox::run($forwardCopyCommand,"noprint");
	    my $reverseCopyCommand="mv $reverseFile ".$ReadGroup."/." if $reverseFile;#reverse file may be absent
	    toolbox::run($reverseCopyCommand,"noprint") if $reverseFile;
	}
    }
    return 1;
}
##################################################################################
#repairing
#From two de-paired files, forward + reverse, will create three files, forward, reverse and single
##################################################################################
sub repairing
{
    toolbox::exportLog("ERROR: pairing::repairing : should get at least two arguments\n",0) if (@_ < 2);

    my($forwardFile,$reverseFile,$directory)=@_;		# recovery of informations
    
    toolbox::existsDir($directory);		# check if the directory exits  
    my $dirOut=defined($directory)?$directory.'/':'./';		#### ????????????????

    #Check Fastq format
    toolbox::checkFormatFastq($forwardFile);		# check if forward file is a fastq one
    toolbox::checkFormatFastq($reverseFile);		# check if reverse file is a fastq one
    
    #Extraction of the name and creation of output names
    my ($forwardTag,$readGroup) = extractName($forwardFile);
    my $forwardFileOut=$dirOut.$forwardTag.".REPAIRING.fastq";
    my ($reverseTag,$ReadGroup_tmp)= extractName($reverseFile);
    my $reverseFileOut=$dirOut.$reverseTag.".REPAIRING.fastq";
   
    my @splitDir = split ("$readGroup", $dirOut);		# recovery of path to create the specific single directory
    my $singleDir = $splitDir[0].$readGroup."_Single/";		# creation of the specific single directory
    toolbox::makeDir($singleDir);				# creation of single directory
    my $singleFileOut=$singleDir.$readGroup."Single.fastq";

    #Opening infiles
    open(my $forwardIn, "<", $forwardFile) or toolbox::exportLog("ERROR: pairing::repairing : Can't open the file $forwardFile $!\n",0);
    open(my $reverseIn, "<", $reverseFile) or toolbox::exportLog("ERROR: pairing::repairing : Can't open the file $reverseFile $!\n",0);
    
    #Opening outfiles
    open(my $mateForward, ">",$forwardFileOut) or toolbox::exportLog("ERROR: pairing::repairing : Can't open the file $forwardFileOut $!\n",0);
    open(my $mateReverse, ">",$reverseFileOut) or toolbox::exportLog("ERROR: pairing::repairing : Can't open the file $reverseFileOut $!\n",0);
    open(my $single,">",$singleFileOut) or toolbox::exportLog("ERROR: pairing::repairing : Can't open the file $singleFileOut $!\n",0);
      
        
    #If input files are gzipped
    if($forwardFile =~ m/\.gz$/)
    {
	#Transforming the output flux in a compressed gz format
	$mateForward = new IO::Compress::Gzip $mateForward or toolbox::exportLog("ERROR: pairing::repairing : Can't create the gz file $mateForward : $GzipError $!\n",0);
	$mateForward->autoflush(1);
	$mateReverse = new IO::Compress::Gzip $mateReverse or toolbox::exportLog("ERROR: pairing::repairing : Can't create the gz file $mateReverse : $GzipError $!\n",0);
	$mateReverse->autoflush(1);
	$single = new IO::Compress::Gzip $single or toolbox::exportLog("ERROR: pairing::repairing : Can't create the gz file $single : $GzipError $!\n",0);
	$single->autoflush(1);
	#Decompressing in input gz flux
	$forwardIn = new IO::Uncompress::Gunzip $forwardIn or toolbox::exportLog("ERROR: pairing::repairing : Can't open the gz file $forwardIn : $GzipError $!\n",0);
	$reverseIn = new IO::Uncompress::Gunzip $reverseIn or toolbox::exportLog("ERROR: pairing::repairing : Can't open the gz file $reverseIn : $GzipError $!\n",0);
    }
    
    #Creating counters
    my $pairedSequences=0;
    my $singleSequences=0;
    
    #Variables
    my %forwardSequences;
    
    #Reading forward input file
    while (<$forwardIn>)		# for each line of the file
    {
	my $line = $_;		# recovery of sequence's name line
	chomp $line;
        next if ($line =~ m/^$/);		# next line if the line is empty
        my $next = $line."\n";
	$line =~ s/\/\d$//;		# Removing the ending /1 or /2 if the file is encoded in Illumina 1.3-1.8
        $line =~ s/\s\d:\w:\d:\w*$//;		#Removing the ending 1:N:1234 or analog if the file is encoded in Illumina 1.9+
	$next .= <$forwardIn>;		# add the sequence line to the "next" variable
	$next .= <$forwardIn>;		# add the informations line to the "next" variable
	$next .= <$forwardIn>;		# add the quality line to the "next" variable
	$forwardSequences{$line}=$next;		# hash listing forward sequences{name of the current sequence}= lines of sequences, infos and quality
    }

    #Comparing with the reverse seq ID
    while (<$reverseIn>)		# for each line of the file
    {
	my $line = $_;		# recovery of sequence's name line
	chomp $line;
	next if ($line =~ m/^$/);		# next line if the line is empty
	my $next = $line."\n";
	$line =~ s/\/\d$//;		# Removing the ending /1 or /2 if the file is encoded in Illumina 1.3-1.8
	$line =~ s/\s\d:\w:\d:\w*$//;		#Removing the ending 1:N:1234 or analog if the file is encoded in Illumina 1.9+
	$next .= <$reverseIn>;		# add the sequence line to the "next" variable
	$next .= <$reverseIn>;		# add the informations line to the "next" variable
	$next .= <$reverseIn>;		# add the quality line to the "next" variable
	if (exists $forwardSequences{$line})		# if it exists a key of the current sequence in the hash of forward sequence
	{
	    my $out = $forwardSequences{$line};
	    print $mateForward $out;		# print in the forward file the coorespondante sequence has the current sequence
	    my $out2 = $next;
	    print $mateReverse $out2;		# print in the reverse file the current sequence
	    delete $forwardSequences{$line}; # To save memory and to conserve only the singles
	    $pairedSequences++; #Increment the number of pairs
	}
	else		# if the current sequence is not listed in the hash of forward sequence, that means it is a single one
	{
	    my $out2 = $next;
	    print $single $out2;		# print in the single file the current sequence
	    $singleSequences++;#Increment the number of single sequences
	}
    }
    
    foreach my $remainingNames (keys %forwardSequences)		# after the comparison of all the reverse sequences, if it remains some forward sequence, that mean it is a single one
    {
        my $out = $forwardSequences{$remainingNames};
        print $single $out;		# print in the single file the current sequence
        $singleSequences++;#Increment the number of single sequences
    }

    if ($singleSequences == 0)		# if there is none single sequence
    {
	my $removeSingle = "rm -Rf $singleFileOut $singleDir";		# remove the single file and the single directory
	toolbox::run($removeSingle);
    }
    
    #Closing files
    close $mateForward;
    close $mateReverse;
    close $single;
    close $forwardIn;
    close $reverseIn;
    
    #Creation of the output
    my $outlog="INFOS: pairing::repairing : ".$pairedSequences." pairs were recovered (".($pairedSequences*2)." sequences), and ".$singleSequences." singles were extracted\n";
    
    #Exporting the output
    toolbox::exportLog($outlog,1);    
}
#########################################
#extractName
#remove fastq from the name, and clean the name
#########################################
sub extractName
{
    my ($file)=@_;		# recovery of informations
    chomp $file; 

    my $name=`basename $file` or toolbox::exportLog("ERROR : $0 : pairing::extractName cannot work on $file name\n",0);		# recovery of the last name present in the path, normally the name of the file

    chomp $name;
    
    my @listName=split /\./, $name; #Separate in a list the filename and its format (ex file1.fastq is separated in file1 and fastq)

    my $shortName = shift @listName ;#pick up what is before the first "."

    my $readGroup=$shortName; #Picking up the readGroup name as the returning of the previous line
    $readGroup =~ s/_.*$//; # Removing what is after any "_" 
       
    #cleaning name
    $shortName=~ s/ /_/g;#removing spaces
    $shortName=~ s/[é|è|ê]/e/g;#removing non utf8 accents from e
    $shortName=~ s/à/a/g;#removing non utf8 accents from a

    return ($shortName,$readGroup);
}
1;

=head1 NAME

    Package I<pairing> 

=head1 SYNOPSIS

        use pairing;
    
        pairing::pairRecognition ($folder);
    
        pairing::createDirPerCouple ($hashOfPairs,$targetDirectory);
    
        pairing::repairing ($forwardFile,$reverseFile,$directory);
	
	pairing::extractName ($file);

=head1 DESCRIPTION

Package pairing will ensure the repairing of Fastq files and the recognition of the pairs

=head2 FUNCTIONS


=head3 pairing::pairRecognition

This module recognize pair files. From a folder which contain the all pairs (and single) file from your analyse, it is able to know which file is mating which other file, even if the names are differents
It takes only one argument, the folder of your analyse.



=head3 pairing::createDirPerCouple

This module create a directory for each pair (and single) listed in a hash you give to it
It takes two arguments: the hash of pair (and single) informations, the name of the directory you want for it



=head3 pairing::repairing

This module recognize pair sequences. After the cleaning stage, some pair files does not contain the same number of sequences anymore, this module extract the new single sequences from the pair files to
get real pair files again.
It takes three arguments: the forward file, the reverse file, the name of the output directory
The FASTQ files can be in plain format or in gz compressed fastq format



=head3 pairing::extractName

This module remove the path and the extention of the file to return only the name of it
It takes only one argument, the file you want the name of



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut

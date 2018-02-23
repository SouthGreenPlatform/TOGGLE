#!/usr/bin/env perl

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
use Data::Dumper;

system("clear");

print"###################################################################################################################################
#
# This script will allow a user-space automatic installation of TOGGLe from the current MASTER version https://github.com/SouthGreenPlatform/TOGGLE
#
###################################################################################################################################\n";

#64bits validation

print "\n\n############################################
##\t Checking if your system is a 64 bits
############################################\n";

my $machineType=`uname -m`;
if ($machineType =~ m/x86_64/)
{
    print "Your system is an intel-like 64 bits:\t$machineType";
}
else
{
	warn "** It seems that your system is not an intel-like 64bits.**
	**Some softwares may not work**
	** $machineType **";
}
    
#Software and lib validation

print "\n\n############################################
##\t Checking installed softwares and libraries
############################################\n";

my $requirements = {"mandatory"=>{  "git" => 0,
                                    "tar" => 0,
                                    "wget"=> 0,
                                    "perl"=>0}};
#print Dumper($requirements);

#Software validation MANDATORY

foreach my $softs (keys %{$requirements->{"mandatory"}})
{
    my $controlCommand = "which $softs 2>/dev/null";
    my $controlRes = `$controlCommand`;
    chomp $controlRes;
    #print "\n$softs --> $controlRes\n";
    die ("The $softs software is not installed and is mandatory for continuing the installation of TOGGLe.\n\n Please contact your administrator for installing it.\n") unless $controlRes;
    $requirements->{"software"}->{$softs} = $controlRes;
}

#Test perl version 5 higher than 5.18 and 5.22
my $perlVersionCom=`perl -v| grep version`;
chomp $perlVersionCom;

#Die if not a perl5
die ("\nPerl 5 is required for TOGGLe to work (5.16 minimum, 5.22 recommended), we cannot continue the installation.\nABORTING...\n") unless $perlVersionCom =~ m/This is perl 5/;

#The output is normally as such:  This is perl 5, version 22, subversion 1 (v5.22.1) built for x86_64-linux-gnu-thread-multi
my @fields = split /version /,$perlVersionCom;
my ($perlVersion) = split /,/, $fields[1];
if ($perlVersion < 16)
{
    die ("Your standard Perl5 version is too old (version 5.$perlVersion) and will block TOGGLe execution. Please contact your administrator to have access to a newer Perl5 version (5.16 minimum, 5.22 recommended)\n");
}
if ($perlVersion < 22)
{
    warn ("Your standard Perl5 version is old (version 5.$perlVersion) and may provoke warnings. You can contact your administrator to have access to a newer Perl5 version (5.22 recommended)\n")
}
    
print "\nAll required softwares (perl, git, wget and tar) are installed.\n";

print "#######################################################
## \tTOGGLe installation
#######################################################";

my $INSTALLPATH;
if (defined $ARGV[0]) 
{
	$INSTALLPATH=$ARGV[0];
}
else
{
	print "\nPlease provide the installation path as ABSOLUTE (e.g. /home/myUserName/TOGGLE):\n";
	$INSTALLPATH = <STDIN>;
	chomp $INSTALLPATH;

	while (1)
	{  
	    if ($INSTALLPATH !~ m/^\//) #Absolute Path ?
	    {
		print "\nThis path is not an absolute one (not starting from '/').\nPlease provide an absolute Path:\n";
		$INSTALLPATH = <STDIN>;
		chomp $INSTALLPATH;
		next;
	    }

	    #Asking for the path to be Ok
	    print "\nIs the correct Install Path is: $INSTALLPATH ? [Y|N]: \n";
	    my $yn;
	    $yn = <STDIN>;
	    chomp $yn;

	    if ($yn =~ m/^Y|^y/)
	    {
		print"\nTOGGLE will be installed in '$INSTALLPATH'\n ";
		last;
	    }
	    elsif ($yn =~ m/^N|^n/ )
	    {
		print "\nPlease, provide the correct installation path:";
		$INSTALLPATH = <STDIN>;
		chomp $INSTALLPATH;
		next;
	    }
	    else
	    {
		print "\nPlease answer 'Y' or 'N'\n";
		next;
	    }
	}

}
mkdir $INSTALLPATH or die ("\nCannot create installation directory $INSTALLPATH:\n$!\n\nAborting installation...\n");

#Cloning current version of TOGGLe

print "\nCloning the current Git MASTER Version";

my $gitCommand = "git clone https://github.com/SouthGreenPlatform/TOGGLE.git $INSTALLPATH";
system("$gitCommand") and die ("\nCannot clone the current version of TOGGLe: $!\n");


chdir $INSTALLPATH or die ("\nCannot go to $INSTALLPATH : $!\nTOGGLe is cloned but not configured\n");

#Parsing automatically the localConfig.pm file for softwares

my $localConfigFile = "modules/localConfig.pm";

open (my $fhConfig, "<", $localConfigFile) or die ("\nCannot read the original localConfig.pm file:\n$!\nTOGGLE is cloned but not configured\n");

while (my $line = <$fhConfig>)
{
    chomp $line;
    next unless $line =~ m/^our \$/;
    my @fieldsType = split /\s|=/, $line;
    my $softwareType = $fieldsType[1];
    $softwareType =~ s/\$//;
    my @fieldsName = split /\//, $line;
    my $softwareName = $fieldsName[-1];
    $softwareName =~ s/";$//;
    $requirements->{"software"}->{$softwareType}=$softwareName;    
}
close $fhConfig;

foreach my $soft (keys %{$requirements->{"software"}})
{
    my $controlCommand;
    my $name =$requirements->{"software"}->{$soft};
    if ($name =~ m/\.jar$/)
    {
        #This soft must be launched using the java -jar command, we will use the find command
        $controlCommand = "find / 2> /dev/null | grep -m 1 $name ";
    }
    elsif ($soft eq "java")
    {
        #This is java
        $controlCommand = "which java 2> /dev/null";
        my $controlJava = `$controlCommand`;
        chomp $controlCommand;
        if ($controlCommand)
        {
            my $javaTypeCom = "$controlJava -version";
            warn("\nYour current standard JAVA is openJDK. Unfortunately, picard-tools require the Oracle JAVA to run, thus you will not be able to use those tools in TOGGLe in the current configuration.\n") if $javaTypeCom =~ m/openjdk/;
        }
    }
    else
    {
        $controlCommand = "which $name 2> /dev/null";
    }
    
    #Launching the command
    my $controlRes = `$controlCommand`;
    chomp $controlRes;
    
    if ($controlRes)
    {
        #The lib is present
        $requirements->{"software"}->{$soft} = $controlRes;
    }
    else
    {
        #The lib is absent
        warn("\nThe $soft software is absent or not in the current PATH. TOGGLe cannot run pipeline using this tool in your current configuration. Other tools will function anyway\n") unless $soft =~ m/toggle/i;
    }
}

#print Dumper($requirements);

print ("\nCONFIGURING YOUR PERSONAL localConfig.pm");

my $configOk = "modules/localConfigTEMP.pm";

open (my $fhConfig2, "<", $localConfigFile) or die ("\nCannot read the original localConfig.pm file:\n$!\nTOGGLe is cloned but not configured\n");
open (my $fhOut, ">", $configOk) or die ("\nCannot create the local configuration file:\n$!\nTOGGLe is cloned but not configured\n");

while (my $line = <$fhConfig2>)
{
    chomp $line;
    if ($line =~ m/^our \$/)
    {
        #the configuration is here
        my @fieldsType = split /\s/, $line;
        my $softwareType = $fieldsType[1];
        $softwareType =~ s/\$//;
        my $localPath=$requirements->{"software"}->{$softwareType};
        
        #If the software is unknwon.
        if ($localPath !~ m/\// && $softwareType !~ m/toggle/i)
        {
        #do not change the line, do nothing
        #print "Cannot provide path for $softwareType, leave undefined\n\n";
        }
        elsif ($softwareType =~ m/toggle/i)
        {
            $line = "our \$toggle = \"$INSTALLPATH\";\n";
        }
        elsif ($softwareType eq "java")
        {
            $line = "our \$java = \"".$localPath." -Xmx12g -jar\";";
        }
        else
        {
            my $newLine = "our \$".$softwareType." = \"";
            #if java requested
            $newLine.="\$java " if $line =~ /\.jar$/;
            $newLine.=$localPath."\";";
            $line = $newLine;
        }
        
    }
   print $fhOut $line;
   print $fhOut "\n";
}
close $fhConfig2;
close $fhOut;

my $cpCommand ="mv $configOk $localConfigFile";
system($cpCommand) and die ("\nCannot copy correct config file as localConfig.pm:\n$!\nTOGGLe is cloned but not configured\n");

print "\nDownloading the various Perl modules";

my $wgetCommand = "wget http://bioinfo-web.mpl.ird.fr/toggle/perlModules.tar.gz";
system("$wgetCommand") and die ("\nCannot download the required Perl libraries: $!\n\nPlease use the CPAN to install them...\n");

# DECLARE VARIABLES WITH PATHS 
my $MODULES=$INSTALLPATH."/modules";

print "\tDecompressing Perl Modules\n";

system ("cd $INSTALLPATH && tar xzf perlModules.tar.gz && rm -Rf perlModules.tar.gz && cp -R perlModules/* $MODULES/. && rm -Rf perlModules && cd") and die ("\nCannot decompress the perl Modules:\n$!\n");

#Adding toggle in the user PERL5LIB path

system ("echo \"\nPERL5LIB=\$PERL5LIB:$MODULES\n\" | cat - >> ~/.bashrc && echo \"\nPATH=\$PATH:$INSTALLPATH\n\" | cat - >> ~/.bashrc && source ~/.bashrc") and die("\nCannot add paths: $!\n");
print "echo \"\nPERL5LIB=\$PERL5LIB:$MODULES\n\" | cat - >> ~/.bashrc && echo \"\nPATH=\$PATH:$INSTALLPATH\n\" | cat - >> ~/.bashrc && source ~/.bashrc";
print "\n The automatic configuration is finished.\n\nPlease use first the test data as recommended on the GitHub https://github.com/SouthGreenPlatform/TOGGLE.\n\nThanks for using TOGGLe\n\n
###########################################################################################################################
#\tCITATION:
#\tTOGGLe, a flexible framework for easily building complex workflows and performing robust large-scale NGS analyses.
#\tChristine Tranchant-Dubreuil, Sebastien Ravel, Cecile Monat, Gautier Sarah, Abdoulaye Diallo, Laura Helou, Alexis Dereeper,
#\tNdomassi Tando, Julie Orjuela-Bouniol, Francois Sabot.
#\tbioRxiv, doi: https://doi.org/10.1101/245480
#\thttps://toggle.southgreen.fr/
###########################################################################################################################";

exit;

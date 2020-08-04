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
use warnings;
use Data::Dumper;
use localConfig;
use Switch;

#This tool will allows developer to add easily a new software in TOGGLe

print "\nWelcome to the tenebrous world of TOGGLe, the best workflow manager you have ever dreamed of...\n";

#Asking for the module name, ie the generic tool
my $module="";
while ($module eq "")
{
    print "\nWhat is the name of your tool (e.g. bwa) ?\n";
    $module = <STDIN>;
    chomp $module;
}

#Asking for the function name, ie the specific tool itself
my $function = "";
while ($function eq "")
{
    print "\nWhat is the name of your function (e.g. bwaAln for bwa aln) ?\n";
    $function = <STDIN>;
    chomp $function;
    $function =~ s/\s//g;
}

# Create the module if needed
$module = lc $module;
my $moduleFile = $module.".pm";

unless (-f "$toggle/modules/$moduleFile")
{   #We create the module if it does not exists
    open (my $fhTemplate, "<", "$toggle/modules/module_template.pm") or die ("\nCannot open the $toggle/modules/module_template.pm module_template file:\n$!\n");
    open (my $fhModule, ">", "$toggle/modules/$moduleFile") or die  ("\nCannot create the $moduleFile file:\n$!\n");

    while (my $line = <$fhTemplate>)
    {
        if ($line =~ m/^package/)
        {
            $line =~ s/module_template/$module/;
        }
        print $fhModule $line;
    }
    close $fhTemplate;
    close $fhModule;
}

#Testing if exists

my $subName = "sub $function";
my $grepCom = "grep -c \"$subName\" $toggle/modules/$moduleFile 2>/dev/null";
my $grep = `$grepCom`;
chomp $grep;
if ($grep)
{   #the function exists...
    die "\nThe function already exists, will quit...\n";
}

#asking for IN, OUT mandatory and version

my $in = "";
my $out = "";
my $mandatory ="";
my $version = "";
my $testParams = "";

while ($in eq "")
{
    print "\nWhat are the different entry formats of your tool, separated by commas (e.g. fastq or fastq,fasta) ?\n";
    $in = <STDIN>;
    chomp $in;
    $in =~ s/\s//g;
}

while ($out eq "")
{
    print "\nWhat are the different output formats of your tool, separated by commas (e.g. fastq or fastq,fasta) ? \n**NOTE if your tool is a dead-end one, such as FASTQC, please provide NA as output format.**\n";
    $out = <STDIN>;
    chomp $out;
    $out =~ s/\s//g;
}

print "\nAre there any mandatory requirement as option for your tool (e.g. reference or gff) ? \n**NOTE if none leave empty.**\n";
$mandatory = <STDIN>;
chomp $mandatory;

while ($version eq "")
{
    print "\nHow do you obtain the version of your tool (e.g. 'bwa 2>&1| grep version' or 'java --version | grep Version') \n";
    $version = <STDIN>;
    chomp $version;
    #Testing if the command version is ok
    `$version`or $version ="";
}

#The testParams value will be used in the fileConfigurator.pm module for block test.
print "\n'Are there any test parameters to include for testing (e.g. '-n 5' for bwa aln tests ) ?\n";
$testParams = <STDIN>;
chomp $testParams;

#Standard command line
my $commandLine="";
while ($commandLine eq "")
{
    print "\nWhat is the standard command line to launch your tool ?\n\t- A file in must be written FILEIN\n\t- a file out must be written FILEOUT\n\t- Options location must be written as [options]\n\t- Reference must be written as REFERENCE\n\t- Gff/Gtf must be written as GFF\n\t- Keyfile must be written as KEYFILE\n\t- Vcf must be written as VCF\n";
    print "\nAn example for bwa aln will be:\n\tbwa aln [options] FILEOUT REFERENCE FILEIN\n\n";
    print "\n**NOTE You may have to adapt and correct this command at this end...\n";
    $commandLine=<STDIN>;
    chomp $commandLine;
}



#Create the function

#Creating the text
my $subText=$subName."\n{\n";

$subText .= "#The standard way to write variables are:\n#REFERENCE = \$reference";
$subText .= '#PLEASE CHECK IF IT IS OK AT THIS POINT!!'."\n\t".'my ($fileIn,$fileOut,';
$subText .= '$reference,' if $commandLine =~ m/REFERENCE/;
$subText .= '$gff,' if $commandLine =~ m/GFF|GTF/;
$subText .= '$keyfile,' if $commandLine =~ m/KEYFILE/;
$subText .= '$vcf,' if $commandLine =~ m/VCF/;
$subText .= '$optionsHachees) = @_;'."\n";


#Testing the type of IN OUT and mandatory to generate the input output files variables

$subText .= "\tmy \$validation = 0;\n\tswitch (1)\n\t{";

my %formatValidator = (
                        fasta =>"\n\t\t".'case ($fileIn =~ m/fasta|fa|fasta\.gz|fa\.gz$/i){$validation = 1 if (checkFormat::checkFormatFasta($fileIn) == 1)}',
                        fastq =>"\n\t\t".'case ($fileIn =~ m/fastq|fq|fastq\.gz|fq\.gz$/i){$validation = 1 if (checkFormat::checkFormatFastq($fileIn) == 1)}',
                        sam => "\n\t\t".'case ($fileIn =~ m/sam$/i){$validation = 1 if (checkFormat::checkFormatSamOrBam($fileIn) == 1)}',
                        bam => "\n\t\t".'case ($fileIn =~ m/bam$/i){$validation = 1 if (checkFormat::checkFormatSamOrBam($fileIn) == 2)}',
                        vcf => "\n\t\t".'case ($fileIn =~ m/vcf|vcf\.gz$/i){$validation = 1 if (checkFormat::checkFormatVcf($fileIn) == 1)}',
                        gff => "\n\t\t".'case ($fileIn =~ m/gff$/i){$validation = 1 if (checkFormat::checkFormatGff($fileIn) == 1)}',
                        bed => "\n\t\t".'case ($fileIn =~ m/bed$/i){$validation = 1 if (checkFormat::checkFormatBed($fileIn) == 1)}'
                        #'phylip'=>"\n\t\t".'case ($fileIn =~ m/phy|phylip$/i){$validation = 1 if (checkFormat::checkFormatPhylip($fileIn) == 1)}'
                        # gtf
                        # nwk, newik, nk
                        # bcf
                        # ped
                        # intervals
                        # txt
                        # sai
                       );

#Format case addition
my @listIn = split/,/,$in;
foreach my $format (@listIn)
{
    $subText.= $formatValidator{$format};
}
#finishing the infos
$subText .= "\n\t\telse {toolbox::exportLog(\"ERROR: $module::$function : The file \$fileIn is not a $in file\\n\",0);}\n\t};\n\tdie (toolbox::exportLog(\"ERROR: $module::$function : The file \$fileIn is not a $in file\\n\",0)) if \$validation == 0;";

#Extract options

$subText .= "\t#Picking up options\n\t".'my $options="";'."\n";
$subText .= "\t".'$options = toolbox::extractOptions($optionsHachees) if $optionsHachees;'."\n\n";

#Adding options in the command

$commandLine =~ s/\[options\]/\$options/;

#Generating command
$subText .= "\t#Execute command\n";
$commandLine =~ s/FILEIN/\$fileIn/;
$commandLine =~ s/FILEOUT/\$fileOut/;
$commandLine =~ s/REFERENCE/\$reference/;
$commandLine =~ s/GFF/\$gff/;
$commandLine =~ s/KEYFILE/\$keyfile/;
$commandLine =~ s/VCF/\$vcf/;
$commandLine =~ s/\[options\]/\$options/;
$subText .= "\tmy \$command = \"\$$commandLine\" ;";

#Executing command and return
$subText .= "\n\treturn 1 if (toolbox::run(\$command));\n";
$subText .= "\ttoolbox::exportLog(\"ERROR: $module::$function : ABORTED\\n\",0);\n}";


#Print in module
my $sedCom = "sed -i 's/^1;\$//' $toggle/modules/$moduleFile";
system("$sedCom") and die ("\nCannot remove the previous 1;:\n$!\n");

open (my $fhModule, ">>", "$toggle/modules/$moduleFile") or die ("\nCannot open for writing the file $moduleFile:\n$!\n");
print $fhModule $subText;
print $fhModule "\n\n1;\n";
close $fhModule;

#localConfig

#check if the software name already exists
$grep ="";
$grepCom = "grep -c \"$module\" $toggle/modules/localConfig.pm 2>/dev/null";
$grep = `$grepCom`;
chomp $grep;
if ($grep)
{   #the function exists...
    warn "\nThe function already exists in localConfig.pm\n";
}
else
{
    #the function is not registered

    #adding the soft name in the export
    my $sedCom = "sed -i 's/^our \@EXPORT=qw(/our \@EXPORT=qw(\$$module /' $toggle/modules/localConfig.pm";
    system("$sedCom") and die ("\nCannot add the software in the \@EXPORT:\n$!\n");

    #adding the path at the end of the file
    $sedCom = "sed -i 's/^1;\$//' $toggle/modules/localConfig.pm";
    system("$sedCom") and die ("\nCannot remove the previous 1;:\n$!\n");

    my $localLine = "#Path to $module\nour \$$module=\"/path/to/$module\";\n\n";
    open (my $fhLocal, ">>", "$toggle/modules/localConfig.pm") or die ("\nCannot open for writing the file localConfig.pm:\n$!\n");
    print $fhLocal $localLine;
    print $fhLocal "\n1;\n";
    close $fhLocal;
}

#function name
$grep ="";
$grepCom = "grep -c \"$function\" $toggle/modules/softwareManagement.pm 2>/dev/null";
$grep = `$grepCom`;
chomp $grep;
if ($grep >= 2)
{   #the function exists...
    warn "\nThe function already exists in softwareManagement.pm\n";
}
else
{
    #the function is not registered

    #adding the soft correct name
    open (my $fhTmp, ">", "/tmp/tempModule.pm") or die ("\nCannot open for writing the temp file /tmp/tempModule.pm:\n$!\n");
    open (my $fhRead, "<", "$toggle/modules/softwareManagement.pm") or die ("\nCannot open for reading the module softwareManagement:\n$!\n");

    while (my $line = <$fhRead>)
    {
        chomp $line;
        if ($line =~ m/#NEW SOFT ADDED AUTOMATICALLY/)
        {
            my $subFunction = $function;
            $subFunction =~ s/$module//;
            my $newName = "\n\t#FOR $function\n\tcase (\$name".' =~ m/^'."$module".'[\s|\.|\-| \/|\\|\|]*'."$subFunction".'/i'."){\$correctedName=\"$function\";} #Correction for $function";
            $line .= "\n";
            $line .= $newName;
        }
        elsif ($line =~ m/#INFOS FOR NEW TOOLS/)
        {
            my $newInfos = "\t'$function'=>{'IN' => '$in',\n\t\t\t'OUT'=>'$out',\n\t\t\t";
            if ($mandatory ne "")
            {
                $newInfos .="'MANDATORY' => '$mandatory',\n\t\t\t";
            }
            $version =~s/"/'/g;
            $newInfos .="'cmdVersion' => \"$version\"},\n";

            $line .= "\n".$newInfos;
        }
        elsif ($line =~ m/#LOG INFO FOR NEW TOOLS/)
        {
            my $subFunction = $function;
            $subFunction =~ s/$module//;
            my $newName = "\n\t#FOR $function\n\tcase (\$softOrder".' =~ m/^'."$module".'.*/i'."){\$softPathVersion{\"$function\"}=\`\$softInfos{\"$function\"}{'cmdVersion'} if not defined \$softPathVersion{\"$function\"};\n\t\t\$softPath{\"$function\"}=\$function if not defined \$softPath{\"$function\"};\n}";
            $line .= "\n";
            $line .= $newName;
        }


        print $fhTmp $line;
        print $fhTmp "\n";
    }
    close $fhTmp;
    close $fhRead;

    my $replaceCom = "cp /tmp/tempModule.pm $toggle/modules/softwareManagement.pm && rm -f /tmp/tempModule.pm";
    system ("$replaceCom") and die ("\nCannot replace the softwareManagement.pm file:\n$!\n");
}

#BLOCK CREATION
my $blockName=$function."Block.txt";
if (-e "$toggle/onTheFly/$blockName") # verifying block does not exist
{
    print "The block exists already\n";
}
else
{
    open (my $fhBlock, ">>", "$toggle/onTheFly/$blockName") or die ("\nCannot open for writing the file $blockName :\n$!\n");
    my $localLine;
    my $format;
    my $formatOut;
    my @formatList;
    my $fileInType;
    my $fileOutName;
    my $functionName=uc($function);

    #populating $format to identify input and output formats
    switch (1)
    {
        case ($in =~ m/fastq/ ) { push @formatList, "fastq\$\|fastq.gz\$\|fq\$\|fq.gz\$"; $fileInType='$fastqForwardIn'}
        case ($in =~ m/fasta/ ) { push @formatList, "fasta\$\|fasta.gz\$\|fa\$\|fa.gz\$"; $fileInType='$fastaFileIn'}
        case ($in =~ m/sam/ ) { push @formatList, "sam\$"; $fileInType='$samFileIn'}
        case ($in =~ m/bam/ ) { push @formatList, "bam\$"; $fileInType='$bamFileIn'}
        case ($in =~ m/vcf/ ) { push @formatList, "vcf\$\|vcf.gz\$"; $fileInType='$vcfFileIn'}
        case ($in =~ m/bed/ ) { push @formatList, "bed\$\|bed.gz\$"; $fileInType='$bedFileIn'}
        case ($in =~ m/ped/ ) { push @formatList, "ped\$\|ped.gz\$"; $fileInType='$pedFileIn'}
        case ($in =~ m/phylip/ ) { push @formatList, "phy\$"; $fileInType='$phylipFileIn'}
        case ($in =~ m/readseq/ ) { push @formatList, "readseq\$"; $fileInType='$readseqFileIn'}
        else {push @formatList, "*"; $fileInType='$fileIn'};
    }
    $format=join('\|',@formatList);
    my $formatPlain = $format;
    $formatPlain =~ s/\$/ /g;
    $formatPlain =~ s/\|/ /g;

     switch (1)
    {
        case ($out =~ m/fastq/ ) {$fileOutName='$fastqForwardOut';$formatOut="fastq"}
        case ($out =~ m/fasta/ ) {$fileOutName='$fastaFileOut';$formatOut="fasta"}
        case ($out =~ m/sam/ ) {$fileOutName='$samFileOut';$formatOut="sam"}
        case ($out =~ m/bam/ ) {$fileOutName='$bamFileOut';$formatOut="bam"}
        case ($out =~ m/vcf/ ) {$fileOutName='$vcfFileOut';$formatOut="vcf"}
        case ($out =~ m/bed/ ) {$fileOutName='$bedFileOut';$formatOut="bed"}
        else {$fileOutName='$fileOut'; $formatOut="txt"};
    }
    #writting general block
    $localLine.="

###########################################
## Block for $function
###########################################

#Correct variable populating

foreach my \$file (\@{\$fileList}) #Checking the type of files that must be $formatPlain
{
    if (\$file =~ m/$format\$/) # the file type is normally $formatPlain
    {
        if ($fileInType ne \"NA\") # Already a $formatPlain recognized
        {
            toolbox::exportLog(\"ERROR : \$0 : there are more than one single $formatPlain file at \$stepName step.\\n\",0);
        }
        else
        {
            $fileInType = \$file;
        }
    }
}

if ($fileInType eq \"NA\") #No $format file found in the previous folder
{
    toolbox::exportLog(\"ERROR : \$0 : No $formatPlain file found in \$previousDir at step \$stepName.\\n\",0);
}

\$softParameters = toolbox::extractHashSoft(\$optionRef,\$stepName);                                # recovery of specific parameters of $function

$fileOutName = \"\$newDir\".\"/\".\"\$readGroup\".\"\.$functionName\.$formatOut\";
$module::$function($fileInType,$fileOutName,";

$localLine .= '$reference,' if $commandLine =~ m/REFERENCE/i;
$localLine .= '$gff,' if $commandLine =~ m/GFF|GTF/i;
$localLine .= '$keyfile,' if $commandLine =~ m/KEYFILE/i;
$localLine .= '$vcf,' if $commandLine =~ m/VCF/i;


$localLine .= "\$softParameters);   # Sending to $function

    ";
    print $fhBlock $localLine;
    close $fhBlock;
}





#TEST MODULE
my $moduleTest=$module."_test.t";
if (-e "$toggle/test/modules/$moduleTest") # verifying module does not exist
{
    print "The module $moduleTest exists already\n";
}
else
{
    open (my $fhModule, ">>", "$toggle/test/modules/$moduleTest") or die ("\nCannot open for writing the file $moduleTest :\n$!\n");
    my $localLine;
    $localLine.="#!/usr/bin/env perl

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
## COMMON MODULE TEST HEADER
######################################################################################################################################

use strict;
use warnings;
use Data::Dumper;

use Test::More \'no_plan\'; #Number of tests, to modify if new tests implemented. Can be changed as \'no_plan\' instead of tests=>11 .
use Test::Deep;

# Load localConfig if primary test is successful
use_ok(\'localConfig\') or exit;
use localConfig;


########################################
# Extract automatically tool name and sub name list
########################################
my (\$toolName,\$tmp) = split /_/ , \$0;
my \$subFile=\$toggle.\"/modules/\".\$toolName.\".pm\";
my \@sub = `grep \"^sub\" \$subFile`or die (\"ERROR: \$0 : Cannot extract automatically sub name list by grep command \\n\$!\\n\");


########################################
#Automatically module test with use_ok and can_ok
########################################

use_ok(\$toolName) or exit;
eval \"use \$toolName\";

foreach my \$subName (\@sub)
{
    chomp (\$subName);
    \$subName =~ s/sub //;
    can_ok(\$toolName,\$subName);
}

#########################################
#Preparing test directory
#########################################
my \$testDir=\"\$toggle/dataTest/\$toolName\".\"TestModule\";
my \$cmd=\"rm -Rf \$testDir ; mkdir -p \$testDir\";
system(\$cmd) and die (\"ERROR: \$0 : Cannot execute the test directory \$testDir (\$toolName) with the following cmd \$cmd\\n\$!\\n\");
chdir \$testDir or die (\"ERROR: \$0 : Cannot go into the test directory \$testDir (\$toolName) with the chdir cmd \\n\$!\\n\");


#########################################
#Creating log file
#########################################
my \$logFile=\$toolName.\"_log.o\";
my \$errorFile=\$toolName.\"_log.e\";
system(\"touch \$testDir/\$logFile \$testDir/\$errorFile\") and die \"\\nERROR: \$0 : cannot create the log files \$logFile and \$errorFile: \$!\\nExiting...\\n\";

######################################################################################################################################
######################################################################################################################################



##########################################
### input output Options
##########################################

my \%optionsHachees = ();                # Hash containing informations
my \$optionHachees = \\%optionsHachees;   # Ref of the hash

##########################################
##### $module::$function
##########################################



is($module::$function(FILEIN,FILEOUT,\$optionHachees),1,'$module::$function  - Test for $function running');

# expected output test
my \$observedOutput = \`ls\`;
my \@observedOutput = split /\\n/,\$observedOutput;
my \@expectedOutput = (\'FILEOUT\',\'".$module."_log.e\',\'".$module."_log.o\');
is_deeply(\\\@observedOutput,\\\@expectedOutput,\'$module::$function - output list\');

################ TODO add test for output content
    ";

    print $fhModule $localLine;
}

close $fhModule;

#TEST BLOCK
my $blockTest=$module."Block.pl";
if (-e "$toggle/test/blocks/$blockTest") # verifying block does not exist
{
    print "The block $blockTest exists already\n";
}
else
{
    open (my $fhBlock, ">>", "$toggle/test/blocks/$blockTest") or die ("\nCannot open for writing the file $blockTest :\n$!\n");
    my $localLine;
    $localLine.="#!/usr/bin/env perl

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
use warnings;
use Test::More 'no_plan';
use Test::Deep;
use fileConfigurator;
use localConfig;
use Data::Dumper;

#####################
## PATH for data test
#####################

# references files
my \$dataRefIrigin = \"\$toggle/data/Bank/referenceIrigin.fasta\";
# input file
my \$dataFastq=\"\$toggle/data/testData/fastq/pairedTwoIndividusIrigin\";


print \"\\n\\n#################################################\\n\";
print \"#### TEST $module \\n\";
print \"#################################################\\n\";

# Remove files and directory created by previous test
my \$testingDir=\"\$toggle/dataTest/$module-noSGE-Blocks\";
my \$cleaningCmd=\"rm -Rf \$testingDir\";
system (\$cleaningCmd) and die (\"ERROR: \$0 : Cannot remove the previous test directory with the command \$cleaningCmd \\n\$!\\n\");

#Creating config file for this test
my \@listSoft = (\"$module\");
fileConfigurator::createFileConf(\\\@listSoft,\"blockTestConfig.txt\");

my \$runCmd = \"toggleGenerator.pl -c blockTestConfig.txt -d \".\$dataFastq.\" -r \".\$dataRefIrigin.\" -o \".\$testingDir;
print \"\\n### Toggle running : \$runCmd\\n\";
system(\"\$runCmd\") and die \"#### ERROR : Can't run TOGGLE for $module\";

# check final results

# expected output content
my \$observedOutput = \`ls \$testingDir/finalResults\`;
my \@observedOutput = split /\\n/,\$observedOutput;
my \@expectedOutput = (\'FILE.EXTENTION\',\'FILE.EXTENTION\');

is_deeply(\\\@observedOutput,\\\@expectedOutput,\'toggleGenerator - Two FILES (no SGE) $module::$function file list \');

# expected output value TO DO

    ";

    print $fhBlock $localLine;
    close $fhBlock;
}


# list of files to check
print "\n#######################################\nFinished...\n\n Please have a look to the following files to check if everything is Ok:\n\n
    - modules/$moduleFile ##NOTE: Please check if the variable is noted as \$bwa and not /usr/bin/bwa !!
    - modules/localConfig.pm
    - modules/softwareManagement.pm ##NOTE: Please check if the variable is noted as \$bwa and not /usr/bin/bwa !!
    - modules/fileConfigurator.pm ##NOTE: Please add line $function =>[\"\"] to add default value for test block
    - onTheFly/$blockName
    - test/modules/$moduleTest
    - test/blocks/$blockTest\n\n#######################################\n";

exit;

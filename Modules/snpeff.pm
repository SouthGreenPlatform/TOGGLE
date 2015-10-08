package snpeff;

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot
#
###################################################################################################################################


use strict;
use warnings;
use localConfig;
use toolbox;
use Data::Dumper;

sub dbCreator { #will prepare a snpEff database 
    my ($genomeFile, $referenceFasta, $name, $optionsHachees)=@_; # $name is not an obligation; If no $name is provided, it will be the same as $genomeFile w/o the extension
    
    #reading the options for snpEff
    my $options;
    if ($optionsHachees)
	{
        $options=toolbox::extractOptions($optionsHachees);		##Get given options
        }
    
    #Giving a name if no name provided
    if (not defined $name)
	{# $name not provided
	#will create a $name
	$name = toolbox::extractName($referenceFasta);	
	}
    
        
    #Checking if the name exists in the config file provided
    #Extraction of config file name from the optionHachees
    my $configFile="snpeff.config";
    if (exists $$optionsHachees{"-c"})
	{
	#a config file different from the standard one has been provided
	$configFile = $$optionsHachees{"-c"};
	}
    #Verifying if existing
    my $testCommand="grep $name $configFile";
    my $commandResult = `$testCommand`;
    chomp $commandResult;
    if ($commandResult)
	{
	#The name is already used in a formatted database signaled in the config file provided
	toolbox::exportLog("WARNING : snpEff::dbCreator : The database name $name is already in use in the current configuration. please dowload/use it or change the name\n",2); # Give a warning but not an error, as we can still work with this already existing db
        return 0;
	}
    
    
    #Detecting the type of genome file based on file extension
    my $fileOption;
    my $fileType;
    
    if ($genomeFile =~ m/\.gtf$/ or $genomeFile =~ m/\.gtf\.gz$/)
	{
	#the file is a gtf file, gunzipped or not

	if(toolbox::gtfFormatValidator($genomeFile)==1)		##The file is a good GTF one
	    {
            toolbox::exportLog("INFOS: snpEff::dbCreator : the file $genomeFile is a GTF file\n",1);
	    }
	else
	    {
            toolbox::exportLog("ERROR: snpEff::dbCreator : The file $genomeFile is not a correct GTF file. ABORTED\n",0); #Not a GTF file
            return 0;
	    }
	
	$fileOption = "-gtf22";
	$fileType = "gtf";
	}
    elsif ($genomeFile =~ m/\.gff$/ or $genomeFile =~ m/\.gff\.gz$/)
	{
	
	if(toolbox::gff3FormatValidator($genomeFile)==1)		##The file is a good GFF one
	    {
            toolbox::exportLog("INFOS: snpEff::dbCreator : the file $genomeFile is a GFF file\n",1);
	    }
	else
	    {
            toolbox::exportLog("ERROR: snpEff::dbCreator : The file $genomeFile is not a correct GFF file. ABORTED\n",0); #Not a GTF file
            return 0;
	    }
	#the file is a GFF file, gunzipped or not
	$fileOption = "-gff3";
	$fileType = "gff";
	}
    
    #Verification of the sequence file as a FASTA one
    if (toolbox::dnaFastaFormatValidator($referenceFasta))
	{
	#The file is a correct FASTA one
	toolbox::exportLog("INFOS : snpEff::dbCreator : the file $referenceFasta is a FASTA file\n",1);
        }
    else
        {
        toolbox::exportLog("ERROR : snpEff::dbCreator : The file $referenceFasta is not a correct FASTA file. ABORTED\n",0); #Not a GTF file
        return 0;
	}
    
    #Preparing the folders for copying the files
    my $grepCom = "grep \"data.dir =\" $configFile";
    my $grepResult = `$grepCom`;
    chomp $grepResult;
    my ($nothing, $dataFolder)=split / = /, $grepResult;
    $dataFolder.="/".$name;
    #Creating the data folder
    my $mkdirLocalCommand ="mkdir -p $dataFolder  && echo \"Dir ok\"";
    print $mkdirLocalCommand,"\n";
    system ($mkdirLocalCommand) and die ("ERROR: snpEff::dbCreator : Cannot create the $dataFolder folder: $!\nABORTED\n!");
    print $dataFolder,"\n";
    my $lsCom=`ls ../DATA-TEST/snpeffTestDir/`;
    my $lsCom2=`ls $dataFolder`;
    print "++../DATA-TEST/snpeffTestDir/\n", $lsCom,"\n**\n",$dataFolder,"\n",$lsCom2;
    
   #BUGGY !!  toolbox::makeDir($dataFolder,0);
    #Copying files
    my $cpGenomeFileCommand = "cp ".$genomeFile." ".$dataFolder."/genes.".$fileType;
    system ("$cpGenomeFileCommand") and die ("ERROR : snpEff::dbCreator : Cannot copy the file $genomeFile in $dataFolder: $!\n. ABORTED\n");
    my $cpSequenceFileCommand = "cp ".$referenceFasta." ".$dataFolder."/sequences.fa";
    system ("$cpSequenceFileCommand") and die ("ERROR : snpEff::dbCreator : Cannot copy the file $referenceFasta in $dataFolder: $!\n. ABORTED\n");
    #Adding the info of the genome to the config file
    my $addline = "#".$name." genome\\n".$name.".genome : ".$name;
    my $addingCommand = "echo \"".$addline."\"| cat - >> ".$configFile;
    print $addingCommand,"\n";
    system($addingCommand) and die ("ERROR : snpEff::dbCreator : Cannot adding the new genome info in the config file $configFile: $!\n. ABORTED\n");
    
    #Preparing the command
    my $command="$snpEff build $options $fileOption -v $name";
    print $command,"\n";
    if(toolbox::run($command)==1)		##The command should be executed correctly (ie return) before exporting the log
	{
            toolbox::exportLog("INFOS: snpEff::dbCreator : correctly done\n",1);		# dbCreator has been correctly done
            return 1;
        }
    else
        {
            toolbox::exportLog("ERROR: snpEff::dbCreator : ABORTED\n",0);		# dbCreator has not been correctly done
            return 0;
        }
   
}

sub snpeffAnnotation { #Will annotate a VCF file based on a already prepared db for SNPeff
    my ($vcf,$database, $outputName, $optionsHachees)=@_;
    # the vcf file $vcf will be annotated using the infos from the already formatted database $database, and the resulting vcf will be outputted in $outputName
    #The options may be the gatk compatibility (-o gatk) and the config file positions
    
    my $options="";

    if ($optionsHachees)
	{
            $options=toolbox::extractOptions($optionsHachees);		##Get given options
        }
    my $command=$snpEff." ".$options." ".$database." ".$vcf." > ".$outputName;
    
    #Execute command
    if(toolbox::run($command)==1)		##The command should be executed correctly (ie return) before exporting the log
	{
            toolbox::exportLog("INFOS: snpEff::snpEffAnnotation : correctly done\n",1);		# snpEffAnnotation has been correctly done
            return 1;
        }
    else
        {
            toolbox::exportLog("ERROR: snpEff::snpEffAnnotation : ABORTED\n",0);		# snpEffAnnotation has not been correctly done
            return 0;
        }
}


1;

=head1 NAME

    Package I<SNPeff> 

=head1 SYNOPSIS

	use snpeff;
    
	snpeff::snpeffAnnotation($vcf,$database,$outputName)
	
=head1 DESCRIPTION

    Package SNPeff (***, http:// ) is a software package for refining SNP/Indel information based on annotation.

=head2 FUNCTIONS

=head3 snpEff::snpeffAnnotation

This module is the core function of SNPeff
It will annotate a given VCF using an annotation database already formatted for SNPeff.
The resulting VCF will have the annotation information (and effects of SNP) within it.


=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
Written by Cecile Monat, Ayite Kougbeadjo, Marilyne Summo, Cedric Farcy, Mawusse Agbessi, Christine Tranchant and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut

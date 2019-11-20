package scheduler;

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

use localConfig;

use toolbox;

use Data::Dumper;

##################################
#
#LAUNCHING
#
##################################

our ($commandLine, $requirement, $sample, $configInfo, $jobList, %jobHash, @errorFileList, $outputDir);

#Here is the core commands for any scheduler: run, acct and queue command.

my %commands =('run' => {'normalRun'=>'sh',
                         'sge'=>'qsub',
						 'slurm'=>'sbatch',
						 'mprun'=>'ccc_msub',
						 'lsf'=>'bsub'},
			   'queue' => {'normalRun'=>'',
                          'sge'=>'qstat -u \$USER',
						  'slurm'=>'squeue',
						  'mprun'=>'ccc_mstat -u \$USER',
						  'lsf'=>'bjobs -u \$USER'});

##DEBUG toolbox::exportLog(Dumper(\%commands),0);

#Here are the infos for parsing data in waiting system

my %parsings = ('JIDposition' => {	'normalRun'=>'',
                                    'sge'=>2,
									'slurm'=>3,
									'mprun'=>3,
									'lsf'=>1},
				'JIDsplitter' => {'normalRun'=>'',
                                  'sge'=>"\\s",
								  'slurm'=>"\\s",
								  'mprun'=>"\\s",
								  'lsf'=>"<\|>"});

##DEBUG toolbox::exportLog(Dumper(\%parsings),0);

sub checkingCapability { #Will test the capacity of launching using various schedulers on the current system

    my $capabilityValue;

    #SGE test
    $capabilityValue = `qsub -help 2>/dev/null | grep usage`; #Will provide a not-empty output if SGE is installed
    chomp $capabilityValue;
	return "sge" if $capabilityValue;

    #SLURM test
    $capabilityValue=`sbatch -V 2>/dev/null | grep slurm`; #Will provide a not-empty output if SLURM is installed
    chomp $capabilityValue;
	return "slurm" if $capabilityValue;


    #mprun test
    $capabilityValue=`ccc_mprun -h 2>&1 | grep usage`; #Will provide a not-empty output if mprun is installed
    chomp $capabilityValue;
    return "mprun" if $capabilityValue;

    #lsf test
    $capabilityValue=`bsub -h 2>&1 | grep -i Synopsis`; #Will provide a not-empty output if lsf is installed
    chomp $capabilityValue;
	return "lsf" if $capabilityValue;

	#None
	return "normaLRun" unless $capabilityValue;

}

#Requirement means the job needs to be achieved for the next steps or it is not a blocking job?
# A zero value (0) means mandatory
# Any non-zero value means not blocking job

sub launcher {

	#Global function for launching, will recover the command to be launch and will re-dispatch to normal or other scheduler

    ($commandLine,$requirement, $sample, $configInfo) = @_;

    #Picking up sample name

    $sample=`basename $sample` or toolbox::exportLog("ERROR : scheduler::launcher : Cannot pickup the basename for $sample: $!\n",0);
    chomp $sample;
    if ($sample =~ m/^\d/)
    {
	#the sample name starts with a number, not a good idea for schedulers
	toolbox::exportLog("WARNING: scheduler::launcher: the sample name for $sample is starting with a number. Its scheduler job will be named s$sample for ensure its successful launching\n",2);
	$sample = "s".$sample;
    } 
     my $schedulerType = &checkingCapability;
	##DEBUG toolbox::exportLog("INFO : scheduler::launcher : Scheduler is $schedulerType\n",0);

    my $runOutput = &schedulerRun($schedulerType);


	if ($runOutput eq "0" && $requirement == 0)
    {
	#The job has to succeed either it will kill all other jobs
        toolbox::exportLog("ERROR: scheduler::launcher on $sample using the scheduler $schedulerType: ".$commandLine."\nThe job cannot be achieved and is mandatory, thus the whole analysis is stop\n",0);
    	return 0;
    }

        #Picking up the output error file
    my @listOne = split /-d\s+/, $commandLine;
    my ($folderOut) = split /\s+/, $listOne[1];
    my $errorLog = `basename $folderOut`;
    chomp $errorLog;
    $errorLog=$folderOut."/".$errorLog."_log.e";
    
    return ($runOutput, $errorLog);
}

#sub normalRun
#{ #For running in normal mode, ie no scheduler
#
#    #BASED ON TOOLBOX::RUN, but will not stop the whole pipeline for an error
#    use Capture::Tiny qw(capture);
#
#    toolbox::exportLog("INFOS: scheduler::normalRun : $commandLine\n",1);
#
#    ##Execute the command
#    my ($result,$stderr)=capture {` $commandLine `};
#
#    ##Log export according to the error
#    if ($?==0) #Success, no error
#    {
#		return 1;
#    }
#    else  #Error, the job cannot be achieved for any reason
#    {
#	##DEBUG
#		toolbox::exportLog("WARNING: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\nThe $sample data set has provoked and error, and was not analyzed anymore.\n",2);
#		return 0;
#}

sub schedulerRun
{ #For any scheduler,will launch a script encapsulating the command line

	my ($schedulerType) = @_;
    #Returning to normal run if no scheduler config provided even if scheduler exists
    $schedulerType = "normalRun" unless toolbox::extractHashSoft($configInfo,$schedulerType);
    
    my $schedulerOptions = "";
    if ($schedulerType ne "normalRun")
    {
        my $schedulerOptionsHash=toolbox::extractHashSoft($configInfo,$schedulerType);
        $schedulerOptions=toolbox::extractOptions($schedulerOptionsHash);
    }

    
	#Picking up ENV variable
    my $envOptionsHash=toolbox::extractHashSoft($configInfo,"env");
    my $envOptions=toolbox::extractOptions($envOptionsHash,"","\n");

	#Picking up location
	my $location = `pwd`;
	chomp $location;

	#Creating a folder for scripts
	my $schedulerFolder = $location."/schedulerJobs";
	unless (-d $schedulerFolder)
	{
		toolbox::makeDir($schedulerFolder,0);
	}

    #Adding scheduler options
    my $launcherCommand = $commands{'run'}{$schedulerType}." ".$schedulerOptions;
    
   
    ##DEBUG    toolbox::exportLog("NORMAL: scheduler::Run: the errorLog is $errorLog;",1);

    #Creating the bash script for slurm to launch the command
    #my $date =`date +%Y_%m_%d_%H_%M_%S`;
    #chomp $date;
    my $scriptName=$schedulerFolder."/".$sample."_schedulerScript.sh";
    my $bashScriptCreationCommand= "echo \"#!/bin/bash\n\n".$envOptions."\n".$commandLine."\n\n\nexit 0;\" | cat - > $scriptName && chmod 777 $scriptName";
    toolbox::run($bashScriptCreationCommand,"noprint");
    $launcherCommand.=" ".$scriptName;
    $launcherCommand =~ s/ +/ /g; #Replace multiple spaces by a single one, to have a better view...

	#launching the job through a bash script
    my $currentJID = `$launcherCommand`;

    if ($!) #There are errors in the launching...
    {
        toolbox::exportLog ("WARNING : $0 : Cannot launch the job for $sample: $!\n",2);
        $currentJID = "";
    }
    
    #Parsing infos and informing logs
    chomp $currentJID;
    $currentJID = "SerialJob" if $schedulerType eq "normalRun"; 
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return -1; #Returning to launcher subprogram the error type
    }


    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in $commands{'run'}{$schedulerType} mode through the command:\n\t$launcherCommand\n",1);

    #Picking job ID
    if ($schedulerType ne "normalRun")
    {
        my @infosList=split ($parsings{'JIDsplitter'}{$schedulerType}, $currentJID);
        $currentJID = $infosList[$parsings{'JIDposition'}{$schedulerType}];
    }

	##DEBUG	toolbox::exportLog($currentJID,2);

    return $currentJID;
}


##################################
#
# WAITING for schedulers ONLY!
#
##################################

sub waiter
{ #Global function for waiting, will recover the jobID to survey and will re-dispatch to scheduler

    ($jobList,my $jobsInfos, $outputDir, my $report) = @_;

    %jobHash = %{$jobsInfos};

    my $schedulerType = &checkingCapability;
    
    my $stopWaiting = &schedulerWait($schedulerType, $outputDir, $report);

    return $stopWaiting;
}


sub schedulerWait
{
    my ($schedulerType, $ouputDir, $report) = @_;
    my $nbRunningJobs = 1;
    my @jobsInError=();

    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {
      last if $schedulerType eq "normalRun";
      #Picking up the number of currently running jobs
      my $statCommand = $commands{'queue'}{$schedulerType}." | egrep -c \"$jobList\"";
      $nbRunningJobs = `$statCommand`;
      chomp $nbRunningJobs;
      sleep 15;
    }

    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 5;
	


	#open file for report job info
	my $openNameFile = $outputDir."/sample_parallel_table.tex";
	
	if (scalar(keys%jobHash) == 1 && defined($jobHash{"global"}))
	{
		 $openNameFile = $outputDir."/sample_global_table.tex";
	}
	
	open (my $fhOut, '>', $openNameFile) or die "\nERROR: $0 : cannot open file $openNameFile. $!\nExiting...\n" if ($report);

	my $outputLineTex = "\\begin{table}[ht]
	\\centering
	\\begin{tabular}{l||r||r}
	Samples & Job ID & Status \\\\\\hline
	" if ($report);

    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tExitStatus
---------------------------------------------------------",1);

    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
        
        my ($currentLine, $outputLine);
        $outputLine = "$individual\t$jobHash{$individual}{output}\t";
		$outputLineTex .= "$individual & $jobHash{$individual}{output} & " if ($report);
        
        my $grepError = `grep "ERROR" $jobHash{$individual}{errorFile}`;
        chomp $grepError;
        
        if ($grepError)
        {
            $currentLine = "Error";
            push @jobsInError, $individual;
		}
        else
        {
            $currentLine = "Normal";
        }
       
    	$outputLine .= $currentLine;
		$outputLineTex .= $currentLine."\\\\" if ($report);
		toolbox::exportLog($outputLine,1);
		
		

    }
	print $fhOut $outputLineTex if ($report);
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
	print $fhOut "\n\\end{tabular}
	\\end{table}" if ($report);

    if (scalar @jobsInError)
    {
		#at least one job has failed
		return \@jobsInError;
    }
    close $fhOut if ($report);
    return 1;
}


1;

=head1 NAME

    Package I<scheduler>

=head1 SYNOPSIS

	use scheduler;

	scheduler::launcher($launcherCommand, $requirement, $sample, $configInfo);

	scheduler::waiter($jobList,$jobInfos);

=head1 DESCRIPTION

    Package scheduler will prepare and launch the different command in a scheduler if possible and requested. Either, it will launch it through the normal way, using a quite identical way to toolbox::run
    !!DIFFERENCE with toolbox::run is that an abnormal job will be a blocking one ONLY if requirement is of a zero value (0). Else, the job error will be warned but it will not block the whole system.

=head2 FUNCTIONS

=head3 scheduler::launcher

This module will prepare and decide to which scheduler sending the command. It takes as arguments the command $launcherCommand, the state of the job (0 mandatory finishing, non-0 not mandatory) $requirement, the sample name $sample and the $configInfo hash for configuration

=head3 scheduler::waiter

This module will allow the job launched through a scheduler to be wait to finish.
It takes as arguments the list of jobs $jobList and the reference of the hash for the informations about the jobs (sample name) $jobInfos



=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Cecile Monat, Christine Tranchant, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://www.southgreen.fr/>

=cut

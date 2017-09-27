package scheduler;

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

use localConfig;

use toolbox;

use Data::Dumper;

##################################
#
#LAUNCHING
#
##################################

our ($commandLine, $requirement, $sample, $configInfo, $jobList, %jobHash);

#Here is the core commands for any scheduler: run, acct and queue command.

my %commands =('run' => {'sge'=>'qsub',
						 'slurm'=>'sbatch',
						 'mprun'=>'ccc_msub',
						 'lsf'=>'bsub'},
			   'acct' => {'sge'=>'qacct -j',
						  'slurm'=>'qacct -j',
						  'mprun'=>'ccc_macct',
						  'lsf'=>'bacct'},
			   'queue' => {'sge'=>'qstat -u \$USER',
						  'slurm'=>'squeue',
						  'mprun'=>'ccc_mstat -u \$USER',
						  'lsf'=>'bjobs -u \$USER'});

##DEBUG toolbox::exportLog(Dumper(\%commands),0);

#Here are the infos for parsing data in waiting system

my %parsings = ('JIDposition' => {	'sge'=>2,
									'slurm'=>3,
									'mprun'=>3,
									'lsf'=>1},
				'JIDsplitter' => {'sge'=>"\\s",
								  'slurm'=>"\\s",
								  'mprun'=>"\\s",
								  'lsf'=>"<\|>"},
				'acctOutSplitter' => {'sge'=>"\\n",
									  'slurm'=>"\\n",
									  'mprun'=>"\\s+",
									  'lsf'=>"\\s+"},
				'acctOutText' => {'sge'=>"exit_status  0",
								  'slurm'=>"COMPLETED",
								  'mprun'=>"COMPLETED",
								  'lsf'=>"Total number of done jobs:\\s*1\\s*Total number of exited jobs:\\s*0\\s*"});

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
    
    my $schedulerType = &checkingCapability;
	##DEBUG toolbox::exportLog("INFO : scheduler::launcher : Scheduler is $schedulerType\n",0);
    
    my $runOutput;
       
    if (defined $$configInfo{$schedulerType})
	{
		$runOutput = &schedulerRun($schedulerType)
	}
	else
	{
		$runOutput = &normalRun
	};
    
    
	if ($runOutput == 0 && $requirement == 0)
    {
	#The job has to succeed either it will kill all other jobs
        toolbox::exportLog("ERROR: scheduler::launcher on $sample using the scheduler $schedulerType: ".$commandLine."\nThe job cannot be achieved and is mandatory, thus the whole analysis is stop\n",0);
    	return 0;
    }
    
    return $runOutput;
}

sub normalRun { #For running in normal mode, ie no scheduler
    
    #BASED ON TOOLBOX::RUN, but will not stop the whole pipeline for an error
    use Capture::Tiny qw(capture);
        
    toolbox::exportLog("INFOS: scheduler::normalRun : $commandLine\n",1);
    
    ##Execute the command
    my ($result,$stderr)=capture {` $commandLine `};
    
    ##Log export according to the error
    if ($?==0) #Success, no error
    {
		return 1;
    }
    else  #Error, the job cannot be achieved for any reason
    {
	##DEBUG
		toolbox::exportLog("WARNING: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\nThe $sample data set has provoked and error, and was not analyzed anymore.\n",2);
		return 0;
    }
    
}

sub schedulerRun
{ #For any scheduler,will launch a script encapsulating the command line
    
	my ($schedulerType) = @_;
    my $schedulerOptionsHash=toolbox::extractHashSoft($configInfo,$schedulerType);
    my $schedulerOptions=toolbox::extractOptions($schedulerOptionsHash);
    
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
	
    #Creating the bash script for slurm to launch the command
    #my $date =`date +%Y_%m_%d_%H_%M_%S`;
    #chomp $date;
    my $scriptName=$schedulerFolder."/".$sample."_schedulerScript.sh";
    my $bashScriptCreationCommand= "echo \"#!/bin/bash\n\n".$envOptions."\n".$commandLine."\n\nexit 0;\" | cat - > $scriptName && chmod 777 $scriptName";
    toolbox::run($bashScriptCreationCommand);
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
    
    unless ($currentJID) #The job has no output in STDOUT, ie there is a problem...
    {
        return 0; #Returning to launcher subprogram the error type
    }
    
    toolbox::exportLog("INFOS: $0 : Correctly launched for $sample in $commands{'run'}{$schedulerType} mode through the command:\n\t$launcherCommand\n",1);
    
    #Picking job ID
    my @infosList=split ($parsings{'JIDsplitter'}{$schedulerType}, $currentJID); 
    $currentJID = $infosList[$parsings{'JIDposition'}{$schedulerType}];
    
	
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
    
    ($jobList,my $jobsInfos) = @_;
    
    %jobHash = %{$jobsInfos};
        
    my $schedulerType = &checkingCapability;
    my $stopWaiting;
          
    if (defined $$configInfo{$schedulerType})
	{
		$stopWaiting = &schedulerWait($schedulerType)
	}

    return $stopWaiting;
}


sub schedulerWait
{
    my ($schedulerType) = @_;
    my $nbRunningJobs = 1;
    my @jobsInError=();
    
    ##Waiting for jobs to finish
    while ($nbRunningJobs)
    {  
      #Picking up the number of currently running jobs
      my $statCommand = $commands{'queue'}{$schedulerType}." | egrep -c \"$jobList\"";
      $nbRunningJobs = `$statCommand`;
      chomp $nbRunningJobs;
      sleep 50;
    }
    
    #Compiling infos about sge jobs: jobID, node number, exit status
    sleep 25;#Needed for acct to register infos...
    toolbox::exportLog("\n#########################################\nJOBS SUMMARY\n#########################################
\n---------------------------------------------------------
Individual\tJobID\tExitStatus
---------------------------------------------------------",1);
    
    foreach my $individual (sort {$a cmp $b} keys %jobHash)
    {
		my $acctCommand = $commands{'acct'}{$schedulerType}." ".$jobHash{$individual}." 2>&1";
		my $acctOutput = `$acctCommand`;
		my $outputLine;
		
		chomp $acctOutput;
		if ($acctOutput =~ "-bash: " or $acctOutput =~ "installed")
		{
		  #IF acct cannot be run on the node
		  $outputLine = "$individual\t$jobHash{$individual}\tNA\tNA\n";
		  toolbox::exportLog($outputLine,1);
		  next;
		}
		
		my @linesQacct = split ($parsings{'acctOutSplitter'}{$schedulerType}, $acctOutput);
		$outputLine = $individual."\t".$jobHash{$individual}."\t";
		my $currentLine ="Error";
		my $parserText = $parsings{'acctOutText'}{$schedulerType};
		  
		##DEBUG	toolbox::exportLog($parserText,2);
		
		while (@linesQacct) #Parsing the qacct output
		{
		  my $line = shift @linesQacct;
		  #Passing the header lines
		  next if $line =~ m/JobID/;
		  next if $line =~ m/^$/;
		  
		  
		  
		  if ($line =~ m/$parserText/) #Picking up exit status
		  {
			  $currentLine = "Normal";
			  last;
		  }
		  
		  next;
		}
		
		#The parsing text has not been found meaning still an error
		if ($currentLine eq "Error")
		{
		  push @jobsInError, $individual;
		}
		$outputLine .= $currentLine;       
		toolbox::exportLog($outputLine,1);
      
    }
    toolbox::exportLog("---------------------------------------------------------\n",1);#To have a better table presentation
  
    if (scalar @jobsInError)
    {
		#at least one job has failed
		return \@jobsInError;
    }
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

package scp;

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
use scheduler;

#Picking up the NFS mount

sub mountPoint {# Based on the mounted volume, will provide a hash with VolName => IP/DNS
	
	my $mounted = `df -h | grep ":"`;
	chomp $mounted;
	my @listVolumes = split /\n/, $mounted;
	my %volumes;
	while (@listVolumes)
	{
		my $currentVolume = shift @listVolumes;
		my @fields = split /\s+/, $currentVolume;
		my ($ip)=split /:/, $fields[0];
		my $name = pop @fields;
		$name =~ s/\///;
		$volumes{$name} = $ip;
	}
	
	return \%volumes;
}

sub transfer2node { #From a list of folder, will perform a rsync over ssh transfer (normally ok in cluster) and provide a list of new name
	#The $folderIn is the original data folder, the $tmpRoot is the basic folder on the node (eg /scratch, or /tmp)
 	my ($folderIn, $tmpRoot,$readgroup) = @_;
	
	my @path = split /\//, $folderIn;
	shift @path; # the form is "/my/path", thus the [0] position is undef 
	my $mount = &mountPoint;
	
	my $origin = $mount->{$path[0]};
	
	#node name
	my $node = `echo \$HOSTNAME`;
	chomp $node;
	
	#User name
	my $user = `echo \$USER`;
	chomp $user;
	
	#Job number acquisition
	my $jobNb;
	my $schedulerType = scheduler::checkingCapability;
	if ($schedulerType eq "sge")
	{
		$jobNb = `echo \$JOB_ID`;
	}
	elsif ($schedulerType eq "slurm" or $schedulerType eq "mprun")
	{
		$jobNb = `echo \$SLURM_JOBID`;
	}
	elsif ($schedulerType eq "lsf")
	{
		$jobNb = `echo \$LSF_JOBID`;
	}
	chomp $jobNb;
	
	#Node folder creation
	my $newFolder = $tmpRoot."/".$user."-".$jobNb."/".$readgroup;
	$newFolder =~ s/\s//g; #Removin extraspaces that hinder the transfer
	system ("mkdir -p $newFolder") and toolbox::exportLog("ERROR: scp::transfer2node: cannot create the destination folder for data $newFolder:\n\t$!\n",0);

	#Transfer
	
	my $rsyncCom = "rsync -vzurL ".$origin.":".$folderIn."/* ".$newFolder."/.";
	if (toolbox::run($rsyncCom)==1)
        {
            toolbox::exportLog("INFOS: scp::transfer2node Ok, data transferred from $origin to $node, in folder $newFolder\n",1);
	    
	    	
		my $refFolder;
		if (-d $origin.":".$folderIn."/../../referenceFiles")
		{
				$refFolder = $tmpRoot."/".$user."-".$jobNb."/referenceFiles";
				system ("mkdir -p $refFolder") and toolbox::exportLog("ERROR: scp::transfer2node: cannot create the destination folder for reference $refFolder:\n\t$!\n",0);
					$refFolder =~ s/\s//g; #Removin extraspaces that hinder the transfer
	
	
			my $rsyncRef = "rsync -vzurL ".$origin.":".$folderIn."/../../referenceFiles/* ".$refFolder."/.";
			toolbox::run($rsyncRef);
		}
	    
            return ($newFolder,$refFolder);
        }
        else
        {
            toolbox::exportLog("ERROR: scp::transfer2node error, data NOT transferred from $origin to $node\n",0);
        }
	
}

sub transfer2origin {#Will transfer data from a node to the original folder
	
	#the $localFolder is the folder on the node, the $folderIn is the original given folder (ie on the network)
	my ($localFolder,$folderIn) = @_;
	
	my @path = split /\//, $folderIn;
	shift @path; # the form is "/my/path", thus the [0] position is undef 
	my $mount = &mountPoint;
	
	my $origin = $mount->{$path[0]};
	
	#node name
	my $node = `echo \$HOSTNAME`;
	chomp $node;
	
	#Transfer
	my $rsyncCom = "rsync -vazur --exclude=0_initial ".$localFolder."/* ".$origin.":".$folderIn."/. && rm -R ".$localFolder;
	if (toolbox::run($rsyncCom)==1)
        {
            toolbox::exportLog("INFOS: scp::transfer2origin Ok, data transferred from $node to $origin, in folder $folderIn\n",1);
            return 1;
        }
        else
        {
            toolbox::exportLog("ERROR: scp::transfer2origin error, data NOT transferred from $node to $origin\n",0);
        }
	
}


1;

=head1 NAME

    Package I<scp> 

=head1 SYNOPSIS

        use scp;
	
	scp::transfer2node($folderIn,tmpRoot);
	
	scp::transfer2origin($localFolder,$folderIn);
    
        
=head1 DESCRIPTION

    this module allows the SCP transfer on nodes when working in scheduler/HPC mode, IF the SCP option is provided in the config file.
    It uses rsync to ensure a better transfer
    
=head2 FUNCTIONS

=head3 scp::transfer2node($folderIn,tmpRoot)

	From a list of folder, will perform a rsync over ssh transfer (normally ok in cluster) and provide a list of new name
	The $folderIn is the original data folder, the $tmpRoot is the basic folder on the node (eg /scratch, or /tmp)
	
=head3  scp::transfer2origin($localFolder,$folderIn)

	Will transfer data from a node to the original folder
	the $localFolder is the folder on the node, the $folderIn is the original given folder (ie on the network)


=head1 AUTHORS

Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
Written by Christine Tranchant, Cecile Monat, Laura Helou, Abdoulaye Diallo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot

=head1 SEE ALSO

L<http://toggle.southgreen.fr/>

=cut
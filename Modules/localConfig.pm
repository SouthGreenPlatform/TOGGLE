package localConfig;



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
use Exporter;

our @ISA=qw(Exporter);
our @EXPORT=qw($bwa $picard $samtools $GATK $cufflinks $pacBioToCA $cutadapt $fastqc $java $snpEff $toggle);

#toggle path
our $toggle="/home/sabotf/scripts/TOGGLE/";


#PATH for Mapping on cluster
our $java = "/usr/local/java/latest/bin/java -Xmx12g -jar";

our $bwa = "/usr/local/bin/bwa";
our $picard = "$java /home/sabotf/sources/picard-tools";


our $samtools = "/usr/local/samtools-0.1.18/samtools-0.1.18/samtools";

our $GATK = "/usr/java/jre1.7.0_51/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK-3.3/GenomeAnalysisTK.jar";
our $fastqc = "/usr/local/FastQC/fastqc";

#PATH for Cufflinks bin on cluster
our $cufflinks = "/usr/local/cufflinks-2.1.1.Linux_x86_64";

#Path for pacBioToCa
our $pacBioToCA = "/home/sabotf/sources/wgs/Linux-amd64/bin/pacBioToCA";

#Path for CutAdapt
our $cutadapt = "/usr/local/cutadapt-1.2.1/bin/cutadapt";

#Path for SNPeff
our $snpEff="$java /home/sabotf/sources/snpEff/snpEff.jar";

1;

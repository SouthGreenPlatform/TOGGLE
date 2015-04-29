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
our $toggle="/home/ortega-abboud/remplacement_marylyne/TOGGLEv02/TOGGLE";

#PATH for Mapping on cluster
our $java = "/usr/local/java/jre7/bin/java";

our $bwa = "/usr/local/bioinfo/bwa/0.7.12/bwa";
#our $picard = "$java -Xmx3g -jar /usr/local/bioinfo/picard-tools/1.119/";
our $picard = "$java -Xmx3g -jar /usr/local/bioinfo/picard-tools/1.130/";

our $samtools = "/usr/local/bioinfo/samtools/1.2/bin/samtools";
our $GATK = "$java -Xmx3g -jar /usr/local/bioinfo/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar";
our $fastqc = "/usr/local/bioinfo/FastQC/0.10.1/fastqc -j /usr/local/java/jre6/bin/java";

#PATH for Cufflinks bin on cluster
our $cufflinks = "/usr/local/bioinfo/cufflinks/cufflinks";

#Path for pacBioToCa
our $pacBioToCA = "/home/sabotf/sources/wgs/Linux-amd64/bin/pacBioToCA";

#Path for CutAdapt
our $cutadapt = "/usr/local/bioinfo/python/3.4.3/bin/cutadapt";

#Path for SNPeff
our $snpEff="$java /home/sabotf/sources/snpEff/snpEff.jar";

1;

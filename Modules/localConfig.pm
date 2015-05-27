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
our @EXPORT=qw($bwa $picard $samtools $GATK $cutadapt $fastqc $java $toggle $fastxTrimmer $tophat2 $bowtie2Build $bowtieBuild $htseqcount);

#toggle path
our $toggle="/path/to/toggle";

#PATH for Mapping on cluster
our $java = "/path/to/java -Xmx12g -jar";

our $bwa = "/path/to/bwa";
our $picard = "$java /path/to/picard-tools";

our $samtools = "/path/to/samtools_O.1.17";
our $GATK = "/path/to/java -Xmx12g -jar /path/to/GenomeAnalysisTK-3.3/GenomeAnalysisTK.jar";
our $fastqc = "/path/to/FastQC/fastqc";

#Path for CutAdapt
our $cutadapt = "/path/to/cutadapt-1.2.1/bin/cutadapt";

##### FOR RNASEQ analysis
#Path for fastq_trimmer
our $fastxTrimmer="/path/to/bin/fastx_trimmer";

#Path for tophat2
our $tophat2="/path/to/tophat2";

#path for bowtie2-build
our $bowtie2Build="/path/to/bowtie2-build";

#path for bowtie-build
our $bowtieBuild="/path/to/bowtie-build";

#path for htseqcount
our $htseqcount = "/path/to/htseq-count";

1;
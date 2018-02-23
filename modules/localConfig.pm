package localConfig;


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
use Exporter;

our @ISA=qw(Exporter);

our @EXPORT=qw($bwa $picard $samtools $GATK $cutadapt $fastqc $java $toggle $fastxTrimmer $tophat2 $bowtie2Build $bowtieBuild $htseqcount $cufflinks $cuffdiff $cuffmerge $tgicl $trinity  $stacks $snpEff $bamutils $crac $cracIndex $bowtie $bowtie2 $atropos $duplicationDetector $plink $bedtools $snmfbin $readseqjar $fastme $abyss $bam2cfg $breakDancer $pindel $fastqStats);

#toggle path
our $toggle="/path/to/toggleFolder";

#PATH for Mapping on cluster
our $java = "/path/to/java -Xmx12g -jar";

our $bwa = "/path/to/bwa";
our $picard = "$java /path/to/picard_tools/picard.jar";

our $samtools = "/path/to/samtools";
our $GATK = "$java -Xmx12g -jar /path/to/GenomeAnalysisTK.jar";
our $fastqc = "/path/to/fastqc";

#Path for CutAdapt
our $cutadapt = "/path/to/cutadapt";

##### FOR RNASEQ analysis
#Path for fastq_trimmer
our $fastxTrimmer="/path/to/fastx_trimmer";

#Path for tophat2
our $tophat2="/path/to/tophat2";

#path for bowtie2-build
our $bowtie2Build="/path/to/bowtie2-build";

#path for bowtie-build
our $bowtieBuild="/path/to/bowtie-build";

#path for htseqcount
our $htseqcount = "/path/to/htseq-count";

#path for Cufflinks
our $cufflinks = "/path/to/cufflinks";
our $cuffdiff = "/path/to/cuffdiff";
our $cuffmerge = "/path/to/cuffmerge";

#path for tgicl
our $tgicl = "/path/to/tgicl";

#path for trinity
our $trinity = "/path/to/trinity";

#path for process_radtags
our $stacks = "/path/to/process_radtags";

#path for snpEff
our $snpEff = "/path/to/snpEff/snpEff.jar";

#path for bamutils
our $bamutils = "/path/to/bamutils";

#path for atropos
our $atropos="/path/to/atropos";

#Path to bowtie
our $bowtie = "/path/to/bowtie";
our $bowtie2 = "/path/to/bowtie2";

#Path to crac
our $crac = "/path/to/crac";
our $cracIndex = "/path/to/crac-index";

#Path to DuplicationDetector
our $duplicationDetector = "/path/to/duplicationDetector";

#Path to BEDtools
our $bedtools = "/path/to/bedtools";

#Path to Abyss/TransAbyss
our $abyss = "/path/to/abyss";
our $transAbyss = "/path/to/transAbyss";

#Path to breakDancer
our $bam2cfg = "/path/to/bam2cfg.pl";
our $breakDancer = "/path/to/breakDancer";

#Path to pindel
our $pindel = "/path/to/pindel";


# path for plink
our $plink="/path/to/plink";

# path to sNMF
our $snmfbin = "/path/to/snmf";

# path to readseq
our $readseqjar = "/path/to//readseq.jar";

#path to FastME
our $fastme= "/path/to/fastme";

#path to fastq-stats
our $fastqStats= "/path/to/fastq-stats";
1;

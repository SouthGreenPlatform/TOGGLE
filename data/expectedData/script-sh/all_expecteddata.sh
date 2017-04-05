#!/bin/sh

echo "
####################################################
######### Suppression du jeu de donnees existant
####################################################
";

cd ..
rm -rf RC3_1_fastqc RC3_2_fastqc RC1_1.ILLUMINA.fastq
rm RC3_1_fastqc.zip RC3_2_fastqc.zip RC3_1.SANGER.fastq RC3_2.SANGER.fastq RC3_1.CUTADAPT.fastq RC3_2.CUTADAPT.fastq RC3_1.REPAIRING.fastq RC3_2.REPAIRING.fastq RC3.REPAIRING.fastq
rm Reference.fasta.amb Reference.fasta.ann Reference.fasta.bwt Reference.fasta.pac Reference.fasta.sa Reference.dict Reference.fasta.fai
rm RC3_1.BWAALN.sai RC3_2.BWAALN.sai RC3.BWAALN.sai RC3.BWASAMPE.sam RC3.BWASAMSE.sam RC3.PICARDTOOLSSORT.bam RC3Single.PICARDTOOLSSORT.bam
rm RC3.PICARDTOOLSSORT.bam.bai RC3Single.PICARDTOOLSSORT.bam.bai RC3.SAMTOOLSVIEW.bam RC3Single.SAMTOOLSVIEW.bam RC3.SAMTOOLSVIEW.bam.bai RC3Single.SAMTOOLSVIEW.bam.bai
rm RC3.GATKREALIGNERTARGETCREATOR.intervals RC3Single.GATKREALIGNERTARGETCREATOR.intervals
rm RC3.GATKINDELREALIGNER.bam RC3.GATKINDELREALIGNER.bai RC3Single.GATKINDELREALIGNER.bam RC3Single.GATKINDELREALIGNER.bai RC3.GATKINDELREALIGNER.bam.bai RC3Single.GATKINDELREALIGNER.bam.bai
rm RC3.PICARDTOOLSMARKDUPLICATES.bam RC3.PICARDTOOLSMARKDUPLICATES.bamDuplicates RC3Single.PICARDTOOLSMARKDUPLICATES.bam RC3Single.PICARDTOOLSMARKDUPLICATES.bamDuplicates
rm GATKHAPLOTYPECALLER.vcf  GATKSELECTVARIANTS.vcf  GATKVARIANTFILTRATION.vcf
rm GATKHAPLOTYPECALLER.vcf.idx  GATKSELECTVARIANTS.vcf.idx  GATKVARIANTFILTRATION.vcf.idx
rm RC3.SAMTOOLSFLAGSTAT.txt
rm RC3_1.FASTXTRIMMER.fastq RC3_2.FASTXTRIMMER.fastq
rm -rf tophat
rm referenceRNASeq*bt2 referenceRNASeq*ebwt

cd script-sh

echo "
####################################################
######## Lancement des scripts de creation des données test
####################################################
";


echo "
################# fastqc
";
sh fastqc.sh;


echo "
################# seqret sanger -> illumina
";
seqret fastq-sanger::../RC1_1.fastq fastq-illumina::../RC1_1.ILLUMINA.fastq
#seqret fastq-illumina::../RC3_1.fastq fastq-sanger::../RC3_1.SANGER.fastq
#seqret fastq-illumina::../RC3_2.fastq fastq-sanger::../RC3_2.SANGER.fastq


echo "
################# cutadapt
";
sh cutadapt.sh;



echo "
################# pairing
";
perl pairing.pl ../RC3_1.CUTADAPT.fastq ../RC3_2.CUTADAPT.fastq ../RC3_1.REPAIRING.fastq ../RC3_2.REPAIRING.fastq ../RC3.REPAIRING.fastq



echo "
################# bwa
";
sh bwa_index.sh;
sh bwa_aln.sh;
sh bwa_sampe.sh;
sh bwa_samse.sh;

echo "
################ picardtools
";
sh picardtools_createSequenceDictionnary.sh;
sh picardtools_sort.sh;

echo "
################ samtools
";
sh samtools_index1.sh;
sh samtools_view.sh;
sh samtools_index2.sh;
sh samtools_faidx.sh;
sh samtools_flagstat.sh;

echo "
################ gatk
";
sh gatk_realignerTargetCreator.sh;
sh gatk_indelRealigner.sh;

echo "
################ samtools
";
sh samtools_index3.sh;

echo "
################ picardtools
";
sh picardtools_markDuplicates.sh;

echo "
################ samtools
";
sh samtools_index4.sh;

echo "
################ gatk
";
sh gatk_haplotypeCaller.sh;
sh gatk_variantFiltration.sh;
sh gatk_selectVariants.sh

echo "
################ fastxTrimmer
";
sh fastxTrimmer.sh;
sh bowtie-build.sh;
sh bowtie2-build.sh;
sh tophat.sh;
sh htseqcount.sh;

exit;

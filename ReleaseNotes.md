
# Versions Notes

## Release 0.3, the 2th of December, 2015

### Modules

- Adding the picardTools AddOrReplaceReadGroup CleanSam ValidateSamFile SamFormatConverter
- Adding the GATK UnifiedGenotyper BaseRecalibrator PrintReads
- Adding the samTools sortSam flagstats depth idxstats
- Adding scheduler.pm module for scheduling system
- Adding the SGE and Slurm schedulers
- Adding cleaning and compression of given steps


### Functions

- *On the fly* creation of pipeline (see MANUAL.md for a detailled explanation)
- Use of **Graphviz** to generate a visual output of the pipeline structure
- Adding the automatic scheduling-aware launching system - For SGE and Slurm
- Adding a cleaning and a compression system to remove/compress chosen intermediate data in order to save hard drive space (see MANUAL.md)
- Use of initial compressed gzipped files
- Possibility of starting from FASTQ, SAM/BAM or VCF (gzipped or not)
- Possibility of working in relative path
- Providing an output folder in order to not modify the input data (working with symbolic links)
- Single creation of index/dictionary for reference, once per pipeline and no more once per sample, only if they do not exist.

## Release 0.2.1, the 15th of May, 2015

### Addition to previous version
* Adding correct SGE version for globalAnalysis.pl
* Adding RNAseq pipeline v1

### Corrections of bugs
* Correction of the ls bug for only a single pair
* Correction of hard coded path
* Validation of current softwares versions
* Modification of the FASTQ format control
* Various minor bugs adjustments

### Software versions tested
* bwa
* Picard-tools 1.124 and Higher
* SAMtools
* GATK
* FastQC
* CutAdapt
* Java


## Release 0.2, 14st of March, 2015

Second version release

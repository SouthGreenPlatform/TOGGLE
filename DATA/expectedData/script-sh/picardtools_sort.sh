/usr/local/java/latest/bin/java -Xmx12g -jar /usr/local/picard-tools-1.130/picard.jar SortSam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT INPUT=../RC3.BWASAMPE.sam OUTPUT=../RC3.PICARDTOOLSSORT.bam
/usr/local/java/latest/bin/java -Xmx12g -jar /usr/local/picard-tools-1.130/picard.jar SortSam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT INPUT=../RC3.BWASAMSE.sam OUTPUT=../RC3Single.PICARDTOOLSSORT.bam
# Extract HEADER
#
#samtools view -H RC3.PICARDTOOLSSORT.bam > headerFile.txt
#echo "@CO	addInfoHeaderTest" | cat - >> headerFile.txt 
#samtools reheader headerFile.txt RC3.PICARDTOOLSSORT.bam > RC3.ADDHEADER.bam
#rm headerFile.txt

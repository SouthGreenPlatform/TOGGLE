/usr/local/java/latest/bin/java -Xmx12g -jar /home/sabotf/sources/picard-tools/SortSam.jar SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT INPUT=../RC3.BWASAMPE.sam OUTPUT=../RC3.PICARDTOOLSSORT.bam
/usr/local/java/latest/bin/java -Xmx12g -jar /home/sabotf/sources/picard-tools/SortSam.jar SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT INPUT=../RC3.BWASAMSE.sam OUTPUT=../RC3Single.PICARDTOOLSSORT.bam
# Extract HEADER
#
#samtools view -H RC3.PICARDTOOLSSORT.bam > headerFile.txt
#echo "@CO	addInfoHeaderTest" | cat - >> headerFile.txt 
#samtools reheader headerFile.txt RC3.PICARDTOOLSSORT.bam > RC3.ADDHEADER.bam
#rm headerFile.txt

#!/bin/sh
#read -p "enter the RNA fastq file:" x


#echo "$INPUT_DIR1"
samtools mpileup -f /1_disk/public_resources/hg19.fa ./$1/accepted_hits.bam | cut -f 1,2,3,4 > ./$1/coverage.pileup

perl coverage.pl $1 100000

perl sort_bed.pl $1 





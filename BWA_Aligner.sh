#!/bin/bash

#“BWA_Aligner.sh” generates aligned/mapped and sorted “sampleX_mapped.bam” file for ‘paired end’ samples

for i in `ls *.fq.gz | cut -f1 -d'.'`
do
        #echo $i
        #echo '@RG\tID:'${i}'\tPL:illumina\tPU:HISEQ200B\tLB:'${i}'\tSM:'${i}'\tDS:'${i}':exome'

bwa mem -R  '@RG\tID:'${i}'\tPL:illumina\tPU:HISEQ200B\tLB:'${i}'\tSM:'${i}'\tDS:'${i}':exome' GenomeRef  ${i}.r1.fq.gz ${i}.r2.fq.gz |samblaster -M | samtools fixmate - - | samtools sort -O bam 
-o "${i}_mapped-bwa.bam"
done

# Variant calling performed on RNA-seq from tumor samples offers a valuable addition to WES DNA-seq to identify new variants,
# which likely comes from higher sequencing coverage of significantly expressed genes.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/
# /usr/bin/GATK/GenomeAnalysisTK.jar

SampleID=$1
Working_Dir=$2

# First adpater trim the reads using CUTADAPT_adapter_trimmer.sh

# Remove rRNA and tRNA contamination in step below
#/usr/bin/bowtie2-2.2.4-linux-x86_64/bowtie2-2.2.4/bowtie2 -x /usr/bin/RESOURCES/Genome/HUMAN/hg19.contamination.fa -p 30 -N 1 -1 R1.AT.fastq -2 R2.AT.fastq --un-conc $wd/$SampleName[$i]\_uncontaminated.fq -S $wd/$SampleName[$i]\.sam

#Run STAR aligner to generate SAM file
#/usr/bin/STAR-STAR_2.4.1d/bin/Linux_x86_64/STAR --runThreadN 30 --runMode alignReads --genomeDir /usr/bin/RESOURCES/Genome/HUMAN/STAR_Genome_hg19  --genomeFastaFiles /usr/bin/RESOURCES/Genome/HUMAN/STAR_Genome_hg19/hg19.fa --sjdbGTFfile /usr/bin/RESOURCES/Genome/HUMAN/STAR_Genome_hg19/NEW_Homo_sapiens.GRCh37.75.gtf --outSAMunmapped Within --outSAMattrRGline ID:$SampleName[$i] PL:illumina PU:$SampleName[$i]\_uncontaminated_1_2 SM:$SampleName[$i] --outFileNamePrefix $wd/$SampleName[$i]\. --sjdbOverhang 100 --limitGenomeGenerateRAM 61000000000 --readFilesIn $wd/$SampleName[$i]\_uncontaminated.1.fq $wd/$SampleName[$i]\_uncontaminated.2.fq

#Running samtools to generate sorted.bam file (next 3 steps)
#/usr/bin/samtools-0.1.19/samtools view -@ 15 -bSh $wd/$SampleName[$i]\.Aligned.out.sam > $wd/$SampleName[$i]\.Aligned.out.bam

#/usr/bin/samtools-0.1.19/samtools sort -@ 10 -m 3G $wd/$SampleName[$i]\.Aligned.out.bam $wd/$SampleName[$i]\.Aligned.out.sorted

#/usr/bin/samtools-0.1.19/samtools index  $wd/$SampleName[$i]\.Aligned.out.sorted.bam

# mark duplicates

#java -Xmx25g -jar /usr/bin/picard-tools-1.140/picard.jar MarkDuplicates I=$Working_Dir/$SampleID\.Aligned.out.sorted.bam O=$Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$Working_Dir/$SampleID\.metrics 

#Split'N'Trim and reassign mapping qualities

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.bam -o $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

# Indel Realignment
java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.bam -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa -o $Working_Dir/$SampleID\.indelrealigner.intervals -known /usr/bin/GATK/Mills.Gold.vcf

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T IndelRealigner -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.bam -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa -targetIntervals $Working_Dir/$SampleID\.indelrealigner.intervals -o $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.bam -known /usr/bin/GATK/Mills.Gold.vcf

/usr/bin/samtools-0.1.19/samtools index $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.bam

#Base Recalibration

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.bam --knownSites /usr/bin/RESOURCES/DBSNP/dbSNP135/dbsnp_135.hg19.vcf -o $Working_Dir/$SampleID\.grp -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa --disable_indel_quals

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T PrintReads -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.bam -BQSR $Working_Dir/$SampleID\.grp -o $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.recal.bam -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa

samtools index $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.recal.bam

#Variant calling

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa -I $Working_Dir/$SampleID\.Aligned.out.sorted.dedupped.split.realigned.recal.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $Working_Dir/$SampleID\.vcf

java -Xmx25g -jar /usr/bin/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /usr/bin/RESOURCES/HG19_GENOME/hg19.fa -V $Working_Dir/$SampleID\.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $Working_Dir/$SampleID\.Filtered.vcf


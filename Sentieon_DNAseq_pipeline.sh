#!/bin/bash

# Sentieon can be run from any server]
# Input for this is adapter trimmed fastq R1 and R2 reads 


if test -n "$1" -a "$2"
then

	# Setting variables:
	THREADS="$2"
	R1=$(basename "$1")
	RN=$(basename "$1" _R1.fastq)
	R2="$RN""_R2.fastq"
	OUT=$(basename "$1" _at_R1.fastq)
	DIRPATH=$(dirname "$1")
	
	if [ ! -d "$OUT" ]
	then

		mkdir "$OUT"
		printf "$OUT\nSTART: "; date
		SECONDS=0

		# Sentieon Pipeline:
		# 1. BWA mapping
		/opt/sentieon-genomics-201611/bin/bwa mem -M -R "@RG\tID:$OUT\tSM:$OUT\tPL:illumina" -t $THREADS -K 10000000 .../RESOURCES/External_data/hg19/hg19.fa "$DIRPATH/$R1" "$DIRPATH/$R2" | /opt/sentieon-genomics-201611/bin/sentieon util sort -o  "$OUT/$OUT.sorted.bam"  -t $THREADS --sam2bam -i -  

		# 2. Calculate metrics
		/opt/sentieon-genomics-201611/bin/sentieon driver -t $THREADS -r .../RESOURCES/External_data/hg19/hg19.fa  -i  "$OUT/$OUT.sorted.bam"  --algo GCBias --summary  "$OUT/$OUT.gc_summary.txt"  "$OUT/$OUT.gc_metrics.txt" --algo MeanQualityByCycle  "$OUT/$OUT.mq_metrics.txt"  --algo QualDistribution  "$OUT/$OUT.qd_metrics.txt" --algo InsertSizeMetricAlgo  "$OUT/$OUT.is_metrics.txt"   --algo  AlignmentStat "$OUT/$OUT.aln_metrics.txt"

		# 3. Remove duplicates
		/opt/sentieon-genomics-201611/bin/sentieon driver -t $THREADS -i "$OUT/$OUT.sorted.bam"   --algo LocusCollector  --fun score_info "$OUT/$OUT.score.txt"
		/opt/sentieon-genomics-201611/bin/sentieon driver -t $THREADS -i "$OUT/$OUT.sorted.bam"  --algo Dedup --rmdup --score_info "$OUT/$OUT.score.txt" --metrics "$OUT/$OUT.dedup_metrics.txt" "$OUT/$OUT.deduped.bam"

		# 4. InDel realignment
		/opt/sentieon-genomics-201611/bin/sentieon driver -r  .../RESOURCES/External_data/hg19/hg19.fa  -t $THREADS  -i "$OUT/$OUT.deduped.bam" --algo Realigner -k .../RESOURCES/External_data/hg19/Mills_and_1000G_gold_standard_indels/all_chr.vcf  "$OUT/$OUT.realigned.bam"

		# 5. Base quality score recalibration (BQSR)
		/opt/sentieon-genomics-201611/bin/sentieon driver -r .../RESOURCES/External_data/hg19/hg19.fa  -t $THREADS -i  "$OUT/$OUT.realigned.bam"  --algo QualCal -k .../RESOURCES/External_data/hg19/dbSNP149/All_20161121.mod.vcf  "$OUT/$OUT.recal_data.table"
		/opt/sentieon-genomics-201611/bin/sentieon driver -r .../RESOURCES/External_data/hg19/hg19.fa   -t $THREADS -i "$OUT/$OUT.realigned.bam" -q "$OUT/$OUT.recal_data.table" --algo QualCal -k .../RESOURCES/External_data/hg19/dbSNP149/All_20161121.mod.vcf   "$OUT/$OUT.recal_data.table.post" 
		/opt/sentieon-genomics-201611/bin/sentieon driver -t $THREADS --algo QualCal --plot --before "$OUT/$OUT.recal_data.table"  --after "$OUT/$OUT.recal_data.table.post" "$OUT/$OUT.recal_data.csv"
		/opt/sentieon-genomics-201611/bin/sentieon plot bqsr -o "$OUT/$OUT.bqsrreport.pdf" "$OUT/$OUT.recal_data.csv"

		# 6. Variant calling
		/opt/sentieon-genomics-201611/bin/sentieon driver -r .../RESOURCES/External_data/hg19/hg19.fa   -t $THREADS -i  "$OUT/$OUT.realigned.bam"   -q "$OUT/$OUT.recal_data.table"   --algo Haplotyper -d .../RESOURCES/External_data/hg19/dbSNP149/All_20161121.mod.vcf --emit_mode gvcf  "$OUT/$OUT.HC.gvcf"

		printf "\nEND: "; date; 
		MIN=$(($SECONDS / 60))
		HRS=$(($MIN / 60))
		DAYS=$(($HRS / 24))
		if [ $DAYS == 1 ]; then d_str="day"; else d_str="days"; fi
		if [ $(($HRS % 24)) == 1 ]; then h_str="hr"; else h_str="hrs"; fi
		printf "RUNTIME: $DAYS $d_str, $(($HRS % 24)) $h_str, $(($MIN % 60)) min & $(($SECONDS % 60)) sec\n\n"

	fi

else
SCRIPT=$(basename "$0")
printf "\nUSAGE: bash $SCRIPT <path/R1.fastq> <INT_#threads>\n\n"
fi

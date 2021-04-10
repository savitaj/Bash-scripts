#!/bin/bash
#MiXCR is a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data
#This program takes R1 and R2 fastq files of TCRseq as input and runs MiXCR to align, assemble and export the TCR clones
#Non-overlaping paired-end reads are always processed in VThenJ mode. JThenV can be used for short reads (~100bp) with full (or nearly full) J gene coverage.


R1=$1
R2=$2
T="${1%_R1.AT.fastq.gz}"	#removing  prefix "_R1.AT.fastq"

P="${T%/*}/"		#Path where input files are available
S="${T##*/}"		#SampleName
O_P=$P$S"/"     #Output path
mkdir -p $P"/"$S	#make sample name directory in the path 

echo "$T	$P	$S	$O_P"
echo " MIXCR ANALYSIS FOR SAMPLE $S STARTS HERE "

mixcr align -s hs -OvParameters.geneFeatureToAlign=VTranscript -OvjAlignmentOrder=JThenV -OsaveOriginalReads=true -r $O_P$S"_alignreports.txt" $R1 $R2 $O_P$S"_alignments.vdjca"

mixcr assemble -r  $O_P$S"_assemble_report.txt" --write-alignments $O_P$S"_alignments.vdjca" $O_P$S"_alignments.clna"

mixcr exportClones -o -t  $O_P$S"_alignments.clna" $O_P$S"_clones.txt"

mixcr exportAlignments $O_P$S"_alignments.clna" $O_P$S"_alignments.txt"

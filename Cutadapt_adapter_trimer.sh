#For Takara, use '.gz' file for R1 & R2; and need to remove first 30 nucleotides from each read in R2.
#Usage: sh CUTADAPT.sh $FASTQR1 $FASTQR2 $SAMPLE_NAME $WORKING_DIR

FASTQR1=$1
FASTQR2=$2
SAMPLE_NAME=$3
WORKING_DIR=$4

#Adapter trimming for R1 reads
/usr/bin/cutadapt -j 30 -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o $WORKING_DIR/$SAMPLE_NAME"_R1_1.fastq.gz" $FASTQR1 
/usr/bin/cutadapt -j 30 -a AGATCGGAAGAGCGTCGTGT -o $WORKING_DIR/$SAMPLE_NAME"_R1_2.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_1.fastq.gz"
/usr/bin/cutadapt -j 30 -a AGATCGGAAGAGCACACGTC -o $WORKING_DIR/$SAMPLE_NAME"_R1_3.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_2.fastq.gz"
/usr/bin/cutadapt -j 30 -a ATCTCGTATGCCGTCTTCTGCTTG -o $WORKING_DIR/$SAMPLE_NAME"_R1_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_3.fastq.gz"

#Adapter trimming for R2 reads
/usr/bin/cutadapt -j 30 -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o $WORKING_DIR/$SAMPLE_NAME"_R2_1.fastq.gz" $FASTQR2
/usr/bin/cutadapt -j 30 -a AGATCGGAAGAGCGTCGTGT -o $WORKING_DIR/$SAMPLE_NAME"_R2_2.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_1.fastq.gz"
/usr/bin/cutadapt -j 30 -a AGATCGGAAGAGCACACGTC -o $WORKING_DIR/$SAMPLE_NAME"_R2_3.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_2.fastq.gz"
/usr/bin/cutadapt -j 30 -a ATCTCGTATGCCGTCTTCTGCTTG -o $WORKING_DIR/$SAMPLE_NAME"_R2_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_3.fastq.gz"

#check paired end reads and min read length 30
perl /usr/scripts/CheckAfterAdapterTrimming.new.pl $WORKING_DIR/$SAMPLE_NAME"_R1_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_AT.cutadapt.fastq" $WORKING_DIR/$SAMPLE_NAME"_R2_AT.cutadapt.fastq"

/usr/scripts/fastq-mcf -x 10 -q 20 -l 30 -S -k 0 -o $WORKING_DIR/$SAMPLE_NAME"_R1.AT.fastq" -o $WORKING_DIR/$SAMPLE_NAME"_R2.AT.fastq" /BXRX_ANALYSIS/AIRR_TCR/AdapterTrim/Adapter_seq.fa $WORKING_DIR/$SAMPLE_NAME"_R1_AT.cutadapt.fastq" $WORKING_DIR/$SAMPLE_NAME"_R2_AT.cutadapt.fastq"

rm -rf $WORKING_DIR/$SAMPLE_NAME"_R1_1.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_2.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_3.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_1.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_2.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_3.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R2_4.fastq.gz" $WORKING_DIR/$SAMPLE_NAME"_R1_AT.cutadapt.fastq" $WORKING_DIR/$SAMPLE_NAME"_R2_AT.cutadapt.fastq"

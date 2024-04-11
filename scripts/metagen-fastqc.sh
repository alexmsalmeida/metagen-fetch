#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

Clean and remove host DNA sequences from FASTQ files.

OPTIONS:
   -h      Show help message
   -t      Number of threads (recommended: 4)
   -c      Referece genome for host decontamination (default: human_hg38)
   -f      Forward or single-end fastq file (.fastq or .fastq.gz) [REQUIRED]
   -r	   Reverse fastq file (.fastq or .fastq.gz) [OPTIONAL]
EOF
}

THREADS=
REF=
FASTQ_R1=
FASTQ_R2=

while getopts “ht:c:r:f:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         t)
             THREADS=${OPTARG}
             ;;
         c)
             REF=${OPTARG}
             ;;
         f)
             FASTQ_R1=$OPTARG
             ;;
	 r)
             FASTQ_R2=$OPTARG
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

timestamp() {
  date +"%H:%M:%S"
}

echo "$(timestamp) [ fastq clean-up ] Parsing command-line"
# check if all required arguments are supplied
if [[ -z ${FASTQ_R1} ]]
then
     echo "ERROR : Please supply input FASTQ files"
     usage
     exit 1
fi

if [[ ! -z ${FASTQ_R2} ]]
then
    N_FWD=$(zcat ${FASTQ_R1} | wc -l)
    N_REV=$(zcat ${FASTQ_R2} | wc -l)
    if [[ ${N_FWD} -eq ${N_REV} ]]
    then
        echo "Number of reads in FWD and REV match, proceeding"
    else
        echo "Number of reads in FWD and REV do not match, there is likely an issue with the files. Exiting..."
        exit 1
    fi
fi


if [[ -z ${REF} ]]
then
    REF="/rfs/project/rfs-42V0sRc5Rk8/databases/bwa/hg38.fa"
fi

if [ ${THREADS} -eq 1 ]
then
    THREADS_SAM=1
else
    THREADS_SAM=$((${THREADS}-1))
fi

if [[ ! -z ${FASTQ_R2} ]]
then
        echo "$(timestamp) [ fastq clean-up ] Cleaning FASTQ files"
        trim_galore --paired ${FASTQ_R1} ${FASTQ_R2} -o $(dirname ${FASTQ_R1}) --j ${THREADS}
        name=${FASTQ_R1%%_1.fastq*}
        echo "$(timestamp) [ fastq clean-up ] Mapping files to host genome: ${REF}"
        bwa mem -M -t ${THREADS} ${REF} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz | samtools view -@ ${THREADS_SAM} -f 12 -F 256 -uS - -o ${name}_both_unmapped.bam 
	samtools sort -@ ${THREADS_SAM} -n ${name}_both_unmapped.bam -o ${name}_both_unmapped_sorted.bam
	bedtools bamtofastq -i ${name}_both_unmapped_sorted.bam -fq ${name}_clean_1.fastq -fq2 ${name}_clean_2.fastq
	echo "$(timestamp) [ fastq clean-up ] Compressing output files"
	pigz ${name}_clean_1.fastq
	pigz ${name}_clean_2.fastq
        echo "$(timestamp) [ fastq clean-up ] Cleaning tmp files"
        rm -rf ${name}_both_unmapped.bam ${name}_both_unmapped_sorted.bam ${name}_*_trimmed.fq.gz ${name}_*_val_*.fq.gz ${name}_*.fastq.gz_*txt
else
        echo "$(timestamp) [ fastq clean-up ] Cleaning FASTQ files"
        trim_galore ${FASTQ_R1} -o $(dirname ${FASTQ_R1}) --j ${THREADS}
        name=${FASTQ_R1%%.fastq*}
        echo "$(timestamp) [ fastq clean-up ] Mapping files to host genome: ${REF}"
        bwa mem -M -t ${THREADS} ${REF} ${name}_trimmed.fq.gz | samtools view -@ ${THREADS_SAM} -f 4 -F 256 -uS - -o ${name}_unmapped.bam
        samtools sort -@ ${THREADS_SAM} -n ${name}_unmapped.bam -o ${name}_unmapped_sorted.bam
        bedtools bamtofastq -i ${name}_unmapped_sorted.bam -fq ${name}_clean.fastq
	echo "$(timestamp) [ fastq clean-up ] Compressing output file"
	pigz ${name}_clean.fastq
        echo "$(timestamp) [ fastq clean-up ] Cleaning tmp files"
        rm -rf ${name}_unmapped.bam ${name}_unmapped_sorted.bam ${name}_trimmed.fq.gz ${FASTQ_R1}_*txt
fi

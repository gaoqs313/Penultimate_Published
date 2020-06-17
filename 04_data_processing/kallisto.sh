#!/bin/bash

SAMPLE=$1
CPU=4
DIR=/scratch_space/qgao/Penultimate/2.data/mouse

mkdir -p $SAMPLE
cd $SAMPLE
zcat $DIR/${SAMPLE}_1.fastq.gz > ${SAMPLE}_1.fastq
zcat $DIR/${SAMPLE}_2.fastq.gz > ${SAMPLE}_2.fastq

#source activate common_tools
INDEX=/scratch_space/qgao/Penultimate/1.reference/mouse/Kallisto/Mus_musculus.GRCm38.98.transcript.kallisto.idx
kallisto quant -i $INDEX -t $CPU -o $SAMPLE ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq

## transcript to gene expression
perl ../kallisto_transcript2gene.pl $SAMPLE/abundance.tsv $SAMPLE.tsv

rm ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq

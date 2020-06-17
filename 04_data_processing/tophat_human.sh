#!/bin/bash

SAMPLE=$1
LEFT=$2
RIGHT=$3

module load python/2.7.15-rhel7
module load tophat/2.1.1

REF=/scratch_space/qgao/Penultimate/1.reference/human/Bowtie2/hg38
GTF=/scratch_space/qgao/Penultimate/1.reference/human/GTF/Homo_sapiens.GRCh38.98.chr.gtf
CPU=4 

tophat -o $SAMPLE -p $CPU -G $GTF $REF $LEFT $RIGHT
samtools index $SAMPLE/accepted_hits.bam

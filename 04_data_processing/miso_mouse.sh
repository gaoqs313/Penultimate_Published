#!/bin/bash

SAMPLE=$1
BAM=$2
SPECIES=$3
CPU=1

REF="/scratch_space/qgao/Penultimate/1.reference/$SPECIES/B_MISO/indexed_SE_events/"

LEN=`samtools view $BAM | head -n1 | cut -f 10 | awk '{print length}'`

MEAN=`head -n1 /scratch_space/qgao/Penultimate/4.insert_len/$SPECIES/$SAMPLE/accepted_hits.bam.insert_len | cut -f 1 -d ',' | cut -f 2 -d '='`
SD=`head -n1 /scratch_space/qgao/Penultimate/4.insert_len/$SPECIES/$SAMPLE/accepted_hits.bam.insert_len | cut -f 2 -d ',' | cut -f 2 -d '='`

source activate miso
miso --run $REF $BAM --output-dir $SAMPLE/ --read-len $LEN --paired-end $MEAN $SD --event-type SE -p $CPU

summarize_miso --summarize-samples $SAMPLE/ $SAMPLE/

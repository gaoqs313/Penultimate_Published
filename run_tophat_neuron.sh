#!/bin/bash

DIR=/scratch_space/qgao/Penultimate/2.data/neuron

cat /scratch_space/qgao/Penultimate/2.data/neuron_SRR_Acc_List.txt | while read sample
do
	lsf_prio 4 4 $sample bash tophat_mouse.sh $sample $DIR/$sample\_1.fastq.gz $DIR/$sample\_2.fastq.gz
done

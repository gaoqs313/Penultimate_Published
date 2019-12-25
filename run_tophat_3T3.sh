#!/bin/bash

DIR=/scratch_space/qgao/Penultimate/2.data/3T3
for sample in ERR498282 ERR498284
do
	lsf_normal 4 4 $sample bash tophat_mouse.sh $sample $DIR/$sample\_1.fastq.gz $DIR/$sample\_2.fastq.gz
done

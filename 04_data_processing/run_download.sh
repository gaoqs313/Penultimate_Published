
## mouse
cat mouse_SRR_Acc_List.txt | while read acc
do
	lsf_normal 2 1 $acc fastq-dump -O mouse --split-3 --gzip $acc
done

cat neuron_SRR_Acc_List.txt | while read acc
do
        lsf_normal 2 1 $acc fastq-dump -O neuron --split-3 --gzip $acc
done

for acc in ERR498282 ERR498284
do
	lsf_normal 2 1 $acc fastq-dump -O 3T3 --split-3 --gzip $acc
done

## human
grep PAIRED human_SraRunTable.txt | cut -f 1 -d ',' > paired_human_SRR_Acc_List.txt
cat paired_human_SRR_Acc_List.txt | while read acc
do
	lsf_normal 2 1 $acc fastq-dump -O human --split-3 --gzip $acc
done



Rscript neuron_stat.R

## get fasta for differential events
perl bed6_to_fasta.pl

## get fasta for background
perl bed6_to_fasta2.pl

findMotifs.pl neuron_diff.fa fasta motifResultsRBP -fastaBg neuron_same.fa -rna -len 6,7,8

## change file name /home/qgao/Tools/Homer/data/knownTFs

## RBP
examine_row.sh /scratch_space/qgao/Penultimate/11.revision/E_Mouse_GO/Neuron/Gene/alltissue.all.tsv 44449 | cut -f 2-3 > Pcbp3.tsv
examine_row.sh /scratch_space/qgao/Penultimate/11.revision/E_Mouse_GO/Neuron/Gene/alltissue.all.tsv 39461 | cut -f 2-3 > Mbnl1.tsv
paste Pcbp3.tsv Mbnl1.tsv  | cut -f 1,2,4 > rbp_expr.tsv

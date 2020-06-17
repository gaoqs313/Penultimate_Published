
## mouse to human
sed 's/^/chr/g' /scratch_space/qgao/Penultimate/1.reference/mouse/F_MISO/Eligibile/mm10.SE.intron.bed | awk -v OFS='\t' '{print $1,$2,$3,$4"__"$1"__"$2"__"$3,$5,$6}' > chr.mm10.SE.intron.bed
liftOver chr.mm10.SE.intron.bed ../mm10ToHg38.over.chain.gz mm10.SE.intron.lift2human.bed mm10.unmapped -minMatch=0.1

## human to mouse
sed 's/^/chr/g' /scratch_space/qgao/Penultimate/1.reference/human/B_MISO/hg38.SE.intron.bed | awk -v OFS='\t' '{print $1,$2,$3,$4"__"$1"__"$2"__"$3,$5,$6}' > chr.hg38.SE.intron.bed
liftOver chr.hg38.SE.intron.bed ../hg38ToMm10.over.chain.gz hg38.SE.intron.lift2mouse.bed hg38.unmapped -minMatch=0.1

## link
perl get_one_to_one_orthlog.pl

## get both
awk '$4==1 && $8==1' liftOver.results | cut -f 7 > both_support_id.txt
 
## format
perl format.pl

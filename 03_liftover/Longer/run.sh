
## mouse to human
sed 's/^/chr/g' /research/rgs01/scratch/users/qgao/Penultimate/1.reference/mouse/TEST/mm10.SE.shorter.intron.bed | awk -v OFS='\t' '{print $1,$2,$3,$4"__"$1"__"$2"__"$3,$5,$6}' > chr.mm10.SE.intron.bed
liftOver chr.mm10.SE.intron.bed ../mm10ToHg38.over.chain.gz mm10.SE.intron.lift2human.bed mm10.unmapped -minMatch=0.1

## human to mouse
sed 's/^/chr/g' /research/rgs01/scratch/users/qgao/Penultimate/1.reference/human/TEST/hg38.SE.shorter.intron.bed | awk -v OFS='\t' '{print $1,$2,$3,$4"__"$1"__"$2"__"$3,$5,$6}' > chr.hg38.SE.intron.bed
liftOver chr.hg38.SE.intron.bed ../hg38ToMm10.over.chain.gz hg38.SE.intron.lift2mouse.bed hg38.unmapped -minMatch=0.1

## link
perl get_one_to_one_orthlog.pl


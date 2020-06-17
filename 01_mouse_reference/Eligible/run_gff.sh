
## convert gtf to gff
cat Mus_musculus.GRCm38.98.chr.gtf | perl gtf_to_gff.pl | sed 's/^chr//g' > Mus_musculus.GRCm38.98.chr.gff

## get long exon
exon_utils --get-const-exons Mus_musculus.GRCm38.98.chr.gff --min-exon-size 1000 --output-dir exons/

## get gene id
awk '$3=="gene"' Mus_musculus.GRCm38.98.chr.gtf | cut -f 2,6 -d '"' | sed 's/"/\t/g' > Mus_musculus.GRCm38.98.geneid.txt


###########################################
############### Longer C termini ##########
###########################################

perl extract_eligible_transcript.pl /scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf /scratch_space/qgao/Penultimate/1.reference/mouse/Genome/mm10.fa

cat StopCodonInLastExon.final.Mus_musculus.GRCm38.98.chr.gtf StopCodonInPenultimateExon.final.Mus_musculus.GRCm38.98.chr.gtf > final.Mus_musculus.GRCm38.98.chr.gtf

cut -f 2,6 -d '"' final.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 2 | uniq | sort -u | wc -l
# 3233
cut -f 2,6 -d '"' final.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 1 | uniq | sort -u | wc -l
#2860

###########################################
############### Shorter C termini #########
###########################################

perl extract_shorter.pl /scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf /scratch_space/qgao/Penultimate/1.reference/mouse/Genome/mm10.fa

cut -f 2,6 -d '"' StopCodonInLastExon.final.shorter.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 2 | uniq | sort -u | wc -l
# 9791
cut -f 2,6 -d '"' StopCodonInLastExon.final.shorter.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 1 | uniq | sort -u | wc -l
# 8472

###########################################
##################### 3N ##################
###########################################

perl extract_3n.pl /scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf

cut -f 2,6 -d '"' StopCodonInLastExon.final.3N.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 2 | uniq | sort -u | wc -l
# 8512
cut -f 2,6 -d '"' StopCodonInLastExon.final.3N.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 1 | uniq | sort -u | wc -l
# 7103

###########################################
################# Internal ################
###########################################

perl extract_internal.pl /scratch_space/qgao/Penultimate/1.reference/mouse/GTF/Mus_musculus.GRCm38.98.chr.gtf

cat StopCodonInLastExon.final.internal.Mus_musculus.GRCm38.98.chr.gtf StopCodonInPenultimateExon.final.internal.Mus_musculus.GRCm38.98.chr.gtf > final.internal.Mus_musculus.GRCm38.98.chr.gtf 

cut -f 2,6 -d '"' final.internal.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 2 | uniq | sort -u | wc -l
# 23193
cut -f 2,6 -d '"' final.internal.Mus_musculus.GRCm38.98.chr.gtf | sed 's/"/\t/g' | cut -f 1 | uniq | sort -u | wc -l
# 16842





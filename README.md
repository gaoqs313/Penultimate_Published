# Penultimate

This is the repository for scripts used in our Genome Biology paper "Splicing-accessible coding 3’UTRs control protein stability and interaction networks"


==========================================================
## Prepare penultimate exon annotation

### Mouse Reference

Version: Mus_musculus.GRCm38.98
Script folder: 01_mouse_reference

Extract eligible transcript
```
perl extract_eligible_transcript.pl Mus_musculus.GRCm38.98.chr.gtf mm10.fa
cat StopCodonInLastExon.final.Mus_musculus.GRCm38.98.chr.gtf StopCodonInPenultimateExon.final.Mus_musculus.GRCm38.98.chr.gtf > final.Mus_musculus.GRCm38.98.chr.gtf
```

Calculate amino acid content
```
perl calculate_AA.pl
Rscript stat.R
```

Convert GTF to MISO GFF3 file
```
perl GenerateMisoSE.pl
```

Build MISO index
```
index_gff --index mm10.SE.gff3 indexed_SE_events/
```

### Human Reference

Version: Homo_sapiens.GRCh38.98
Script folder: 02_human_reference

Extract eligible transcript
```
perl extract_eligible_transcript.pl Homo_sapiens.GRCh38.98.chr.gtf hg38.fa
cat StopCodonInLastExon.final.Homo_sapiens.GRCh38.98.chr.gtf StopCodonInPenultimateExon.final.Homo_sapiens.GRCh38.98.chr.gtf > final.Homo_sapiens.GRCh38.98.chr.gtf
```

Calculate amino acid content
```
perl calculate_AA.pl
Rscript stat.R
```

Convert GTF to MISO GFF3 file
```
perl GenerateMisoSE.pl
```

Build MISO index
```
index_gff --index hg38.SE.gff3 indexed_SE_events/
```

### Liftover

```
bash run_liftOver.sh
```

==========================================================
## Process RNA-seq data

Script folder: 04_data_processing

### 1. download
```
bash run_download.sh
```

### 2. align using tophat2
```
bash tophat.sh
```

### 3. run miso
```
bash insert_len.sh
bash miso.sh
```
### 4. run kallisto
```
bash kallisto.sh
```

==========================================================
## Do statistics and plotting   



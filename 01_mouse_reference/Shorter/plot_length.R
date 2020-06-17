
x <- read.table("short_AA_content.tsv.length", row.names=1)

pdf("length_distribution.pdf", height = 5)
boxplot(x[,1], outline=F, las=1, ylab="Amino acid length")
dev.off()

sum(x[,1]>=20)
# 2862


df <- read.table("AA_content.tsv", stringsAsFactor=F)
colnames(df) <- c("id", "aa", "type", "perc")
df$perc <- 100*df$perc

fs <- df[df$type=='fs',]
fl <- df[df$type=='fl',]
cds <- df[df$type=='cds',]
utr <- df[df$type=='utr',]

dif1 <- fs
dif1$perc <- dif1$perc - fl$perc
dif2 <- utr
dif2$perc <- dif2$perc - cds$perc

data <- rbind(dif1, dif2)

dif1s <- split(dif1, dif1$aa)
dif1sm <- lapply(dif1s, function(x) mean(x$perc))
dif1smb <- do.call(cbind, dif1sm)

data$aa <- factor(data$aa, levels=colnames(dif1smb)[order(dif1smb)], ordered=T)

library(ggplot2)

pdf("Human_AAcontent_comparison.pdf", height=5)
p <- ggplot(data,aes(x = aa, y = perc, fill = type)) +
  geom_boxplot(outlier.shape=NA) +
  theme_bw() +
  scale_fill_manual(name=NULL, labels=c("Frameshift (fs) - Full_length (fl)", "UTR - CDS"), values=c("grey40", "white")) +
  xlab("") +
  ylim(-20, 20) +
  geom_hline(yintercept=0, color="grey70") +
  ylab("Amino acid content difference") +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position=c(0.28, 0.91)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()




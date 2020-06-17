library(ggplot2)
library(VennDiagram)

pdf("Venn.pdf")
pb <- venn.diagram(list(Mouse=1:687,Human=162:799), fill = c("#8DA0CB", "#A6D854"),  alpha = c(0.5, 0.5), cex = 2, filename = NULL)
grid.draw(pb)
dev.off()



df <- read.table("go.results", sep="\t")
colnames(df) <- c("GO", "Enrichment")

pdf("GO_overlap.pdf")
ggplot(df, aes(x=reorder(GO,Enrichment), y=Enrichment)) +
  geom_col(fill="grey30") +
  coord_flip() +
  labs(x="Pathway", y="Enrichment") +
  theme(legend.position = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=10, family = "Helvetica")) +
  theme(axis.title=element_text(size=12))
dev.off()



library(ggplot2)

df <- read.table("enrichr.tsv", sep="\t")
colnames(df) <- c("GO", "Enrichment")
df$Enrichment <- -log10(df$Enrichment)

pdf("GO_overlap2.pdf")
ggplot(df, aes(x=reorder(GO,Enrichment), y=Enrichment)) +
  geom_col(fill="grey30") +
  coord_flip() +
  labs(x="Pathway", y="-log10(adjusted p value)") +
  theme(legend.position = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=15, family = "Helvetica")) +
  theme(axis.text.y=element_text(size=10)) +
  scale_y_continuous(limits = c(0,6), breaks=c(0,3,6)) +
  theme(plot.margin=unit(c(0.5,1,0.5,0.5), "cm")) +
  theme(axis.title.x = element_text(hjust = 1))
dev.off()

pdf("UV_GSEA.pdf")
p <- ggplot(x, aes(reorder(NAME, NES), NES)) +
  geom_col(fill = "red3") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15, family = "Helvetica")) +
  theme(axis.text.y=element_text(size=8)) +
  scale_y_continuous(limits = c(0,2.5), breaks=c(0,1,2)) +
  theme(plot.margin=unit(c(0.5,1,0.5,0.5), "cm")) +
  theme(axis.title.x = element_text(hjust = 1))
print(p)
dev.off()


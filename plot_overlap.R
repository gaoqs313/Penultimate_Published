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

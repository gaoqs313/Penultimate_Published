
library(ggplot2)

bar_data <- as.data.frame(matrix(c(
"fp", "53.0",
"short", "51.1",
"long", "37.2",
"long10_Ulti", "32.0",
"long10", "29.4",
"long10_Penulti", "22.4"),byrow=T, ncol=2))
colnames(bar_data) <- c("Position", "Percentage")
bar_data$Percentage <- as.numeric(as.character(bar_data$Percentage))
bar_data$Position <- factor(bar_data$Position, levels=c("fp", "short", "long", "long10_Ulti", "long10", "long10_Penulti"), ordered=T)

pdf("conservation_comparison.pdf", width=11)
p <- ggplot(bar_data, aes(x=Position, y=Percentage)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("") +
  ylab("Percentage of exons") +
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white"))+ 
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position = "top")
print(p)
dev.off()



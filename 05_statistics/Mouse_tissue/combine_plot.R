
library(ggplot2)

bar_data <- as.data.frame(matrix(c(
 "PSI<90", "Internal_fs", "No", "36226", "41.6",
 "PSI<90", "Internal_fs", "Yes", "50867", "58.4",
 "PSI<90", "Internal_fp", "No", "25771", "42.9",
 "PSI<90", "Internal_fp", "Yes", "34236", "57.1",
 "PSI<90", "Penulti_fp", "No", "2763", "32.5",
 "PSI<90", "Penulti_fp", "Yes", "5749", "67.5",
 "PSI<90", "Penulti_short", "No", "2974", "30.4",
 "PSI<90", "Penulti_short", "Yes", "6817", "69.6",
 "PSI<90", "Penulti_long", "No", "1081", "33.4",
 "PSI<90", "Penulti_long", "Yes", "2152", "66.6",
 "PSI<85", "Internal_fs", "No", "59554", "68.4",
 "PSI<85", "Internal_fs", "Yes", "27539", "31.6",
 "PSI<85", "Internal_fp", "No", "41402", "69.0",
 "PSI<85", "Internal_fp", "Yes", "18605", "31.0",
 "PSI<85", "Penulti_fp", "No", "3741", "43.9",
 "PSI<85", "Penulti_fp", "Yes", "4771", "56.1",
 "PSI<85", "Penulti_short", "No", "4079", "41.7",
 "PSI<85", "Penulti_short", "Yes", "5712", "58.3",
 "PSI<85", "Penulti_long", "No", "1487", "46.0",
 "PSI<85", "Penulti_long", "Yes", "1746", "54.0",
 "PSI<80", "Internal_fs", "No", "70509", "81.0",
 "PSI<80", "Internal_fs", "Yes", "16584", "19.0",
 "PSI<80", "Internal_fp", "No", "48673", "81.1",
 "PSI<80", "Internal_fp", "Yes", "11334", "18.9",
 "PSI<80", "Penulti_fp", "No", "4402", "51.7",
 "PSI<80", "Penulti_fp", "Yes", "4110", "48.3",
 "PSI<80", "Penulti_short", "No", "4836", "49.4",
 "PSI<80", "Penulti_short", "Yes", "4955", "50.6",
 "PSI<80", "Penulti_long", "No", "1708", "52.8",
 "PSI<80", "Penulti_long", "Yes", "1525", "47.2"),byrow=T, ncol=5))
colnames(bar_data) <- c("PSI", "Position", "Support", "Number", "Percentage")
bar_data$Percentage <- as.numeric(as.character(bar_data$Percentage))
bar_data$Position <- factor(bar_data$Position, levels=c("Internal_fp", "Internal_fs", "Penulti_fp", "Penulti_short", "Penulti_long"), ordered=T)

## Fig S1F ##
pdf("all_internal_fs_fp_penultimate_different_cutoffs.pdf", width=11)
p <- ggplot(bar_data, aes(x=Position, y=Percentage, fill=Support)) +
  geom_bar(stat = "identity", colour="black") +
  theme_bw() +
  facet_wrap(.~PSI) +
  scale_fill_manual(name="RNA-seq support", values=c("#999999", "#E31A1C")) +
  xlab("") +
  ylab("Percentage of exons") +
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "top")
print(p)
dev.off()



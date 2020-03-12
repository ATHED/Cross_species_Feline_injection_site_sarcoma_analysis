cellline_data <- read.table("DEGs_cellline_cancer_vs_normal.txt",
                            sep="\t",
                            header=TRUE)

primary_data <- read.table("DEG_primary_with_muscle.txt",
                           sep="\t",
                           header=TRUE)

merged_data <- merge(cellline_data, primary_data, by.x=0, by.y=0)

merged_data_sig <- subset(merged_data, padj.x <= 0.05 & padj.y <= 0.05)

temp <- merged_data_sig[,c(1,2,4)]

cancer_VS_normal_matrix <- read.table("counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt",
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE,
                                      quote="",
                                      comment.char="#")

gene_name_matrix <- cancer_VS_normal_matrix[, c(1,28)]

temp2 <- gene_name_matrix[which(gene_name_matrix$Geneid %in% temp$Row.names),]
target <- temp$Row.names
temp3 <- temp2[match(target, temp2$Geneid),]

temp.full <- data.frame(temp3, temp[,2:3])

temp.yup.xdown <- temp.full[which(temp.full$log2FoldChange.x < 0 & temp.full$log2FoldChange.y > 0),]
temp.yup.xup <- temp.full[which(temp.full$log2FoldChange.x > 0 & temp.full$log2FoldChange.y > 0),]
temp.ydown.xdown <- temp.full[which(temp.full$log2FoldChange.x < 0 & temp.full$log2FoldChange.y < 0),]
temp.ydown.xup <- temp.full[which(temp.full$log2FoldChange.x > 0 & temp.full$log2FoldChange.y < 0),]

write.table(temp.yup.xdown, file="GeneList(PrimaryUPandCelllineDOWN).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(temp.yup.xup, file="GeneList(PrimaryUPandCelllineUP).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(temp.ydown.xdown, file="GeneList(PrimaryDownandCelllineDOWN).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(temp.ydown.xup, file="GeneList(PrimaryDOWNandCelllineUP).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)


library(ggplot2)
ggplot(merged_data_sig,
       aes(log2FoldChange.x, log2FoldChange.y)) +
    theme_classic(base_size=18) +
    theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
    geom_point(alpha=0.3) +
    xlab("log2(cancer/normal), cell line") +
    ylab("log2(cancer/normal), primary") +
    ggsave("plot_cellline_vs_primary.pdf")

signs <- data.frame(cbind(sign(merged_data_sig$log2FoldChange.x),
                          sign(merged_data_sig$log2FoldChange.y)))

library(plyr)
ctable <- matrix(count(signs, vars=c("X1","X2"))$freq[c(3,4,1,2)], ncol=2, byrow=TRUE)

print(sprintf("Number of DEGs in primary data: %d", length(which(primary_data$padj <= 0.05))))
print(sprintf("Number of DEGs in cell line data: %d", length(which(cellline_data$padj <= 0.05))))
print(sprintf("Number of genes that are differentially expressed in both: %d", nrow(merged_data_sig)))


print(fisher.test(ctable))


####################################
#clean
#######################
#rm(list=ls(all=TRUE))

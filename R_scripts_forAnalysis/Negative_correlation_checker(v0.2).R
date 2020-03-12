################################
#data import
##################################
#cancer VS. normal
cancer_VS_normal_matrix <- read.table("counts_annot_phase1AND2_norm_withLogRatio(v0.2).txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
cancer_VS_normal_matrix <- cancer_VS_normal_matrix[, c(1,28, 51,55)]

#sulforophane VS. untreated
sulforophane_VS_untreated_matrix <- read.table("counts_annot_sulforapha ne_norm(vCatname).txt",
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE,
                                      quote="",
                                      comment.char="#")
sulforophane_VS_untreated_matrix <- sulforophane_VS_untreated_matrix[, c(1,28,39,43)]

####################################################################
#cross match to see genes differentially expressed (<0.05) in both 
###################################################################
cancer_VS_normal_compressed <- cancer_VS_normal_matrix[which(cancer_VS_normal_matrix$padj_cellLine < 0.05),]

sulforophane_VS_untreated_compressed <- sulforophane_VS_untreated_matrix[which(sulforophane_VS_untreated_matrix$padj < 0.05),]

#reduce to intersection of genes
common_genid <- intersect(cancer_VS_normal_compressed$Geneid, sulforophane_VS_untreated_compressed$Geneid)

cancer_VS_normal_compressed_c <- cancer_VS_normal_compressed[cancer_VS_normal_compressed$Geneid %in% common_genid, ]
sulforophane_VS_untreated_compressed_c <- sulforophane_VS_untreated_compressed[sulforophane_VS_untreated_compressed$Geneid %in% common_genid, ]

############################
#output 
###########################
write.table(cancer_VS_normal_compressed_c,
            file="cancer_VS_normal_intersect(withGeneName).txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

write.table(sulforophane_VS_untreated_compressed_c,
            file="sulforophane_VS_untreated_intersect(withGeneName).txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)


##########################
#correlation
############################

correlation_CVN_SVU <- cor(cancer_VS_normal_compressed_c$log2FoldChange_cellLine, sulforophane_VS_untreated_compressed_c$log2FoldChange,
                           method = "pearson")

#############################
#scatter plot
###############################
#g.cell.line <- plot(cancer_VS_normal_compressed_c$log2FoldChange_cellLine, sulforophane_VS_untreated_compressed_c$log2FoldChange,
#                    xlab="cancer_VS_normal",
#                    ylab="sulforophane_VS_untreated")

dat <- data.frame(cancer_VS_normal_compressed_c$log2FoldChange_cellLine, sulforophane_VS_untreated_compressed_c$log2FoldChange)
rownames(dat) <- cancer_VS_normal_compressed_c$Geneid
colnames(dat) <- c("log2(cancer/normal)","log2(sulforophane/untreated)")


library(ggplot2)
g.cell.line <- ggplot(dat, aes(x=dat$`log2(cancer/normal)`, y=dat$`log2(sulforophane/untreated)`)) +
                      geom_point() +
                      geom_vline(xintercept = 0) +
                      geom_hline(yintercept = 0) +
                      xlab(expression(log[2](cancer/normal))) +
                      ylab(expression(log[2](sulforophane/untreated))) +
                      geom_text(label=cancer_VS_normal_compressed_c$Associated.Gene.Name, hjust=0, vjust=0, angle= 45, size=8) +
                      theme(text = element_text(size=20) )

pdf("CVN_VS_SVU(vLF).pdf", width = 16, height =12)
print(g.cell.line)
dev.off()

###########################
#fisher test
############################

fisher.df <- matrix(c(12, 2, 1, 4),nrow = 2)
fisher.test(fisher.df)


####################################
#clean
#######################
#rm(list=ls(all=TRUE))










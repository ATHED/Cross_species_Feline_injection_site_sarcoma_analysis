################################
#data import
##################################
#cancer VS. normal
cancer_VS_normal_matrix <- read.table("../datafiles/counts_annot_phase1AND2_norm_withLogRatio(v0.2).txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")


data_tissue_vs_cells <- subset(cancer_VS_normal_matrix,
                               padj <= 0.05 & padj_cellLine <= 0.05)

print(table(sign(data_tissue_vs_cells[, c("log2FoldChange", "log2FoldChange_cellLine")])))

g.cell.line <- ggplot(data_tissue_vs_cells,
                      aes(x=log2FoldChange, y=log2FoldChange_cellLine)) +
    geom_point(colour="#505050", size=2) +
    theme_minimal(base_size=18) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(expression(log[2](sarcoma/skin))) +
    ylab(expression(log[2](sarcoma/fibroblast))) +
        xlim(-8, 8) +
        ylim(-6, 6) +
    ggsave("scatter_plot_primary_tisue_vs_cell_lines.pdf", width=4.5, height=4.5)



cancer_VS_normal_matrix_slim <- cancer_VS_normal_matrix[, c("Geneid", "Associated.Gene.Name",
                                                       "log2FoldChange_cellLine",
                                                       "padj_cellLine")]

#sulforophane VS. untreated
sulforophane_VS_untreated_matrix <- read.table("../datafiles/counts_annot_sulforaphane_norm(vCatname).txt",
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE,
                                      quote="",
                                      comment.char="#")
sulforophane_VS_untreated_matrix <- sulforophane_VS_untreated_matrix[, c("Geneid", "Associated.Gene.Name",
                                                                         "log2FoldChange", "padj")]

####################################################################
#cross match to see genes differentially expressed (<0.05) in both 
###################################################################
cancer_VS_normal_compressed <- cancer_VS_normal_matrix_slim[which(cancer_VS_normal_matrix_slim$padj_cellLine <= 0.05),]

sulforophane_VS_untreated_compressed <- sulforophane_VS_untreated_matrix[which(sulforophane_VS_untreated_matrix$padj <= 0.05),]

#reduce to intersection of genes
common_geneid <- intersect(cancer_VS_normal_compressed$Geneid, sulforophane_VS_untreated_compressed$Geneid)

cancer_VS_normal_compressed_c <- cancer_VS_normal_compressed[cancer_VS_normal_compressed$Geneid %in% common_geneid, ]
sulforophane_VS_untreated_compressed_c <- sulforophane_VS_untreated_compressed[sulforophane_VS_untreated_compressed$Geneid %in% common_geneid, ]

############################
#output 
###########################
write.table(cancer_VS_normal_compressed_c,
            file="../datafiles/cancer_VS_normal_intersect(withGeneName)(v20170412).txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

write.table(sulforophane_VS_untreated_compressed_c,
            file="../datafiles/sulforophane_VS_untreated_intersect(withGeneName)(v20170412).txt",
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


gene_names <- cancer_VS_normal_compressed_c$Associated.Gene.Name
gene_names[gene_names == ""] <- rownames(dat)[gene_names==""]

library(ggplot2)
g.cell.line <- ggplot(dat, aes(x=dat$`log2(cancer/normal)`, y=dat$`log2(sulforophane/untreated)`)) +
    geom_point(colour="gray") +
    theme_minimal(base_size=18) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(expression(log[2](sarcoma/fibroblast))) +
    ylab(expression(log[2](sulforaphane/untreated))) +
    xlim(-6, 6) +
    geom_text(label=gene_names, vjust=0,
              colour="#505050",
              size=c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.5, 3, 3, 3, 2, 2, 2, 3),
              angle=c(0, 0, 0, 0, 0, 0, -30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -30, 0),
              hjust=c(0, 0, 0, 0, 0.4, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 01),
              nudge_x=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
              nudge_y=c(0, 0, -0.1, 0.05, -0.1, 0, 0.0, 0.05, -0.05, 0, 0, -0.03, 0, 0, 0, 0, 0, 0, -0.07)) +
    ggsave("scatter_plot_sulforaphane_vs_cancereffect.pdf", width=4.5, height=4.5)

g.cell.line <- ggplot(dat, aes(x=dat$`log2(cancer/normal)`, y=dat$`log2(sulforophane/untreated)`)) +
    geom_point(colour="#101010", size=2) +
    theme_minimal(base_size=18) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(expression(log[2](sarcoma/fibroblast))) +
    ylab(expression(log[2](sulforaphane/untreated))) +
    xlim(-6, 6) +
#    geom_text(label=gene_names, vjust=0,
#              colour="#505050",
#              size=c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.5, 3, 3, 3, 2, 2, 2),
#              angle=c(0, 0, 0, 0, 0, 0, -30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -30),
#              hjust=c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1),
#              nudge_x=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#              nudge_y=c(0, 0, -0.1, 0.05, -0.12, 0, 0.0, 0.05, -0.05, 0, 0, -0.03, 0, 0, 0, 0, 0, 0)) +
    ggsave("scatter_plot_sulforaphane_vs_cancereffect_unannot.pdf", width=4.5, height=4.5)

###########################
#fisher test
############################

fisher.df <- matrix(c(12, 2, 1, 4),nrow = 2)
print(fisher.test(fisher.df))


####################################
#clean
#######################
#rm(list=ls(all=TRUE))










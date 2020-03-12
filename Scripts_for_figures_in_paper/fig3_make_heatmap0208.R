################################
#data import
##################################
#cancer VS. normal
cancer_VS_normal_matrix <- read.table("./datafiles/counts_annot_phase1AND2_norm_withLogRatio(v20170412)_JUNE.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")

labels <- ifelse(
    cancer_VS_normal_matrix$Associated.Gene.Name != "",
    paste(" (", cancer_VS_normal_matrix$Associated.Gene.Name, ")", sep=""),
    cancer_VS_normal_matrix$Associated.Gene.Name)

rownames(cancer_VS_normal_matrix) <- paste(cancer_VS_normal_matrix$Geneid, labels, sep="")

data_tissue_vs_cells <- subset(cancer_VS_normal_matrix,
                               padj <= 0.05 & padj_cellLine <= 0.05)[, c("X1009_normal",
                                                                         "TigSK_S39.1",
                                                                         "LB.Skin_S40",
                                                                         "X1009_sarcoma",
                                                                         "X1093_sarcoma",
                                                                         "X1033_sarcoma",
                                                                         "Tig.Skin.Fib.C_S35",
                                                                         "Tig.Skin.Fib.SFN_S36",
                                                                         "LB.Skin.Fib.C_S37",
                                                                         "LB.Skin.Fib.SFN_S38",
                                                                         "Ela_C_S31.1",
                                                                         "Ela_SFN_S32.1",
                                                                         "Kaiser_C_S33.1",
                                                                         "Kaiser_SFN_S34.1",
                                                                         "log2FoldChange",
                                                                         "log2FoldChange_cellLine")]

names(data_tissue_vs_cells) <- c("skin.1",
                                 "skin.2",
                                 "skin.3",
                                 "sarcoma.1",
                                 "sarcoma.2",
                                 "sarcoma.3",
                                 "fibrobls.1a",
                                 "fibrobls.1b",
                                 "fibrobls.2a",
                                 "fibrobls.2b",
                                 "fissline.1a",
                                 "fissline.1b",
                                 "fissline.2a",
                                 "fissline.2b",
                                 "lfc_tissue",
                                 "lfc_cells")

row_order <- order(-sign(data_tissue_vs_cells$lfc_tissue),
                   -sign(data_tissue_vs_cells$lfc_cells),
                   -data_tissue_vs_cells$lfc_tissue - data_tissue_vs_cells$lfc_cells)

reordered_df <- data.frame(t(scale(t(data_tissue_vs_cells[row_order, 1:14]), center=TRUE, scale=FALSE)))

reordered_df$label <- factor(1:nrow(reordered_df), labels=rownames(reordered_df))

library(reshape2)
melted_df <- melt(reordered_df)

library(ggplot2)
ggplot(melted_df, aes(label, variable)) +
    geom_tile(aes(fill=value), colour="transparent") +
    scale_fill_gradient2(low="blue", high="yellow") +
    scale_x_discrete(position="top") +
    theme_minimal(base_size=9) +
    theme(axis.title.y=element_blank(),
          legend.position="bottom",
          legend.text=element_text(size=9),
            legend.title=element_blank(),
          axis.title.x=element_blank(),
            axis.text.x=element_text(size=2,angle=90,hjust=0,vjust=0.5)) + ggsave("heatmap.pdf", width=12, height=3)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
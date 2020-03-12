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

label <- read.table("label.txt", header=TRUE)
rownames(label) <- label[,1]

common_genid <- intersect(cancer_VS_normal_matrix$Geneid, label$Ensembl_Gene_ID)

cancer_VS_normal_compressed_c <- cancer_VS_normal_matrix[cancer_VS_normal_matrix$Geneid %in% common_genid, ]
cancer_VS_normal_compressed_c <- cancer_VS_normal_compressed_c[,c(1,51,54,55)]
rownames(cancer_VS_normal_compressed_c) <- cancer_VS_normal_compressed_c[,1]

cancer_VS_normal_compressed_o <- cancer_VS_normal_compressed_c[rownames(label),,drop=F]

write.table(cancer_VS_normal_compressed_o,
            file="attach.txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
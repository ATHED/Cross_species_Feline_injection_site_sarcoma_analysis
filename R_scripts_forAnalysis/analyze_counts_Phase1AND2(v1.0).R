##########################
#import necessary library
#############################
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

##################################
#import datasets
################################
count_matrix <- read.table("phase1AND2_counts.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")

#cat_gene_annot <- read.table("Ensembl_85_FelisCatus62_gene_annot.txt",
#                             header=TRUE,
#                             sep="/t",
#                             quote="",
#                             stringsAsFactors=FALSE)
cat_gene_annot <- read.csv("Ensembl_85_FelisCatus62_gene_annot.txt")


## apparently there is no cytogenetic band information for the cat genome, so this just removes a column that would be all "NA" anyhow
cat_gene_annot <- cat_gene_annot[, setdiff(1:ncol(cat_gene_annot),which(names(cat_gene_annot) == "Band"))]

#####################################################
#bind the annotation data with count data
##############################################
count_matrix_annot <- cbind(count_matrix,
                            cat_gene_annot[match(count_matrix$Geneid,
                                                 cat_gene_annot$Gene.ID),])

sample_names <- c("1099_normal",
                  "1099_sarcoma",
                  "1033_sarcoma",
                  "1093_sarcoma",
                  "Ela_C_S31",
                  "Ela_SFN_S32",
                  "Kaiser_C_S33",
                  "Kaiser_SFN_S34",
                  "Tig-Skin-Fib-C_S35",
                  "Tig-Skin-Fib-SFN_S36",
                  "LB-Skin-Fib-C_S37",
                  "LB-Skin-Fib-SFN_S38",
                  "TigSK_S39",
                  "LB-Skin_S40")

names(count_matrix_annot)[7:20] <- sample_names

#column_data <- data.frame(sarcoma=factor(c(0,1,1,1)),
#                          type=rep("paired-read",4),
#                          row.names=sample_names)

#dds <- DESeqDataSetFromMatrix(countData=count_matrix_annot[,7:20],
#                               colData=column_data,
#                               design=~sarcoma)


##################################################
#Deseq pipeline
#####################################################
condition_2 <- factor(c("normal","cancer","cancer","cancer","cancer","cancer","cancer","cancer","normal","normal","normal","normal","normal","normal"))
condition_1 <- factor(c("primary tissue", "primary tissue", "primary tissue", "primary tissue", "cell line", "cell line",
                         "cell line", "cell line", "cell line", "cell line", "cell line", "cell line", "primary tissue", "primary tissue"))

column_data <- data.frame(condition_1, condition_2, row.names=sample_names)

##########################
#Deseq for primary tissue
####################################
column_data.primary.tissue <- column_data[which(column_data$condition_1 == "primary tissue"),]
column.data.primary.tissue <- data.frame(design=factor(column_data.primary.tissue[,2]), row.names=rownames(column_data.primary.tissue))

match(colnames(count_matrix_annot), rownames(column_data.primary.tissue))
countData.primary.tissue <- count_matrix_annot[ ,c(7:10, 19,  20)]


primary.tissue.dds <- DESeqDataSetFromMatrix(countData=countData.primary.tissue, colData=column.data.primary.tissue
                                             ,design=~design)

primary.tissue.dds <- estimateSizeFactors(primary.tissue.dds)

norm_counts.primary.tissue <- counts(primary.tissue.dds, normalized=TRUE)

#log transformation
log2_norm_counts.primary.tissue <- log2(1+norm_counts.primary.tissue)
#sarcoma_vs_normal_log2_ratio <- apply(log2_norm_counts[,2:4],1,mean) - log2_norm_counts[,1]
    
#count_matrix_annot_norm <- cbind(count_matrix_annot, log2_norm_counts, sarcoma_vs_normal_log2_ratio)



#write.table(count_matrix_annot_norm[order(sarcoma_vs_normal_log2_ratio, decreasing=TRUE),],
#            file="counts_annot_norm.txt",
#            sep="\t",
#            row.names=FALSE,
#            col.names=TRUE,
#            quote=FALSE)

####################################################
#compute an intraspecies log2 expression ratio and p‑value 
#(for the test of equal means of the log2 expression levels of the two sample groups) 
#based on a negative binomial distribution-based count model
#################################################
df_primary.tissue <- results(DESeq(primary.tissue.dds), contrast=c("design","cancer", "normal"))

head(df_primary.tissue)

#count_matrix_annot_norm <- cbind(count_matrix_annot, log2_norm_counts.primary.tissue, df_primary.tissue)

##########################
#Deseq for cell line
####################################
column_data.cell.line <- column_data[which(column_data$condition_1 == "cell line"),]
column.data.cell.line <- data.frame(design=factor(column_data.cell.line[,2]), row.names=rownames(column_data.cell.line))


match(colnames(count_matrix_annot), rownames(column_data.cell.line))
countData.cell.line <- count_matrix_annot[ ,c(11:18)]

cell.line.dds <- DESeqDataSetFromMatrix(countData=countData.cell.line, colData=column.data.cell.line
                                             ,design=~design)

cell.line.dds <- estimateSizeFactors(cell.line.dds)

norm_counts.cell.line <- counts(cell.line.dds, normalized=TRUE)

#log transformation
log2_norm_counts.cell.line <- log2(1+norm_counts.cell.line)
#sarcoma_vs_normal_log2_ratio <- apply(log2_norm_counts[,2:4],1,mean) - log2_norm_counts[,1]

#count_matrix_annot_norm <- cbind(count_matrix_annot, log2_norm_counts, sarcoma_vs_normal_log2_ratio)

####################################################
#compute an intraspecies log2 expression ratio and p‑value 
#(for the test of equal means of the log2 expression levels of the two sample groups) 
#based on a negative binomial distribution-based count model
#################################################
df_cell.line <- results(DESeq(cell.line.dds), contrast=c("design","cancer", "normal"))

head(df_cell.line)

out_names.cell.line <- c("baseMean_cellLine", "log2FoldChange_cellLine", "lfcSE_cellLine", "stat_cellLine"
                         , "pvalue_cellLine", "padj_cellLine")
names(df_cell.line) <- out_names.cell.line

count_matrix_annot_norm <- cbind(count_matrix_annot, log2_norm_counts.primary.tissue, df_primary.tissue 
                                 , log2_norm_counts.cell.line, df_cell.line)

#########################################
#data output
############################################
 




##############################################################
fix_hugo_symbol <- function(allcaps_symbol) {
    paste(c(substring(allcaps_symbol, 1, 1), tolower(substring(allcaps_symbol, 2))), collapse="")
}

#############################################
#primary tissue
#############################################
cell.line.data_for_gsea <- count_matrix_annot_norm[ ,c(21, 24, 28, 42:55)]
cell.line.data_for_gsea <- cell.line.data_for_gsea[which(cell.line.data_for_gsea$Associated.Gene.Name != ""),]
cell.line.data_for_gsea[,1] <- sapply(cell.line.data_for_gsea[,1], fix_hugo_symbol)
cell.line.data_for_gsea <- cell.line.data_for_gsea[ ,c(3,15)]
write.table(cell.line.data_for_gsea[order(cell.line.data_for_gsea$stat,decreasing=TRUE),],
            file="cell_line_gsea.rnk",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

#############################################
#cell line
#############################################
primary.tissue.data_for_gsea <- count_matrix_annot_norm[ ,c(21, 24, 28, 30:41)]
primary.tissue.data_for_gsea <- primary.tissue.data_for_gsea[which(primary.tissue.data_for_gsea$Associated.Gene.Name != ""),]
primary.tissue.data_for_gsea[,1] <- sapply(primary.tissue.data_for_gsea[,1], fix_hugo_symbol)
primary.tissue.data_for_gsea <- primary.tissue.data_for_gsea[ ,c(3,13)]

write.table(primary.tissue.data_for_gsea[order(primary.tissue.data_for_gsea$stat,decreasing=TRUE),],
            file="primary_tissue_gsea.rnk",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

#data_for_gsea <- count_matrix_annot_norm[,c(17, 27)]
#data_for_gsea <- data_for_gsea[which(data_for_gsea$Associated.Gene.Name != ""),]
#data_for_gsea[,1] <- sapply(data_for_gsea[,1], fix_hugo_symbol)
#write.table(data_for_gsea[order(data_for_gsea[,2],decreasing=TRUE),],
#            file="data_for_gsea.rnk",
#            sep="\t",
#            quote=FALSE,
#            row.names=FALSE,
#           col.names=FALSE)


####################################
#clean
#######################
#rm(list=ls(all=TRUE))

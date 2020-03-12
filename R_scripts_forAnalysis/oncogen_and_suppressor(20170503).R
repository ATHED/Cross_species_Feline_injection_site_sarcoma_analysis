##################################
#import datasets
################################
log2_matrix <- read.table("counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
#filter out based on both primary tissue and cell line 
#log2_matrix <- subset(log2_matrix, log2_matrix$padj <= 0.05 & log2_matrix$padj_cellLine <= 0.05)

#filter out only based on primary tissue
log2_matrix <- subset(log2_matrix, log2_matrix$padj <= 0.05)

#filter out only based on cell line
#log2_matrix <- subset(log2_matrix, log2_matrix$padj_cellLine <= 0.05)

log2_matrix_compress <- log2_matrix[, c(1, 28, 37, 51)]
colnames(log2_matrix_compress) <- c("gene_id", "gene_name", "log2_primary_tissue", "log2_cellline")


##########################
#import list of oncogene and TSG
#############################
oncogene_list <- read.table("(minzhao)Oncogene_list.txt",
                           sep="\t",
                           header=F,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")

TSG_list <- read.table("TSG_list.txt",
                            sep="\t",
                            header=F,
                            stringsAsFactors=FALSE,
                            quote="",
                            comment.char="#")

num_intersection_2 <- intersect(oncogene_list$V1, TSG_list$V1)


temp_oncogene <- log2_matrix_compress[which(log2_matrix_compress$gene_name %in% oncogene_list$V1),]

temp_oncogene$type <- "oncogene"

temp_TSG <- log2_matrix_compress[which(log2_matrix_compress$gene_name %in% TSG_list$V1),]

temp_TSG$type <- "TSG"

num_intersection <- intersect(temp_oncogene$gene_name, temp_TSG$gene_name)

TSG_oncogen_list <- rbind(oncogene_list, TSG_list)

temp_rest <- log2_matrix_compress[-which(log2_matrix_compress$gene_name %in% TSG_oncogen_list$V1),]

temp_rest$type <- "other genes"

temp_oncogene <- temp_oncogene[-which(temp_oncogene$gene_name %in% num_intersection),]
temp_TSG <- temp_TSG[-which(temp_TSG$gene_name %in% num_intersection),]

log2.cvn.full <- rbind(temp_oncogene, temp_TSG, temp_rest)

log2.cvn.primary <- log2.cvn.full[,c(3,5)]
log2.cvn.cellline <- log2.cvn.full[,c(4,5)]

library(ggplot2)

primary.tissue.graph <- ggplot(log2.cvn.primary, aes(x=log2_primary_tissue)) +xlab(expression(paste(log[2](cancer/normal), " in cat"))) + 
  geom_density(aes(group=type, colour=type, fill=type) , alpha=0.3)
primary.tissue.graph <- primary.tissue.graph + theme_bw() + theme(axis.title = element_text(size = 14)
                                                        , axis.text = element_text(size=14), legend.text = element_text(size = 14))

pdf("log2ration_density_cancerVSnormal_primary_tissue(20170503).pdf")
print(primary.tissue.graph)
dev.off()

library(ggplot2)

#cell.line.graph <- ggplot(log2.cvn.cellline, aes(x=log2_cellline)) +xlab(expression(paste(log[2](cancer/normal), " in cat"))) + 
#                                                                    geom_density(aes(group=type, colour=type, fill=type) , alpha=0.3)
#cell.line.graph <- cell.line.graph + theme_bw() + theme(axis.title = element_text(size = 14)
#                                                        , axis.text = element_text(size=14), legend.text = element_text(size = 14))

#pdf("log2ration_density_cancerVSnormal_cell_line(20170503).pdf")
#print(cell.line.graph)
#dev.off()

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
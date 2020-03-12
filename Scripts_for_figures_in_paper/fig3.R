##################################
#import datasets
################################
log2_matrix <- read.table("../datafiles/counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")


oncogene_list <- read.table("../datafiles/(minzhao)Oncogene_list.txt",
                           sep="\t",
                           header=F,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")$V1

TSG_list <- read.table("../datafiles/TSG_list.txt",
                            sep="\t",
                            header=F,
                            stringsAsFactors=FALSE,
                            quote="",
                            comment.char="#")$V1

log2_matrix$type <- rep("none", nrow(log2_matrix))
log2_matrix$type[na.omit(match(oncogene_list, log2_matrix$Associated.Gene.Name))] <- "oncogene"
log2_matrix$type[na.omit(match(TSG_list, log2_matrix$Associated.Gene.Name))] <- "TSG"

library(ggplot2)

data_tissue <- subset(log2_matrix, padj <= 0.05)

ggplot(data_tissue, aes(x=log2FoldChange)) +
    xlab(expression(paste(log[2](cancer/normal)))) +
    geom_line(aes(group=type, colour=type), stat="density") +
    theme_minimal(base_size=10) +
    theme(legend.title=element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size=10),
          legend.text = element_text(size = 10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggsave("fig3_primary_tissue.pdf", width=3, height=2)

counts <- table(data.frame(sgn=sign(data_tissue$log2FoldChange),
                           type=data_tissue$type))
colnames(counts) <- c("TSG", "neither", "oncogene")
props_df <- data.frame(100 * counts / matrix(rep(apply(counts, 2, sum), each=2), ncol=3))
props_df$sgn <- factor(props_df$sgn, labels=c("downregulated", "upregulated"))
ggplot(props_df, aes(fill=sgn, x=type, y=Freq)) +
    geom_bar(position="dodge",
             stat="identity") +
    theme_minimal(base_size=12) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1)) +
    ylab("percent") +
    labs(fill="") +
    ggsave("fig3_tissue_props.pdf", width=3, height=2.5)

data_cells <- subset(log2_matrix, padj_cellLine <= 0.05)

ggplot(data_cells, aes(x=log2FoldChange_cellLine)) +
    xlab(expression(paste(log[2](cancer/normal)))) + 
    geom_line(aes(group=type, colour=type), stat="density") +
    theme_minimal(base_size=10) +
    theme(legend.title=element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size=10),
          legend.text = element_text(size = 10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggsave("fig3_cell_lines.pdf", width=3, height=2)

counts <- table(data.frame(sgn=sign(data_cells$log2FoldChange_cellLine),
                           type=data_cells$type))
colnames(counts) <- c("TSG", "neither", "oncogene")
props_df <- data.frame(100 * counts / matrix(rep(apply(counts, 2, sum), each=2), ncol=3))
props_df$sgn <- factor(props_df$sgn, labels=c("downregulated", "upregulated"))
ggplot(props_df, aes(fill=sgn, x=type, y=Freq)) +
    geom_bar(position="dodge",
             stat="identity") +
    theme_minimal(base_size=12) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1)) +
    ylab("percent") +
    labs(fill="") +
    ggsave("fig3_cells_props.pdf", width=3, height=2.5)


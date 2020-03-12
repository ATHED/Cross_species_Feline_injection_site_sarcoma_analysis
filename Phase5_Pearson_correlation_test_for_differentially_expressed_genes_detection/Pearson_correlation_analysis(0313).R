##########################
#import necessary library
#############################
#install.packages("openxlsx")
#library(openxlsx)

##################################
#import datasets
################################
mRNA.expr.data.og <- read.table("matProdMSE.csv",
                           sep=",",
                           header=T,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
gene.expr.data.og <- read.table("Differentially_Expressed_Genes_primary_tissues_cancer_vs_normal.csv",
                           sep=",",
                           header=T,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#",
                           fill = T)

##############################################
#data selection
##################################################
colnames(mRNA.expr.data.og)
mRNA.expr.data.beforeLog <- mRNA.expr.data.og[,c(1,2,8,11,14,17,20,23,26)]
#mRNA.expr.data.beforeLog$fca_cocoa_STS_raw <- as.numeric(as.character(mRNA.expr.data.beforeLog$fca_cocoa_STS_raw))
mRNA.expr.data.beforeLog[,c(3:9)] <- lapply(mRNA.expr.data.beforeLog[,c(3:9)], function(x) as.numeric(as.character(x)))
mRNA.expr.data <- cbind(mRNA.expr.data.beforeLog[,c(1:2)], log(mRNA.expr.data.beforeLog[,c(3:9)]+1))

#list.interested.mRNA <- c("rr74363-2a-3p-miR", "rr2413226-5a", "rr378630-2a", "rr1381690-4a-5p-miR","rr3035894-2a", "rr840812-2b-5p-miR", "rr46761-1a",
#                          "rr1175091-1b-5p-miR", "rr1175088-3a-3p-miR", "rr2711500-2a-5p-miR", "rr2054804-4a")
list.interested.mRNA <- c("mir-1","mir-21","mir-133","mir-142","mir-143","mir-182","mir-205","mir-221","let-7a")
mRNA.expr.data.sel <- mRNA.expr.data[which(mRNA.expr.data$mir_Family %in% list.interested.mRNA),]

#
colnames(gene.expr.data.og)
gene.expr.data <- gene.expr.data.og[,c(1,2,4,7:34)]
colnames(gene.expr.data)
list.interested.gene <- c("X1033_sarcoma.1", "X1093_sarcoma.1", "Tig.Skin.Fib.C_S35.1", "Tig.Skin.Fib.SFN_S36.1", "TigSK_S39.1")
gene.expr.data.sel <- gene.expr.data[,c(1,2,3,18,19,20,21,22)]

#adjust the sequence and name of samples to be the same
mRNA.expr.data.cal <- mRNA.expr.data.sel[,c(1, 4, 6, 7,8,9)]
rownames(mRNA.expr.data.cal) <- mRNA.expr.data.cal[,1]
mRNA.expr.data.cal <- mRNA.expr.data.cal[,-1]
colnames(mRNA.expr.data.cal) <- c("Cocoa_sarcoma", "Mojo_sarcoma","Tigger_muscle","Tiger_skin","Tigger_sarcoma")

#
gene.expr.data.cal <- gene.expr.data.sel[,c(1,7,6,4,8,5)]
rownames(gene.expr.data.cal) <- gene.expr.data.cal[,1]
gene.expr.data.cal <- gene.expr.data.cal[,-1]
colnames(gene.expr.data.cal) <- c("Cocoa_sarcoma", "Mojo_sarcoma","Tigger_muscle","Tiger_skin","Tigger_sarcoma")

######################################################
#build an empty martix
############################################################
corr.matrix <- matrix(, nrow = nrow(mRNA.expr.data.cal), ncol = nrow(gene.expr.data.cal))
rownames(corr.matrix) <- rownames(mRNA.expr.data.cal)
colnames(corr.matrix) <- rownames(gene.expr.data.cal)

################################################
#conversion to numeric
#######################################################
mRNA.expr.data.num <- data.matrix(mRNA.expr.data.cal)
gene.expr.data.num <- data.matrix(gene.expr.data.cal)

################################################################
#calculate pearson correlation
#######################################################
for (i in 1:nrow(mRNA.expr.data.cal)){
  for (j in 1:nrow(gene.expr.data.cal))
    corr.matrix[i,j] <- cor(mRNA.expr.data.num[i,], gene.expr.data.num[j,], method = "pearson")
}

###############################################
#output matrix
###########################################
write.table(corr.matrix, file = "Known_9_family_mRNA_pearson_correlation(0313).tsv", append = FALSE, quote = F, sep = "\t",
          eol = "\n", na = "NA", dec = ".", row.names = TRUE,
          col.names = TRUE, qmethod = c("escape", "double"),
          fileEncoding = "")



####################################
#clean
#######################
#rm(list=ls(all=TRUE))

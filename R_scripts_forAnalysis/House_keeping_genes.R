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

cancer_VS_normal_matrix <- cancer_VS_normal_matrix[, c(1,28, 30:35, 42:49)]

cancer_VS_normal_compressed <- cancer_VS_normal_matrix[which(cancer_VS_normal_matrix$Associated.Gene.Name %in% c("GAPDH","HMBS", "YWHAZ", "PPIA") 
                                                             | cancer_VS_normal_matrix$Geneid %in% "ENSFCAG00000032017"),]

cancer_VS_normal_compressed[4, 2] <- "YWHAZ"

dat <- cancer_VS_normal_compressed[,3:16]
rownames(dat) <- cancer_VS_normal_compressed[,2]

##############
#PPIA
###############
ppia.mean <- mean(as.numeric(dat[1,]))
GAPDH.mean <- mean(as.numeric(dat[2,]))
HMBS.mean <- mean(as.numeric(dat[3,]))
YWHAZ.mean <- mean(as.numeric(dat[4,]))

###############
ppia.sd <- sd(as.numeric(dat[1,]))
GAPDH.sd <- sd(as.numeric(dat[2,]))
HMBS.sd <- sd(as.numeric(dat[3,]))
YWHAZ.sd <- sd(as.numeric(dat[4,]))

ppia.cv <- ppia.sd/ppia.mean
GAPDH.cv <- GAPDH.sd/GAPDH.mean
HMBS.cv <- HMBS.sd/HMBS.mean
YWHAZ.cv <- YWHAZ.sd/YWHAZ.mean

cat("coefficients of variation of PPIA:  " , ppia.cv)
cat("coefficients of variation of GAPDH: " , GAPDH.cv)
cat("coefficients of variation of HMBS:  " , HMBS.cv)
cat("coefficients of variation of YWHAZ:  " , YWHAZ.cv)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
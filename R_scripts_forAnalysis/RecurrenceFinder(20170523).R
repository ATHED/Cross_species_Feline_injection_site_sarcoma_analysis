##################################
#import datasets
################################
normal.1009 <- read.table("1009_normal_fusion_candidates.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
sarcoma.1009 <- read.table("1009_sarcoma_fusion_candidates.txt",
                          sep="\t",
                          header=TRUE,
                          stringsAsFactors=FALSE,
                          quote="",
                          comment.char="#")
sarcoma.1033 <- read.table("1033_sarcoma_fusion_candidates.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
sarcoma.1093 <- read.table("1093_sarcoma_fusion_candidates.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
skin.s039 <- read.table("s039_skin_fusion_candidates.txt",
                           sep="\t",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
skin.s040 <- read.table("s040_skin_fusion_candidates.txt",
                        sep="\t",
                        header=TRUE,
                        stringsAsFactors=FALSE,
                        quote="",
                        comment.char="#")

####################################
#full dataframe of fusion genes
######################################
normal.1009.compress <- normal.1009[,c(1:3, 6, 8)]
normal.1009.compress[, "ResultFrom"] <- "normal.1009"

sarcoma.1009.compress <- sarcoma.1009[,c(1:3, 6, 8)]
sarcoma.1009.compress[, "ResultFrom"] <- "sarcoma.1009"

sarcoma.1033.compress <- sarcoma.1033[,c(1:3, 6, 8)]
sarcoma.1033.compress[, "ResultFrom"] <- "sarcoma.1033"

sarcoma.1093.compress <- sarcoma.1093[,c(1:3, 6, 8)]
sarcoma.1093.compress[, "ResultFrom"] <- "sarcoma.1093"

skin.s039.compress <- skin.s039[,c(1:3, 6, 8)]
skin.s039.compress[, "ResultFrom"] <- "skin.s039"

skin.s040.compress <- skin.s040[,c(1:3, 6, 8)]
skin.s040.compress[, "ResultFrom"] <- "skin.s040"

full.compress <- rbind(sarcoma.1009.compress, sarcoma.1033.compress, sarcoma.1093.compress)
length(which(full.compress$ResultFrom == c('normal.1009')))
length(which(full.compress$ResultFrom == c('skin.s039')))
length(which(full.compress$ResultFrom == c('skin.s040')))

search.name <- "RIOK1--ENSFCAG00000027290"
result <- which(full.compress$FusionName == search.name)
length(result)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
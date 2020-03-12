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

##############################
#unique list of fusion genes
##############################
normal.1009.names <- normal.1009$FusionName
sarcoma.1009.names <- sarcoma.1009$FusionName
sarcoma.1033.names <- sarcoma.1033$FusionName
sarcoma.1093.names <- sarcoma.1093$FusionName
skin.s039.names <- skin.s039$FusionName
skin.s040.names <- skin.s040$FusionName

full.names <- c(normal.1009.names, sarcoma.1009.names, sarcoma.1033.names, sarcoma.1093.names, skin.s039.names, skin.s040.names)
unique.full.names <- unique(full.names)

full.names.occur <- table(full.names)

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

full.compress <- rbind(normal.1009.compress,skin.s039.compress, skin.s040.compress,
                       sarcoma.1009.compress, sarcoma.1033.compress, sarcoma.1093.compress)
length(which(full.compress$ResultFrom == c('normal.1009')))
length(which(full.compress$ResultFrom == c('skin.s039')))
length(which(full.compress$ResultFrom == c('skin.s040')))

#full.compress.sorted <- full.compress[order(full.compress$LeftBreakpoint, full.compress$RightBreakpoint),]

full.compress.unique <- full.compress[!duplicated(full.compress[4:5]),]

length(which(full.compress.unique$ResultFrom == 'normal.1009'))
length(which(full.compress$ResultFrom == c('skin.s039')))
length(which(full.compress$ResultFrom == c('skin.s040')))

full.compress.NoNormal <- full.compress.unique[-which(full.compress.unique$ResultFrom == 'normal.1009'),]
full.compress.NoNormal <- full.compress.NoNormal[-which(full.compress.NoNormal$ResultFrom == 'skin.s039'),]
full.compress.NoNormal <- full.compress.NoNormal[-which(full.compress.NoNormal$ResultFrom == 'skin.s040'),]

full.compress.NoNormal.sorted <- full.compress.NoNormal[order(full.compress.NoNormal$LeftBreakpoint, full.compress.NoNormal$RightBreakpoint),]
full.compress.NoNormal.sortedbyJR <- full.compress.NoNormal.sorted[order(full.compress.NoNormal.sorted$JunctionReadCount
                                                                         , full.compress.NoNormal.sorted$SpanningFragCount, decreasing = TRUE),]

write.table(full.compress.NoNormal.sortedbyJR, file="SacromaOnly_fusion_genes_RankedbyJR(20170515).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=F)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
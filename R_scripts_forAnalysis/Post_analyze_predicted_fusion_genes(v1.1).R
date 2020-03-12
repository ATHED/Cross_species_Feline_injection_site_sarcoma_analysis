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

full.compress <- rbind(normal.1009.compress, sarcoma.1009.compress, sarcoma.1033.compress, sarcoma.1093.compress,
                       skin.s039.compress, skin.s040.compress)
full.compress.sorted <- full.compress[order(full.compress$LeftBreakpoint, full.compress$RightBreakpoint),]

#full.compress.repeat <- full.compress.sorted[duplicated(full.compress.sorted["LeftBreakpoint"]) 
#                                             & duplicated(full.compress.sorted["RightBreakpoint"]),]

full.compress.repeat <- full.compress.sorted[duplicated(full.compress.sorted[,c(4:5)]) 
                                             | duplicated(full.compress.sorted[,c(4:5)], fromLast=TRUE) ,]
#full.compress.repeat.table <- table(full.compress.repeat$FusionName)
full.compress.repeat.table <- table(full.compress.repeat$LeftBreakpoint, full.compress.repeat$RightBreakpoint)
full.compress.repeat.table <- data.frame(full.compress.repeat.table)
full.compress.repeat.table <- full.compress.repeat.table[full.compress.repeat.table$Freq !=0, ]
  
colnames(full.compress.repeat.table) <- c("LeftBreakpoint", "RightBreakpoint", "Frequency")

temp.ResultFrom <- data.frame()
temp.ResultFrom <- aggregate(full.compress.repeat$ResultFrom ~ full.compress.repeat$LeftBreakpoint + full.compress.repeat$RightBreakpoint
                                    , data=full.compress.repeat, FUN=paste)
colnames(temp.ResultFrom) <- c("LeftBreakpoint", "RightBreakpoint", "ResultFrom")

temp.FusionName <- data.frame()
temp.FusionName <- aggregate(full.compress.repeat$FusionName ~ full.compress.repeat$LeftBreakpoint + full.compress.repeat$RightBreakpoint
                             , data=full.compress.repeat, FUN=paste)
colnames(temp.FusionName) <- c("LeftBreakpoint", "RightBreakpoint", "FusionName")

temp <- data.frame()
temp <- aggregate(cbind(full.compress.repeat$JunctionReadCount, full.compress.repeat$SpanningFragCount) ~ full.compress.repeat$LeftBreakpoint + full.compress.repeat$RightBreakpoint
                             , data=full.compress.repeat, FUN=sum)
colnames(temp) <- c("LeftBreakpoint", "RightBreakpoint", "TotalJunctionReadCount","TotalSpanningFragCount")

new.temp <-data.frame()
new.temp <- merge(full.compress.repeat.table, temp, by= c("LeftBreakpoint","RightBreakpoint"))
new.temp <- merge(temp.FusionName, new.temp, by= c("LeftBreakpoint","RightBreakpoint"))
new.temp <- merge(new.temp, temp.ResultFrom, by= c("LeftBreakpoint","RightBreakpoint"))
typeof(new.temp$ResultFrom)
typeof(new.temp$FusionName)

string.temp <- vector()
name.temp <- vector()
for(i in 1:nrow(new.temp)){
  string.temp[i]<- paste( unlist(new.temp[i,]$ResultFrom), collapse=', ')
  name.temp[i]<- paste( unlist(new.temp[i,]$FusionName), collapse=', ')
}
new.temp$ResultFrom <- string.temp
new.temp$FusionName <- name.temp
typeof(new.temp$ResultFrom)
typeof(new.temp$FusionName)

#unique.full.compress <- full.compress.repeat[,c(1,4,5)]
#unique.temp <- unique(unique.full.compress)

write.table(new.temp, file="unique_fusion_genes(20170510).txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=F)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
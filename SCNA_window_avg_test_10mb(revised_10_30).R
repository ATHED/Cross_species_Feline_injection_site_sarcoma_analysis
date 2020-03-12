########################################
#import count data
###########################################

count.full <- read.table("counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt",
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE,
                         quote="",
                         comment.char="#")
colnames(count.full)
count.compress <- count.full[,c(1,24,25,28,37,40,51,54)]
colnames(count.compress) <- c("id", "chr", "start_bp", "symbol", "log2_primary_cvn", "p_value_primary_tissue"
                              , "log2_celline_cvn", "p_value_celline")

count.compress$avg_log2_primary <- 0
count.compress$avg_log2_celline <- 0

#####################################
#count the levels
##############################################
count.compress.factor <- factor(count.compress$chr)
num.obj <- data.frame(table(count.compress$chr))

index.only.one <- which(num.obj$Freq == 1)
index.need.avg <- which(num.obj$Freq != 1)

###################################
#split into various data frames
##########################################
count.data.X <- split(count.compress, count.compress$chr)

#####################################
#assign "avg" to all the sub data frame with only one obj
###################################################
#test
#count.data.X[[110]]

for(i in index.only.one){
  #notice: use this to get an element from lists
  count.data.X[[i]][9] <- count.data.X[[i]][5]
  count.data.X[[i]][10] <- count.data.X[[i]][7]
}

########################################
#assign "avg" to all within 1mb window
####################################################
#test
#count.data.X[[1]]
#count.data.X[[1]][1,3]

#install "zoo" package for roll apply
#install.packages("zoo")
library(zoo)

#set window size to be 10mb
window.size <- 10^7

#width.list <- as.vector(rep(0, nrow(count.data.X[[1]])))
#width.list[1]

#build a big list to record all the window size
width.big.list <- list()
for(i in index.need.avg){
  tmp <- as.vector(rep(0, nrow(count.data.X[[i]])))
  width.big.list[[i]] <- tmp
}

#assign value to all lists

#part 1: test for one chromosome, verified same result
#for(i in 1:nrow(count.data.X[[1]])){
#  for(j in i:nrow(count.data.X[[1]])){
#    if(count.data.X[[1]][i,3] + window.size >= count.data.X[[1]][j,3]){
#      width.list[i] <- width.list[i] + 1
#    }
#  }
#}

#calculate all window width for all the rest chromosome, which only takes 5 mins
for(x in index.need.avg){
  for(i in 1:nrow(count.data.X[[x]])){
    for(j in i:nrow(count.data.X[[x]])){
      if(count.data.X[[x]][i,3] + window.size >= count.data.X[[x]][j,3]){
        width.big.list[[x]][i] <- width.big.list[[x]][i] + 1
      }
    }
  }
}

#save the width big list

lapply(width.big.list, function(x) write.table( data.frame(x), 'width_window_list_10mb_revised_10_30.txt'  , append= T, sep='/t' ))

#save which chr needs avg
write.table(index.need.avg, "index_for_avg_10mb.txt")

#test the list 
list.of.levels <- as.list(levels(count.compress.factor))
not.need.chr <- list.of.levels[index.only.one]

#test on the first list, verified
#nrow( count.data.X[[1]] )
#test <- rollapply(count.data.X[[1]]$log2_primary_cvn, width = width.big.list[[1]], FUN = mean, align = "left", na.rm = T)

#calculate average based on all window sizes
for(x in index.need.avg){
  count.data.X[[x]][9] <- rollapply(count.data.X[[x]]$log2_primary_cvn, width = width.big.list[[x]], FUN = mean, align = "left", na.rm = T)
  count.data.X[[x]][10] <- rollapply(count.data.X[[x]]$log2_celline_cvn, width = width.big.list[[x]], FUN = mean, align = "left", na.rm = T)
}


#######################################
#write output
##########################################

#lapply(count.data.X, function(i){  capture.output( print(  summary(i) ) , 
#                                           file="SCNA_window_10mb_avg.txt", append=TRUE)})

lapply(count.data.X, function(x) write.table( data.frame(x), 'SCNA_window_10mb_avg_revised_10_30.txt'  , append= T, sep='\t' ))


############################################
#out put top 5 and bottom 5 entries
#################################################
#install.packages("rlist")
#library(rlist)

#outstanding.up.list <- list()
#outstanding.down.list <- list()
#for(x in index.need.avg){
#  outstanding.up.list[[x]] <- head(sort(count.data.X[[x]]$avg_log2_primary, decreasing = T, index.return = T), 5)
#  outstanding.down.list[[x]] <- head(sort(count.data.X[[x]]$avg_log2_primary, decreasing = F, index.return = T), 5)
#
#}

#initiate empty lists
tmp.primary.list <- list()
top.primary.list <- list()
down.primary.list <- list()

tmp.celline.list <- list()
top.celline.list <- list()
down.celline.list <- list()

for(x in index.need.avg){
  tmp.primary.list[[x]] <- na.omit(count.data.X[[x]][order(count.data.X[[x]]$avg_log2_primary), ])
  top.primary.list[[x]] <- head(tmp.primary.list[[x]], 10)
  down.primary.list[[x]] <- tail(tmp.primary.list[[x]], 10)
  
  tmp.celline.list[[x]] <- na.omit(count.data.X[[x]][order(count.data.X[[x]]$avg_log2_celline), ])
  top.celline.list[[x]] <- head(tmp.celline.list[[x]], 10)
  down.celline.list[[x]] <- tail(tmp.celline.list[[x]], 10)
}

##########################################
#out put
#############################################
lapply(top.primary.list, function(x) write.table( data.frame(x), 'top10_primary_avg.txt'  , append= T, sep="\t" ))
lapply(top.celline.list, function(x) write.table( data.frame(x), 'top10_celline_avg.txt'  , append= T, sep="\t" ))

lapply(down.primary.list, function(x) write.table( data.frame(x), 'bottom10_primary_avg.txt'  , append= T, sep="\t" ))
lapply(down.celline.list, function(x) write.table( data.frame(x), 'bottom10_celline_avg.txt'  , append= T, sep="\t" ))




####################################
#clean
#######################
#rm(list=ls(all=TRUE))
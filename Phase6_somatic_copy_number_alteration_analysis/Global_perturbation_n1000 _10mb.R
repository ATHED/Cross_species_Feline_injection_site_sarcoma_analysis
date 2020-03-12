##################################
#data import
#######################################
data.full <- read.table("beforePerturbation_SCNA_data_10mb.txt",
                         sep=",",
                         header=TRUE,
                         stringsAsFactors=FALSE,
                         quote="",
                         comment.char="#")

#####################################
#count the levels
##############################################
count.compress.factor <- factor(data.full$chr)
num.obj <- data.frame(table(data.full$chr))

index.less.four <- which(num.obj$Freq < 4)

###################################
#split into various data frames
##########################################
count.data.X <- split(data.full, data.full$chr)

########################################
#assign "avg" to all within 1mb window
####################################################

#install "zoo" package for roll apply
#install.packages("zoo")
library(zoo)

#set window size to be 10mb
window.size <- 10^7

#build a big list to record all the window size
width.big.list <- list()
for(i in 1:30){
  tmp <- as.vector(rep(0, nrow(count.data.X[[i]])))
  width.big.list[[i]] <- tmp
}

#calculate all window width for all the rest chromosome, which only takes 5 mins
#the size of loop decided by the number of levels of count.compress.factor
for(x in 1:30){
  for(i in 1:nrow(count.data.X[[x]])){
    for(j in i:nrow(count.data.X[[x]])){
      if(count.data.X[[x]][i,3] + window.size >= count.data.X[[x]][j,3]){
        width.big.list[[x]][i] <- width.big.list[[x]][i] + 1
      }
    }
  }
}

#add the window size on to the data
for(x in 1:30){
  count.data.X[[x]]$window_size <- width.big.list[[x]]
}

#######################################
#write output
#########################################

lapply(count.data.X, function(x) write.table( data.frame(x), 'stage1_SCNA_data_10mb.txt'  , append= T, sep=',' ))

############################################
#read as a data frame
#################################################
data.compress <- read.table("stage2_SCNA_data_10mb.txt",
                        sep=",",
                        header=TRUE,
                        stringsAsFactors=FALSE,
                        quote="",
                        comment.char="#")

data.clean <- data.compress[!(data.compress$window_size < 4), ]

#store the vector of actual avg readings
actual.avg.log2.primary <- as.vector(data.clean$avg_log2_primary)
actual.avg.log2.celline <- as.vector(data.clean$avg_log2_celline)

#####################################
#perturbation part
######################################################
num.col <- nrow(data.clean)
perturbation.avg.primary.df <- data.frame(matrix(nrow = 1000, ncol = num.col))
perturbation.avg.celline.df <- data.frame(matrix(nrow = 1000, ncol = num.col))

for (i in 1:1000){
  perturbation.avg.primary.df[i, ] <- sample(actual.avg.log2.primary, replace = F)
  perturbation.avg.celline.df[i, ] <- sample(actual.avg.log2.celline, replace = F)
}

empirical.p.value.df <- data.frame(matrix(nrow = num.col, ncol = 3))
empirical.p.value.df[, 1] <- data.clean$id
colnames(empirical.p.value.df) <- c("id", "empirical_p_value_primary", "empirical_p_value_celline")

#define n and r
n <- 1000
r.primary <- 0
r.celline <- 0
empircal.value.primary.temp <-0
empircal.value.celline.temp <-0

for(i in 1:num.col){
  r.primary <- sum(abs(perturbation.avg.primary.df[,i]) >= abs(actual.avg.log2.primary[i]), na.rm = T)
  r.celline <- sum(abs(perturbation.avg.celline.df[,i]) >= abs(actual.avg.log2.celline[i]), na.rm = T)
  empircal.value.primary.temp <- (r.primary + 1)/ (n + 1)
  empircal.value.celline.temp <- (r.celline + 1)/ (n + 1)
  
  empirical.p.value.df[i, 2] <- empircal.value.primary.temp
  empirical.p.value.df[i, 3] <- empircal.value.celline.temp
}

###########################################
#adjust p-values for multiple comparisons
##################################################
#empirical.p.value.df$adjust_p_value_primary <- p.adjust(empirical.p.value.df$empirical_p_value_primary, method = "BH", n = length(empirical.p.value.df$empirical_p_value_primary))
#empirical.p.value.df$adjust_p_value_celline <- p.adjust(empirical.p.value.df$empirical_p_value_celline, method = "BH", n = length(empirical.p.value.df$empirical_p_value_celline))


only.primary.p.value.df.compress <- empirical.p.value.df[which(empirical.p.value.df$empirical_p_value_primary <= 0.001), ]

data.clean.compress <- data.clean[,c(1, 5:6)]
data.id.list <- data.clean.compress[which(abs(data.clean.compress$log2_primary_cvn) >= 1),]

result.df <- only.primary.p.value.df.compress[which(only.primary.p.value.df.compress$id %in% data.id.list$id),]

#OR.empirical.p.value.df.compress <- empirical.p.value.df[which(empirical.p.value.df$empirical_p_value_primary <= 0.01 | empirical.p.value.df$empirical_p_value_celline <= 0.01), ]
#AND.empirical.p.value.df.compress <- empirical.p.value.df[which(empirical.p.value.df$empirical_p_value_primary <= 0.01 & empirical.p.value.df$empirical_p_value_celline <= 0.01), ]



###################################################
#out put
############################################################
#write.table(OR.empirical.p.value.df.compress ,
#            file="empirical_p_value_less_0.01_primaryORcelline_10mb.txt",
#            sep="\t",
#            row.names=FALSE,
#            col.names=TRUE,
#            quote=FALSE)

#write.table(AND.empirical.p.value.df.compress ,
#            file="empirical_p_value_less_0.01_primaryANDcelline_10mb.txt",
#            sep="\t",
#            row.names=FALSE,
#            col.names=TRUE,
#            quote=FALSE)

write.table(result.df,
            file="empirical_p_value_less_0.001_OnlyPrimaryAbove1_10mb.txt",
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

############################################
#visualization
#############################################
#install.packages("ggplot2")

library(ggplot2)

g.p.primary.all <- ggplot(empirical.p.value.df, aes(x=empirical_p_value_primary)) + geom_histogram(binwidth=.05, colour="black", fill="white")
pdf("p_value_primary_all.pdf")
print(g.p.primary.all)
dev.off()

g.p.celline.all <- ggplot(empirical.p.value.df, aes(x=empirical_p_value_celline)) + geom_histogram(binwidth=.05, colour="black", fill="white")
pdf("p_value_celline_all.pdf")
print(g.p.celline.all)
dev.off()

g.p.pANDc.celline <- ggplot(AND.empirical.p.value.df.compress, aes(x=empirical_p_value_celline)) + geom_density()
pdf("p_value_celline_both_less_0.01.pdf")
print(g.p.pANDc.celline)
dev.off()

g.p.pANDc.primary <- ggplot(AND.empirical.p.value.df.compress, aes(x=empirical_p_value_primary)) + geom_density()
pdf("p_value_primary_both_less_0.01.pdf")
print(g.p.pANDc.primary)
dev.off()



####################################
#clean
#######################
#rm(list=ls(all=TRUE))
##################################
#data import
#######################################
data.full <- read.table("../datafiles/Matched_10mb_windows_withBreen2009.tsv",
                        sep="\t",
                        header=TRUE,
                        stringsAsFactors=FALSE,
                        quote="",
                        comment.char="#")

########################################################
#separate into primary tissue and celline dataset
######################################################
#primary
p.data <- data.full[, c(2,6)]
#celline
c.data <- data.full[, c(2,8)]

###########################################################
#label Breen's up and down with addtional column
##############################################################
p.data <- na.omit(p.data)
p.data$Breen <- NA
c.data <- na.omit(c.data)
c.data$Breen <- NA

for (i in 1:nrow(p.data)){
  if (p.data[i,2] > 0){
    p.data[i,3] <- "Up"
  }else{
    p.data[i,3] <- "Down"
  }
}


for (i in 1:nrow(c.data)){
  if (c.data[i,2] > 0){
    c.data[i,3] <- "Up"
  }else{
    c.data[i,3] <- "Down"
  }
}

##################################
#print out
##################################
#######################
#librarys and install from github
##############################
library(devtools)
#install_github("ggbiplot", "vqv")

library(ggbiplot)


# Density plots
ggplot(p.data, aes(x=log2_primary_cvn, colour=Breen)) +
    geom_line(stat="density") +
    theme_minimal(base_size=10) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(colour="") +
    xlab(expression(log[2](cancer/normal))) +
    ggsave("Dist_vs_SCNA_primary_tissue.pdf", width=2.75, height=2)


ggplot(c.data, aes(x=log2_celline_cvn, colour=Breen)) +
    geom_line(stat="density") +
    theme_minimal(base_size=10) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(colour="") +
    xlab(expression(log[2](cancer/normal))) +
    ggsave("Dist_vs_SCNA_cell_lines.pdf", width=2.75, height=2)
    
####################################
#clean
#######################
#rm(list=ls(all=TRUE))

##################################
#import datasets
################################
og_matrix <- read.csv("Merged_FCH_Sarc_SR_data.csv",
                           sep=",",
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")


up_matrix <- og_matrix[which(og_matrix$Type == "UP"),]
mix_matrix <- og_matrix[which(og_matrix$Type == "MIX"),]
down_matrix <- og_matrix[which(og_matrix$Type == "DOWN"),]

#######################################################
#install.packages("tidyverse")
#install.packages("ggplot2", dependencies = TRUE)

#devtools::install_github("tidyverse/rlang", build_vignettes = TRUE, force = T)
#devtools::install_github("tidyverse/ggplot2",force = T)

library(ggplot2)
packageVersion("ggplot2")

####################################################
#build a new data matix
up_r_matrix <- matrix(, nrow = 3*nrow(up_matrix), ncol = 3)

for(i in 1:nrow(up_matrix)){
  t = 3*i - 2
  up_r_matrix[t:(t+2),1] <- up_matrix[i,1]
  up_r_matrix[t:(t+2),2] <- c("Feline", "Canine", "Human")
  up_r_matrix[t,3] <- log(up_matrix[i,2])
  up_r_matrix[t+1,3] <- log(up_matrix[i,3])
  up_r_matrix[t+2,3] <- log(up_matrix[i,4])
}
up.df <- data.frame(up_r_matrix)
colnames(up.df) <- c("Gene", "Animal", "Value")

####################################################################################
up.df$Value <- as.numeric(as.character(up.df$Value))
up.df$Animal <- factor(up.df$Animal, labels=c("Feline", "Canine", "Human"))
ggplot(up.df, aes(fill=Animal, x=Gene, y=Value)) +
    scale_fill_manual(values=c("white", "lightgrey","#404040")) +
    geom_bar(position="dodge",
             stat="identity", colour = "black") +
    coord_flip() +
    theme_minimal(base_size=12) +
    theme(axis.title.x=element_text(size=10),
          axis.text.x=element_text(size=10), 
          axis.ticks.x=element_blank()) +
    #ylab("percent") +
    ylab(expression(log[2](cancer/normal))) +
    labs(fill="") +
    ggsave("fig7_UP_group.pdf", width=9, height=8)

#######################################################

####################################################
#build a new data matix
down_r_matrix <- matrix(, nrow = 3*nrow(down_matrix), ncol = 3)

for(i in 1:nrow(down_matrix)){
  t = 3*i - 2
  down_r_matrix[t:(t+2),1] <- down_matrix[i,1]
  down_r_matrix[t:(t+2),2] <- c("Feline", "Canine", "Human")
  down_r_matrix[t,3] <- log(down_matrix[i,2])
  down_r_matrix[t+1,3] <- log(down_matrix[i,3])
  down_r_matrix[t+2,3] <- log(down_matrix[i,4])
}
down.df <- data.frame(down_r_matrix)
colnames(down.df) <- c("Gene", "Animal", "Value")

####################################################################################
down.df$Value <- as.numeric(as.character(down.df$Value))
down.df$Animal <- factor(down.df$Animal, labels=c("Feline", "Canine", "Human"))
ggplot(down.df, aes(fill=Animal, x=Gene, y=Value)) +
  scale_fill_manual(values=c("white", "lightgrey","#404040")) +
  geom_bar(position="dodge",
           stat="identity", colour = "black") +
  coord_flip() +
  theme_minimal(base_size=12) +
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=10), 
        axis.ticks.x=element_blank()) +
  #ylab("percent") +
  ylab(expression(log[2](cancer/normal))) +
  labs(fill="") +
  ggsave("fig7_DOWN_group.pdf", width=9, height=8)

####################################################
#build a new data matix
mix_r_matrix <- matrix(, nrow = 3*nrow(mix_matrix), ncol = 3)

for(i in 1:nrow(mix_matrix)){
  t = 3*i - 2
  mix_r_matrix[t:(t+2),1] <- mix_matrix[i,1]
  mix_r_matrix[t:(t+2),2] <- c("Feline", "Canine", "Human")
  mix_r_matrix[t,3] <- log(mix_matrix[i,2])
  mix_r_matrix[t+1,3] <- log(mix_matrix[i,3])
  mix_r_matrix[t+2,3] <- log(mix_matrix[i,4])
}
mix.df <- data.frame(mix_r_matrix)
colnames(mix.df) <- c("Gene", "Animal", "Value")

####################################################################################
mix.df$Value <- as.numeric(as.character(mix.df$Value))
mix.df$Animal <- factor(mix.df$Animal, labels=c("Feline", "Canine", "Human"))

mix.dat1 <- subset(mix.df, Value >= 0)
mix.dat2 <- subset(mix.df, Value < 0)

#ggplot() +
#  scale_fill_manual(values=c("white", "lightgrey","#404040")) +
#  geom_bar(data = mix.dat1, aes(fill=Animal, x=Gene, y=Value), position = position_dodge(preserve = "single"), 
#           stat="identity", colour = "black") +
#  geom_bar(data = mix.dat2, aes(fill=Animal, x=Gene, y=Value), position = position_dodge(preserve = "single"),
#           stat="identity", colour = "black") +
#  #ylab(expression(paste(log[2](cancer/normal)))) +
#  coord_flip() +
#  theme_minimal(base_size=12) +
#  theme(axis.title.x=element_text(size=10),
#        axis.text.x=element_text(size=10), 
#        axis.ticks.x=element_blank()) +
#  #ylab("percent") +
#  ylab(expression(log[2](cancer/normal))) +
#  labs(fill="") +
#  ggsave("fig7_MIX_group.pdf", width=9, height=8)

ggplot(data=mix.df, aes(y=Value, x=Gene)) +
  geom_bar(aes(fill=Animal), position=position_dodge(), stat="identity", colour="black") +
  scale_fill_manual(values=c("white", "lightgrey","#404040")) +
  theme_minimal(base_size=12) +
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=10), 
        axis.ticks.x=element_blank()) +    
  ylab(expression(log[2](cancer/normal))) +
  labs(fill="") +
  coord_flip() +
  ggsave("fig7_MIX_group.pdf", width=3.75, height=8)

####################################
#clean
#######################
#rm(list=ls(all=TRUE))

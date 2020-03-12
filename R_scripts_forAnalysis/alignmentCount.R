############################
#setup the working directory
###########################
getwd()
setwd("/Users/Weiqi0/Desktop/RNAseq_Cat")

###########################################
#read in the data
##########################################
phase1.data <- read.table ("phase1_counts.txt", sep="\t", header=T)
phase2.data <- read.table ("phase2_counts.txt", sep="\t", header=T)

################################################################
#data slice into useful part(remain 5th column for sanity check)
#################################################################
phase1.data.temp <- phase1.data[,6:10]
phase2.data.temp <- phase2.data[,6:16]

##############################################################
#count num of nonzero columns
############################################################
phase1.nonzero <- colSums(phase1.data.temp != 0)
phase1.nonzero
phase2.nonzero <- colSums(phase2.data.temp != 0)
phase2.nonzero

##################################################
#transform into precentage of nonzero columns
###################################################
phase1.nonzero.percent <- 100*phase1.nonzero / 21890
phase1.nonzero.percent
phase2.nonzero.percent <- 100*phase2.nonzero / 21890
phase2.nonzero.percent

####################################
#clean
#######################
#rm(list=ls(all=TRUE))
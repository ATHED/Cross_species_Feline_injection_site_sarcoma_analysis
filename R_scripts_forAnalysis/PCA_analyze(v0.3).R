#################################################
#data import
############################################

##################################
#primary.tissue
#################################
phase1AND2.data <- read.table ("counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt", sep="\t", header=T)
primary.tissue.phase1AND2.num.data <- phase1AND2.data[,30:35]

rownames(primary.tissue.phase1AND2.num.data) <- phase1AND2.data[,1]

condition_1 <- c("normal","cancer","cancer","cancer","normal","normal")

#ND2.num.data) <- c("normal","cancer","cancer","cancer","cancer","cancer","cancer","cancer","normal","normal","normal","normal","normal","normal")

primary.tissue.train.data <- t(as.matrix(primary.tissue.phase1AND2.num.data))

primary.tissue.df_f <- primary.tissue.train.data[,apply(primary.tissue.train.data, 2, var, na.rm=TRUE) > 2]

primary.tissue.pca_res_cat <- prcomp(primary.tissue.df_f,
                        center=TRUE)
names(primary.tissue.pca_res_cat)
print(primary.tissue.pca_res_cat)

###############################################
#cell line
#################################################
cell.line.phase1AND2.num.data <- phase1AND2.data[,42:49]

rownames(cell.line.phase1AND2.num.data) <- phase1AND2.data[,1]

condition_2 <- c("cancer","cancer","cancer","cancer","normal","normal","normal","normal")

#ND2.num.data) <- c("normal","cancer","cancer","cancer","cancer","cancer","cancer","cancer","normal","normal","normal","normal","normal","normal")

cell.line.train.data <- t(as.matrix(cell.line.phase1AND2.num.data))

cell.line.df_f <- cell.line.train.data[,apply(cell.line.train.data, 2, var, na.rm=TRUE) > 2]

cell.line.pca_res_cat <- prcomp(cell.line.df_f,
                                     center=TRUE)
#names(cell.line.pca_res_cat)
#print(cell.line.pca_res_cat)


##################################
#scatter graph
#####################################
#biplot(pca_res_cat, scale = 0)

#######################
#librarys and install from github
##############################
library(devtools)
#install_github("ggbiplot", "vqv")

library(ggbiplot)

############################
#primary.tissue plot
###################################
g.primary.tissue <- ggbiplot(primary.tissue.pca_res_cat, obs.scale = 1, var.scale = 1, 
              groups = condition_1, ellipse = F, 
              circle = T, var.axes = F, pc.biplot = TRUE, labels = rownames(primary.tissue.df_f), labels.size = 5)
g.primary.tissue <- g.primary.tissue + scale_color_discrete(name = '')
g.primary.tissue <- g.primary.tissue + theme(legend.direction = 'horizontal', 
               legend.position = 'top',legend.text = element_text(size=10), plot.title = element_text(size=10))

pdf("PCA_primary.tissue_Felis6.2(vLF).pdf")
print(g.primary.tissue)
dev.off()

###########################################
#cell line plot
########################################
g.cell.line <- ggbiplot(cell.line.pca_res_cat, obs.scale = 1, var.scale = 1, 
                             groups = condition_2, ellipse = F, 
                             circle = TRUE, var.axes = F, pc.biplot = TRUE, labels = rownames(cell.line.df_f), labels.size = 5)
g.cell.line <- g.cell.line + scale_color_discrete(name = '')
g.cell.line <- g.cell.line + theme(legend.direction = 'horizontal', 
                                             legend.position = 'top', legend.text = element_text(size=10), plot.title = element_text(size=10))

pdf("PCA_cell.line_Felis6.2(vLF).pdf")
print(g.cell.line)
dev.off()


####################################
#clean
#######################
#rm(list=ls(all=TRUE))
#################################################
#data import
############################################

##################################
#primary.tissue
#################################
phase1AND2.data <- read.table ("../datafiles/counts_annot_phase1AND2_norm_withLogRatio(v20170412).txt", sep="\t", header=T)
primary.tissue.phase1AND2.num.data <- phase1AND2.data[,30:35]

rownames(primary.tissue.phase1AND2.num.data) <- phase1AND2.data[,1]

condition_1 <- c("muscle","FISS","FISS","FISS","skin","skin")

#ND2.num.data) <- c("normal","cancer","cancer","cancer","cancer","cancer","cancer","cancer","normal","normal","normal","normal","normal","normal")

primary.tissue.train.data <- t(as.matrix(primary.tissue.phase1AND2.num.data))

primary.tissue.df_f <- primary.tissue.train.data[,apply(primary.tissue.train.data, 2, var, na.rm=TRUE) > 2]

primary.tissue.pca_res_cat <- prcomp(primary.tissue.df_f,
                        center=TRUE)
explained_var_primary_tissue <- summary(primary.tissue.pca_res_cat)$importance[2,]

library(ggplot2)
exp_var_prim_df <- data.frame(component=colnames(primary.tissue.pca_res_cat$x),
                         explained_variance=100*explained_var_primary_tissue)
ggplot(data=exp_var_prim_df, aes(component, explained_variance)) +
    geom_bar(stat="identity") +
    theme_minimal(base_size=10) +
    ylab("% explained variance") +
    ggsave("explained_variance_primary.pdf", width=2.0, height=2.0)

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

explained_var_cell_lines <- summary(cell.line.pca_res_cat)$importance[2,]

library(ggplot2)
exp_var_cells_df <- data.frame(component=colnames(cell.line.pca_res_cat$x),
                         explained_variance=100*explained_var_cell_lines)
ggplot(data=exp_var_cells_df[1:7,], aes(component, explained_variance)) +
    geom_bar(stat="identity") +
    theme_minimal(base_size=10) +
    ylab("% explained variance") +
    ggsave("explained_variance_celllines.pdf", width=2.25, height=2.0)

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
ggbiplot(primary.tissue.pca_res_cat, obs.scale = 1, var.scale = 1, 
                             groups = condition_1, ellipse = TRUE, 
                             circle = TRUE, var.axes = FALSE, pc.biplot = TRUE) +
    geom_point(aes(colour=condition_1), size=1) +
    scale_color_discrete(name = '') +
    theme_bw(base_size=10) +
    xlab(paste('PC1 (', as.character(round(100*explained_var_primary_tissue["PC1"], 0)), "% var. exp.)", sep='')) +
    ylab(paste('PC2 (', as.character(round(100*explained_var_primary_tissue["PC2"], 0)), "% var. exp.)", sep='')) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.key.width=unit(0.1, "inch"),
          legend.direction = 'horizontal', 
          legend.position = 'top',
          legend.text = element_text(size=9),
          plot.title = element_text(size=10)) +
    ggsave("pca-primarytissue.pdf", width=2, height=2)

###########################################
#cell line plot
########################################

ggbiplot(cell.line.pca_res_cat, obs.scale = 1, var.scale = 1, 
         groups = condition_2, ellipse = TRUE, 
         circle = TRUE, var.axes = FALSE, pc.biplot = TRUE) +
    geom_point(aes(colour=condition_2), size=1) +
    scale_color_discrete(name = '') +
    theme_bw(base_size=10) +
    xlab(paste('PC1 (', as.character(round(100*explained_var_cell_lines["PC1"], 0)), "% var. exp.)", sep='')) +
    ylab(paste('PC2 (', as.character(round(100*explained_var_cell_lines["PC2"], 0)), "% var. exp.)", sep='')) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.key.width=unit(0.1, "inch"),
          legend.direction = 'horizontal', 
          legend.position = 'top',
          legend.text = element_text(size=9)) +
    ggsave("pca-celllines.pdf", width=2, height=2)



####################################
#clean
#######################
#rm(list=ls(all=TRUE))

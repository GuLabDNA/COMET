#SuppleFig.7

library(dplyr) 
library(openxlsx)
library(tidyverse)
library(reshape2)
library(data.table) #detDT,fread
library(RSpectra)
library(Rtsne)
library(ggrepel)
library(ggpubr)

#do parallel
library(doParallel)  
library(foreach)
library(impute) #impute
#create folder
createFolder <- function(pathtemp,nameFolder){
  folder <- paste0(pathtemp, "/", nameFolder)
  
  if (file.exists(folder)) {
    cat("The folder already exists")
  } else {dir.create(folder)}
}

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

## load R codes
load(file=file.path(pathbf2,"Draft/Suppl_Data_Coding", "25-08-04-TCGA_xr.betas_60k_exGroups.Rdata"))

anno_train$TCGA_Project_mcf <- ifelse(anno_train$TCGA_Project_mcf=="CTRL (CSF)",
                                      "CTRL",(ifelse(anno_train$TCGA_Project_mcf=="STAD",
                                                     "MCF UGI-AD",anno_train$TCGA_Project_mcf)))
anno_train2 <- anno_train[!anno_train$TCGA_ID %in% c("BF3974XR26C","BF3615XR26C",
                                                     "BF3274XR25B","BF3311XR2326",
                                                     "BF3431XR26B","BF3676XR26C",
                                                     "BF1152LW5XR25A","BF3957XR26C",
                                                     "BF3488XR26B","BF3489XR26B","BF3388XR25B",
                                                     "BF1033vXR22B","BF1153LW6XR25A",
                                                     "BF3302XR25B","BF3905XR25B","BF3683XR21B",
                                                     "BF3400XR26B"),]
table(anno_train2$TCGA_Project_mcf)

p_filter=20000
beta_train_or <- beta_train_or[!rownames(beta_train_or) %in% 
                                 c("BF3974XR26C","BF3615XR26C","BF3274XR25B",
                                   "BF3311XR2326","BF3431XR26B","BF3676XR26C",
                                   "BF1152LW5XR25A","BF3957XR26C","BF3488XR26B",
                                   "BF3489XR26B","BF3388XR25B","BF1033vXR22B",
                                   "BF1153LW6XR25A","BF3302XR25B","BF3905XR25B",
                                   "BF3683XR21B","BF3400XR26B"),]

imp_df_filter <- beta_train_or[,order(-apply(beta_train_or,2,sd))[1:p_filter]]
rownames(imp_df_filter) <- gsub("_beta","",rownames(imp_df_filter))

or_index <- anno_train2$TCGA_ID
imp_df_filter_or <- imp_df_filter[or_index, ]

## Remove any NA
betas_trans <- imp_df_filter_or[,!colSums(!is.finite(imp_df_filter_or))]


rm(.Random.seed, envir=globalenv())  #remove seed

pca <- prcomp_svds(betas_trans,k=30)
#perplex_vec <- c(3,4,5,6,8,10)
perplex_vec=5;i=1

  rm(.Random.seed, envir=globalenv())  #remove seed
  set.seed(91235)
  
  # calculate tSNE
  res <- Rtsne(pca$x, pca=FALSE, max_iter=2500,theta=0.0,verbose=T, perplexity=perplex_vec[i])
  
  meth_class <- as.factor(anno_train2$TCGA_Project_mcf)
  levels(meth_class) #6
  # shape_class <- factor(anno_train$Sample_tp.2cat, order=T,
  #                       levels=c("BF","FFPE"), labels=c("BF","FFPE"))
  # levels(shape_class)
  
  tsne_plot.df <- 
    data.frame(x = res$Y[,1], y = res$Y[,2],
               cluster = factor(meth_class, labels=as.vector(levels(meth_class))),
               #shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
               Sentrix_ID = rownames(pca$x))
  
  #different color pattern        
  meth3 <- factor(tsne_plot.df$cluster, order=T,
                  levels=c("BRCA","CTRL","DLBC","LUAD","MCF LGI-AD","PAAD"))
  hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
  tsne_plot.df$hclust = meth3
  hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% dplyr::select(x, y) %>% 
    summarise(x=median(x), y=median(y))
  
  col_code <- c("BRCA"="#6e1e33","CTRL"="#77DD77","MCF LGI-AD"="#F0C507",
                "LUAD"="#b890e0","DLBC"="#e4703a",
                "MCF LGI-AD"="#F0C507","PAAD"="#27780D")
  
labeled_tsne1 <-
          ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=TRUE))+
          xlab("t-SNE 1") + ylab("t-SNE 2") +
          geom_point()+
          scale_x_continuous(breaks=seq(-300, 300, by = 100),
                             labels=seq(-300, 300, by = 100), limits=c(-301, 301), expand = c(0, 0)) +
          scale_y_continuous(breaks=seq(-300, 300, by = 100),
                             labels=seq(-300, 300, by = 100), limits=c(-301,301), expand = c(0, 0)) + 
          theme_bw() +
          geom_label_repel(aes(label = hclust), data=hc.norm.cent, 
                           label.size = 0.02, max.overlaps = Inf) + 
          scale_colour_manual(values = col_code, labels=meth3) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position="right")
labeled_tsne1

#######################################
#color by sample types
rm(.Random.seed, envir=globalenv())  #remove seed
set.seed(91235)
perpn=5
# calculate tSNE
res <- Rtsne(pca$x, pca=FALSE, max_iter=2500,theta=0.0,verbose=T, perplexity= perpn)

# meth_class <- as.factor(anno_train$TCGA_Project_subgrp)
# levels(meth_class) #3
meth_class <- factor(anno_train2$Sample_tp.3cat, order=T,
                     levels=c("ABDO","CSF","FFPE","FNA","PLEU"), 
                     labels=c("ABDO","CSF","FFPE","FNA","PLEU"))
levels(meth_class)

tsne_plot.df <- 
  data.frame(x = res$Y[,1], y = res$Y[,2],
             cluster = factor(meth_class, labels=as.vector(levels(meth_class))),
             #shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
             Sentrix_ID = rownames(pca$x))

#different color pattern        
meth3 <- factor(tsne_plot.df$cluster, order=T,
                levels=c("ABDO","CSF","FFPE","FNA","PLEU"))
hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
tsne_plot.df$hclust = meth3
hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% dplyr::select(x, y) %>% 
  summarise(x=median(x), y=median(y))

col_code <- c("ABDO"="#CC79A7","CSF"="#D55E00", "FFPE"="#0072B2",
              "FNA"="#F0E442", "PLEU"="#009E73")

labeled_tsne2 <-
        ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=F))+
        labs(x="t-SNE 1", y="t-SNE 2")+
        geom_point() + 
        #scale_shape_manual(values=c(19,17)) + 
        scale_x_continuous(breaks=seq(-300, 300, by = 100),
                           labels=seq(-300, 300, by = 100), limits=c(-301,301), expand = c(0, 0)) +
        scale_y_continuous(breaks=seq(-300, 300, by = 100),
                           labels=seq(-300, 300, by = 100), limits=c(-301,301), expand = c(0, 0)) + 
        theme_bw() +
        scale_colour_manual(values = col_code) +
        theme(legend.position="right")
labeled_tsne2

#######################################
#color by sex
anno_train2_sex <- anno_train2[!anno_train2$TCGA_ID %in% 
                                 c("BF1159LW15XR25A","BF1123XR25A"),]

imp_df_filter_sex <- imp_df_filter[!rownames(imp_df_filter) %in% 
                                     c("BF1159LW15XR25A","BF1123XR25A"),]

or_index <- anno_train2_sex$TCGA_ID
imp_df_filter_sex.or <- imp_df_filter_sex[or_index, ]

## Remove any NA
betas_trans <- imp_df_filter_or[,!colSums(!is.finite(imp_df_filter_or))]

rm(.Random.seed, envir=globalenv())  #remove seed
pca <- prcomp_svds(betas_trans,k=30)

rm(.Random.seed, envir=globalenv())  #remove seed
set.seed(91235)
perpn=5
# calculate tSNE
res <- Rtsne(pca$x, pca=FALSE, max_iter=2500,theta=0.0,verbose=T, perplexity= perpn)

# levels(meth_class) #3
meth_class <- factor(anno_train2$Sex, order=T,
                     levels=c("F","M"), 
                     labels=c("Female","Male"))
levels(meth_class)

tsne_plot.df <- 
  data.frame(x = res$Y[,1], y = res$Y[,2],
             cluster = factor(meth_class, labels=as.vector(levels(meth_class))),
             #shape.f = factor(shape_class, labels=as.vector(levels(shape_class))),
             Sentrix_ID = rownames(pca$x))

meth3 <- factor(tsne_plot.df$cluster, order=T, levels=c("Female","Male"))
hc.norm = hclust(dist(tsne_plot.df[,c(1,2)]))
tsne_plot.df$hclust = meth3
hc.norm.cent = tsne_plot.df %>% group_by(hclust) %>% dplyr::select(x, y) %>% 
  summarise(x=median(x), y=median(y))

col_code <- c("Female"="#E47F7B","Male"="#2072B3")

labeled_tsne3 <-
        ggplot(data=tsne_plot.df, aes(x=x, y=y, color = hclust, label=F))+
        labs(x="t-SNE 1", y="t-SNE 2")+
        geom_point() + 
        #scale_shape_manual(values=c(19,17)) + 
        scale_x_continuous(breaks=seq(-300, 300, by = 100),
                           labels=seq(-300, 300, by = 100), limits=c(-301,301), expand = c(0, 0)) +
        scale_y_continuous(breaks=seq(-300, 300, by = 100),
                           labels=seq(-300, 300, by = 100), limits=c(-301,301), expand = c(0, 0)) + 
        theme_bw() +
        scale_colour_manual(values = col_code) +
        theme(legend.position="right")
labeled_tsne3
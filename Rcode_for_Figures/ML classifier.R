# ML classifier
  
# 1.cross-validation on your training set: (no need here)
# This process will involve splitting your training set into 'k' subsets, training your model on 'k-1' subsets, and then testing the model on the left out '1' subset. This gives us an unbiased estimate of model evaluation metric that helps in understanding how well the model is likely to perform on unseen data.
# 
# 2.Train the final model on the full training set: 
#   You take the model settings (like hyperparameters) that performed the best on average during the cross-validation process and use those settings to train a model on the full training set.
# 
# 3.Validate the final model on your external test set: 
#   You take the final model from step 2, use it to predict the target variable on your external test set, and then compare the predictions with the actual values in your test set to calculate performance metrics.


# load packages
library(openxlsx)
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table) # setDT
library(tibble) #tibble::rownames_to_column
library(randomForest)
library(glmnet)
library(ranger)
library(writexl)
#install.packages("Matrix")

#do parallel
library(parallel)
library(doParallel)
library(doMC)
library(HandTill2001) # Multiple Class Area under ROC Curve

### function cbind.fill, combine cols with different lengths
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# Function to get column names with highest values
getTopThree <- function(df) {
  
  # Apply a function to each row of the dataframe
  topThree <- apply(df, 1, function(row) {
    # Sort the row in decreasing order to get the top 3 values
    top_index <- order(-row)[1:3]
    # Extract the top 3 values
    values <- row[top_index]
    # Get the names of the columns containing the top 3 values
    names <- names(df)[top_index]
    # Return a list containing the names and values
    data.frame(top1_cal_label=names[1],top1_cal_score=values[1],
               top2_cal_label=names[2],top2_cal_score=values[2],
               top3_cal_label=names[3],top3_cal_score=values[3])
    
  })
  return(topThree)
}

# load functions
pathmeth <- "D:/Jingru/methylation"
pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

source(file.path(pathmeth,"mnp_training-master/R","makefolds.R"))
source(file.path(pathmeth,"mnp_training-master/R","train.R"))
source(file.path(pathmeth,"mnp_training-master/R","calculateCVfold.R"))
source(file.path(pathmeth,"mnp_training-master/R","batchadjust.R"))
source(file.path(pathmeth,"mnp_training-master/R","MNPprocessIDAT_functions.R"))

#external validation set, our real samples
options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999999999)

folds <- 5
ntrees <- 500  
cores <- 20
seed <- 76803149
p <- 5000 
p_filter1 <- 15000 

# Train the final model on the full training set
rm(.Random.seed, envir=globalenv())

#1. load betas of training set
load(file=file.path(pathbf2,"Draft/Suppl_Data_Coding", 
                    "24-10-10-TCGA.sub.array_hg38_betas_60k_t250_TCGA.MESO.SarCtrl.Rdata"))
#TCGA_hg38_betas.array, anno_TCGA.array

#intersect with XR covered CpGs
load(file=file.path(pathref,"diffmeth",
                    "array.450k.hg38_XRfixed.interBlocks_250328.RData")) #261071

#filter probes
CG_df <- as.data.frame(colnames(TCGA_hg38_betas.array))
CG_df$flag="df"
colnames(CG_df)[1]="CG"
CG_df_m <- merge(CG_df,array.450k.hg38_XRfixed.interBlocks,by="CG") #21529

TCGA_hg38_betas.array_filter <- TCGA_hg38_betas.array[,CG_df_m$CG]

anno_TCGA.array$TCGA_Project_mcf <-anno_TCGA.array$TCGA_Project
table(anno_TCGA.array$TCGA_Project_mcf)

anno_TCGA.array <- within(anno_TCGA.array,{
  TCGA_Project_mcf <- NA
  TCGA_Project_mcf[TCGA_Project %in% c("COAD","READ")] <- "MCF LGI-AD"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("STAD","ESCA_AD")] <- "MCF UGI-AD"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("ESCA_SCC","HNSC","LUSC")] <- "MCF SCC"
  TCGA_Project_mcf[TCGA_Project %in% c("CTRL (BF)","CTRL (BLOOD)","CTRL (MUS)",
                                       "CTRL (REA)")] <- "Control"
  TCGA_Project_mcf[TCGA_Project_subgrp == "Low_fraction"] <- "Low_fraction"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("OV","UCEC_AD","UCEC_Other","UCS")] <- "MCF GYN"
  TCGA_Project_mcf[TCGA_Project_subgrp == "PNET"] <- "pNET"
  TCGA_Project_mcf[TCGA_Project_subgrp == "CESC_AD"] <- "CESC_AD"
  TCGA_Project_mcf[TCGA_Project_subgrp == "CESC_SCC"] <- "CESC_SCC"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("LGG","GBM")] <- "MCF GBM"
  TCGA_Project_mcf[TCGA_Project_subgrp == "TGCT_S"] <- "TGCT_S"
  TCGA_Project_mcf[TCGA_Project_subgrp == "TGCT_NS"] <- "TGCT_NS"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("BRCA_L","BRCA_B","BRCA")] <- "BRCA"
  
})

anno_TCGA.array$TCGA_Project_mcf <- ifelse(is.na(anno_TCGA.array$TCGA_Project_mcf),
                                           anno_TCGA.array$TCGA_Project,
                                           anno_TCGA.array$TCGA_Project_mcf)
table(anno_TCGA.array$TCGA_Project_mcf,useNA="ifany")


#2. load betas of test set
load(file=file.path(pathbf2,"Draft/Suppl_Data_Coding", 
                    "24-10-10-TCGA.sub.array.test_hg38_betas_60k_t250_TCGA.MESO.SarCtrl.Rdata"))
#TCGA_hg38_betas.sub.test, anno_TCGA_sub.test
TCGA_betas.test <- TCGA_hg38_betas.sub.test[,CG_df_m$CG]

TCGA_test_trans <- as.data.frame(t(TCGA_betas.test))
TCGA_test_trans$CG=colnames(TCGA_betas.test)

anno_TCGA_sub.test <- within(anno_TCGA_sub.test,{
  TCGA_Project_mcf <- NA
  TCGA_Project_mcf[TCGA_Project %in% c("COAD","READ")] <- "MCF LGI-AD"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("STAD","ESCA_AD")] <- "MCF UGI-AD"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("ESCA_SCC","HNSC","LUSC","SCC")] <- "MCF SCC"
  TCGA_Project_mcf[TCGA_Project %in% c("CTRL (BF)","CTRL (BLOOD)","CTRL (MUS)",
                                       "CTRL (REA)")] <- "Control"
  TCGA_Project_mcf[TCGA_Project_subgrp == "Low_fraction"] <- "Low_fraction"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("OV","UCEC_AD","UCEC_Other","UCS")] <- "MCF GYN"
  TCGA_Project_mcf[TCGA_Project_subgrp == "PNET"] <- "pNET"
  TCGA_Project_mcf[TCGA_Project_subgrp == "CESC_AD"] <- "CESC_AD"
  TCGA_Project_mcf[TCGA_Project_subgrp == "CESC_SCC"] <- "CESC_SCC"
  TCGA_Project_mcf[TCGA_Project_subgrp == "LUAD"] <- "LUAD"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("LGG","GBM")] <- "MCF GBM"
  TCGA_Project_mcf[TCGA_Project_subgrp == "TGCT_S"] <- "TGCT_S"
  TCGA_Project_mcf[TCGA_Project_subgrp == "TGCT_NS"] <- "TGCT_NS"
  TCGA_Project_mcf[TCGA_Project_subgrp %in% c("BRCA_L","BRCA_B","BRCA")] <- "BRCA"
})

anno_TCGA_sub.test$TCGA_Project_mcf <- ifelse(is.na(anno_TCGA_sub.test$TCGA_Project_mcf),
                                              anno_TCGA_sub.test$TCGA_Project,
                                              anno_TCGA_sub.test$TCGA_Project_mcf)
table(anno_TCGA_sub.test$TCGA_Project_mcf,useNA="ifany")

#======================
#3. combine samples
load(file=file.path(pathbf2,"Draft/Suppl_Data_Coding", 
                    "25-07-23-allsamp_TCGA_betas.ImpByTest.k5_trans.RData")) #use all BFs
#allsamp_betas_TCGA_trans

samp1 <- openxlsx::read.xlsx(file.path(pathbf2,"Draft",
                                       "25-1202-Supple.Tables.xlsx"),
                             sheet=2, skipEmptyRows=FALSE, colNames = TRUE)
samp1 <- samp1 %>% 
  rename(sample_tp = `Sample.type(PLEU,ABDO,FNA,CSF)`,
         composite_score2 = composite_score,
         composite_top1match = `composite_top1match(comment(9,.unclear.diag.or.no.refs.or.na;.2,.indeterminate;.3.organ.not.subtypes)`)
samp1 <- samp1[samp1$sample_tp!='CSF',]

samp2 <- samp1[samp1$Tumor_Purity_ichor>=0.495 & samp1$TCGA_Project %in%
                 c("ACC","BLCA","BRCA","CHOL","DLBC","GBM","KIRC","KIRP",
                   "LAML","LIHC","MCF GYN","LUAD","MCF LGI-AD","MCF UGI-AD",
                   "MESO","PAAD","PAAD/CHOL","pNET","PRAD","SCC","SKCM","TGCT_NS",
                   "TGCT_S","THCA","UVM") 
               & samp1$Total.Reads>=10000000
               & !is.na(samp1$`tSNE(Reference.TCGA.Tumors)`),]
samp_co <- samp1[samp1$TCGA_Project=="CTRL"& samp1$Total.Reads>=10000000
                 & !is.na(samp1$`tSNE(Reference.TCGA.Tumors)`),]
samp_coca <- unique(rbind(samp2,samp_co))

samp2_tcga <- 
  samp_coca[,c("TCGA_ID","sample_tp","TCGA_Project","Tumor_Purity_ichor")] %>%  
  mutate(Sentrix_ID=paste0(TCGA_ID,"-5cov.TCGA60k_original"),
         Diag_TCGA_mcf=ifelse(TCGA_Project=="CTRL","Control",TCGA_Project))

colnames(samp2_tcga) <- c("TCGA_ID","material","TCGA_Project_mcf","purity.abs",
                          "Sentrix_ID","TCGA_Project")
samp2_tcga$Sample_type <- samp2_tcga$material
samp2_tcga$material <- 'cfDNA'

anno_add <- read.xlsx(file.path(pathbf2,"Draft/Suppl_Data_Coding", 
                                     "MockSamples_for_MLclassifier.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
anno_add$TCGA_Project_mcf <- anno_add$TCGA_Project
anno_add$material <- "gDNA"
anno_addm <- rbind(samp2_tcga,anno_add)

############# load Sample Set 2 #############

# load(file=file.path(pathbf2,"Draft/Suppl_Data_Coding", 
#                      "25-10-17-valid.allsamp.106CSFs_TCGA_betas.ImpByTest.k5_trans.RData")) # use CSF
# 
# samp_csf <- read.xlsx(file.path(pathbf2,"sample_info",
#                               "test_validationCSF_list_improve.xlsx"),
#                       sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
# samp_csf$Sentrix_ID <- paste0(samp_csf$TCGA_ID,"-10cov.TCGA60k_original")
# samp_csf2 <- samp_csf[c("TCGA_ID","TCGA_Project","Diag_TCGA_mcf",
#                         "Tumor_Purity_ichor","Sentrix_ID","sample_tp")]
# 
# samp_csf.co <- read.xlsx(file.path(pathbf2,"sample_info",
#                                 "BF_ctrls.xlsx"),
#                       sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
# 
# samp_csf3 <- samp_csf2[(samp_csf2$Tumor_Purity_ichor>0.495 |
#                         samp_csf2$Tumor_Purity_ichor<0.05)
#                        & !samp_csf2$TCGA_ID %in% samp_csf.co$Sample_ID, ] %>%
#   mutate(material="cfDNA") %>% 
#   rename(purity.abs = Tumor_Purity_ichor,
#          TCGA_Project_mcf = Diag_TCGA_mcf,
#          Sample_type=sample_tp)
# 
# anno_add <- read.xlsx(file.path(pathbf2,"Draft/Suppl_Data_Coding", 
#                                 "MockSamples_for_MLclassifier.xlsx"),
#                       sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
# anno_add$TCGA_Project_mcf <- anno_add$TCGA_Project
# anno_add$material <- "gDNA"
# anno_addm <- rbind(samp_csf3,anno_add)

#######################################
#4. combine our sample with the test set, including annotation file

anno_test2 <- anno_TCGA_sub.test %>% 
  select(TCGA_ID,Sentrix_ID,material,TCGA_Project,TCGA_Project_mcf,purity.abs)
table(anno_test2$TCGA_Project_mcf)

anno_test2$Sample_type<-"Tissue"

anno_test_final <- rbind(anno_test2,anno_addm)
anno_test_final <- 
  anno_test_final[,-which(names(anno_test_final) %in% "Sample_type")]

#5. get the whole dataset, and extract training and test set
anno_TCGA.array$or<- 1:nrow(anno_TCGA.array)
anno_test_final$or <- 
  (nrow(anno_TCGA.array)+1):(nrow(anno_TCGA.array)+nrow(anno_test_final))

anno_train.test <- 
  rbind(anno_TCGA.array[,-which(names(anno_TCGA.array) %in% "TCGA_Project_subgrp")],
        anno_test_final)
anno_train.test <- anno_train.test[order(anno_train.test$or),]
y <- anno_train.test$TCGA_Project_mcf

row_keep <- as.vector(anno_addm$Sentrix_ID)
addm_selected <- allsamp_betas_TCGA_trans[row_keep, ]

df <- rbind(TCGA_hg38_betas.array_filter,TCGA_betas.test,addm_selected) %>%
  mutate_if(is.character, as.numeric) # Remove NAs if any

# sd pre filtering to 20k probes, to speed up the example, by rows instead of columns
df_filter <- df[,order(-apply(df,2,sd))[1:p_filter1]]

or_index <- match(anno_train.test$Sentrix_ID, rownames(df_filter))

# Order the data frame by the index
df_filter <- df_filter[or_index, ]

indextrain <- 1:nrow(TCGA_hg38_betas.array_filter)
x.train <- df_filter[indextrain, ]
x.test <- df_filter[-indextrain, ]
y.train <- as.factor(y[indextrain])
y.test <- as.factor(y[-indextrain])

message("performing variable selection ...",Sys.time())
message("cores: ",cores)
message("ntrees: ",ntrees)  
message("n: ",nrow(x.train))
message("p: ",ncol(x.train))  

load(file.path(pathmeth,"TCGA/training/CV_seq_addBFCtrl_mcf","nfolds.RData"))

scores <- list()
idx <- list()
for(i in 1:length(nfolds)){
  fname <- paste0("CVfold.",i,".",0,".RData")
  load(file.path(pathmeth,"TCGA/training/CV_seq_addBFCtrl_mcf",fname))
  scores[[i]] <- rf.scores
  idx[[i]] <- nfolds[[i]][[1]][[1]]$test
}
scores <- do.call(rbind,scores)

idx <- unlist(idx)
y.cv <- anno_TCGA.array$TCGA_Project_mcf[idx]   


i="T9.413BFs";seed=62838125

# FUNCTION RFres_TCGA
RFres_TCGA <- function(i, seed){
  rm(.Random.seed, envir=globalenv())
  
  #5. Train the random forest model on the training set
  rf.varsel <- rfp(x=x.train, y=as.factor(y.train),
                   mc=cores, ntree=ntrees,
                   sampsize=rep(min(table(as.factor(y.train))),
                                length(table(as.factor(y.train)))),
                   importance=TRUE)
  
  # get permutation variable importance
  imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
  
  # reduce data matrix
  or <- order(imp.meandecrease,decreasing=T)
  
  rf.pred <- randomForest(y.train ~ ., data=x.train[,or[1:p]],
                          ntree=ntrees,
                          strata=y.train,
                          mtry=sqrt(ncol(x.train)),
                          sampsize=rep(min(table(y.train)),length(table(y.train))),
                          proximity=TRUE,
                          oob.prox=TRUE,
                          importance=TRUE,
                          keep.inbag=TRUE,
                          do.trace=FALSE,
                          seed=seed)

  #6. Validate the final model on your external test set
  # Calibration, use the RF-cv scores of training set
  rm(.Random.seed, envir=globalenv())
  set.seed(seed,"L'Ecuyer")
  
  cl <- makeCluster( 20 )  
  registerDoParallel(cl) 
  
  # fit multinomial logistic ridge regression model
  suppressWarnings(cv.calfit <- 
                     cv.glmnet(y=y.cv, x=scores, family="multinomial",
                               type.measure="mse", alpha=0, nlambda=100,
                               lambda.min.ratio=10^-6, parallel=TRUE))
  stopCluster(cl)    
  
  #7. Predict on the test set using the trained model
  rf.scores <- predict(rf.pred, newdata=x.test[,colnames(x.train[,or[1:p]])],type="prob")
  #type = c("response", "terms"), 
  ys <- colnames(rf.scores)[apply(rf.scores,1,which.max)]
  ys.score <- as.data.frame(apply(rf.scores, 1, function(x) max(x, na.rm = TRUE)))
  
  err <- sum(ys!=y.test)/nrow(x.test)
  message("misclassification error: ",err)
  
  #8. refit the cv.calfit model to get the calibrated score
  message("calibrating raw scores",Sys.time())
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response",
                   s=cv.calfit$lambda.1se)[,,1] 
  # use lambda estimated by 10fold CVlambda
  
  # Compute accuracy
  yp <- colnames(probs)[apply(probs,1,which.max)]
  yp.score <- as.data.frame(apply(probs, 1, function(x) max(x, na.rm = TRUE)))
  
  err2 <- sum(yp!=y.test)/nrow(anno_test_final)
  print(paste("err2: ", round(err2, 4)))
  
  tops <- getTopThree(as.data.frame(probs))
  
  tops_df <- do.call(rbind, lapply(tops, data.frame, stringsAsFactors = FALSE))
  tops_df <- data.frame(tops_df)
  tops_df$Sentrix_ID <- rownames(tops_df)
  
  m2 <- merge(anno_train.test,tops_df,by="Sentrix_ID",all.y=T)
  score_compare1<- cbind(rownames(probs),ys,ys.score)
  colnames(score_compare1) <- c("Sentrix_ID","pred_orig_label","pred_orig_score")
  score_compare1$grp=c(rep("Array", times=nrow(anno_TCGA_sub.test)),
                       c(rep("FLEXseq", times=nrow(anno_addm))))
  #score_compare1$grp=c(rep("XR-methylSeq", times=nrow(anno_add_fil)))
  score_compare2 <- merge(m2,score_compare1,by="Sentrix_ID")
  ca_cal <- subset(score_compare2,grp=="FLEXseq")
  ca_cal$ref="TCGA"
  ca_cal$test=i
  
  errs <- sum(y.test!=ys)/nrow(anno_test_final)
  errp <- sum(y.test!=yp)/nrow(anno_test_final)
  message("overall misclassification error scores: ",errs) 
  message("overall misclassification error calibrated: ",errp) 

  #9. Save results
  save(rf.scores,probs,y,ys,yp,score_compare2,errs,errp,
       file=file.path(pathbf2,"results_df",
                      paste0(i,"-validation.cv.score.k5_136BFs_TCGA.XR.RData")))
  
  write_xlsx(score_compare2, file.path(pathbf2,"results_df",
                                       paste0(i,"-score.comp.k5_TCGA_136BFs_TCGA.XR.xlsx")),
             col_names = TRUE,format_headers = TRUE)
  
  
  write_xlsx(ca_cal, file.path(pathbf2,"results_df",
                               paste0(i,"-score.comp.k5_TCGA_onlySamples.136BFs_TCGA.XR.xlsx")),
             col_names = TRUE,format_headers = TRUE)
  
  pred_mat <- as.data.frame(probs)
  pred_mat$Sentrix_ID <- rownames(pred_mat)
  write_xlsx(pred_mat, file.path(pathbf2,"results_df",
                                 paste0(i,"-score.comp.k5_TCGA_136BFs_matrix.xlsx")),
             col_names = TRUE,format_headers = TRUE)
  #score_compare has the calibrated and original score
  #yp is the calibrated predicted groups, ys is the original predicted groups
  return(ca_cal)
}

#Fig.6, SupplFig.9, SupplFig.10

######* load packages *######
library(tidyr)
library(dplyr)
library(reshape2)
library(stringr)
library(openxlsx)
library(writexl) #save x
library(readr)
library(data.table) #setDT
library(scales)
library(forcats) 
library(purrr) # map function
library(rstatix)

library(ggplot2)
library(ggbreak) #scale_y_break
library(ggrepel) #geom_text_repel
library(ggpubr) #ggarrange
library(RColorBrewer)
library(ggbeeswarm)
library(gridExtra)
library(grid)
library(pROC)

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"


# Make confusion matrix

res_diag <- read.xlsx(file.path(pathbf2,"Draft","25-1202-Supple.Tables.xlsx"),
                      sheet=2, skipEmptyRows=FALSE, colNames = TRUE)
res_diag <- res_diag %>% 
  rename(sample_tp = `Sample.type(PLEU,ABDO,FNA,CSF)`,
         composite_score2 = composite_score,
         composite_top1match = `composite_top1match(comment(9,.unclear.diag.or.no.refs.or.na;.2,.indeterminate;.3.organ.not.subtypes)`)

res_diag_comb1 <- within(res_diag, {
  ML_if_match.f <- NA
  ML_if_match.f[`ML_if_match(score0.5)` == 0] <- "Not matched"
  ML_if_match.f[`ML_if_match(score0.5)`  == 1] <- "Matched"
  ML_if_match.f[`ML_if_match(score0.5)`  == 2] <- "Indeterminate"
  ML_if_match.f[`ML_if_match(score0.5)`  == 9] <- "N/a"
  ML_if_match.f[`ML_if_match(score0.5)`  == 8] <- "Reads<15M"
  
  deconv_if_match.f <- NA
  deconv_if_match.f[deconv_if_match == 0] <- "Not matched"
  deconv_if_match.f[deconv_if_match == 1] <- "Matched"
  deconv_if_match.f[deconv_if_match == 2] <- "Indeterminate"
  deconv_if_match.f[deconv_if_match == 3] <- "Organ Matched"
  deconv_if_match.f[deconv_if_match == 9] <- "N/a"
  
  tSNE_match.f <- NA
  tSNE_match.f[tSNE_match == 0] <- "Not matched"
  tSNE_match.f[tSNE_match == 1] <- "Matched"
  tSNE_match.f[tSNE_match == 9] <- "Indeterminate"
  
  composite_top1match.f<- NA
  composite_top1match.f[composite_top1match == 0] <- "Not matched"
  composite_top1match.f[composite_top1match == 1] <- "Matched"
  composite_top1match.f[composite_top1match == 2] <- "Indeterminate"
  composite_top1match.f[composite_top1match == 3] <- "Organ Matched"
  composite_top1match.f[composite_top1match == 9] <- "N/a"
  
  composite_top2match.f<- NA
  composite_top2match.f[`composite_top2match`  == 0] <- "Not matched"
  composite_top2match.f[`composite_top2match`  == 1] <- "Matched"
  composite_top2match.f[`composite_top2match`  == 2] <- "Indeterminate"
  composite_top2match.f[`composite_top2match`  == 3] <- "Organ Matched"
  composite_top2match.f[`composite_top2match`  == 9] <- "N/a"
  
  cyto_catx2 <- NA
  cyto_catx2[cyto_catx == "Atypical"] <- "Atypical"
  cyto_catx2[cyto_catx == "Suspicious"] <- "Suspicious"
  cyto_catx2[cyto_catx == "Negative"] <- "Negative"
  cyto_catx2[cyto_catx == "Positive"] <- "Malignant"
  cyto_catx2[cyto_catx == "Unknown"] <- "Non-diagnostic"
  
  diag_caco <- NA
  diag_caco[TCGA_Project %in% c("Unknown","Unclear")] <- "Unclear"
  diag_caco[TCGA_Project %in% c("CTRL","Control")] <- "Control"  
  diag_caco[TCGA_Project == "Benign"] <- "Benign"  
  
  analy_grp <- NA
  analy_grp[Tumor_Purity_ichor >= 0.495 &
              Total.Reads>=15000000 & !is.na(top1_cal_label.MCF)] <- "analy_ca_ML"
  analy_grp[TCGA_Project=="CTRL"] <- "analy_co_ML"
  
  agen <- as.numeric(Age)
})

res_diag_comb1$analy_grp <- 
  ifelse(is.na(res_diag_comb1$analy_grp),"analy_deconv",res_diag_comb1$analy_grp)
res_diag_comb1$diag_caco <- 
  ifelse(is.na(res_diag_comb1$diag_caco),"Cancers",res_diag_comb1$diag_caco)

res_diag_comb1 <- within(res_diag_comb1, {
  purity_catn <- NA
  purity_catn[Tumor_Purity_ichor >=0.495] <- 1
  purity_catn[Tumor_Purity_ichor>=0.05 & Tumor_Purity_ichor <0.495] <- 2
  purity_catn[Tumor_Purity_ichor<0.05 & !is.na(Tumor_Purity_ichor)] <- 3
  
  quant_dna <- as.numeric(`Quant.DNA(ng/uL)`)
  
  dx_cat.f = factor(as.factor(diag_caco),order=T,
                    levels=c("Benign","Cancers", "Control","Unclear"),
                    labels=c("Benign","Cancers", "Control","Unclear")) 
  
  idn = as.numeric(gsub('BF', '', BF_id))
  idn2 <- idn
  idn2[BF_id=="BF1152LW5"] <- 1152
  idn2[BF_id=="BF1153LW6"] <- 1153
  idn2[BF_id=="BF1158LW14"] <- 1158
  idn2[BF_id=="BF1159LW15"] <- 1159
  idn2[BF_id=="BF1168LW25"] <- 1168
  idn2[BF_id=="BF1169LW26"] <- 1169
  idn2[BF_id=="BF3713dil10ng"] <- 3713
})

res_diag_comb1$idn3 <- 
  as.numeric(ifelse(str_sub(res_diag_comb1$BF_id, -1, -1)=='A',
                    gsub('BF', '', substr(res_diag_comb1$BF_id,1,
                                          nchar(res_diag_comb1$BF_id) - 1)), res_diag_comb1$idn2))
table(res_diag_comb1$`Sample.type(PLEU,ABDO,FNA,CSF)`)

res_diag3_filter <- res_diag_comb1[res_diag_comb1$Total.Reads>=10000000,] #501

res_diag3_filter <- 
  res_diag3_filter[!(res_diag3_filter$TCGA_Project %in% 
    c("ACC","Adeno/carcinoma", "AFX/PDS", "AS", "Benign", "CHOL", 
      "FMS","GIST", "Haem", "LAML", "LMO",
      "LMS", "MALIGNANT", "MECA", "MESO", "Myeloma", 
      "NET", "PRAD", "SARC", "SARC (MPNST-like)", "SKCM",
      "SWN", "SYSA", "TGCT_NS", "TGCT_S", "THCA/SARC","MDS",
      "ThoraxAD", "Unclear","Unknown", "UVM", "YST",
      "NET in esophagus","NET in lung","NET in small bowel","pNET",
      "MCF IDH GLM", "MCF MB G3G4") & 
       (res_diag3_filter$Tumor_Purity_ichor < 0.495 & 
       res_diag3_filter$Tumor_Purity_ichor>=0.05)),]
#NET in esophagus has two results: SCC/MCG UGI
#NET in lung has two results: SCC/LUAD

res_diag3_filter <- 
  res_diag3_filter[!(res_diag3_filter$TCGA_Project %in% 
    c("Adeno/carcinoma","Haem","NET","NET in esophagus","NET in lung",
      "NET in small bowel","SBAD","SYSA","ThoraxAD","Unclear","Unknown",
      "MCF IDH GLM", "MCF MB G3G4") & 
       res_diag3_filter$Tumor_Purity_ichor >= 0.495),]

#remove FNA of tumor fraction <50% and all types <5%
res_diagm3_filter_rmFNA <- res_diag3_filter %>% 
  mutate(composite_top2match.f=ifelse(composite_top2match.f=="Organ Matched",
                                      "Matched", composite_top2match.f),
         composite_class2=ifelse(`composite_class(score0.5)` %in% c("BRCA_B", "BRCA_L"),"BRCA",
                                 (ifelse(`composite_class(score0.5)` %in% c("Head-Neck","LUSC","SCC"),"MCF SCC",
                                         (ifelse(`composite_class(score0.5)`=="SBAD","MCF LGI-AD",
                                                 `composite_class(score0.5)`))))),
         TCGA_Project=ifelse(TCGA_Project=="LUAD/LUSC","LUAD",
                             (ifelse(TCGA_Project=="SBAD","MCF LGI-AD",
                                     (ifelse(TCGA_Project=="NET in esophagus","Head-Neck",
                                             (ifelse(TCGA_Project=="KIRC","Kidney",
                                                     (ifelse(TCGA_Project=="SCC","MCF SCC",TCGA_Project)))))))))) %>% 
  filter(sample_tp %in% c("ABDO","PLEU","CSF") | 
           (sample_tp=="FNA" & Tumor_Purity_ichor>0.495) | Tumor_Purity_ichor<0.05) 

######* whole prediction matrix for ROC curve *######

#load data
comp_mat <- read.xlsx(file.path(pathbf2,"Draft/Suppl_Data_Coding","allBF.CSFs_MLmatrix_dfGT05.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

# comp_mat <- read.xlsx(file.path(pathbf2,""Draft/Suppl_Data_Coding","allBF.CSFs_MLmatrix_dfall.xlsx"),
#                       sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

# Assuming res_matrix and classes are already defined
classes <- c("BRCA", "DLBC", "LUAD", "MCF.GYN", "MCF.LGI-AD", "MCF.UGI-AD", "MCF.SCC", "PAAD")
# Create empty vectors to store AUC values and colors for the legend
auc_values <- numeric(length(classes))
auc_ci_lower <- numeric(length(classes))
auc_ci_upper <- numeric(length(classes))
line_colors <- c("#6e1e33", "#e4703a", "#b890e0", "#0581be", "#F0C507", 
                 "#63632f", "#27780d", "#734ba9")

pdf(file=file.path(pathbf2,"Draft/Suppl_Data_Coding",
                   "allBF.CSFs_CompositeClass.roc_8tumors.GT05.pdf"),
    width = 5.5, height = 5.5, useDingbats = FALSE)

# pdf(file=file.path(pathbf2,"Draft/Suppl_Data_Coding",
#     "allBF.CSFs_CompositeClass.roc_8tumors.all.pdf")),
#     width = 5.5, height = 5.5, useDingbats = FALSE)

# Loop through each class to create and plot ROC curves
for (i in seq_along(classes)) {
  class_index <- classes[i]
  
  respon <- as.numeric(ifelse(comp_mat$TCGA_Project2 == class_index,1,0))
  roc_object <- pROC::roc(response = respon, predictor = comp_mat[, class_index], ci = TRUE)
  
  # Store the AUC value and the confidence interval bounds
  auc_values[i] <- pROC::auc(roc_object)
  ci_values <- pROC::ci(roc_object)
  auc_ci_lower[i] <- ci_values[1]
  auc_ci_upper[i] <- ci_values[3]
  
  # Plot the first curve or add subsequent curves
  if (i == 1) {
    # Plot the first curve
    plot(roc_object, 
         xlab = "Specificity", ylab = "Sensitivity",
         col = line_colors[i], lty = 1, lwd = 3.6,
         print.auc = FALSE)
  } else {
    # Add subsequent curves with the lines() function
    lines(roc_object, col = line_colors[i], lwd = 3.6)
  }
}

# Create a legend with class names, AUC values, and confidence intervals
legend_text <- paste0(classes, ", AUC = ", round(auc_values, 2), " (", 
                      round(auc_ci_lower, 2), "-", round(auc_ci_upper, 2), ")")

legend("bottomright",
       legend = legend_text,col = line_colors,
       lwd = 3,cex=0.8,
       title = "Diagnoses")
dev.off()

######* Make confusion matrix *######
#Matrix GT05
res_gt05 <- res_diagm3_filter_rmFNA[
  res_diagm3_filter_rmFNA$composite_top2match.f %in% c("Matched","Not matched") &
    res_diagm3_filter_rmFNA$Tumor_Purity_ichor>=0.05,]

con_matrix <- table(True = res_gt05$TCGA_Project, 
                    Predicted = res_gt05$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

a3 <- data.frame(True=c(lab_pred,"LIHC","THCA"),
                 Predicted=rep("THCA",times=16),
                 Freq=rep(0,times=16))
a4 <- data.frame(True=c(lab_pred,"LIHC","THCA"),
                 Predicted=rep("LIHC",times=16),
                 Freq=rep(0,times=16))

con_df_m <- rbind(con_df,a3,a4)
con_df_m <- unique(con_df_m)

con_df_mat.gt05 <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BLCA","BRCA","DLBC","Kidney","LIHC",
                                "LUAD","MCF GBM","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","MESO","PAAD","pNET","PRAD","THCA"),
                       labels=c("BLCA","BRCA","DLBC","Kidney","LIHC",
                                "LUAD","MCF GBM","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","MESO","PAAD","pNET","PRAD","THCA")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("THCA","PRAD","pNET","PAAD","MESO","MCF UGI-AD",
                                "MCF SCC","MCF LGI-AD","MCF GYN","MCF GBM","LUAD","LIHC",
                                "Kidney","DLBC","BRCA","BLCA"),
                       labels=c("THCA","PRAD","pNET","PAAD","MESO","MCF UGI-AD",
                                "MCF SCC","MCF LGI-AD","MCF GYN","MCF GBM","LUAD","LIHC",
                                "Kidney","DLBC","BRCA","BLCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_gt05 <- ggplot(con_df_mat.gt05, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "Composite classification", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_gt05


# Matrix1
######## composite classification of BF and FNAgt50 ########
res_bf.FNAgt50_fracGT50_1 <- 
  res_diagm3_filter_rmFNA[(res_diagm3_filter_rmFNA$Tumor_Purity_ichor>0.495 &
                            res_diagm3_filter_rmFNA$Total.Reads>=15000000),]

# make matrix
res_bf.FNAgt50_fracGT50 <-
  res_bf.FNAgt50_fracGT50_1[res_bf.FNAgt50_fracGT50_1$composite_top2match.f 
                            %in% c("Matched","Not matched"),]
con_matrix <- table(True = res_bf.FNAgt50_fracGT50$TCGA_Project, 
                    Predicted = res_bf.FNAgt50_fracGT50$composite_class2)
table(res_bf.FNAgt50_fracGT50$TCGA_Project)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

a3 <- data.frame(True=c(lab_pred,"THCA"),
                 Predicted=rep("THCA",times=14),
                 Freq=rep(0,times=14))

con_df_m <- rbind(con_df,a3)
con_df_m <- unique(con_df_m)

con_df_mat1 <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BLCA","BRCA","DLBC","Kidney",
                                "LUAD","MCF GBM","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","MESO","pNET","PRAD","THCA"),
                       labels=c("BLCA","BRCA","DLBC","Kidney",
                                "LUAD","MCF GBM","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","MESO","pNET","PRAD","THCA")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("THCA","PRAD","pNET","MESO","MCF UGI-AD",
                                "MCF SCC","MCF LGI-AD","MCF GYN","MCF GBM","LUAD",
                                "Kidney","DLBC","BRCA","BLCA"),
                       labels=c("THCA","PRAD","pNET","MESO","MCF UGI-AD",
                                "MCF SCC","MCF LGI-AD","MCF GYN","MCF GBM","LUAD",
                                "Kidney","DLBC","BRCA","BLCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_gt50 <- ggplot(con_df_mat1, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "COMET of tumor fraction >=50%", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_gt50

# Matrix2
# get sampels with frac  0.2-0.5
res_bf.FNAgt50_frac2050_1 <- 
  res_diagm3_filter_rmFNA[(res_diagm3_filter_rmFNA$Tumor_Purity_ichor<=0.495 &
                            res_diagm3_filter_rmFNA$Tumor_Purity_ichor>0.195) |
                           (res_diagm3_filter_rmFNA$Tumor_Purity_ichor>0.495 & 
                              res_diagm3_filter_rmFNA$Total.Reads<15000000),]
#no samples with Tumor_Purity_ichor>0.495 & TotalReads<15000000
table(res_diagm3_filter_rmFNA$TCGA_Project)
res_bf.FNAgt50_frac2050 <-
  res_bf.FNAgt50_frac2050_1[res_bf.FNAgt50_frac2050_1$composite_top2match.f 
                            %in% c("Matched","Not matched"),]

freq1 <- as.data.frame(table(res_bf.FNAgt50_frac2050$TCGA_Project))
freq2 <- as.data.frame(table(res_bf.FNAgt50_frac2050$composite_class2))
freqm=merge(freq1,freq2,by="Var1",all=T)

con_matrix <- table(True = res_bf.FNAgt50_frac2050$TCGA_Project, 
                    Predicted = res_bf.FNAgt50_frac2050$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_true
lab_pred <- names(table(con_df$Predicted))
lab_pred 

a2 <- data.frame(True=rep("Kidney",times=11),
                 Predicted=lab_pred,
                 Freq=rep(0,times=11))


con_df_m <- rbind(con_df,a2)
con_df_m <- unique(con_df_m)

con_df_mat2 <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BLCA","BRCA","DLBC","Kidney","LUAD","MCF GBM",
                                "MCF GYN","MCF LGI-AD","MCF SCC","MCF UGI-AD","PAAD"),
                       labels=c("BLCA","BRCA","DLBC","Kidney","LUAD","MCF GBM",
                                "MCF GYN","MCF LGI-AD","MCF SCC","MCF UGI-AD","PAAD")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("PAAD","MCF UGI-AD","MCF SCC","MCF LGI-AD","MCF GYN",
                                "MCF GBM","LUAD","Kidney","DLBC","BRCA","BLCA"),
                       labels=c("PAAD","MCF UGI-AD","MCF SCC","MCF LGI-AD","MCF GYN",
                                "MCF GBM","LUAD","Kidney","DLBC","BRCA","BLCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_frac2050 <- ggplot(con_df_mat2, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "COMET of tumor fraction 20%-50%", y = "Path. diagnosis", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  #scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_frac2050

# Matrix 3
# get sampels with frac  0.05-0.2
res_bf.FNAgt50_frac0520_1 <- 
  res_diagm3_filter_rmFNA[(res_diagm3_filter_rmFNA$Tumor_Purity_ichor<0.195 &
                             res_diagm3_filter_rmFNA$Tumor_Purity_ichor>=0.05),]

res_bf.FNAgt50_frac0520 <-
  res_bf.FNAgt50_frac0520_1[res_bf.FNAgt50_frac0520_1$composite_top2match.f 
                            %in% c("Matched","Not matched"),]

con_matrix <- table(True = res_bf.FNAgt50_frac0520$TCGA_Project, 
                    Predicted = res_bf.FNAgt50_frac0520$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_true
lab_pred <- names(table(con_df$Predicted))
lab_pred 

a2 <- data.frame(True=rep("MCF SCC",times=10),
                 Predicted=c(lab_pred,"LIHC"),
                 Freq=rep(0,times=10))
a3 <- data.frame(True=c(lab_pred,"LIHC"),
                 Predicted=rep("LIHC",times=10),
                 Freq=rep(0,times=10))

con_df_m <- rbind(con_df,a2,a3)
con_df_m <- unique(con_df_m)

con_df_mat3 <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BRCA","DLBC","Kidney","LIHC","LUAD","MCF GYN",
                                "MCF LGI-AD","MCF SCC","MCF UGI-AD","PAAD"),
                       labels=c("BRCA","DLBC","Kidney","LIHC","LUAD","MCF GYN",
                                "MCF LGI-AD","MCF SCC","MCF UGI-AD","PAAD")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("PAAD","MCF UGI-AD","MCF SCC",
                                "MCF LGI-AD","MCF GYN","LUAD","LIHC",
                                "Kidney","DLBC","BRCA"),
                       labels=c("PAAD","MCF UGI-AD","MCF SCC",
                                "MCF LGI-AD","MCF GYN","LUAD","LIHC",
                                "Kidney","DLBC","BRCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_frac0520 <- ggplot(con_df_mat3, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "COMET of tumor fraction 5%-20%", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_frac0520

# Matrix of ABDO
res_abdo <- res_diagm3_filter_rmFNA[
  res_diagm3_filter_rmFNA$composite_top2match.f %in% c("Matched","Not matched") &
    res_diagm3_filter_rmFNA$sample_tp=="ABDO" &
    res_diagm3_filter_rmFNA$Tumor_Purity_ichor>=0.05,]

con_matrix <- table(True = res_abdo$TCGA_Project, 
                    Predicted = res_abdo$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

a3 <- data.frame(True=c(lab_pred,"LIHC","Kidney"),
                 Predicted=rep("Kidney",times=10),
                 Freq=rep(0,times=10))
a4 <- data.frame(True=c(lab_pred,"LIHC","Kidney"),
                 Predicted=rep("LIHC",times=10),
                 Freq=rep(0,times=10))

con_df_m <- rbind(con_df,a3,a4)
con_df_m <- unique(con_df_m)

con_df_mat.abdo <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BRCA","DLBC","Kidney","LIHC",
                                "LUAD","MCF GYN","MCF LGI-AD",
                                "MCF UGI-AD","MESO","PAAD"),
                       labels=c("BRCA","DLBC","Kidney","LIHC",
                                "LUAD","MCF GYN","MCF LGI-AD",
                                "MCF UGI-AD","MESO","PAAD")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("PAAD","MESO","MCF UGI-AD",
                                "MCF LGI-AD","MCF GYN","LUAD","LIHC",
                                "Kidney","DLBC","BRCA"),
                       labels=c("PAAD","MESO","MCF UGI-AD",
                                "MCF LGI-AD","MCF GYN","LUAD","LIHC",
                                "Kidney","DLBC","BRCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_abdo <- ggplot(con_df_mat.abdo, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "Composite classification", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  #scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_abdo

# Matrix of PLEU
res_pleu <- res_diagm3_filter_rmFNA[
  res_diagm3_filter_rmFNA$composite_top2match.f %in% c("Matched","Not matched") &
    res_diagm3_filter_rmFNA$sample_tp=="PLEU" &
    res_diagm3_filter_rmFNA$Tumor_Purity_ichor>=0.05,]

con_matrix <- table(True = res_pleu$TCGA_Project, 
                    Predicted = res_pleu$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

a3 <- data.frame(True=c(lab_pred,"THCA"),
                 Predicted=rep("THCA",times=10),
                 Freq=rep(0,times=10))

con_df_m <- rbind(con_df,a3)
con_df_m <- unique(con_df_m)

con_df_mat.pleu <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BLCA","BRCA","DLBC","Kidney",
                                "LUAD","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","THCA"),
                       labels=c("BLCA","BRCA","DLBC","Kidney",
                                "LUAD","MCF GYN","MCF LGI-AD","MCF SCC",
                                "MCF UGI-AD","THCA")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("THCA","MCF UGI-AD","MCF SCC",
                                "MCF LGI-AD","MCF GYN","LUAD",
                                "Kidney","DLBC","BRCA","BLCA"),
                       labels=c("THCA","MCF UGI-AD","MCF SCC",
                                "MCF LGI-AD","MCF GYN","LUAD",
                                "Kidney","DLBC","BRCA","BLCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_pleu <- ggplot(con_df_mat.pleu, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "Composite classification", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  #scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_pleu

# Matrix of FNA
res_fna <- res_diagm3_filter_rmFNA[
  res_diagm3_filter_rmFNA$composite_top2match.f %in% c("Matched","Not matched") &
    res_diagm3_filter_rmFNA$sample_tp=="FNA" &
    res_diagm3_filter_rmFNA$Tumor_Purity_ichor>0.495,]

con_matrix <- table(True = res_fna$TCGA_Project, 
                    Predicted = res_fna$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

con_df_m <- con_df
con_df_m <- unique(con_df_m)

con_df_mat.fna <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BLCA","BRCA","LUAD","MCF LGI-AD","MCF SCC",
                                "pNET","PRAD"),
                       labels=c("BLCA","BRCA","LUAD","MCF LGI-AD","MCF SCC",
                                "pNET","PRAD")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("PRAD","pNET","MCF SCC",
                                "MCF LGI-AD","LUAD","BRCA","BLCA"),
                       labels=c("PRAD","pNET","MCF SCC",
                                "MCF LGI-AD","LUAD","BRCA","BLCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_fna <- ggplot(con_df_mat.fna, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "Composite classification", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  #scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_manual(values = c("#2979FF",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_fna

# Matrix of CSF
res_csf <- res_diagm3_filter_rmFNA[
  res_diagm3_filter_rmFNA$composite_top2match.f %in% c("Matched","Not matched") &
    res_diagm3_filter_rmFNA$sample_tp=="CSF" &
    res_diagm3_filter_rmFNA$Tumor_Purity_ichor>=0.05,]

con_matrix <- table(True = res_csf$TCGA_Project, 
                    Predicted = res_csf$composite_class2)
con_df <- as.data.frame(con_matrix)
colnames(con_df)

lab_true <- names(table(con_df$True))
lab_pred <- names(table(con_df$Predicted))
lab_true
lab_pred

a1 <- data.frame(True=rep("MCF SCC",times=8),
                 Predicted=c(lab_pred,"MCF UGI-AD"),
                 Freq=rep(0,times=8))
a2 <- data.frame(True=rep("Kidney",times=8),
                 Predicted=c(lab_pred,"MCF UGI-AD"),
                 Freq=rep(0,times=8))

a3 <- data.frame(True=c(lab_pred,"MCF UGI-AD"),
                 Predicted=rep("MCF UGI-AD",times=8),
                 Freq=rep(0,times=8))

con_df_m <- rbind(con_df,a1,a2,a3)
con_df_m <- unique(con_df_m)

con_df_mat.csf <- con_df_m %>%
  mutate(True_char=as.character(True),
         Pred_char=as.character(Predicted),
         Pred_f=factor(as.factor(Pred_char),order=T,
                       levels=c("BRCA","DLBC","Kidney","LUAD",
                                "MCF GBM","MCF GYN","MCF SCC","MCF UGI-AD"),
                       labels=c("BRCA","DLBC","Kidney","LUAD",
                                "MCF GBM","MCF GYN","MCF SCC","MCF UGI-AD")),
         True_f=factor(as.factor(True_char),order=T,
                       levels=c("MCF UGI-AD","MCF SCC","MCF GYN","MCF GBM","LUAD",
                                "Kidney","DLBC","BRCA"),
                       labels=c("MCF UGI-AD","MCF SCC","MCF GYN","MCF GBM","LUAD",
                                "Kidney","DLBC","BRCA")),
         col_grp=ifelse(True_char==Pred_char,1,ifelse(True_char!=Pred_char & Freq!=0, 2,0)),
         col_grp_f=factor(as.factor(col_grp),
                          levels=c(1,2,0),
                          labels=c("Correct classification","Misclassified","No cases")))

#plot
pmat_csf <- ggplot(con_df_mat.csf, aes(x = Pred_f, y = True_f, fill = col_grp_f)) +
  geom_tile(color = "black") +
  labs(x = "Composite classification", y = "Gold standard", title = NULL) +
  geom_text(aes(label = Freq), vjust = 0.5, hjust = 0.5, color="black") +
  #scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_manual(values = c("#2979FF","#cc0000",'transparent')) + #"#0202d390"
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none")
pmat_csf

######* Accuracy by tumor fraction *######
mat_BFCSF_comb <- 
  read.xlsx(file.path(pathbf2,"Output","25-11-05-matrix.accuracy.byFraction_allBFCSFs.xlsx"),
            sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

mat_BFCSF_comb$frac_grp.f <- factor(as.factor(mat_BFCSF_comb$frac_grp), order=T,
                                    levels=c(1,2,3,4,5,6),
                                    labels=c(">=50%","20%-50%","5-20%","<5%","Total","GT05"))

p_accu.line <- ggplot(data = mat_BFCSF_comb[mat_BFCSF_comb$col_grp==1 &
                                              mat_BFCSF_comb$frac_grp %in% c(1,2,3),], 
                      aes(x = frac_grp.f, y = accu, group=1)) +
  geom_line(color="blue",linewidth=1.3) +
  geom_point(size=2) +
  scale_y_continuous(breaks=seq(0,100, by = 20),
                     labels=seq(0,100, by = 20),
                     limits=c(0,101), expand=c(0,0)) +
  labs(x="Tumor fraction",y="Accuracy (%)") +
  theme(legend.position="none",
        legend.title=element_blank())
p_accu.line    

######* Accuracy by TF trend *###########
mat_comb_all <- read.xlsx(file.path(pathbf2, "Output",
                                    "25-11-05-matrix.accuracy.byFraction_allBFCSFs.xlsx"),
                          sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
# 1. Define the number of successes (k) and total trials (n) for each group
successes <- c(mat_comb_all[1,2],mat_comb_all[3,2],mat_comb_all[5,2])
totals <- c(mat_comb_all[1,6],mat_comb_all[3,6],mat_comb_all[5,6])
dose_scores <- 1:3

# 2. Create a data frame for use with logistic regression
df_trend <- data.frame(
  score = dose_scores,
  success = successes,
  failure = totals - successes)

# 3. Convert to a single row-per-trial format (optional, but good for GLM)
df_long <- data.frame(
  dose = rep(df_trend$score, df_trend$success + df_trend$failure),
  outcome = rep(c(rep(1, df_trend$success[1]), rep(0, df_trend$failure[1]),
                  rep(1, df_trend$success[2]), rep(0, df_trend$failure[2]),
                  rep(1, df_trend$success[3]), rep(0, df_trend$failure[3]))
  ))

library(stats) 
contingency_table <- rbind(successes, totals - successes)

# scores = dose_scores tells R which numeric weights to use for the trend test.
trend_test_ca <- 
  prop.trend.test(x = successes,n = totals,score = dose_scores)
trend_test_ca$p.value #X-squared = 1.4456, df = 1, p-value = 0.2292

######* Accuracy by sample type *######

mat_BFCSF_comb2 <- read.xlsx(file.path(pathbf2, "Output",
                                       "25-11-05-matrix.accuracy.bySampleType_allBFCSFs.xlsx"),
                             sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

mat_BFCSF_comb2$frac_grp.f <- factor(as.factor(mat_BFCSF_comb2$frac_grp),order=T,
                                     levels=c(2,1,3,4),
                                     labels=c("PLEU","ABDO","FNA","CSF"))

p_accu.line2 <- ggplot(data = mat_BFCSF_comb2[mat_BFCSF_comb2$col_grp==1,], 
                       aes(x = frac_grp.f, y = accu, group=1)) +
  geom_line(color="blue",linewidth=1.3) +
  geom_point(size=2) +
  scale_y_continuous(breaks=seq(0,100, by = 20),
                     labels=seq(0,100, by =20),
                     limits=c(0,101), expand=c(0,0)) +
  labs(x="Sample type",y="Accuracy (%)") +
  theme(legend.position="none",
        legend.title=element_blank())
p_accu.line2    

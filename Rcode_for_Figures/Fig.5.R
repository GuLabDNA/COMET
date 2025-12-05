#Fig.5
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
library(ComplexHeatmap)
library(circlize) # For the color ramp

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

######* Make heatmap *######

# write_xlsx(samp_block_meth2,
#            file.path(pathbf2,"Draft","SupplData_markers_heatmap.xlsx"),
#            col_names = TRUE, format_headers = TRUE)
samp_block_meth2<- read.xlsx(file.path(pathbf2,"Draft",
                                "SupplData_markers_heatmap.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

samp_ls1 <- sapply(str_split(colnames(samp_block_meth2)[4:83], "_",  n = 2), `[`, 1)
samp_ls.se <- samp_ls[samp_ls$TCGA_ID %in% samp_ls1,]
# prepare the annotation info
anno_df <- data.frame(
  TCGA_ID = samp_ls.se$TCGA_ID,
  `Sample type` = samp_ls.se$sample_tp,
  Diagnosis = factor(as.factor(samp_ls.se$TCGA_Project),
                     levels = c("BRCA", "DLBC", "LUAD", "MCF GYN", "PAAD","CTRL"),
                     labels = c("BRCA", "DLBC", "LUAD", "MCF GYN", "PAAD","CTRL")),
  status = factor(as.factor(samp_ls.se$status),
                  levels = c("Cancer", "Non-cancer"),
                  labels = c("Cancer", "Non-cancer"))
)
rownames(anno_df) <- anno_df$TCGA_ID
anno_df <- anno_df[anno_df$Diagnosis!="PAAD",] #!anno_df$TCGA_ID %in% c("BF3201XR25B","BF3095XR25A") &
anno_df <- anno_df[order(anno_df$Diagnosis,anno_df$Sample.type),]
col_names <- paste0(anno_df$TCGA_ID,"_beta")

samp_block_meth3 <- as.matrix(samp_block_meth2[,col_names])
dim(samp_block_meth3) #198  74

# Create the top annotation object
ha_top <- HeatmapAnnotation(
  Diagnosis = anno_df$Diagnosis,
  `Sample type` = anno_df$Sample.type,
  col = list(
    # Define colors for the 'site' annotation
    Diagnosis = c("BRCA"="#6e1e33","DLBC"="#e4703a","LUAD"="#b890e0",
                  "MCF GYN"="#0581be","CTRL"="#27780d"),
    # Define colors for the 'status' annotation
    `Sample type` = c("ABDO" = "#D25D5D", "PLEU" = "#B2CD9C","FNA" = "#F0C507")),
  border = TRUE)

# --- Prepare the Heatmap Colors ---
# Create a color gradient from blue to white to red
# You can use a divergent palette for z-scores
color_fun <- colorRamp2(c( 0, 0.5,1), c("blue", "white", "red"))

# --- Create the Heatmap ---
ht <- Heatmap(
  samp_block_meth3,
  name = "Methylation",
  col = color_fun,
  # Row and column settings
  cluster_rows = TRUE,
  cluster_columns = FALSE, # Set to FALSE to use manual ordering
  show_row_names = FALSE,
  show_column_names = FALSE,
  # Split the rows by their site
  row_split = as.vector(samp_block_meth2$target),
  # Add space between the row clusters (e.g., 5mm)
  row_gap = unit(2, "mm"),
  column_title = "Tissue-specific methylation regions",
  row_title = NULL,
  # Add the annotations
  top_annotation = ha_top
)
ht

######* plots of Z-scores *######

#import data
res_diag <- read.xlsx(file.path(pathbf2,"Draft",
                                "25-1202-Supple.Tables.xlsx"),
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

res_diag3_filter <- res_diag_comb1[res_diag_comb1$Total.Reads>=10000000
                              & res_diag_comb1$sample_tp!='CSF',]
res_diag3_deconv.noref <- 
  res_diag3_filter[!res_diag3_filter$TCGA_Project %in% 
    c("ACC","Adeno/carcinoma", "AFX/PDS", "AS", "Benign", "CHOL", 
      "FMS","GIST", "Haem", "LAML", "LMO",
      "LMS", "MALIGNANT", "MECA", "MESO", "Myeloma", 
      "NET", "PRAD", "SARC", "SARC (MPNST-like)", "SKCM",
      "SWN", "SYSA", "TGCT_NS", "TGCT_S", "THCA/SARC","MDS",
      "ThoraxAD", "Unclear","Unknown", "UVM", "YST",
      "NET in esophagus","NET in lung","NET in small bowel","pNET"),] #299
#NET in esophagus has two results: SCC/MCG UGI
#NET in lung has two results: SCC/LUAD

res_diag3_deconv.noref2 <- res_diag3_deconv.noref %>% 
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
  filter(`Sample.Type(Sub2)` %in% c("ABDO","PLEU") &
           (Tumor_Purity_ichor>=0.05)) #128


# write_xlsx(ca_long_Zabdo,
#            file.path(pathbf2,"Draft","SupplData_deconv_boxplot.xlsx"),
#            col_names = TRUE, format_headers = TRUE)
ca_long_Zabdo <- read.xlsx(file.path(pathbf2,"Draft",
                                "2SupplData_deconv_boxplot.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

ca_df = res_diag3_deconv.noref2[,c("BF_id","TCGA_Project",
                                   "sample_tp","Tumor_Purity_ichor")]
case_freq <- as.data.frame(table(ca_df$TCGA_Project))
luad_id <- ca_df[ca_df$TCGA_Project=="LUAD","BF_id"]
brca_id <- ca_df[ca_df$TCGA_Project=="BRCA","BF_id"]
lympho_id <- ca_df[ca_df$TCGA_Project=="DLBC","BF_id"]
gyn_id <- ca_df[ca_df$TCGA_Project=="MCF GYN","BF_id"]
paad_id <- ca_df[ca_df$TCGA_Project=="PAAD","BF_id"]
scc_id <- ca_df[ca_df$TCGA_Project %in% c("SCC","MCF SCC"),"BF_id"]
ugi_id <- ca_df[ca_df$TCGA_Project=="MCF UGI-AD","BF_id"]
lgi_id <- ca_df[ca_df$TCGA_Project=="MCF LGI-AD","BF_id"]

ca_long_Zcomb <- rbind(ca_long_Zabdo,ca_long_Zpleu)
ca_long_Zcomb <- merge(ca_long_Zcomb,ca_df[ca_df$sample_tp!="FNA",
                                           c("BF_id","Tumor_Purity_ichor")],
                       by="BF_id",all.x =T) #3996 #2634
ca_long_Zcomb <- within(ca_long_Zcomb, {
  purity_catn <- NA
  purity_catn[Tumor_Purity_ichor >=0.495] <- 1
  purity_catn[Tumor_Purity_ichor>=0.05 & Tumor_Purity_ichor <0.495] <- 2
})

ca_long_Zcomb <- ca_long_Zcomb %>% 
  mutate(grp_luad =ifelse(BF_id %in% luad_id, yes="LUAD", no="Non-LUAD"),
         grp_brca =ifelse(BF_id %in% brca_id, yes="BRCA", no="Non-BRCA"),
         grp_lympho =ifelse(BF_id %in% lympho_id, yes="DLBC", no="Non-DLBC"), 
         grp_gyn =ifelse(BF_id %in% gyn_id, yes="MCF GYN", no="Non-MCF.GYN"), 
         grp_paad =ifelse(BF_id %in% paad_id, yes="PAAD", no="Non-PAAD"), 
         grp_scc =ifelse(BF_id %in% scc_id, yes="MCF SCC", no="Non-MCF.SCC"), 
         grp_ugi =ifelse(BF_id %in% ugi_id, yes="MCF UGI-AD", no="Non-MCF.UGI.AD"), 
         grp_lgi =ifelse(BF_id %in% lgi_id, yes="MCF LGI-AD", no="Non-MCF.LGI.AD"), 
         
         grp_luad.f=factor(as.factor(grp_luad), order=T,
                           levels=c("LUAD", "Non-LUAD"),
                           labels=c("LUAD", "Non-LUAD")),
         grp_brca.f=factor(as.factor(grp_brca), order=T,
                           levels=c("BRCA", "Non-BRCA"),
                           labels=c("BRCA", "Non-BRCA")),
         grp_lympho.f=factor(as.factor(grp_lympho), order=T,
                             levels=c("DLBC", "Non-DLBC"),
                             labels=c("DLBC", "Non-DLBC")),
         grp_gyn.f=factor(as.factor(grp_gyn), order=T,
                          levels=c("MCF GYN", "Non-MCF.GYN"),
                          labels=c("MCF GYN", "Non-MCF.GYN")),
         grp_paad.f=factor(as.factor(grp_paad), order=T,
                           levels=c("PAAD", "Non-PAAD"),
                           labels=c("PAAD", "Non-PAAD")),
         grp_scc.f=factor(as.factor(grp_scc), order=T,
                          levels=c("MCF SCC", "Non-MCF.SCC"),
                          labels=c("MCF SCC", "Non-MCF.SCC")),
         grp_ugi.f=factor(as.factor(grp_ugi), order=T,
                          levels=c("MCF UGI-AD", "Non-MCF.UGI.AD"),
                          labels=c("MCF UGI-AD", "Non-MCF.UGI.AD")),
         grp_lgi.f=factor(as.factor(grp_lgi), order=T,
                          levels=c("MCF LGI-AD", "Non-MCF.LGI.AD"),
                          labels=c("MCF LGI-AD", "Non-MCF.LGI.AD")),
         purity_catn.f = factor(as.factor(purity_catn),order=T,
                                levels=c("1","2"),labels=c(">=0.5","0.05-0.5"))
  )

ca_df_grp <- within(ca_df, {
  purity_catn <- NA
  purity_catn[Tumor_Purity_ichor >=0.495] <- 1
  purity_catn[Tumor_Purity_ichor>=0.05 & Tumor_Purity_ichor <0.495] <- 2
  purity_catn[Tumor_Purity_ichor<0.05 & !is.na(Tumor_Purity_ichor)] <- 3
})

ca_df_grp <- ca_df_grp %>% 
  mutate(grp_luad =ifelse(BF_id %in% luad_id, yes="LUAD", no="Non-LUAD"),
         grp_brca =ifelse(BF_id %in% brca_id, yes="BRCA", no="Non-BRCA"),
         grp_lympho =ifelse(BF_id %in% lympho_id, yes="DLBC", no="Non-DLBC"), 
         grp_gyn =ifelse(BF_id %in% gyn_id, yes="MCF GYN", no="Non-MCF.GYN"), 
         grp_paad =ifelse(BF_id %in% paad_id, yes="PAAD", no="Non-PAAD"), 
         grp_scc =ifelse(BF_id %in% scc_id, yes="MCF SCC", no="Non-MCF.SCC"), 
         grp_ugi =ifelse(BF_id %in% ugi_id, yes="MCF UGI-AD", no="Non-MCF.UGI.AD"), 
         grp_lgi =ifelse(BF_id %in% lgi_id, yes="MCF LGI-AD", no="Non-MCF.LGI.AD"), 
         
         grp_luad.f=factor(as.factor(grp_luad), order=T,
                           levels=c("LUAD", "Non-LUAD"),
                           labels=c("LUAD", "Non-LUAD")),
         grp_brca.f=factor(as.factor(grp_brca), order=T,
                           levels=c("BRCA", "Non-BRCA"),
                           labels=c("BRCA", "Non-BRCA")),
         grp_lympho.f=factor(as.factor(grp_lympho), order=T,
                             levels=c("DLBC", "Non-DLBC"),
                             labels=c("DLBC", "Non-DLBC")),
         grp_gyn.f=factor(as.factor(grp_gyn), order=T,
                          levels=c("MCF GYN", "Non-MCF.GYN"),
                          labels=c("MCF GYN", "Non-MCF.GYN")),
         grp_paad.f=factor(as.factor(grp_paad), order=T,
                           levels=c("PAAD", "Non-PAAD"),
                           labels=c("PAAD", "Non-PAAD")),
         grp_scc.f=factor(as.factor(grp_scc), order=T,
                          levels=c("MCF SCC", "Non-MCF.SCC"),
                          labels=c("MCF SCC", "Non-MCF.SCC")),
         grp_ugi.f=factor(as.factor(grp_ugi), order=T,
                          levels=c("MCF UGI-AD", "Non-MCF.UGI.AD"),
                          labels=c("MCF UGI-AD", "Non-MCF.UGI.AD")),
         grp_lgi.f=factor(as.factor(grp_lgi), order=T,
                          levels=c("MCF LGI-AD", "Non-MCF.LGI.AD"),
                          labels=c("MCF LGI-AD", "Non-MCF.LGI.AD")),
         purity_catn.f = factor(as.factor(purity_catn),order=T,
                                levels=c("1","2","3"),labels=c(">=0.5","0.05-0.5","<0.05"))
  )

###################################
#LUAD
lung1 <- ca_long_Zcomb %>% filter(CellType=="Lung-Ep-Alveo")  

lung1 %>% wilcox_test(z.score ~ grp_luad.f, paired = F)
#z.score LUAD   Non-LUAD    26   102      2329 0.0000000029
lung1_can <- lung1[lung1$grp_luad.f=="LUAD",]
table(lung1_can$or_zscore)
#  1  2  4  6  9 15 16 17 18 
# 15  3  1  1  1  2  1  1  1  
lung1_top <- lung1[lung1$grp_luad.f=="Non-LUAD",]

p_luad <- ggplot(data=lung1) + 
  geom_boxplot(aes(x = grp_luad.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_luad.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") + 
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Lung alveolar epi.")+
  theme(
        legend.position="none",
        legend.title=element_blank())
p_luad

nref=21
p_luad_or <- ggplot(data=lung1, aes(x = or_zscore, y = z.score, 
                                    color = grp_luad.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#b890e0","grey65"))+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5,nref, by = 5),nref),
                     labels=c(1,2,3,seq(5,nref, by = 5),nref),
                     limits=c(0.5,nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Lung alveolar epi.")+
  theme(
        legend.position="none",
        legend.title=element_blank()) 
p_luad_or

###################################
#BRCA
brca1 <- ca_long_Zcomb %>% filter(CellType %in% c("Breast-Luminal-Ep"))
table(brca1$grp_brca.f)

brca1 %>% wilcox_test(z.score ~ grp_brca.f, paired = F)
#z.score BRCA   Non-BRCA    20   108      1897 0.000000084
brca1_can <- brca1[brca1$grp_brca.f=="BRCA",]
table(brca1_can$or_zscore)
#  1  9 13 17 
# 17  1  1  1 
brca1_top <- brca1[brca1$grp_brca.f=="Non-BRCA",]

# plot
p_brca <- ggplot(data=brca1) + 
  geom_boxplot(aes(x = grp_brca.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_brca.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +  
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Breast luminal/basal epi.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_brca

p_brca_or <- ggplot(data=brca1, aes(x = or_zscore, y = z.score, 
                                    color = grp_brca.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#6e1e33","grey65"))+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5,nref, by = 5),nref),
                     labels=c(1,2,3,seq(5,nref, by = 5),nref),
                     limits=c(0.5,nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Breast luminal/basal epi.")+
  theme(
        legend.position="none",
        legend.title=element_blank()) 
p_brca_or

##############################
#DLBC
lympho1 <- ca_long_Zcomb %>% filter(CellType=="Blood-B") 
table(lympho1$grp_lympho.f)

lympho1 %>% wilcox_test(z.score ~ grp_lympho.f, paired = F)
#z.score DLBC   Non-DLBC    12   116      1336 0.000000172
lympho1_can <- lympho1[lympho1$grp_lympho.f=="DLBC",]
table(lympho1_can$or_zscore)
#  1  6 14 
# 10  1  1 
lympho1_top <- lympho1[lympho1$grp_lympho.f=="Non-DLBC",]

p_lympho <- ggplot(data=lympho1) +
  geom_boxplot(aes(x = grp_lympho.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_lympho.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="B-cell")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_lympho

p_lympho_or <- ggplot(data=lympho1, aes(x = or_zscore, y = z.score, 
                                        color = grp_lympho.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#e4703a","grey65"))+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5, nref, by = 5), nref),
                     labels=c(1,2,3,seq(5, nref, by = 5), nref),
                     limits=c(0.5, nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="B-cell")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_lympho_or

#########################
#MCF GYN
gyn1 <- ca_long_Zcomb %>% filter(CellType=="Ovary-Ep") 
table(gyn1$grp_gyn.f)

gyn1 %>% wilcox_test(z.score ~ grp_gyn.f, paired = F)
#z.score MCF GYN Non-MCF.GYN    34    94      3128 1.56e-16
gyn1_can <- gyn1[gyn1$grp_gyn.f=="MCF GYN",]
table(gyn1_can$or_zscore)
#  1  2  3  4 
# 31  1  1  1
gyn1_top <- gyn1[gyn1$grp_gyn.f=="Non-MCF.GYN",]

p_gyn <- ggplot(data=gyn1) +
  geom_boxplot(aes(x = grp_gyn.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_gyn.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-15,0,30,60,90),
                     labels=c(-15,0,30,60,90),
                     limits=c(-15,92), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Ovarian ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_gyn

p_gyn_or <- ggplot(data=gyn1, aes(x = or_zscore, y = z.score, 
                                  color = grp_gyn.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#0581be","grey65"))+
  scale_y_continuous(breaks=c(-15,0,30,60,90),
                     labels=c(-15,0,30,60,90),
                     limits=c(-15,92), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5, nref, by = 5), nref),
                     labels=c(1,2,3,seq(5, nref, by = 5), nref),
                     limits=c(0.5, nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Ovarian ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_gyn_or

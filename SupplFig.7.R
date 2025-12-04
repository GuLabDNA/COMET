#SupplFig. 7

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

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

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


ca_long_Zabdo <- read.xlsx(file.path(pathbf2,"Draft",
                                     "SupplData_deconv_boxplot.xlsx"),
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

nref=21
#########################
#MCF LGI-AD
lgi1 <- ca_long_Zcomb %>% filter(CellType %in% c("Colon-Ep")) 
lgi1 %>% wilcox_test(z.score ~ grp_lgi.f, paired = F)
#z.score MCF LGI-AD Non-MCF.LGI.AD     7   121       823 0.000029
lgi1_can <- lgi1[lgi1$grp_lgi.f=="MCF LGI-AD",]
lgi1_top <- lgi1[lgi1$grp_lgi.f=="Non-MCF.LGI.AD",]

coad1 <- lgi1 %>% 
  mutate(grp_coad=ifelse(BF_id %in% c("BF1159LW15","BF3012","BF1123","BF3063","BF3164"),
                         yes="COAD", no="Non-COAD"),
         grp_coad.f=factor(as.factor(grp_coad), order=T,
                           levels=c("COAD", "Non-COAD"),
                           labels=c("COAD", "Non-COAD")))
coad1_can <- coad1[coad1$BF_id %in% c("BF1159LW15","BF3012","BF1123","BF3063","BF3164"),]

coad1_top <- coad1[coad1$grp_coad.f=="Non-COAD",]

p_coad <- ggplot(data=coad1) +
  geom_boxplot(aes(x = grp_coad.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_coad.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Colon ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_coad

p_coad_or <- ggplot(data=coad1, aes(x = or_zscore, y = z.score, 
                                    color = grp_coad.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#f0c507","grey65"))+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5, nref, by = 5), nref),
                     labels=c(1,2,3,seq(5, nref, by = 5), nref),
                     limits=c(0.5, nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Colon ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_coad_or

#########################
#MCF UGI-AD
ugi1 <- ca_long_Zcomb %>% filter(CellType %in% c("Gastric-Ep")) 
ugi1 %>% wilcox_test(z.score ~ grp_ugi.f, paired = F)
#z.score MCF UGI-AD Non-MCF.UGI.AD    11   117      1086 0.000171

ugi1_can <- ugi1[ugi1$grp_ugi.f=="MCF UGI-AD",]          
ugi1_top <- ugi1[ugi1$grp_ugi.f=="Non-MCF.UGI.AD",]

stad1 <- ugi1 %>% 
  mutate(grp_stad=ifelse(BF_id %in% c("BF3102","BF1033","BF3112","BF3905","BF3266",
                                      "BF3379","BF1153LW6","BF3388"),
                         yes="STAD", no="Non-STAD"),
         grp_stad.f=factor(as.factor(grp_stad), order=T,
                           levels=c("STAD", "Non-STAD"),
                           labels=c("STAD", "Non-STAD")))
table(stad1_can$or_zscore)
# 1  2 12 
# 5  2  1
stad1_top <- stad1[stad1$grp_stad.f=="Non-STAD",]

p_stad <- ggplot(data=stad1) +
  geom_boxplot(aes(x = grp_stad.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_stad.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Gastric ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_stad

p_stad_or <- ggplot(data=stad1, aes(x = or_zscore, y = z.score, 
                                    color = grp_stad.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#F79246","grey65"))+
  scale_y_continuous(breaks=c(-10,0,25,50),
                     labels=c(-10,0,25,50),
                     limits=c(-10,52), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5, nref, by = 5), nref),
                     labels=c(1,2,3,seq(5, nref, by = 5), nref),
                     limits=c(0.5, nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Gastric ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_stad_or

#########################
#PAAD
paad1 <- ca_long_Zcomb %>% filter(CellType %in% c("Pancreas-Acinar"))
paad1 %>% wilcox_test(z.score ~ grp_paad.f, paired = F)
#z.score PAAD   Non-PAAD     6   122       466 0.262
paad1_can <- paad1[paad1$grp_paad.f=="PAAD",]
table(paad1_can$or_zscore)
# 1  7 17 20 
# 3  1  1  1
paad1_top <- paad1[paad1$grp_paad.f=="Non-PAAD",]

p_paad <- ggplot(data=paad1) +
  geom_boxplot(aes(x = grp_paad.f, y = z.score),width=0.4) + 
  geom_beeswarm(aes(x = grp_paad.f, y = z.score, shape = purity_catn.f),
                cex=4, priority = "density", corral = "wrap") +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_y_continuous(breaks=c(-15,0,15,30),
                     labels=c(-15,0,15,30),
                     limits=c(-15,31), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Pancreatic acinar/duct ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_paad

p_paad_or <- ggplot(data=paad1, aes(x = or_zscore, y = z.score, 
                                    color = grp_paad.f)) + 
  geom_point(size=1.8) +
  geom_hline(yintercept=-3,color = "grey80", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=3,color = "grey80", linetype="dashed", linewidth=0.5)+
  scale_color_manual(values=c("#734ba9","grey65"))+
  scale_y_continuous(breaks=c(-15,0,15,30),
                     labels=c(-15,0,15,30),
                     limits=c(-15,31), expand=c(0,0)) +
  scale_x_continuous(breaks=c(1,2,3,seq(5, nref, by = 5), nref),
                     labels=c(1,2,3,seq(5, nref, by = 5), nref),
                     limits=c(0.5, nref+0.5), expand=c(0,0)) +
  labs(y="Z-score",x=NULL,title="Pancreatic acinar/duct ep.")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_paad_or

##########################
# TF vs. zscore
library(ggpmisc)
require(broom)

res_diag3 <- res_diag_comb1 %>% 
  rename(top1_deconv2 = `Top1.Deconvolution.(CellType-zscore)`,
         top2_deconv2 = `Top2.Deconvolution.(CellType-zscore)`,
         top3_deconv2 = `Top3.Deconvolution.(CellType-zscore)`) %>% 
  filter(sample_tp!='CSF')

zcomb <- ca_long_Zcomb[ca_long_Zcomb$or_zscore==1,]
zcomb_diag <- merge(zcomb,res_diag3[res_diag3$Tumor_Purity_ichor>=0.05 & 
                                      res_diag3$sample_tp %in% c("ABDO","PLEU"), 
                                    c("BF_id","Final.Diagnosis","TCGA_Project","Tumor_Purity_ichor",
                                      "purity_catn",
                                      "top1_deconv2","top2_deconv2","top3_deconv2")])

zcomb_diag <- zcomb_diag %>% 
  mutate(top1_zscore = as.numeric(sapply(str_split(top1_deconv2, "_",  n = 2), `[`, 2)),
         purity_catn.f = factor(as.factor(purity_catn),order=T,
                                levels=c("1","2","3"),labels=c(">=0.5","0.05-0.5","<0.05")))

#Line chart
pz2 <- ggplot(zcomb_diag, aes(x=Tumor_Purity_ichor, y=top1_zscore)) +
  geom_point(color="#005b96") + 
  geom_smooth(method=lm,na.rm = TRUE,color="#005b96") +
  stat_poly_eq(use_label(c("eq", "R2","P")),
               label.x = "right",
               formula = y ~ x,
               parse = TRUE) +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "gray80", linewidth=0.5) +
  #geom_hline(yintercept = 2, linetype="dashed", color = "gray20", linewidth=1) +
  geom_hline(yintercept = 3, linetype="dashed", color = "gray70", linewidth=1) +
  labs(x = "Tumor fraction", y = "Z-score") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.2),
                     labels=seq(0, 1, by = 0.2), limits=c(-0.005,1), expand = c(0, 0)) + 
  scale_y_continuous(breaks=seq(0, 100, by = 20),
                     labels=seq(0, 100, by = 20), limits=c(-0.01,100.1), expand = c(0, 0)) + 
  theme(panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(linewidth=1.1,color="black",fill=NA),
        text=element_text(size=14, face="bold"),
        axis.text.x = element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.position="none",
        legend.title=element_blank()) 
pz2

fit = lm(top1_zscore ~ Tumor_Purity_ichor, data=zcomb_diag)
stat_all <- broom::glance(fit)
print(stat_all)

##############################
# 4 cancer types
brca1 <- ca_long_Zcomb %>% filter(CellType %in% c("Breast-Luminal-Ep"))
brca_uniq <- brca1[!duplicated(brca1$BF_id) & brca1$grp_brca=="BRCA",] #20

lung1 <- ca_long_Zcomb %>% filter(CellType=="Lung-Ep-Alveo")  
lung_uniq <- lung1[lung1$grp_luad=="LUAD",] #26

gyn1 <- ca_long_Zcomb %>% filter(CellType=="Ovary-Ep") 
gyn_uniq <- gyn1[gyn1$grp_gyn=="MCF GYN",] #34

lympho1 <- ca_long_Zcomb %>% filter(CellType=="Blood-B") 
lympho_uniq <- lympho1[lympho1$grp_lympho=="DLBC",] #12

can4_comb <- rbind(brca_uniq,lung_uniq,gyn_uniq,lympho_uniq)
can4_comb <- within(can4_comb, {
  diag <- NA
  diag[grp_brca == "BRCA"] <- "BRCA"
  diag[grp_luad=="LUAD"] <- "LUAD"
  diag[grp_gyn=="MCF GYN"] <- "MCF GYN"
  diag[grp_lympho=="DLBC"] <- "DLBC"})

can4_comb2 <- merge(can4_comb[,c(1:10,27,28)],
                    res_diag3[,c("BF_id","tSNE(Reference.TCGA.Tumors)","Final.Diagnosis",
                    "top1_deconv2","top2_deconv2","top3_deconv2")],all.x=T)
can4_comb2$top1_zscore <- as.numeric(sapply(str_split(can4_comb2$top1_deconv2, "_",  n = 2), `[`, 2))

table(can4_comb2$purity_catn.f,can4_comb2$diag)
#          BRCA DLBC LUAD MCF GYN
# >=0.5       7    5    2      14
# 0.05-0.5   13    7   24      20
# 14/34=41.2%
# 7/20=35%
# 5/12=41.7%
# 2/26=7.7%

pz_can <- ggplot(can4_comb2, aes(x=Tumor_Purity_ichor, y=top1_zscore, color=diag)) +
  geom_point() + 
  scale_color_manual(values=line_colors) +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "gray70", linewidth=0.5) +
  geom_hline(yintercept = 3, linetype="dashed", color = "gray70", linewidth=1) +
  labs(x = "Tumor purity", y = "Z-score") +
  facet_wrap(~ diag, ncol=2) +
  scale_x_continuous(breaks=seq(0, 1, by = 0.2),
                     labels=seq(0, 1, by = 0.2), limits=c(-0.005,1), expand = c(0, 0)) + 
  scale_y_continuous(breaks=seq(0, 100, by = 25),
                     labels=seq(0, 100, by = 25), limits=c(-0.01,100), expand = c(0, 0)) + 
  theme(legend.position="none",
        legend.title=element_blank()) 
pz_can


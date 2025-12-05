# #Supplementary Fig. 4
pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

res_diag <- read.xlsx(file.path(pathbf2,"Draft",
                                "25-1202-Supple.Tables.xlsx"),
                      sheet=2, skipEmptyRows=FALSE, colNames = TRUE)
res_diag <- res_diag %>% 
  rename(composite_score2 = composite_score,
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

res_diag3 <- res_diag_comb1[res_diag_comb1$`Sample.type(PLEU,ABDO,FNA,.CSF)`!="CSF",] #413
res_diag3 <- res_diag3 %>% 
  rename(sample_tp = `Sample.type(PLEU,ABDO,FNA,.CSF)`)
         
########## DNA concentration by sample type ##########
re_order <- c("ABDO", "PLEU", "FNA")

res_diag3_or <- res_diag3[order(match(res_diag3$sample_tp, re_order)), ]
res_diag3_or <- res_diag3_or %>% 
  mutate(
    sample_tp.f = factor(as.factor(sample_tp),order=T,levels = re_order))

res_diag3_or <- res_diag3_or[res_diag3_or$quant_dna>=0.1,] #3 samples<0.1

p_inputAll <- ggplot(data=res_diag3_or, 
                     aes(x = sample_tp.f, y = quant_dna, color=sample_tp.f)) +
  geom_boxplot(width=0.7) +
  geom_beeswarm(priority="density",cex=1.3,size=1.7)+
  labs(y="DNA concentration (ng/uL)",x=NULL)+ 
  scale_color_manual(values = c("#009392","#9CCB86","#EEB479")) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100, 1000,10000),labels = c(0.1, 1, 10, 100, 1000,10000),
    limits=c(0.1,1e5)) + #labels = scales::number_format(),
  theme() 
p_inputAll  

stat_c <- res_diag3_or %>% 
  group_by(sample_tp.f) %>%
  reframe(
    count = n(),
    mean = mean(quant_dna, na.rm = TRUE),
    sd = sd(quant_dna, na.rm = TRUE),
    median=median(quant_dna, na.rm = TRUE),
    iqr=quantile(quant_dna, na.rm = TRUE)
  )

########## DNA concentration by diag_caco ##########

res_diag3_rm <- res_diag3[res_diag3$quant_dna>=0.1,] #3 samples <0.1
p_dxcat <- 
  ggplot(data=res_diag3_rm, aes(x = dx_cat.f, y = quant_dna, color=dx_cat.f)) +
  geom_boxplot(width=0.7) +
  geom_beeswarm(priority="density",cex=1.3,size=1.5)+
  labs(y="DNA concentration (ng/uL)",x=NULL)+
  scale_y_log10(labels = c(0.1, 1, 10, 100, 1000,10000),
                breaks = c(0.1, 1, 10, 100, 1000,10000), limits=c(0.1,1e5))+
  scale_color_manual(values = c("#EEB479","#F8766D","#00BFC4","#00A9FF"))+ 
  theme()
p_dxcat

stat_df <- res_diag3_rm %>% 
  group_by(dx_cat.f) %>%
  reframe(
    count = n(),
    mean = mean(quant_dna, na.rm = TRUE),
    sd = sd(quant_dna, na.rm = TRUE),
    median=median(quant_dna, na.rm = TRUE),
    iqr=quantile(quant_dna, na.rm = TRUE))

########## DNA concentration by cytology ##########
res_diag3_rm <- res_diag3[res_diag3$cyto_catx2!="Non-diagnostic",] #3 samples <0.1

ex_cyto.neg <- res_diag3_rm[res_diag3_rm$cyto_catx2=="Negative",]
ex_cyto.neg.gt1 <- ex_cyto.neg[ex_cyto.neg$quant_dna>=1,]
#91/138=65.9%
#boxplot: input by dx cat

res_diag3_rm <- res_diag3_rm[res_diag3_rm$quant_dna>=0.1,] #3 samples <0.1
res_diag3_rm$cyto_catx2.f <- factor(as.factor(res_diag3_rm$cyto_catx2),order=T,
                                    levels=c("Atypical","Malignant","Negative","Suspicious"),
                                    labels=c("Atypical","Malignant","Benign","Suspicious"))
table(res_diag3_rm$cyto_catx2.f)
p_cytocat <- 
  ggplot(data=res_diag3_rm, aes(x = cyto_catx2.f, y = quant_dna, color=cyto_catx2.f)) +
  geom_boxplot(width=0.7) +
  geom_beeswarm(priority="density",cex=1.3,size=1.3)+
  labs(y="DNA concentration (ng/uL)",x="Cytology test")+
  scale_y_log10(labels = c(0.1, 1, 10, 100, 1000,10000),
                breaks = c(0.1, 1, 10, 100, 1000,10000), limits=c(0.1,1e6))+
  scale_color_manual(values = c("#ED775A","#B2CD9C","#FAD691","#5BC0DE"))+ 
  theme()
p_cytocat

########## tumor fraction by sample type ##########

p_tf_ca <- 
  ggplot(data=res_diag3_or[res_diag3_or$diag_caco=="Cancers",], 
         aes(x = sample_tp.f, y = Tumor_Purity_ichor, color=sample_tp.f)) +
  geom_boxplot(width=0.7) +
  #geom_beeswarm(priority="density",cex = 0.5, size = 0.8, groupOnX = TRUE)+
  geom_beeswarm(priority="density",cex = 1.3, size = 1.5)+
  labs(y="Tumor fraction",x=NULL)+
  scale_color_manual(values = c("#009392","#9CCB86","#EEB479"))+ 
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.3), expand=c(0,0)) +
  theme()  
p_tf_ca

###### Tumor purity with diff references ######

df1 <- read.xlsx(file.path(pathbf2, "Draft","TPextract_ManualCheck.xlsx"),
                 sheet=1, skipEmptyRows=FALSE, colNames = TRUE)

df2 <- df1[!df1$Name %in% c("BF3119XR25A_comb","BF3120XR25A_comb",
                            "BF3121XR25A_comb","BF3194XR26A_comb",
                            "BF3201XR25B_comb","BF3203XR26A_comb",
                            "BF3206XR26A_comb","BF3593XR26A_comb",
                            "BF3376XR26A_comb","BF3003XR25A_comb","BF3056XR25A_comb"),1:3] %>% 
  mutate(
    temp1=sapply(str_split(REF, "\\.",  n = 2), `[`, 1),
    sample_tp=sapply(str_split(temp1, "ref",  n = 2), `[`, 1),
    ref_tp=sapply(str_split(temp1, "ref",  n = 2), `[`, 2),
    TCGA_ID=gsub("_comb","",Name))

df2 <- df2[,c("TCGA_ID","sample_tp","ref_tp","Tumor_Purity")]  

df2m <- merge(df2, res_diag3[,c("TCGA_ID","Sample.Type(Sub1)","Sample.Type(Sub2)")],
              by="TCGA_ID",all.y=T)

df2m$Tumor_Purity2 <- ifelse(df2m$Tumor_Purity<0.05, 0.025,df2m$Tumor_Purity)

dfm_trans <- data.table::dcast(setDT(df2m), TCGA_ID ~ ref_tp,value.var = "Tumor_Purity")

dfm_trans <- dfm_trans %>%
  mutate(ref_new = coalesce(ABDO,CSF,LN,PLEU),
         ref_new = ifelse(is.na(ref_new),Plasma,ref_new),
         ref_purity_diff=Plasma-ref_new,
         cna_ref_plasma=ifelse(Plasma>=0.05,1,0),
         cna_ref_bf=ifelse(ref_new>=0.05,1,0))

dfm_trans <- dfm_trans[dfm_trans$ref_purity_diff!=0.000000,] #rm 20 samples
table(dfm_trans$cna_ref_plasma,dfm_trans$cna_ref_bf)

p_scatter <-ggplot(data=dfm_trans, aes(x=Plasma,y=ref_new)) +
  geom_point()
p_scatter

mean_diff <- mean(dfm_trans$ref_purity_diff)
lower <- mean_diff - 1.96*sd(dfm_trans$ref_purity_diff)
upper <- mean_diff + 1.96*sd(dfm_trans$ref_purity_diff)

#create Bland-Altman plot
p_altman <- ggplot(data=dfm_trans, aes(x = Plasma, y = ref_purity_diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff, linetype="dashed") +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  labs(y="Difference Between Measurements",x="Average Measurement",title=NULL) +
  scale_y_continuous(breaks=seq(-0.2,0.4, by = 0.2),
                     labels=seq(-0.2,0.4, by = 0.2),
                     limits=c(-0.21, 0.41), expand=c(0,0)) +
  theme() 
p_altman

######* Tumor fraction of matched pairs *######

############ df_abdo #################
df_abdo <- df2m[df2m$sample_tp=="ABDO",]
df_abdo$ref_grp <- ifelse(df_abdo$ref_tp=="ABDO",1,2)
df_abdo$paired <- rep(1:(nrow(df_abdo)/2), each = 2)

df_abdo %>% wilcox_test(Tumor_Purity2 ~ ref_tp, paired = TRUE)
#1 Tumor_Purity2 ABDO   Plasma    78    78      738. 0.625

p_abdo <- ggplot(data=df_abdo, aes(x = ref_tp, y = Tumor_Purity, 
                                   color = ref_tp)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#009392","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65",  alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.625", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",title="Peritoneal fluid/ascites")+
  theme(legend.position="none",
        legend.title=element_blank())
p_abdo

############ df_pleu #################
df_pleu <- df2m[df2m$sample_tp=="PLEU",]
df_pleu$ref_grp <- ifelse(df_pleu$ref_tp=="PLEU",1,2)
df_pleu$paired <- rep(1:(nrow(df_pleu)/2), each = 2)
df_pleu$ref_tp.f <- factor(as.factor(df_pleu$ref_tp),order=T,
                           levels=c("PLEU","Plasma"),
                           labels=c("PLEU","Plasma")) 
df_pleu_or <- df_pleu[order(df_pleu$ref_tp.f),]

df_pleu %>% wilcox_test(Tumor_Purity2 ~ ref_tp.f, paired = TRUE)
#1 Tumor_Purity2 PLEU   Plasma   173   173     3070. 0.827

p_pleu <- ggplot(data=df_pleu_or, aes(x = ref_tp.f, y = Tumor_Purity, 
                                      color = ref_tp.f)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#9CCB86","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65",  alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.827", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",title="Pleural fluid")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_pleu

############ df_fna #################
df_fna <- df2m[df2m$sample_tp=="FNA",]

df_fna_or <- df_fna[order(df_fna[,"Sample.Type(Sub1)"]),]
df_fna_or <- df_fna_or[!df_fna_or$TCGA_ID %in% c(
  "BF3139XR26B","BF3439XR26B","BF3463XR26B","BF3596XR26C", #FNA-BV=4
  "BF3380XR26B","BF3398XR26B","BF3453XR26B") #FNA-Parotid=3
  & !df_fna_or$`Sample.Type(Sub1)`=="FNA-Others",]
table(df_fna_or$`Sample.Type(Sub1)`)
# FNA-ABDO FNA-Plasma   FNA-PLEU  FNA-Spine 
#      124         78         72          8
df_fna_or$paired <- c(rep(1:(124/2), each = 2),rep(1:(78/2), each = 2),
                      rep(1:(72/2), each = 2),rep(1:(8/2), each = 2))

df_fna_abdo <- df_fna_or[df_fna_or$`Sample.Type(Sub1)`=="FNA-ABDO",]
df_fna_abdo$ref_grp <- ifelse(df_fna_abdo$ref_tp=="ABDO",1,2)

df_fna_abdo %>% wilcox_test(Tumor_Purity2 ~ ref_tp, paired = TRUE)
#1 Tumor_Purity2 ABDO   Plasma    62    62       358 0.245

p_fna_abdo <- ggplot(data=df_fna_abdo, 
                     aes(x = ref_tp, y = Tumor_Purity, color = ref_tp)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#00BFC4","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65", alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.245", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",
       title="FNA-liver/pancreas/stomach/peritoneum/pelvis samples")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_fna_abdo

#####################
df_fna_pleu <- df_fna_or[df_fna_or$`Sample.Type(Sub1)`=="FNA-PLEU",]
df_fna_pleu$ref_grp <- ifelse(df_fna_pleu$ref_tp=="PLEU",1,2)
df_fna_pleu$ref_tp.f <- factor(as.factor(df_fna_pleu$ref_tp),order=T,
                               levels=c("PLEU","Plasma"),
                               labels=c("PLEU","Plasma")) 
df_fna_pleu_or <- df_fna_pleu[order(df_fna_pleu$ref_tp.f),]

df_fna_pleu_or %>% wilcox_test(Tumor_Purity2 ~ ref_tp.f, paired = TRUE)
#1 Tumor_Purity2 PLEU   Plasma    36    36       142  0.59

p_fna_pleu <- ggplot(data=df_fna_pleu_or, 
                     aes(x = ref_tp.f, y = Tumor_Purity, color = ref_tp.f)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#b890e0","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65", alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.590", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",title="FNA-lung/mediastinum/chest wall samples")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_fna_pleu

#####################

df_fna_plasma <- df_fna_or[df_fna_or$`Sample.Type(Sub1)`=="FNA-Plasma",]
df_fna_plasma_or <- df_fna_plasma[order(df_fna_plasma$TCGA_ID),]
df_fna_plasma$ref_grp <- ifelse(df_fna_plasma$ref_tp=="LN",1,2)


df_fna_plasma %>% wilcox_test(Tumor_Purity2 ~ ref_tp, paired = TRUE)
#1 Tumor_Purity2 LN     Plasma    39    39       208   0.1

p_fna_plasma <- ggplot(data=df_fna_plasma, 
                       aes(x = ref_tp, y = Tumor_Purity, color = ref_tp)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#e98c61","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65", alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.1", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",title="FNA-LN")+
  theme(
        legend.position="none",
        legend.title=element_blank()) 
p_fna_plasma

#####################
df_fna_csf <- df_fna_or[df_fna_or$`Sample.Type(Sub1)`=="FNA-Spine",]
df_fna_csf$ref_grp <- ifelse(df_fna_csf$ref_tp=="CSF",1,2)

df_fna_csf %>% wilcox_test(Tumor_Purity2 ~ ref_tp, paired = TRUE)
#1 Tumor_Purity2 CSF    Plasma     4     4         0 0.181

p_fna_csf <- ggplot(data=df_fna_csf, 
                    aes(x = ref_tp, y = Tumor_Purity, color = ref_tp)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#00A9FF","#f71010")) +
  geom_line(aes(group = paired), linetype=1, linewidth=1, color="gray65", alpha=0.4)+ 
  geom_point(aes(fill=ref_tp,group=paired),size=2,shape=16) +
  annotate("text", x = 1.5, y = 1.0, 
           label = "p = 0.181", size = 4, color = "black") +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels=seq(0,1, by = 0.25),
                     limits=c(0,1.15), expand=c(0,0)) +
  labs(y="Tumor fraction",x="Reference",title="FNA-spine tissue samples")+
  theme(legend.position="none",
        legend.title=element_blank()) 
p_fna_csf


pline_m <- ggarrange(p_abdo,p_pleu,p_fna_abdo,p_fna_pleu,p_fna_plasma,p_fna_csf,
                     common.legend = F, legend= NULL,
                     ncol = 3, nrow = 2)
pline_m

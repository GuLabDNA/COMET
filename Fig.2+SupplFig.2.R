#Fig. 2 and Suppl Fig. 2

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
library(ggforce) # for 'geom_arc_bar'
library(ggpubr) #ggarrange
library(beeswarm)
library(RColorBrewer)
library(ggbeeswarm)
library(gridExtra)
library(grid)
library(ggalluvial)

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"

######* load data *######
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

######* Stacked bar chart of distributions *######
###### Distribution of Sample Set 1 ######
res_diag3 <- res_diag_comb1[res_diag_comb1$sample_tp!="CSF",] #413

###### Distribution of longitudinal samples of Sample set 1 ######
# res_diag3_longi <- res_diag3[res_diag3$idn3>=3000 & res_diag3$idn3<=3215 &
#                                res_diag3$sample_tp!="FNA" & !is.na(res_diag3$BF_id),]
# res_diag3 <- res_diag3_longi

###### Distribution of Sample Set 2 ######
# res_diag3_csf <- res_diag_comb1[res_diag_comb1$sample_tp=="CSF",] #106
# res_diag3 <- res_diag3_csf


df_cyto <- as.data.frame(table(res_diag3$cyto_catx2, res_diag3$diag_caco))
df_sum <- df_cyto %>% 
  group_by(Var2) %>% 
  summarize(sum = sum(Freq))
df_cyto$Method <- "Cytology test"

df_cyto2 <- df_cyto %>% 
  mutate(sum = rep(as.vector(df_sum$sum), each=5),
         perc=round(Freq/sum*100, 1),
         n_perc=paste0(Freq," (",perc,")"))
df_cyto2_trans <- data.table::dcast(setDT(df_cyto2), Var1 ~ Var2,value.var = "n_perc")

df_cnv <- as.data.frame(table(res_diag3$CNA, res_diag3$diag_caco))
df_cnv$Method <- "CNA analysis"

df_cnv2 <- df_cnv %>% 
  mutate(sum = rep(as.vector(df_sum$sum), each=2),
         perc=round(Freq/sum*100, 1),
         n_perc=paste0(Freq," (",perc,")"))
df_cnv2_trans <- data.table::dcast(setDT(df_cnv2), Var1 ~ Var2,value.var = "n_perc")
a1<- data.frame(Var1="Cytology test", Cancers="", Control="")
a2<- data.frame(Var1="CNA test", Cancers="", Control="")

#dfm2_trans_all <- rbind(a1,df_cyto2_trans,a2,df_cnv2_trans)

dfall_vis <- rbind(df_cyto,df_cnv)
colnames(dfall_vis) <- c("Results","Diagnosis","Count","Method")
dfall_vis <- dfall_vis %>%
  mutate(dx_cat.f = factor(as.factor(Diagnosis),order=T,
                           levels=c("Cancers", "Control"),
                           labels=c("Cancer", "Control")), 
         Method.f = factor(as.factor(Method),order=T,
                           levels=c("Cytology test","CNA analysis"),
                           labels=c("Cytology test","CNA analysis")))

# Create a bar plot visualizing the counts
p_bar2 <- ggplot(dfall_vis, aes(x = Results, y = Count, fill = dx_cat.f)) +
  geom_bar(stat = "identity", position = "stack",width=0.6) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("#DC3C22","#DBE4C9"))+ 
  facet_wrap(~ Method.f, scales = "free_x", ncol = 2) +
  labs(x = NULL, y = "Count", fill = "Diagnosis") +
  scale_y_continuous(breaks=seq(0,100,by=20),
                     labels=seq(0,100,by=20),
                     limits=c(-0,101), expand=c(0,0)) +
  theme_minimal()
p_bar2


##### boxplot of tumor fraction vs. cyto ######
res_diag3 <- res_diag_comb1[res_diag_comb1$sample_tp!="CSF",] #413

res_diag3_rm <- res_diag3[res_diag3$cyto_catx2!="Non-diagnostic",] #1
res_diag3_rm$cyto_catx2 <- ifelse(res_diag3_rm$cyto_catx2=="Negative","Benign",res_diag3_rm$cyto_catx2)
res_diag3_rm$cyto_catx2.f <- factor(as.factor(res_diag3_rm$cyto_catx2),order=T,
                                    levels=c("Malignant","Benign","Atypical","Suspicious"),
                                    labels=c("Malignant","Benign","Atypical","Suspicious"))

p_cytocat <- 
  ggplot(data=res_diag3_rm, aes(x = cyto_catx2.f, y = Tumor_Purity_ichor, color=cyto_catx2.f)) +
  geom_boxplot() +
  geom_beeswarm(priority="density",cex=1.5,size=1.7,corral = "wrap")+
  labs(y="Tumor fraction",x="Cytology test")+
  scale_color_manual(values = c("#ED775A","#B2CD9C","#FAD691","#5BC0DE"))+ 
  scale_y_continuous(breaks=seq(0,1,by=0.25),
                     labels=seq(0,1,by=0.25),
                     limits=c(0,1.7), expand=c(0,0)) +
  theme() 
p_cytocat

######* Sankey plot of Cancers *#######

###### cyto vs. CNA for Sample Set 1 ######
res_diag3 <- res_diag_comb1[res_diag_comb1$sample_tp!="CSF",] #413

###### cyto vs. CNA for longitudinal samples of Sample set 1 ######
# res_diag3_longi <- res_diag3[res_diag3$idn3>=3000 & res_diag3$idn3<=3215 &
#                                res_diag3$sample_tp!="FNA" & !is.na(res_diag3$BF_id),]
# res_diag3 <- res_diag3_longi

###### cyto vs. CNA for Sample Set 2 ######
# res_diag3_csf <- res_diag_comb1[res_diag_comb1$sample_tp=="CSF",] #106
# res_diag3 <- res_diag3_csf

df_allu <- res_diag3[,c("BF_id","CNV","cyto_catx2","dx_cat.f")]

df_allu2 <- df_allu[df_allu$dx_cat.f %in% c("Cancers"),] %>% 
  mutate(cnv=ifelse(CNV=="Pos","Positive","Negative"),
         cyto_catx2=ifelse(cyto_catx2=="Negative","Benign",cyto_catx2)) %>% select(-CNV)

df_cnv <- data.frame(BF_id=df_allu2$BF_id,response=df_allu2$cnv,
                     group="CNA analysis")
df_cyto <- data.frame(BF_id=df_allu2$BF_id,response=df_allu2$cyto_catx2,
                      group="Cytology test")

df_comb <- rbind(df_cyto,df_cnv)
df_comb$freq=1

response_order <- c("Malignant","Benign","Atypical",
                    "Suspicious","Non-diagnostic","Positive", "Negative") # Example order

df_comb2 <- df_comb %>%
  mutate(
    group = factor(group, levels = c("Cytology test", "CNA analysis")),
    response = factor(response, levels = response_order) # Order 'response' factor
  ) %>% 
  group_by(BF_id, group, response) %>%
  summarise(freq = n(), .groups = 'drop')

# Plotting with counts on strata
p_allu2 <- ggplot(df_comb2,
                  aes(x = group, stratum = response, alluvium = BF_id, y = freq,
                      fill = response)) +
  geom_stratum(alpha = 0.8) +
  geom_flow(alpha = 0.7)+
  geom_flow()+
  # Add counts to the strata
  geom_text(stat = "stratum", aes(label = after_stat(count)),
            size = 3, color = "black", nudge_x = 0 # Adjust to move labels horizontally
  ) +
  geom_text(stat = "flow", aes(label = after_stat(count)), # Label with frequency
            nudge_x = 0.2, size = 3, color = "black") + # Set label size and color
  scale_fill_manual(values = c(
    "Atypical" =  "#FAD691", "Benign" = "#B2CD9C",
    "Malignant" = "#D6A99D",  "Suspicious" = "#5BC0DE", "Non-diagnostic"="grey70",
    "Positive" = "#D25D5D", "Negative" = "#5E936C")) + #F3A26D "#D6A99D" #"#F08B51"
  scale_x_discrete(expand = c(.15, .05)) +
  labs(x = NULL,y = NULL,fill = "Response Status") +
  theme_minimal() 
p_allu2

######* Sankey plot in Controls *######
df_allu3 <- df_allu[df_allu$dx_cat.f %in% c("Control"),] %>% 
  mutate(cnv=ifelse(CNV=="Pos","Positive","Negative"),
         cyto_catx2=ifelse(cyto_catx2=="Negative","Benign",cyto_catx2)) %>% select(-CNV)

df_cnv <- data.frame(BF_id=df_allu3$BF_id,response=df_allu3$cnv,
                     group="CNA analysis")
df_cyto <- data.frame(BF_id=df_allu3$BF_id,response=df_allu3$cyto_catx2,
                      group="Cytology test")

df_comb <- rbind(df_cyto,df_cnv)
df_comb$freq=1

df_comb3 <- df_comb %>%
  mutate(
    group = factor(group, levels = c("Cytology test", "CNA analysis")),
    response = factor(response, levels = response_order) # Order 'response' factor
  ) %>% 
  group_by(BF_id, group, response) %>%
  summarise(freq = n(), .groups = 'drop')

# Plotting with counts on strata
p_allu3 <- ggplot(df_comb3,
                  aes(x = group, stratum = response, alluvium = BF_id, y = freq,
                      fill = response)) +
  geom_stratum(alpha = 0.8) +
  geom_flow(alpha = 0.7)+
  geom_flow()+
  # Add counts to the strata
  geom_text(stat = "stratum", aes(label = after_stat(count)),
            size = 3, color = "black", nudge_x = 0 # Adjust to move labels horizontally
  ) +
  geom_text(stat = "flow", aes(label = after_stat(count)), # Label with frequency
            nudge_x = 0.2, size = 3, color = "black") + # Set label size and color
  scale_fill_manual(values = c(
    "Atypical" =  "#FAD691", "Benign" = "#B2CD9C",
    "Malignant" = "#D6A99D",  "Suspicious" = "#5BC0DE", "Non-diagnostic"="grey70",
    "Positive" = "#D25D5D", "Negative" = "#5E936C")) + 
  scale_x_discrete(expand = c(.15, .05)) +
  labs(x = NULL, y = NULL, fill = "Response Status") +
  theme_minimal()
p_allu3

p_allu_allm <- ggarrange(p_allu2,p_allu3,
                         common.legend = T, legend="right",
                         ncol = 1, nrow = 2, heights = c(3, 1))
p_allu_allm

############## Suppl Fig. 2 #################
res_cnv_chk <- within(res_diag_comb1, {
  purity_cut03 <- NA
  purity_cut03[Tumor_Purity_ichor >=0.03] <- "Positive"
  purity_cut03[Tumor_Purity_ichor<0.03] <- "Negative"
  
  purity_cut05 <- NA
  purity_cut05[Tumor_Purity_ichor >=0.05] <- "Positive"
  purity_cut05[Tumor_Purity_ichor<0.05] <- "Negative"
})

df_cnv03 <- data.frame(BF_id=res_cnv_chk$BF_id,response=res_cnv_chk$purity_cut03,
                     group="CNA analysis (cutoff 3%)")
df_cnv05 <- data.frame(BF_id=res_cnv_chk$BF_id,response=res_cnv_chk$purity_cut05,
                      group="CNA analysis (cutoff 5%)")

df_standard <- data.frame(BF_id=res_cnv_chk$BF_id,response=res_cnv_chk$dx_cat.f,
                          group="Diagnosis")
df_comb <- rbind(df_cnv03,df_standard,df_cnv05)
df_comb$freq=1

response_order <- c("Cancers","Control","Benign","Unclear","Positive", "Negative") # Example order

df_comb2 <- df_comb %>%
  mutate(
    group = factor(group, levels = c("CNA analysis (cutoff 3%)","Diagnosis", "CNA analysis (cutoff 5%)")),
    response = factor(response, levels = response_order) # Order 'response' factor
  ) %>% 
  group_by(BF_id, group, response) %>%
  summarise(freq = n(), .groups = 'drop')

# Plotting with counts on strata
p_allu1 <- ggplot(df_comb2,
                  aes(x = group, stratum = response, alluvium = BF_id, y = freq,
                      fill = response)) +
  geom_stratum(alpha = 0.8) +
  geom_flow(alpha = 0.7)+
  geom_flow()+
  # Add counts to the strata
  geom_text(stat = "stratum", aes(label = after_stat(count)),
            size = 3, color = "black", nudge_x = 0 # Adjust to move labels horizontally
  ) +
  geom_text(stat = "flow", aes(label = after_stat(count)), # Label with frequency
            nudge_x = 0.2, size = 3, color = "black") + # Set label size and color
  scale_fill_manual(values = c(
    "Cancers"="#D1A980","Control"="#DEE8CE","Benign" = "#B2CD9C",
    "Unclear"="grey70",
    "Positive" = "#D25D5D", "Negative" = "#5E936C")) + #F3A26D "#D6A99D" #"#F08B51"
  scale_x_discrete(expand = c(.15, .05)) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Response Status") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
p_allu1

df_comb <- rbind(df_cnv03,df_cnv05)
df_comb$freq=1

df_comb3 <- df_comb %>%
  mutate(
    group = factor(group, levels = c("CNA analysis (cutoff 3%)", "CNA analysis (cutoff 5%)")),
    response = factor(response, levels = response_order) # Order 'response' factor
  ) %>% 
  group_by(BF_id, group, response) %>%
  summarise(freq = n(), .groups = 'drop')

# Plotting with counts on strata
p_allu3 <- ggplot(df_comb3,
                  aes(x = group, stratum = response, alluvium = BF_id, y = freq,
                      fill = response)) +
  geom_stratum(alpha = 0.8) +
  geom_flow(alpha = 0.7)+
  geom_flow()+
  # Add counts to the strata
  geom_text(stat = "stratum", aes(label = after_stat(count)),
            size = 3, color = "black", nudge_x = 0 # Adjust to move labels horizontally
  ) +
  geom_text(stat = "flow", aes(label = after_stat(count)), # Label with frequency
            nudge_x = 0.2, size = 3, color = "black") + # Set label size and color
  scale_fill_manual(values = c(
    "Atypical" =  "#FAD691", "Benign" = "#B2CD9C",
    "Malignant" = "#D6A99D",  "Suspicious" = "#5BC0DE", "Non-diagnostic"="grey70",
    "Positive" = "#D25D5D", "Negative" = "#5E936C")) + 
  scale_x_discrete(expand = c(.15, .05)) +
  labs(x = NULL, y = NULL, fill = "Response Status") +
  theme_minimal()
p_allu3

res_cnv_chk$pur_perc <- round(res_cnv_chk$Tumor_Purity_ichor*100,1)

p_pur <- ggplot(data=res_cnv_chk) +
  geom_boxplot(aes(x = dx_cat.f, y = pur_perc),width=0.4) + 
  geom_beeswarm(aes(x = dx_cat.f, y = pur_perc),
                cex=6, priority = "density", corral = "wrap") +
  geom_hline(yintercept=3,color = "grey70", linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept=5,color = "grey70", linetype="dashed", linewidth=0.5)+
  scale_y_log10(
    breaks = c(1, 10, 100), # Define major tick marks
    labels = c("1", "10", "100"), # Define labels for major ticks
    guide = guide_axis_logticks() )+# Add minor log ticks
  labs(y="Tumor fraction (%)",x="Group")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(linewidth=1.1,color="black",fill=NA),
        #text=element_text(size=15, face="bold"),
        axis.text.x = element_text(face="bold",angle=40,hjust=1),
        axis.text.y = element_text(face="bold"),
        axis.title=element_text(face="bold"),
        legend.position="none",
        legend.title=element_blank(),
        panel.spacing = unit(0.2, "lines")) 
p_pur

# Add a square area using annotate()
p_pur2 <- p_pur + annotate("rect", xmin = 0, xmax = 5, ymin = 3, ymax = 5,
             alpha = 0.3, fill = "red")
p_pur2

p_pur_all <- ggarrange(p_allu1,p_pur2,
                         common.legend = F,
                         ncol = 2, nrow = 1, widths = c(1.5, 1))
p_pur_all

pdf(file=file.path(pathbf2,"Output",
                   paste0(str_sub(Sys.Date(), start = 3),
                          "-all519samples_CNAcutoff.pdf")), 
    width = 7, height = 4, useDingbats = FALSE)
p_pur_all
dev.off()


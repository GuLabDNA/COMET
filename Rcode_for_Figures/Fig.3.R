# Make plots for controlled FNA samples (supernatant, pellet, cytology)
# Date: 2025-11-11

pathbf2 <- "/Users/jingruyu/Documents/Studies/BF_deconv"
library(dplyr)
library(tidyverse)
library(stringr)
library(openxlsx)
library(data.table) #detDT
library(scales)
library("writexl") 
library(matrixStats) #colMeans, colMedians
library(rstatix) #wilx

library(ggplot2)
library(ggpubr) #ggarrange
library(ggrepel) #geom_text_repel
library(forcats)
library(UpSetR)

############################################
# FNA controls
cofna <- read.xlsx(file.path(pathbf2, "Draft",
                             "25-1202-Supple.Tables.xlsx"),
                   sheet=8, skipEmptyRows=FALSE, colNames = TRUE)
colnames(cofna)
cofna1 <- within(cofna, {
cyto <- NA
cyto[Cytology.from.ACL == "Atypical"] <- "Atypical"
cyto[Cytology.from.ACL == "Suspicious"] <- "Suspicious"
cyto[Cytology.from.ACL == "Negative"] <- "Benign"
cyto[Cytology.from.ACL %in% c("Positive","Adenocarcinoma")] <- "Malignant"
cyto[Cytology.from.ACL %in% c("Insufficient")] <- "Insufficient cells/non-diagnostic"
})

cofna1_ex <- cofna1[!is.na(cofna1$cyto),]
cofna.cyto_ex <- cofna1_ex[,c("BF","CNV","Tumor.fraction","FNA.name","Sample.type",
                              "Diagnosis","Sample.group","Distance.from.margin",
                              "Cytology.from.ACL","cyto")]
cofna.cyto_ex$Sample.group="Cytology"
cofna.cyto_ex$CNV <- cofna1_ex$cyto

cofna.super_ex <- cofna1[cofna1$Sample.group=="Supernatant",
                         c("BF","CNV","Tumor.fraction","FNA.name","Sample.type",
                           "Diagnosis","Sample.group","Distance.from.margin",
                           "Cytology.from.ACL","cyto")]
cofna.pellet_ex <- cofna1[cofna1$Sample.group=="Pellet",
                          c("BF","CNV","Tumor.fraction","FNA.name","Sample.type",
                            "Diagnosis","Sample.group","Distance.from.margin",
                            "Cytology.from.ACL","cyto")]
cofna.pellet_ex$BF <- cofna.super_ex$BF
cofna2 <- rbind(cofna.cyto_ex,cofna.super_ex,cofna.pellet_ex)
table(cofna2$CNV)
cofna3 <- within(cofna2, {
  res_cyto.cna <- NA
  res_cyto.cna[CNV == "Atypical"] <- "Atypical"
  res_cyto.cna[CNV == "Suspicious"] <- "Suspicious"
  res_cyto.cna[CNV == "Benign"] <- "Benign"
  res_cyto.cna[CNV == "Malignant"] <- "Malignant"
  res_cyto.cna[CNV == "Insufficient cells/non-diagnostic"] <- 
    "Insufficient cells/non-diagnostic"
  res_cyto.cna[CNV == "N/a"] <- "Insufficient DNA"
  res_cyto.cna[CNV == "Pos"] <- "Postive"
  res_cyto.cna[CNV == "Neg"] <- "Negative"
  })

cofna3 <- cofna3 %>%
  mutate(BF_id=paste0("BF", BF),
         BF_id.f=factor(BF_id, levels = rev(unique(BF_id))),
         FNA_id.f=factor(FNA.name, order=T,
         levels = c("1_Pancreas_PAAD","2_Pancreas_pNET","3_Rectum_CRC"),
         labels=c("1_Pancreas_PAAD","2_Pancreas_pNET","3_Rectum_CRC")),
         FNA_id.grp=paste(FNA.name, Sample.group, sep = "-"),
         FNA_id.grp.f=factor(FNA_id.grp, levels = rev(unique(FNA_id.grp))),
         #Distance=ifelse(Distance==0.1,0.12,ifelse(Distance==0,0.1,Distance)))
         Sample.group.f=factor(as.factor(Sample.group), order=T,
                             levels=c("Pellet","Supernatant","Cytology"),
                             labels=c("Pellet CNA","Supernatant CNA","Cytology")),
         Distance_t = factor(as.factor(Distance.from.margin), order=T,
                             levels=c(0,0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25),
                             labels=c("0","Adjacent","0.5","1","2","3","4","5",
                                      "6","7","8","10","10","15","20","25")),
         tumor_frac=round(Tumor.fraction*100,1)
         )


#heatmap
cofna4 <- cofna3[cofna3$Distance.from.margin!=9,]

p_heat <- ggplot(cofna4, aes(x = Distance_t, y = Sample.group.f, fill = res_cyto.cna)) +
  geom_tile(color = "white") +
  scale_fill_manual(values=c("Malignant"="#DC3C22","Benign"="#B2CD9C",
                              "Atypical"="#FAD691", "Suspicious"="#5BC0DE",
                              "Insufficient cells/non-diagnostic"="grey75", 
                              "Insufficient DNA"="grey75",
                              "Postive"="#DC3C22","Negative"="#B2CD9C"))+
  facet_wrap(~ FNA_id.f, ncol = 1) +
  labs(fill = "Test results", y = "Test type", x = "Distance from tumor (cm)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color="black",fill=NA),
        axis.text.x = element_text(size=13,face="bold",angle = 45, hjust = 1),
        axis.text.y = element_text(size=13, face="bold"),
        axis.title=element_text(size=13,face="bold"),
        legend.position="bottom",
        legend.title=element_text(face = "bold"))
p_heat

p_line <- ggplot(cofna4, aes(x = Distance_t, y = tumor_frac, 
                          group = Sample.group.f, color = Sample.group.f)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 5, linetype = 2,color = "gray", linewidth = 0.5) +
  scale_y_continuous(breaks = seq(0,100,by=25),
                     labels = seq(0,100,by=25),limits=c(0,100),expand=c(0,0)) +
  scale_color_manual(values=c("Supernatant CNA"="black","Pellet CNA"="#5E936C"))+
  labs(x="Distance from tumor (cm)", y="Tumor fraction (%)",color = "Source") +
  facet_wrap(~ FNA_id.f, ncol = 1) +
  theme(
    legend.position="bottom",
    legend.title=element_text(face = "bold")
  )
p_line

library(patchwork)
comb_plot <- p_heat + p_line
#line chart 1/3 the width of the heatmap:
final_plot <- comb_plot + plot_layout(widths = c(3, 2.8))
final_plot


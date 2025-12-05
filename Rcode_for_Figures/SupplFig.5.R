#SupplFig.5

library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(openxlsx)
library(data.table) #detDT
library(scales)

library(ggpubr) #ggarrange
library(RColorBrewer)

# EPIC
load(file=file.path(pathref,"RefData/K562_EPIC", "k562_epic_hg38.RData")) 
#k562_hg38_betas

#0-based data, with CpGs from top and bottom strand as an entity
k562_epic1 <- k562_hg38_betas[,c(2,3,4,5)]
colnames(k562_epic1)[4] <- "beta"

k562_epic1_chr <- k562_epic1 %>%
  mutate(chr_cpg_loc = paste0(chr,"_",start,"_",end))

# FLEXseq
f="K562LG50ngXR28_comb"
k562_lg50 <- fread(file.path(pathK562,
             paste0(f,"-merged.cov.SNPfiltered_clean.txt")),quote="",sep="\t")
colnames(k562_lg50) <- c("chr","start","end","meth.cnt","depth")

k562_lg50_2 <- as.data.frame(k562_lg50) %>% filter(depth>=15) %>%
  mutate(chr_cpg_loc = paste0(chr,"_",start,"_",end),
         beta=round(meth.cnt/depth, digits=5))

# FLEXseq intersected with EPIC
beta_inter <- as.data.frame(fread(
  (file.path(pathK562,"K562LG50ng_interBlocks_15maxcov.txt")),
  header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(beta_inter) <- c("chr","start","end","chr_marker_loc",
                          "meth.cnt.x","total.x","chr_cpg_loc","beta")
k562int_15cov <- 
  beta_inter[!is.na(beta_inter$beta),c("chr","start","end","meth.cnt.x","total.x",
                                       "chr_cpg_loc","beta")]

colnames(k562int_15cov)[4]="meth.cnt"
colnames(k562int_15cov)[5]="total"

#make heatmap
p_flex.epic <- CorrelatePlot3(samp1=k562_lg50_2, samp2=k562_epic1_chr, 
                              assay_xr='FLEXseq (original)', assay="EPIC", covnum=15)
p_flex.epic_1 <- p_flex.epic + ggtitle("Pearson's r = 0.959")+
  theme(plot.title = element_text(hjust = 0.5))

p_flexInter.epic <- CorrelatePlot3(samp1=k562int_15cov, samp2=k562_epic1_chr, 
                                   assay_xr='FLEXseq (intersected)', assay="EPIC", covnum=15)
p_flexInter.epic_1 <- p_flexInter.epic + ggtitle("Pearson's r = 0.922")+
  theme(plot.title = element_text(hjust = 0.5))

p_flexInter.or <- CorrelatePlot3(samp1=k562int_15cov, samp2=k562_lg50_2, 
                                 assay_xr='FLEXseq (intersected)', assay="FLEXseq (orginal)", covnum=15)
p_flexInter.or_1 <- p_flexInter.or + ggtitle("Pearson's r = 0.965")+
  theme(plot.title = element_text(hjust = 0.5))

pm_k562 <- ggarrange(p_flex.epic_1,p_flexInter.epic_1,p_flexInter.or_1,
                     common.legend = TRUE, legend= "bottom",
                     nrow = 1)
pm_k562 

gp1 <- 
  read.xlsx(xlsxFile = file.path(pathref,"Output/coverage", 
                                 "K562LG50ngXR28_comb-intBlocks_graph.xlsx"),
            sheet = 3, skipEmptyRows = TRUE)
#stack by grp
gp2 <- gp1 %>%
  filter(Pack_grp %in% c("Original","Total")) %>% filter(Group!="15X") %>%
  mutate(Pack_grp.f = factor(as.factor(Pack_grp), order=T,
                             levels=c("Total","Original"),
                             labels=c("Increased by intersection","Original")),
         Group.f = factor(as.factor(Group), order=T,
                          levels=c("1X","5X","10X","20X","30X"),
                          labels=c("1X","5X","10X","20X","30X"))
  ) %>%
  arrange(desc(Value))

# The transformation factor
transf_fact <- max(gp2$Value)/max(gp2$Perc)

p1 <- ggplot(data=gp2,  mapping = aes(x=Group.f, y=Value, 
                                      fill = Pack_grp.f, color = Pack_grp.f)) + 
  geom_bar(stat = "identity", width = 0.5, position = position_identity(), color = NA) + 
  geom_line(aes(y = transf_fact*Perc, group=Pack_grp.f)) + 
  geom_point(aes(y = transf_fact*Perc, group=Pack_grp.f), shape=16, size=1.6, color="black") + 
  # Add second OY axis; note the transformation back (division)
  scale_y_continuous(name=expression("CpG Frequency (10e6)"),
                     sec.axis = sec_axis(trans = ~ . / transf_fact,
                                         name = "Percentage (%)"), limits=c(0,35)) +
  #geom_hline(yintercept = c(5,10,15,20), color="grey") +
  labs(x="Coverage") +
  scale_fill_manual(values=c("red","#0042F9"))+
  scale_color_manual(values = c("red","#0042F9"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color="black",fill=NA),
        axis.text.x = element_text(size=13, angle=0, hjust = 1, face="bold"),
        axis.text.y = element_text(size=13, face="bold"),
        axis.title.x = element_text(size=13, face="bold"),
        axis.title.y = element_text(size=13, face="bold"),
        axis.title.y.right = element_text(size=13, face="bold"),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=10)
  ) 
p1

# dev.new()  
# pdf(file=file.path(pathref,"Output/coverage", "K562LG50ngXR28_comb-intBlocks_graph.pdf"), 
#     width = 3.5, height = 3.8, useDingbats = FALSE)
# p1
# dev.off()

gp3 <- 
  read.xlsx(xlsxFile = file.path(pathref,"Output/coverage", 
                                 "K562LG50ngXR28_comb-intBlocks_graph.xlsx"),
            sheet = 4, skipEmptyRows = TRUE)
#stack by grp
gp4 <- gp3 %>%
  filter(Pack_grp %in% c("Original","Total")) %>%
  mutate(Pack_grp.f = factor(as.factor(Pack_grp), order=T,
                             levels=c("Total","Original"),
                             labels=c("Increased by intersection","Original")),
         Cov.f = factor(as.factor(Cov), order=T,
                        levels=c("5","10","15"),
                        labels=c("5X","10X","15X"))) %>%
  arrange(desc(Value))

# The transformation factor
transf_fact <- max(gp4$Value)/max(gp4$Perc)

p2 <- ggplot(data=gp4,  mapping = aes(
  x=Cov.f, y=Value, fill = Pack_grp.f, color = Pack_grp.f)) + 
  geom_bar(stat = "identity", width = 0.5, 
           position = position_identity(), color = NA) + 
  geom_line(aes(y = transf_fact*Perc, group=Pack_grp.f)) + 
  geom_point(aes(y = transf_fact*Perc, group=Pack_grp.f), 
             shape=16, size=1.6, color="black") + 
  # Add second OY axis; note the transformation back (division)
  scale_y_continuous(name=expression("CpG Frequency (10e3)"),
                     sec.axis = sec_axis(transform = ~ . / transf_fact,
                                         name = "Percentage (%)"), limits=c(0,40)) +
  #geom_hline(yintercept = c(5,10,15,20), color="grey") +
  labs(x="Coverage") +
  scale_fill_manual(values=c("red","#0042F9"))+
  scale_color_manual(values = c("red","#0042F9"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color="black",fill=NA),
        axis.text.x = element_text(size=13, face="bold"),
        axis.text.y = element_text(size=13, face="bold"),
        axis.title.x = element_text(size=13, face="bold"),
        axis.title.y = element_text(size=13, face="bold"),
        axis.title.y.right = element_text(size=13, face="bold"),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=10)
  ) 
p2

dev.new()
pdf(file=file.path(pathref,"Output/coverage", "K562LG50ngXR28_comb-intBlocks.vsTCGA_graph.pdf"),
    width = 3.5, height = 3.8, useDingbats = FALSE)
p2
dev.off()

pm_k562_cov <- ggarrange(p1,p2, common.legend = T, 
                         legend="bottom", nrow = 1, widths=c(1,0.7))

dev.new()  
pdf(file=file.path(pathref,"Output/coverage", "K562LG50ngXR28_comb-intBlocks_coverage_graph.pdf"), 
    width = 8, height = 4, useDingbats = FALSE)
pm_k562_cov
dev.off()

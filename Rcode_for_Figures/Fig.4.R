#Fig.4

res_mat1 <- read.xlsx(file.path(pathbf2,"clinical_info",
                                "T9.413BFs-score.comp.k5_TCGA_136BFs_matrix.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
res_mat1 <- res_mat1[994:1129,]
res_mat2 <- read.xlsx(file.path(pathbf2,"Clinical_info",
                                "T9.413BFs-score.comp.k5_TCGA_onlySamples.136BFs_TCGA.XR.xlsx"),
                      sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
res_mat2 <- res_mat2[1:136,]

res_matrix <- merge(res_mat1,res_mat2[,c("Sentrix_ID","TCGA_Project","purity.abs")],
                    by="Sentrix_ID")


res_matrix <- within(res_matrix, {
  TCGA_Project2 <- NA
  TCGA_Project2[TCGA_Project == "MCF GYN"] <- "MCF.GYN"
  TCGA_Project2[TCGA_Project == "MCF UGI-AD"] <- "MCF.UGI-AD"
  TCGA_Project2[TCGA_Project == "MCF LGI-AD"] <- "MCF.LGI-AD"
  TCGA_Project2[TCGA_Project %in% c("MCF SCC","SCC")] <- "MCF.SCC"
  TCGA_Project2[TCGA_Project == "MCF GBM"] <- "MCF.GBM"
})

res_matrix$TCGA_Project2 <- 
  ifelse(is.na(res_matrix$TCGA_Project2),res_matrix$TCGA_Project,res_matrix$TCGA_Project2)

classes <- c("BRCA","DLBC","LUAD","MCF.GYN","MCF.SCC","MCF.LGI-AD","MCF.UGI-AD","PAAD")
res_df <- read.xlsx(file.path(pathbf2,"Clinical_info",
                              "T9.413BFs-score.comp.k5_TCGA_onlySamples.136BFs_TCGA.XR.xlsx"),
                    sheet=1, skipEmptyRows=FALSE, colNames = TRUE)
res_df <- res_df[,c("TCGA_ID","Sentrix_ID","top1_cal_label","top1_cal_score")]
res_df2 <- merge(res_df,res_matrix,by="Sentrix_ID")
res_matrix_8can <- res_df2[res_df2$TCGA_Project2 %in% classes,]
res_matrix_8can <- within(res_matrix_8can, {
  purity_catn2 <- NA
  purity_catn2[purity.abs >=0.7] <- 1
  purity_catn2[purity.abs>=0.5 & purity.abs<0.7] <- 2
})
res_matrix_8can <- res_matrix_8can %>% 
  mutate(TCGA_Project2.f = factor(as.factor(TCGA_Project2), order=T,
                                  levels=c("BRCA", "DLBC", "LUAD", "MCF.GYN", "MCF.UGI-AD", "MCF.LGI-AD","MCF.SCC", "PAAD"),
                                  labels=c("BRCA", "DLBC", "LUAD", "MCF GYN", "MCF UGI-AD", "MCF LGI-AD","MCF SCC", "PAAD")),
         purity_catn2.f = factor(as.factor(purity_catn2), order=T,
                                 levels=c(1,2),
                                 labels=c(">=0.7", "0.5-0.7")))
res_matrix_8can_rm.co <- res_matrix_8can[!res_matrix_8can$top1_cal_label=="Control",]

######* boxplot of ML classifier scores *######
p_score <- 
  ggplot(data=res_matrix_8can_rm.co, aes(x = TCGA_Project2.f, y = top1_cal_score,fill=purity_catn2.f)) +
  geom_boxplot(width=0.7) +
  geom_jitter(aes(color = purity_catn2.f), 
              position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.7), 
              size = 2.6,alpha=0.8) +# Adds jittered points with a stroke
  scale_fill_manual(values = c("#00BFC4","#B2CD9C")) + 
  scale_color_manual(values = c("black","black")) +
  labs(y="Top ML classifier score",x="Group",fill = "Tumor fraction")+
  scale_y_continuous(breaks=seq(0,1, by = 0.2),
                     labels=seq(0,1, by = 0.2),
                     limits=c(0, 1.05), expand=c(0,0)) +
  theme(legend.position = "bottom")
p_score


#compare scores
classes <- c("BRCA", "DLBC", "LUAD", "MCF.GYN", "MCF.UGI-AD", "MCF.LGI-AD","MCF.SCC", "PAAD")
res_matrix_8can_rm.co %>% 
  filter(TCGA_Project2==classes[i]) %>% 
  wilcox_test(top1_cal_score ~ purity_catn2.f, paired = F)


######* ROC curve *######
# Assuming res_matrix and classes are already defined
classes <- c("BRCA", "DLBC", "LUAD", "MCF.GYN", "MCF.UGI-AD", "MCF.LGI-AD","MCF.SCC", "PAAD")
# Create empty vectors to store AUC values and colors for the legend
auc_values <- numeric(length(classes))
auc_ci_lower <- numeric(length(classes))
auc_ci_upper <- numeric(length(classes))
line_colors <- c("#6e1e33", "#e4703a", "#b890e0", "#0581be", "#63632f", "#F0C507","#27780d", "#734ba9")

pdf("MLclass.roc_BFs_8tumors.pdf",
    width = 5.5, height = 5.5, useDingbats = FALSE)
# Loop through each class to create and plot ROC curves
for (i in seq_along(classes)) {
  class_index <- classes[i]
  
  respon <- as.numeric(ifelse(res_matrix$TCGA_Project2 == class_index,1,0))
  roc_object <- pROC::roc(response = respon, predictor = res_matrix[, class_index], ci = TRUE)
  
  # Store the AUC value and the confidence interval bounds
  auc_values[i] <- pROC::auc(roc_object)
  ci_values <- pROC::ci(roc_object)
  auc_ci_lower[i] <- ci_values[1]
  auc_ci_upper[i] <- ci_values[3]
  
  # Plot the first curve or add subsequent curves
  if (i == 1) {
    # Plot the first curve
    plot(roc_object, 
         main = "ROC Curves for One-vs-All Classification",
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

######* threshould of score *######
res_matrix_pred <- res_matrix %>%
  select(ACC:UVM) %>%
  mutate(predicted_class = names(.)[apply(., 1, which.max)]) %>%
  mutate(max_score = apply(res_matrix %>% select(ACC:UVM), 1, max)) %>%
  bind_cols(res_matrix %>% select(TCGA_Project2), .) # Bring back the gold standard

# Create a sequence of score thresholds to evaluate
cuts <- seq(0, 0.6, by = 0.1)

# Calculate accuracy and coverage for each threshold
df_plot_data <- data.frame(
  threshold = cuts,accuracy = NA_real_,coverage = NA_real_)

for (i in 1:length(cuts)) {
  use_cuts <- cuts[i]
  # Find samples that meet the prediction score threshold
  df_filtered <- res_matrix_pred %>%
    filter(max_score >= use_cuts)
  # Calculate coverage (proportion of samples above the threshold)
  coverage <- nrow(df_filtered) / nrow(res_matrix_pred)
  # Calculate top-1 accuracy for the covered samples
  if (nrow(df_filtered) > 0) {
    accuracy <- sum(df_filtered$predicted_class == df_filtered$TCGA_Project2) / nrow(df_filtered)
  } else {
    accuracy <- NA_real_ # No samples, so accuracy is not applicable
  }
  # Store the results
  df_plot_data[i, "accuracy"] <- accuracy
  df_plot_data[i, "coverage"] <- coverage
}

# Reshape the data to a long format for easier plotting with ggplot2
df_plot_data_long <- df_plot_data %>%
  pivot_longer(cols = c(accuracy, coverage),
               names_to = "metric",values_to = "value")

# Create the plot
p_cut <- ggplot(df_plot_data_long, aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = 0.5, linetype = 2,color = "gray", linewidth = 0.5) +
  scale_color_manual(values = c("accuracy" = "blue", "coverage" = "#00BFC4"),
                     labels = c("Accuracy", "Samples above threshould")) +
  labs(x = "Prediction score threshold", y = NULL, color = "Metric") +
  theme(legend.position = "bottom")
p_cut

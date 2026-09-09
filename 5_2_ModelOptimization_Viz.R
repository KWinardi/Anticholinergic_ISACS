# Optimisation 2
## Assessing sensitivity, specificity, PPV, and NPV

library(caret)
library(ggplot2)
library(ggupset)
library(tidyverse)
library(pROC)
library(clipr)
library(cowplot)

mytheme <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.key = element_blank(),
                 legend.box.background = element_blank(),
                 legend.title=element_text(size=11.5), 
                 legend.text=element_text(size=10.5),
                 plot.title = element_text(face="bold",  size=14, hjust = 0.5),
                 plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
                 axis.title = element_text(size = 12, color = "black"),
                 axis.text = element_text(size = 11, color = "black"),
                 strip.text = element_text(size = 11, color = "black"),
                 panel.border = element_rect(fill = NA), 
                 strip.background = element_rect(fill = NA), 
                 axis.ticks.length=unit(.15, "cm"))

# Defining models
folder_path <- "Results/Optimization/"
model.desc <- read.csv(paste0(folder_path, "/Model_desc.csv"))

# Load model and calculate diagnostic accuracies
file_list <- list.files(path = folder_path, pattern = "_res\\.csv$", full.names = T)
results_list1 <- list()
results_list2 <- list()
for (file in file_list) {
  file_name <- basename(file)
  parts <- str_split(file_name, "_", simplify = TRUE)
  paramA <- parts[1]
  paramB <- str_remove(parts[2], "_res\\.csv")
  df <- read.csv(file)
  tests <- unique(df$distance.metric)
  
  results_list1[[length(results_list1) + 1]] <- df %>%
    mutate(Cholinergic.module = paramB,
           Weights = paramA)
  
  for (t in tests) {
    DF_DIAG <- df %>% 
      dplyr::filter(distance.metric %in% t) %>%
      dplyr::filter(Ground_truth %in% c("True positive", "True negative")) %>%
      mutate(Ground_truth = case_when(Ground_truth == "True positive"~"True anticholinergic",
                                      Ground_truth == "True negative"~"Not an anticholinergic"),
             Diag_accuracy = case_when(Ground_truth == "True anticholinergic" & zscore_sig == "Significant"~"TP",
                                       Ground_truth == "True anticholinergic" & zscore_sig == "NS"~"FN",
                                       Ground_truth == "Not a anticholinergic" & zscore_sig == "Significant"~"FP",
                                       Ground_truth == "Not a anticholinergic" & zscore_sig == "NS"~"TN")
      ) %>% group_by(Diag_accuracy) %>% count
    
    TP <- DF_DIAG %>% dplyr::filter(Diag_accuracy == "TP") %>% pull(n) %>%
      {\(x) if (length(x) == 0) 0 else x}()
    TN <- DF_DIAG %>% dplyr::filter(Diag_accuracy == "TN") %>% pull(n) %>%
      {\(x) if (length(x) == 0) 0 else x}()
    FP <- DF_DIAG %>% dplyr::filter(Diag_accuracy == "FP") %>% pull(n) %>%
      {\(x) if (length(x) == 0) 0 else x}()
    FN <- DF_DIAG %>% dplyr::filter(Diag_accuracy == "FN") %>% pull(n) %>%
      {\(x) if (length(x) == 0) 0 else x}()
    
    # Sensitivity (Recall)
    sensitivity <- TP / (TP + FN)
    
    # Specificity
    specificity <- TN / (TN + FP)
    
    # Positive Predictive Value (PPV/Precision)
    ppv <- TP / (TP + FP)
    
    # Negative Predictive Value (NPV)
    npv <- TN / (TN + FN)
    
    # Accuracy
    accuracy <- (TP + TN)/(TP + TN + FP + FN)
    
    # Store as list
    results_list2[[length(results_list2) + 1]] <- data.frame(Cholinergic.module = paramB,
                                                             Weights = paramA,
                                                             Test = t,
                                                             Sensitivity = round(sensitivity, 3),
                                                             Specificity = round(specificity, 3),
                                                             PPV = round(ppv, 3),
                                                             NPV = round(npv, 3),
                                                             Accuracy = round(accuracy, 3))
    }
}

colours = c( "#A54657",  "#582630", "#F7EE7F", 
             "#D94F70",  "#A8C686",
             "#4DAA57","#F1A66A","#F26157", 
             "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8", "#4C3B4D")

# ROC/AUC curve
roc_points <- right_join(model.desc, bind_rows(results_list1), 
                         by = c("Cholinergic.module", "Weights")) %>%
  dplyr::filter(Ground_truth %in% c("True positive", "True negative")) %>%
  mutate(Ground_truth = case_when(Ground_truth == "True positive"~1,
                                  Ground_truth == "True negative"~0),
         zscore_sig = case_when(zscore_sig == "Significant"~1,
                                zscore_sig == "NS"~0)) %>%
  group_by(Model, Cholinergic.module, Weights, distance.metric) %>%
  group_modify(~{
    roc_obj <- pROC::roc(.x$Ground_truth, .x$zscore, quiet = T)
    coords_df <- pROC::coords(roc_obj, "all", ret = "all", 
                              transpose = FALSE) %>%
      group_by(fpr) %>%
      slice_max(order_by = tpr, n = 1) %>%
      bind_rows(tibble(fpr = 0, tpr = 0, threshold = NA)) %>%
      arrange(fpr, tpr) %>%
      ungroup()
    auc_value <- pROC::auc(roc_obj)[1]
    coords_df$AUC <- auc_value 
    as_tibble(coords_df)
  }) %>% mutate(Model = as.factor(Model))

plot.A <- roc_points %>%
  dplyr::filter(distance.metric %in% c("Closest_b")) %>% 
  ggplot(aes(x = fpr, y = tpr, color = Model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "black") +  
  labs(title = "ROC curves",
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)",
       color = "Model") +
  # facet_wrap(~distance.metric, nrow = 3) +
  scale_color_manual(values = c("#A54657",  "#582630", 
                                "#D94F70",  "#A8C686",
                                "#4DAA57","#F1A66A",
                                "#F26157", "#4C3B4D",
                                "#679289", "#33658A",
                                "#F6AE2D","#86BBD8")) +
  mytheme +
  theme(panel.grid.major = element_line(color = "white", linewidth = 1),
        panel.spacing.x = unit(1.1, "lines"))

# Summary result
results_df <- right_join(model.desc, bind_rows(results_list2), 
                         by = c("Cholinergic.module", "Weights")) %>%
  full_join(., roc_points %>% dplyr::select(Cholinergic.module, Weights, 
                                            distance.metric, AUC, Model) %>% 
              distinct() %>% mutate(Model = as.numeric(Model)),
            by = c("Cholinergic.module", "Weights", "Model", 
                   c("Test" = "distance.metric"))) %>%
  mutate(Cholinergic.module = case_when(Cholinergic.module == "mAChR + enzyme + transporter + delirium"~"CAP + Downstream",
                                        
                                        T~Cholinergic.module)) %>%
  mutate(Cholinergic.module = factor(Cholinergic.module,
                                     levels = c("CAP",
                                                "CAP + Downstream",
                                                "KEGG-defined")),
         Weights = factor(Weights, 
                          levels = c("None", "Linear", "Linear piecewise", 
                                     "Hyperbola"))
  ) %>%
  pivot_longer(cols = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "AUC"),
               names_to = "Diagnostic",
               values_to = "Value")

custom_order <- list(c("None", "CAP"),
                     c("None", "CAP", "Downstream"),
                     c("None", "KEGG-defined"),
                     c("Linear", "CAP"),
                     c("Linear", "CAP", "Downstream"),
                     c("Linear", "KEGG-defined"),
                     c("Linear piecewise", "CAP"),
                     c("Linear piecewise", "CAP", "Downstream"),
                     c("Linear piecewise", "KEGG-defined"),
                     c("Hyperbola", "CAP"),
                     c("Hyperbola", "CAP", "Downstream"),
                     c("Hyperbola", "KEGG-defined")
)

# For plotting order
burner <- results_df %>% 
  rowwise() %>%
  mutate(dup_times = 13 - Model) %>%
  ungroup() %>%
  uncount(weights = dup_times) %>%
  mutate(x = 0)

## Multiple distance
plot.B <- results_df %>% 
  mutate(x = 1) %>%
  bind_rows(., burner) %>%
  mutate(Label = pmap(list(Weights, Cholinergic.module), function(Weights, Cholinergic.module){
    c(paste0(Weights), 
      if(Cholinergic.module == "CAP + Downstream") 
        c("CAP", "Downstream") else paste0(Cholinergic.module))
  })) %>%
  mutate(Cholinergic.module = case_when(x == 1~Cholinergic.module,
                                        T~NA),
         Diagnostic = factor(Diagnostic, levels = c("Sensitivity",
                                                    "Specificity", 
                                                    "PPV", "NPV",
                                                    "Accuracy", "AUC"))) %>% 
  dplyr::filter(Test %in% c("Closest_b")) %>% 
  ggplot(aes(x=Label, y=Value)) +
  geom_rect(xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_rect(xmin = 9.5, xmax = 12.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_line(#position = position_dodge(width = 0.5),
    aes(group = Test)) +
  geom_point(size = 3.5, pch = 21,  color = "#F26157", fill = "#F26157"
             #position = position_dodge(width = 0.5),
             # aes(fill = Test, color = Test)
             ) +
  scale_x_upset(sets = c("CAP", "Downstream", "KEGG-defined", 
                         "None", "Linear", "Linear piecewise", "Hyperbola"),
                position="top", name = "Model Optimization",
                order_by = "freq") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_color_manual(values = c(rep("black", 6)), 
                     na.value = "transparent") +
  # scale_fill_manual(values = c("#F7EE7F", "#A54657", "#F26157",
  #                              "#33658A","#86BBD8", "#4DAA57"), 
  #                   na.value = "transparent", na.translate = F) +
  facet_grid(Diagnostic~., scales = "free_y") + 
  labs(y = "") + # , fill = "Distance metric"
  guides(color = "none") +
  theme_combmatrix(combmatrix.label.text = element_text(size=12),
                   combmatrix.label.extra_spacing = 5) +
  mytheme

plot.B2 <- results_df %>% 
  mutate(x = 1) %>%
  bind_rows(., burner) %>%
  mutate(Label = pmap(list(Weights, Cholinergic.module), function(Weights, Cholinergic.module){
    c(paste0(Weights), 
      if(Cholinergic.module == "CAP + Downstream") 
        c("CAP", "Downstream") else paste0(Cholinergic.module))
  })) %>%
  mutate(Cholinergic.module = case_when(x == 1~Cholinergic.module,
                                        T~NA),
         Diagnostic = factor(Diagnostic, levels = c("Sensitivity",
                                                    "Specificity", 
                                                    "PPV", "NPV",
                                                    "Accuracy", "AUC"))) %>% 
  dplyr::filter(Test %in% c("Closest_b", "Closest_ab")) %>%
  mutate(line.grp = paste0(Test, Weights)) %>%
  ggplot(aes(x=Label, y=Value)) +
  geom_line(#position = position_dodge(width = 0.5),
            aes(group = line.grp)) +
  geom_point(size = 3.5, pch = 21,  
             #position = position_dodge(width = 0.5),
             aes(fill = Test, color = Test)) +
  scale_x_upset(sets = c("CAP", "Downstream", "KEGG-defined", 
                         "None", "Linear", "Linear piecewise", "Hyperbola"),
                position="top", name = "Model Optimization",
                order_by = "freq") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_color_manual(values = c(rep("black", 6)), 
                     na.value = "transparent") +
  scale_fill_manual(values = c("#A54657", "#F26157"), 
                    na.value = "transparent", na.translate = F) +
  facet_grid(Diagnostic~., scales = "free_y") + 
  labs(y = "", fill = "Distance metric") +
  guides(color = "none") +
  theme_combmatrix(combmatrix.label.text = element_text(size=12),
                   combmatrix.label.extra_spacing = 5) +
  mytheme


## Just one distance
plot.C <- results_df %>% 
  mutate(x = 1) %>%
  bind_rows(., burner) %>%
  dplyr::filter(Test %in% "Closest_ab") %>%
  mutate(Label = pmap(list(Weights, Cholinergic.module), function(Weights, Cholinergic.module){
    c(paste0(Weights), 
      if(Cholinergic.module == "CAP + Downstream") 
        c("CAP", "Downstream") else paste0(Cholinergic.module))
  })) %>%
  mutate(Cholinergic.module = case_when(x == 1~Cholinergic.module,
                                        T~NA),
         Diagnostic = factor(Diagnostic, levels = c("Sensitivity",
                                                    "Specificity", 
                                                    "PPV", "NPV",
                                                    "Accuracy", "AUC"))) %>% 
  ggplot(aes(x=Label, y=Value)) +
  geom_line(aes(group = Weights)) +
  geom_point(size = 3.5, pch = 21,  
             aes(fill = Cholinergic.module, color = Cholinergic.module)) +
  scale_x_upset(sets = c("CAP", "Downstream", "KEGG-defined", 
                         "None", "Linear", "Linear piecewise", "Hyperbola"),
                position="top", name = "",
                order_by = "freq") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_color_manual(values = c(rep("black", 3)), 
                     na.value = "transparent") +
  scale_fill_manual(values = c("#A54657",  "#582630", "#F7EE7F"), 
                    na.value = "transparent", na.translate = F) +
  facet_grid(Diagnostic~., scales = "free_y") + 
  labs(y = "", fill = "Cholinergic module",
       title = "Closest (ab)") +
  guides(color = "none") + 
  theme_combmatrix(combmatrix.label.text = element_text(size=12),
                   combmatrix.label.extra_spacing = 5) +
  mytheme

png("Results/Accuracy_1.png", width = 14, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.A, plot.B, nrow = 1, rel_widths = c(1,1.25),
                   axis = "b", align = "h")
dev.off()

png("Results/Accuracy_2.png", width = 8.75, height = 8.75, units = "in", res = 300)
plot.B2
dev.off()

write_clip(results_df %>% pivot_wider(names_from = "Diagnostic", values_from = "Value"))


### For Plotting Only ###
results_df <- read.csv(paste0(folder_path, "/PerformanceSummary.csv")) %>%
  mutate(Cholinergic.module = factor(Cholinergic.module,
                                     levels = c("CAP",
                                                "CAP + downstream",
                                                "CAP + effector",
                                                "CAP + downstream + effector"
                                                )),
         Weights = factor(Weights, 
                          levels = c("None", "Linear", "Linear piecewise", 
                                     "Hyperbola"))
  ) %>%
  pivot_longer(cols = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "AUC"),
               names_to = "Diagnostic",
               values_to = "Value")

custom_order <- list(c("None", "CAP"),
                     c("None", "CAP", "Downstream"),
                     c("None", "CAP", "Effector"),
                     c("None", "CAP", "Downstream", "Effector"),
                     
                     c("Linear", "CAP"),
                     c("Linear", "CAP", "Downstream"),
                     c("Linear", "CAP", "Effector"),
                     c("Linear", "CAP", "Downstream", "Effector"),
                     
                     c("Linear piecewise", "CAP"),
                     c("Linear piecewise", "CAP", "Downstream"),
                     c("Linear piecewise", "CAP", "Effector"),
                     c("Linear piecewise", "CAP", "Downstream", "Effector"),
                     
                     c("Hyperbola", "CAP"),
                     c("Hyperbola", "CAP", "Downstream"),
                     c("Hyperbola", "CAP", "Effector"),
                     c("Hyperbola", "CAP", "Downstream", "Effector")
)

# For plotting order
burner <- results_df %>% 
  rowwise() %>%
  mutate(dup_times = 17 - Model) %>%
  ungroup() %>%
  uncount(weights = dup_times) %>%
  mutate(x = 0)

result.p <- results_df %>% 
  mutate(x = 1) %>%
  bind_rows(., burner) %>%
  mutate(Label = pmap(list(Weights, Cholinergic.module), function(Weights, Cholinergic.module){
    c(paste0(Weights), 
      str_split(Cholinergic.module, " \\+ ")[[1]] %>% str_to_title() %>% str_replace("Cap", "CAP"))
  })) %>%
  mutate(Cholinergic.module = case_when(x == 1~Cholinergic.module,
                                        T~NA),
         Diagnostic = factor(Diagnostic, levels = c("Sensitivity",
                                                    "Specificity", 
                                                    "PPV", "NPV",
                                                    "Accuracy", "AUC"))) %>% 
  dplyr::filter(Distance %in% c("Closest")) %>%
  ggplot(aes(x=Label, y=Value)) +
  geom_rect(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_rect(xmin = 12.5, xmax = 16.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_line(#position = position_dodge(width = 0.5),
    aes(group = Distance)) +
  geom_point(size = 3.5, pch = 21,  color = "#F26157", fill = "#F26157"
             #position = position_dodge(width = 0.5),
             # aes(fill = Test, color = Test)
  ) +
  scale_x_upset(sets = c("CAP", "Downstream", "Effector", 
                         "None", "Linear", "Linear piecewise", "Hyperbola"),
                position="top", name = "Model Optimization",
                order_by = "freq") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_color_manual(values = c(rep("black", 6)), 
                     na.value = "transparent") +
  # scale_fill_manual(values = c("#F7EE7F", "#A54657", "#F26157",
  #                              "#33658A","#86BBD8", "#4DAA57"), 
  #                   na.value = "transparent", na.translate = F) +
  facet_grid(Diagnostic~., scales = "free_y") + 
  labs(y = "") + # , fill = "Distance metric"
  guides(color = "none") +
  theme_combmatrix(combmatrix.label.text = element_text(size=12),
                   combmatrix.label.extra_spacing = 5) +
  mytheme
# dev.off()


# Linear: wij = 1 + (20 - bij)/20
df1 <- data.frame(bij = seq(0, 40, length.out = 200))
df1$wij <- 1 + (20 - df1$bij) / 20
p1 <- ggplot(df1, aes(bij, wij)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = expression(atop("Linear", w[ij] == 1 + frac(20 - b[ij], 20))),
       x = expression(b[ij]), y = expression(w[ij])) +
  mytheme

# Linear piecewise: wij = 1 + (5-bij)/5 for bij<5, else 1
df2 <- data.frame(bij = seq(0, 20, length.out = 200))
df2$wij <- ifelse(df2$bij < 5, 1 + (5 - df2$bij) / 5, 1)
p2 <- ggplot(df2, aes(bij, wij)) +
  geom_line(color = "darkorange", linewidth = 1) +
  labs(title = expression(atop("Linear piecewise",
                               w[ij] == group("{", atop(1 + frac(5 - b[ij], 5) ~ (b[ij] < 5), 1 ~ (b[ij] >= 5)), ""))),
       x = expression(b[ij]), y = expression(w[ij])) +
  mytheme

# Hyperbola: wij = 1 + 1/bij
df3 <- data.frame(bij = seq(0.1, 20, length.out = 200))
df3$wij <- 1 + 1 / df3$bij
p3 <- ggplot(df3, aes(bij, wij)) +
  geom_line(color = "forestgreen", linewidth = 1) +
  labs(title = expression(atop("Hyperbola", w[ij] == 1 + frac(1, b[ij]))),
       x = expression(b[ij]), y = expression(w[ij])) +
  mytheme

weight.p <- plot_grid(p1, p2, p3, ncol = 1)

BR.file <- ggdraw() + 
  draw_image("Biorender/SignalingCascade.png") + 
  theme(plot.margin = margin(t=0.15, l=0.1, r=0.1, b=0.15, unit = "cm"))

tiff(paste0("Results/SFigure1.tiff"), units="in", height=7.5, width=14, res=300)
plot_grid(BR.file, weight.p, result.p, nrow = 1, 
          rel_widths = c(0.85, 0.5, 0.9),
          labels = c("(a)", "(b)", "(c)"))
dev.off()


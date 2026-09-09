# Decoding Anticholinergic Effects
## Only analyse those drugs with complete network and pKi results
source("6_1_RedefiningAnticholinergics_Setup.R")

# Merge all
result.w.anno.full <- distance.df.long %>%
  # Select proximity measure of interest
  dplyr::filter(distance.metric %in% "Closest_b") %>%
  dplyr::select(-distance.metric) %>% drop_na(zscore) %>%
  # Anticholinergic scale
  left_join(., anticholinergics.df, c("Drug" = "Generic.Name")) %>%
  # ORCA
  left_join(., ORCA.df, "Drug") %>%
  mutate(ORCA = factor(ORCA,
                       levels = c("High", "Moderate", "Low", "No"))) %>%
  # ACSBC
  left_join(., ACSBC.df, "Drug") %>%
  mutate(ACSBC.Anticholinergic.Activity = factor(ACSBC.Anticholinergic.Activity,
                                                 levels = c("High", "Moderate", "Low", "No")),
         ACSBC.Sedative.Activity = factor(ACSBC.Sedative.Activity,
                                          levels = c("High", "Moderate", "Low", "No")),
         Concordance.ACSBC.ORCA = factor(case_when(ORCA == ACSBC.Anticholinergic.Activity~ORCA,
                                                   T~NA),
                                         levels = c("High", "Moderate", "Low", "No"))) %>%
  # SAA
  left_join(., SAA.df.average, 
            c("Drug" = "Medication")) %>% 
  # ATS
  left_join(., ATS.df, c("Drug"="Medication")) %>% 
  # ABS 
  left_join(., ABS.df, "Drug") %>%
  # ATS
  left_join(., AEC.df, "Drug") %>% 
  # Ground truth
  mutate(Ground_truth = case_when(Drug %in% true_pos~"True positive",
                                  Drug %in% true_neg~"True negative",
                                  Drug %in% false_pos~"False positive",
                                  T~NA),
         Drug_category = drug_category_map[Drug],
         Drug_category = ifelse(is.na(Drug_category), "Other", Drug_category),
         Drug_category = factor(Drug_category,
                                levels = c("Prototypical anticholinergics",
                                           "Urinary incontinence agent",
                                           "TCA", "Sedating antihistamine",
                                           "Bisphosphonate", "Statin", 
                                           "ACE inhibitor", "DOAC",
                                           "Acetylcholinesterase inhibitor",
                                           "Organophosphate", 
                                           "Cholinergic agonist")),
         Ground_truth = factor(Ground_truth,
                               levels = c("True positive", "True negative",
                                          "False positive"))) %>%
  # Binding affinity
  left_join(., BA.mAChR.formatted, c("Drug" = "DRUGNAME")) %>%
  # Predicted ADMET
  left_join(., drug.ADMET, c("Drug"="Name")) %>%
  # QSAR
  left_join(., PSICHIC.df, "Drug") %>%
  rowwise() %>%
  mutate(QSAR_score = mean(c_across(all_of(Outcomes[1:5])), na.rm = T)) 

result.w.anno.full %>% drop_na(Ground_truth, Drug_category) %>% 
  group_by(Ground_truth, Drug_category) %>% 
  dplyr::summarize(CHRM1_QSAR = mean(CHRM1_QSAR, na.rm = T),
                   CHRM2_QSAR = mean(CHRM2_QSAR, na.rm = T),
                   CHRM3_QSAR = mean(CHRM3_QSAR, na.rm = T),
                   CHRM4_QSAR = mean(CHRM4_QSAR, na.rm = T),
                   CHRM5_QSAR = mean(CHRM5_QSAR, na.rm = T))

### Clustering ###
result.w.anno.full <- result.w.anno.full %>% drop_na(QSAR_score)
dim(result.w.anno.full)

## Looking at pKi of Ground truth
roc_points <- result.w.anno.full %>%
  dplyr::filter(Ground_truth %in% c("True positive", "True negative", "False positive")) %>%
  mutate(Ground_truth.num = case_when(Ground_truth == "True positive"~1,
                                      Ground_truth == "True negative"~0,
                                      Ground_truth == "False positive"~0),
         zscore_sig = case_when(zscore_sig == "Significant"~1,
                                zscore_sig == "NS"~0)) %>% ungroup 

plot.A <- result.w.anno.full %>% 
  dplyr::select(-ACHE_QSAR, -CHAT_QSAR, -SLC5A7_QSAR, -SLC18A3_QSAR) %>%
  pivot_longer(cols = c("QSAR_score", ends_with("_QSAR"))) %>%
  mutate(name = case_when(name == "QSAR_score"~"Mean averaged",
                          T~gsub("_QSAR", "", name))) %>%
  ggplot(aes(x = value, color = name)) +
  geom_density(linewidth = 1.25) + 
  scale_color_manual(values = c( wes_palette("FantasticFox1"), "grey20")) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(x = "PSICHIC score", y = "Density", color = "Protein") +
  mytheme

roc_obj <- roc_points %>%
  dplyr::select(-ACHE_QSAR, -CHAT_QSAR, -SLC5A7_QSAR, -SLC18A3_QSAR) %>%
  pivot_longer(cols = c("QSAR_score", ends_with("_QSAR"))) %>%
  mutate(name = case_when(name == "QSAR_score"~"Mean averaged",
                          T~gsub("_QSAR", "", name))) %>%
  group_by(name) %>%
  group_modify(~{
    roc_obj <- pROC::roc(.x$Ground_truth.num, .x$value, quiet = T)
    coords_df <- pROC::coords(roc_obj, "all", ret = "all", 
                              transpose = FALSE) %>%
      group_by(fpr) %>%
      slice_max(order_by = tpr, n = 1) %>%
      bind_rows(tibble(fpr = 0, tpr = 0, threshold = NA)) %>%
      arrange(fpr, tpr) %>%
      ungroup()
    coords_df$AUROC <- pROC::auc(roc_obj)[1]
    coords_df$Best_threshold <- pROC::coords(roc_obj,
                                             "best",
                                             best.method = "youden",
                                             ret = "threshold")
    as_tibble(coords_df)
  }) %>% mutate(Model = paste0(name, " (AUC-ROC = ", round(AUROC, 3), ")"))

KW.test <- result.w.anno.full %>% 
  drop_na(Ground_truth) %>% 
  kruskal_test(QSAR_score ~ Ground_truth)    

dunn.test <- result.w.anno.full %>% 
  drop_na(Ground_truth) %>% 
  dunn_test(QSAR_score ~ Ground_truth, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Ground_truth") %>%
  mutate(y.position = y.position + 0.25) %>% dplyr::select(-groups)

plot.D <- roc_points %>% 
  ggplot(aes(x = Ground_truth, y = QSAR_score)) +
  geom_violin() +
  geom_boxplot(outliers = F, alpha = 0.3, width = 0.15, fill = "grey70") +
  geom_point(size = 3.5, alpha = 0.65, color = "white", pch =21,
             position = position_jitterdodge(jitter.width = 0.4,
                                             dodge.width = 0.5),
             aes(fill = Drug_category)) + 
  geom_hline(yintercept = roc_obj$Best_threshold[1,1], linetype = "dashed") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.08) +
  labs(x = "", y = "PSICHIC score", 
       subtitle = get_test_label(KW.test, detailed = TRUE),
       title = "Ground Truth", fill = "Drug class") +
  scale_fill_manual(values = c("Prototypical anticholinergics" = "#E1C16E",
                               "Urinary incontinence agent" = "#F4A261",
                               "TCA" = "#E76F51", 
                               "Sedating antihistamine" = "#D45D5D",
                               "Bisphosphonate" = "#617a6c", 
                               "Statin" = "#3bff91", 
                               "ACE inhibitor" = "#5c8745", 
                               "DOAC" = "#075c2d",
                               "Acetylcholinesterase inhibitor" = "grey10",
                               "Organophosphate" = "grey30", 
                               "Cholinergic agonist" = "grey50")) + 
  guides(fill = guide_legend(override.aes = list(color = "black", pch = c(21))),
         shape = guide_legend(override.aes = list(color = "black", pch = c(19,17)))) +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))

plot.E <- ggplot(roc_obj, aes(x = fpr, y = tpr, color = Model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "black") +
  scale_color_manual(values = c( wes_palette("FantasticFox1"), "grey20")) +
  labs(title = "ROC of Ground Truth Drugs",
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)",
       color = "") +
  mytheme +
  theme(panel.grid.major = element_line(color = "white", linewidth = 1),
        panel.spacing.x = unit(1.1, "lines"),
        legend.position = c(0.975, 0.025),
        legend.background = element_blank(),
        legend.justification = c(1, 0))

drugList <- result.w.anno.full$QSAR_score
names(drugList) <- result.w.anno.full$Drug
drugList <- sort(drugList, decreasing = T)
DSEA.result <- GSEA(drugList, minGSSize = 5,
                    TERM2GENE = ATC.df.enrich)

plot.C <- DSEA.result@result %>%
  left_join(., ATC.df %>% mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")"),
                                 name_1st = str_to_sentence(name_1st)) %>%
              distinct(name_1st, ATC_3rd),
            c("ID" = "ATC_3rd")) %>%
  dplyr::filter(p.adjust < 0.01) %>%
  mutate(log.p = -log10(p.adjust)) %>%
  ggdotchart(., x = "Description", y = "NES",
             add = "segments", 
             add.params = list(color = "name_1st", size = 1.25),
             color = "name_1st", 
             sorting = "descending", rotate = T,
             dot.size = "log.p",
             ggtheme = theme_pubr()) +
  scale_size(range = c(3,6.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  guides(color = "none") +
  ATC.col +
  labs(x = "", y = "NES",
       title = "ATC Level 3\nSet Enrichment Analysis",
       size = "-log10(FDR-adjusted p-value)") +
  theme(legend.position = "bottom",
        plot.title = element_text(face="bold",size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",size=13, hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        strip.text = element_text(size = 11, color = "black"),
        axis.ticks.length=unit(.15, "cm")) 

tiff(paste0("Results/SFigure5.tiff"), units="in", height=10, width=18.5, res=300)
plot_grid(plot_grid(plot.A, plot.D, plot.E, ncol = 1, #plot.B, 
                    rel_heights = c(0.7,1,1), align = "v", axis = "lr",
                    labels = c("(a)", "(b)", "(c)")),
          plot.C, nrow = 1, labels = c("", "(d)"), rel_widths = c(0.7,1))
dev.off()

# The other CAPs
roc_obj <- roc_points %>%
  pivot_longer(cols = c(ACHE_QSAR, CHAT_QSAR, SLC5A7_QSAR, SLC18A3_QSAR)) %>%
  mutate(name = gsub("_QSAR", "", name)) %>%
  group_by(name) %>%
  group_modify(~{
    roc_obj <- pROC::roc(.x$Ground_truth.num, .x$value, quiet = T)
    coords_df <- pROC::coords(roc_obj, "all", ret = "all",
                              transpose = FALSE) %>%
      group_by(fpr) %>%
      slice_max(order_by = tpr, n = 1) %>%
      bind_rows(tibble(fpr = 0, tpr = 0, threshold = NA)) %>%
      arrange(fpr, tpr) %>%
      ungroup()
    coords_df$AUROC <- pROC::auc(roc_obj)[1]
    coords_df$Best_threshold <- pROC::coords(roc_obj,
                                             "best",
                                             best.method = "youden",
                                             ret = "threshold")
    as_tibble(coords_df)
  }) %>% mutate(Model = paste0(name, " (AUC-ROC = ", round(AUROC, 3), ")"))

tiff(paste0("Results/SFigure6.tiff"), units="in", height=5.5, width=12, res=300)
KW.test <- roc_points %>%
  pivot_longer(cols = c(ACHE_QSAR, CHAT_QSAR, SLC5A7_QSAR, SLC18A3_QSAR)) %>%
  mutate(name = gsub("_QSAR", "", name)) %>% 
  group_by(name) %>%
  kruskal_test(value ~ Ground_truth)    

dunn.test <- roc_points %>%
  pivot_longer(cols = c(ACHE_QSAR, CHAT_QSAR, SLC5A7_QSAR, SLC18A3_QSAR)) %>%
  mutate(name = gsub("_QSAR", "", name)) %>% 
  group_by(name) %>% 
  dunn_test(value ~ Ground_truth, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Ground_truth") %>%
  mutate(y.position = y.position + 0.25) %>% dplyr::select(-groups)

plot_grid(roc_points %>%
            pivot_longer(cols = c(ACHE_QSAR, CHAT_QSAR, SLC5A7_QSAR, SLC18A3_QSAR)) %>%
            mutate(name = gsub("_QSAR", "", name)) %>%
            ggplot(aes(x = Ground_truth, y = value)) +
            geom_violin() +
            geom_boxplot(outliers = F, alpha = 0.3, width = 0.15, fill = "grey70") +
            geom_point(size = 2.5, alpha = 0.65, color = "white", pch =21,
                       position = position_jitterdodge(jitter.width = 0.4,
                                                       dodge.width = 0.5),
                       aes(fill = Drug_category)) +
            scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
            facet_wrap(.~name, ncol = 2, scales = "free_y") +
            stat_pvalue_manual(dunn.test,
                               label = "p.adj.signif") +
            labs(x = "", y = "PSICHIC score", 
                 title = "Ground Truth", fill = "Drug class") +
            scale_fill_manual(values = c("Prototypical anticholinergics" = "#E1C16E",
                                         "Urinary incontinence agent" = "#F4A261",
                                         "TCA" = "#E76F51",
                                         "Sedating antihistamine" = "#D45D5D",
                                         "Bisphosphonate" = "#617a6c",
                                         "Statin" = "#3bff91",
                                         "ACE inhibitor" = "#5c8745",
                                         "DOAC" = "#075c2d",
                                         "Acetylcholinesterase inhibitor" = "grey10",
                                         "Organophosphate" = "grey30",
                                         "Cholinergic agonist" = "grey50")) +
            guides(fill = guide_legend(override.aes = list(color = "black", pch = c(21))),
                   shape = guide_legend(override.aes = list(color = "black", pch = c(19,17)))) +
            mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1)),
          ggplot(roc_obj, aes(x = fpr, y = tpr, color = Model)) +
            geom_line(size = 1) +
            geom_abline(linetype = "dashed", color = "black") +
            scale_color_manual(values = c(wes_palette("GrandBudapest2")[1],
                                          wes_palette("GrandBudapest1")[2],
                                          wes_palette("GrandBudapest2")[3:4]
            )) +
            labs(title = "ROC of Ground Truth Drugs",
                 x = "False Positive Rate (1 - Specificity)",
                 y = "True Positive Rate (Sensitivity)",
                 color = "") +
            mytheme +
            theme(panel.grid.major = element_line(color = "white", linewidth = 1),
                  panel.spacing.x = unit(1.1, "lines"),
                  legend.position = c(0.975, 0.025),
                  legend.background = element_blank(),
                  legend.justification = c(1, 0)),
          nrow = 1, rel_widths = c(1, 0.6),
          align = "h", axis = "tb")
dev.off()

library(ClusterR)
clust.vars <- c("distance.value", "zscore", 
                "CHRM1_QSAR", "CHRM2_QSAR","CHRM3_QSAR", "CHRM4_QSAR","CHRM5_QSAR"
                ,"CHAT_QSAR", "ACHE_QSAR", "SLC5A7_QSAR", "SLC18A3_QSAR"
)
cluster.df <- result.w.anno.full %>%
  column_to_rownames("Drug") %>%
  dplyr::select(all_of(clust.vars)) 
colnames(cluster.df) <- c("Distance", "Z-score",
                          "CHRM1 pKi", "CHRM2 pKi", "CHRM3 pKi", "CHRM4 pKi", "CHRM5 pKi"
                          ,"CHAT pKi", "ACHE pKi", "SLC5A7 pKi", "SLC18A3 pKi"
)
boxplot(cluster.df)
## Hemicholinium-3 is an indirect anticholinergic as it inhibits SLC5A7(CHT1)
## Vesamicol is also an indirect anticholinergic by inhibiting SLC18A3 (VAChT)

cluster.df <- scale(cluster.df)
boxplot(cluster.df)

####################
# Including receptor expression levels

## CHRM expression across organ (HPA)
CHRM.exp <- read_tsv("6_HPA/CHRM_EXP.tsv") %>% 
  dplyr::filter(Gene != "VLDLR") %>%
  rbind(., read_tsv("6_HPA/gene_name_ACHE.tsv")) %>%
  rbind(., read_tsv("6_HPA/gene_name_CHAT.tsv")) %>%
  rbind(., read_tsv("6_HPA/gene_name_SLC5A7.tsv")) %>%
  rbind(., read_tsv("6_HPA/gene_name_SLC18A3.tsv")) %>%
  dplyr::select(-"Gene description") %>%
  column_to_rownames("Gene")

organLabels <- read.csv("6_HPA/organLabels.csv") %>%
  mutate(labels = str_to_sentence(gsub("_", " ", Organ)),
         System = factor(System))
map <- setNames(organLabels$labels, organLabels$Data_col_names)
CHRM.exp.tissue <- CHRM.exp[,organLabels[organLabels$Dataset == "Tissue expression", "Data_col_names"]]
colnames(CHRM.exp.tissue) <- map[colnames(CHRM.exp.tissue)]

CHRM.exp.brain <- CHRM.exp[,organLabels[organLabels$Dataset == "Brain expression", "Data_col_names"]]
colnames(CHRM.exp.brain) <- map[colnames(CHRM.exp.brain)]

col_fun = colorRamp2(c(0, 0.5, 5, 10, 25, 50, 70), 
                     c("white","#EDD8DA", "#d6a9ab", "#b57d7f", "#824446", "#6e2b2d", "#4a1712")) 

ha <- HeatmapAnnotation(df = organLabels %>% dplyr::filter(Dataset == "Tissue expression") %>%
                          column_to_rownames("labels") %>% dplyr::select(System),
                        show_annotation_name = F,
                        col = list(System = setNames(pal_rickandmorty("schwifty")(length(levels(organLabels$System))),
                                                     levels(organLabels$System))))

ht1 <- ComplexHeatmap::Heatmap(CHRM.exp.tissue, col = col_fun,
                               border = T, border_gp = gpar(col = "gray20"),
                               show_heatmap_legend = F, 
                               cluster_rows = F, top_annotation = ha,
                               column_title = "Tissue\nexpression")
ht2 <- ComplexHeatmap::Heatmap(CHRM.exp.brain, col = col_fun,
                               border = T, border_gp = gpar(col = "gray20"),
                               name = "nTPM",
                               cluster_rows = F,
                               column_title = "Brain\nexpression")
# tiff("Results/CHRM_EXP.tiff", units="in", height=4, width=14, res=300)
draw(ht1 + ht2, auto_adjust = T)
# dev.off()

## Compare tissue and brain datasets
CHRM.exp %>% rownames_to_column("Gene") %>% pivot_longer(-Gene, names_to = "Data_col_names") %>% 
  left_join(., organLabels, "Data_col_names") %>% dplyr::select(-Data_col_names, -gganatogram_organ) %>%
  pivot_wider(values_from = "value", names_from = "Dataset") %>% drop_na() %>%
  ggplot(aes(x = `Tissue expression`, y = `Brain expression`)) +
  geom_point(aes(color = Gene)) + facet_wrap(~labels, scales = "free")

## Since there are some overlapping tissues between the 2 datasets, use brain one
temp <- intersect(colnames(CHRM.exp.brain), colnames(CHRM.exp.tissue))
CHRM.exp.tissue.final <- CHRM.exp.tissue %>% dplyr::select(-all_of(temp))
CHRM.exp.tissue.final <- cbind(CHRM.exp.tissue.final, CHRM.exp.brain)

col_fun = circlize::colorRamp2(c(0, log1p(c(0.5, 10, 50, 100, 150, 250))),
                               c("white","#EDD8DA", "#d6a9ab", "#b57d7f", "#824446", "#6e2b2d", "#4a1712")) 
col_fun2 = circlize::colorRamp2(c(-1,0,1), 
                                c("#002c58", "white","#580000"))

ha <- HeatmapAnnotation(df = organLabels %>% 
                          distinct(labels, System) %>% column_to_rownames("labels"),
                        show_annotation_name = F,
                        col = list(System = setNames(pal_rickandmorty("schwifty")(length(levels(organLabels$System))),
                                                     levels(organLabels$System))))
lgd_scale <- Legend(col_fun = col_fun, title = "nTPM", 
                    at = c(0, log1p(c(0.5, 10, 50, 100, 150, 250))),
                    labels =c(0, 0.5, 10, 50, 100, 150, 250),
                    direction = "vertical", title_position = "topleft",
                    legend_height = unit(3, "cm"), legend_width = unit(6, "cm"))

## Here at least 1 organ with nTPM > 0.5 
protein_tissue <- CHRM.exp.tissue.final[, colSums(CHRM.exp.tissue.final >= 0.5, na.rm = F) > 0] %>% as.matrix
print(paste0("Originally, there were ", ncol(CHRM.exp.tissue.final), " tissues"))
print(paste0("After filtering 0.5 TPM in at least 1 tissue, there were ", ncol(protein_tissue), " tissues left"))
tissue_cols <- colnames(protein_tissue)

drug_receptor <- result.w.anno.full %>%
  dplyr::select(Drug, ends_with("_QSAR")) %>% 
  column_to_rownames("Drug") 
colnames(drug_receptor) <- gsub("_QSAR", "", colnames(drug_receptor))

all(colnames(drug_receptor) == rownames(protein_tissue))

boxplot(drug_receptor)
boxplot(scale(drug_receptor))
boxplot(protein_tissue)
boxplot(log1p(protein_tissue))

drug_tissue_score <- scale(drug_receptor) %*% log1p(protein_tissue) 
cluster.df <- cbind(as.data.frame(cluster.df), drug_tissue_score) %>% dplyr::select(-ends_with("pKi"))
boxplot(cluster.df)

write_clip(log1p(protein_tissue)) # Supplementary Table 6

ha.anno <- organLabels %>% 
  dplyr::filter(labels %in% colnames(drug_tissue_score)) %>%
  distinct(labels, System) %>% column_to_rownames("labels")

ha <- HeatmapAnnotation(System = ha.anno[colnames(protein_tissue),],
                        show_annotation_name = F,
                        col = list(System = setNames(pal_rickandmorty("schwifty")(length(levels(organLabels$System))),
                                                     levels(organLabels$System))))

drug_tissue_score_plot.df <- drug_tissue_score %>% as.data.frame %>%
  rownames_to_column("Drug") %>% 
  pivot_longer(-Drug) %>% 
  inner_join(., result.w.anno.full %>% distinct(Drug, Ground_truth, Drug_category) %>% drop_na(),
             "Drug") %>%
  mutate(Truth_simple = case_when(Ground_truth == "True positive"~"Group1",
                                  T~"Group2"),
         Tissue = case_when(name %in% colnames(CHRM.exp.brain)~"Brain tissue",
                            T~"Other organs"))

region_order <- drug_tissue_score_plot.df %>%
  group_by(name, Truth_simple) %>%
  summarise(mean_val = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Truth_simple, values_from = mean_val) %>%
  mutate(diff = abs(Group1 - Group2)) %>%
  arrange(diff)

drug_tissue_score_plot.1 <- drug_tissue_score_plot.df %>%
  mutate(name = factor(name, levels = region_order$name)) %>%
  ggplot(aes(x = value, y = fct_rev(name), fill = Drug_category)) +
  geom_density_ridges(alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("Prototypical anticholinergics" = "#E1C16E",
                               "Urinary incontinence agent" = "#F4A261",
                               "TCA" = "#E76F51", 
                               "Sedating antihistamine" = "#D45D5D",
                               "Bisphosphonate" = "#617a6c", 
                               "Statin" = "#3bff91", 
                               "ACE inhibitor" = "#5c8745", 
                               "DOAC" = "#075c2d",
                               "Acetylcholinesterase inhibitor" = "grey10",
                               "Organophosphate" = "grey30", 
                               "Cholinergic agonist" = "grey50")) +
  labs(x = "Tissue-specific drug exposure score", y = "", fill = "") +
  # facet_grid(~Tissue, scales = "free_y") +
  theme_ridges(grid = F, center_axis_labels = T) +
  theme(legend.position = "right",
        plot.title = element_text(face="bold",size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",size=13, hjust = 0.5),
        strip.text = element_text(size = 11, color = "black"),
        panel.border = element_rect(fill = NA), 
        strip.background = element_rect(fill = NA))

BR.tissue.score <- ggdraw() + 
  draw_image("Biorender/TissueSpecificDrugExposure.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

tiff("Results/SFigure7.tiff", units="in", height=17, width=12.5, res=300)
plot_grid(BR.tissue.score,
          plot_grid(NULL, grid.grabExpr(draw(ComplexHeatmap::Heatmap(log1p(protein_tissue), col = col_fun,
                                                                     border = T, border_gp = gpar(col = "gray20"),
                                                                     heatmap_legend_param = list(title = "ln(1+nTPM)"),
                                                                     show_heatmap_legend = T, 
                                                                     cluster_rows = T, top_annotation = ha,
                                                                     column_title = "Expression across tissues"))), NULL,
                    nrow = 1, rel_widths = c(0.02, 1, 0.02)),
          plot_grid(NULL,drug_tissue_score_plot.1, NULL,
                    nrow = 1, rel_widths = c(0.12, 1, 0.12)),
          ncol = 1, labels = c("(a)", "(b)", "(c)"),
          rel_heights = c(0.45, 0.65, 1.35))
dev.off()

# K-means
Elbow1 <- ClusterR::Optimal_Clusters_KMeans(cluster.df[,c(1:2)],
                                            max_clusters = 10, num_init = 10, seed = 123,
                                            criterion = "WCSSE",
                                            initializer = "kmeans++")
nk <- 4
fit <- ClusterR::KMeans_rcpp(cluster.df[,c(1:2)], 
                             clusters = nk, num_init = 10, max_iters = 100, seed = 123,
                             initializer = "kmeans++")
ord <- order(rowMeans(fit$centroids)) %>% rev()
drug_class <- factor(match(fit$cluster, ord),
                     levels = 1:nk,
                     labels = c("No", "Low", "Moderate","High"),
                     ordered = TRUE)
result.w.anno.full$distanceCluster <- drug_class

Elbow2 <- ClusterR::Optimal_Clusters_KMeans(cluster.df[,-c(1:2)],
                                            max_clusters = 10, num_init = 10, seed = 123,
                                            criterion = "WCSSE",
                                            initializer = "kmeans++")
nk <- 4
fit <- ClusterR::KMeans_rcpp(cluster.df[,-c(1:2)],
                             clusters = nk, num_init = 10, max_iters = 100, seed = 123,
                             initializer = "kmeans++")
ord <- order(rowMeans(fit$centroids))
drug_class <- factor(match(fit$cluster, ord),
                     levels = 1:nk,
                     labels = c("No", "Low", "Moderate","High"),
                     ordered = TRUE)
result.w.anno.full$affinityCluster <- drug_class

# Try Latent Profile Analysis 
# library(tidyLPA)
# test <- cluster.df[,c(1:2)] %>%
#   estimate_profiles(4,
#                     variances = c("equal", "varying"),
#                     covariances = c("zero", "varying")) 
# 
# test %>% compare_solutions(statistics = c("AIC", "BIC"))

# Preliminary view
Dend <- factoextra::fviz_dend(hclust(dist(t(scale(cluster.df)), method = "euclidean"), "complete"), 
                              k = 2, horiz = T, color_labels_by_k = F, k_colors = "jco",
                              cex = 0.65, lwd = 0.1,
                              rect = T, rect_border = "jco", rect_fill = T) + 
  theme_void() + 
  labs(title = "") +
  theme(legend.position = "none")
Dend$layers[[1]]$aes_params$linewidth <- 1 # Overwrite thickness

Elbow.df <- data.frame(Elbow1, Elbow2) 
colnames(Elbow.df) <- c("Network", "Tissue-specific drug exposure")

Elbow.df <- Elbow.df %>% rownames_to_column("K") %>%
  pivot_longer(cols = -K, names_to = "Data", values_to = "WCSSE") %>%
  mutate(K = as.numeric(K))

cluster.plot.1 <- Elbow.df %>%
  ggplot(aes(x = K, y = WCSSE, group = Data)) +
  geom_line() + geom_point(size = 2, color = "grey20") +
  geom_vline(data=filter(Elbow.df, Data=="Network"), aes(xintercept=4), colour="red", linetype = "dashed") + 
  geom_vline(data=filter(Elbow.df, Data=="Tissue-specific drug exposure"), aes(xintercept=4), colour="red", linetype = "dashed") + 
  labs(x = "Cluster (k)") +
  scale_x_continuous(breaks = seq(1,10,1)) +
  facet_wrap(.~Data, scales = "free") +
  mytheme + theme(legend.position = "none")

cluster.plot.2 <- result.w.anno.full %>%
  ggplot(aes(x = distance.value, y = zscore, 
             color = distanceCluster)) +
  stat_density_2d(geom = "polygon", bins = 75, fill = NA,
                  aes(color = distanceCluster)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = rev(pal_jama()(4))) + # Reversed
  geom_rug(aes(color = distanceCluster)) +
  labs(x = "Network distance", y = "Network Z-score",
       title = "K-mean clustering of\nnetwork results") +
  scale_color_manual(values = rev(pal_jama()(4))) + # Reversed
  mytheme + theme(legend.position = "none")

pca.dat <- prcomp(drug_tissue_score)

pc_eigenvalues <- pca.dat$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

pc_eigenvalues <- pc_eigenvalues[c(1:2),]

cluster.plot.3 <- pca.dat$x %>% as.data.frame %>%
  rownames_to_column("Drug") %>%
  left_join(., result.w.anno.full, "Drug") %>%
  ggplot(aes(x = PC1, y = PC2, color = affinityCluster)) +
  stat_density_2d(geom = "polygon", bins = 75, fill = NA,
                  aes(color = affinityCluster)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = rev(pal_jama()(4))) + # Reversed
  geom_rug(aes(color = affinityCluster)) +
  labs(x = paste0("PC1 (", round(pc_eigenvalues$pct[1], digits = 2), "%)"),
       y = paste0("PC2 (", round(pc_eigenvalues$pct[2], digits = 2), "%)"),
       title = "K-mean clustering of\ntissue-specific drug exposure score") +
  scale_color_manual(values = rev(pal_jama()(4))) + # Reversed
  mytheme + theme(legend.position = "none")

table(result.w.anno.full$distanceCluster, result.w.anno.full$affinityCluster)

tiff("Results/SFigure8.tiff", units="in", height=12.5, width=8.5, res=300)
plot_grid(plot_grid(NULL, Dend, NULL, rel_widths = c(0.05, 1, 0.1), nrow = 1), 
          cluster.plot.1,
          plot_grid(cluster.plot.2, cluster.plot.3, nrow = 1),
          ncol = 1, labels = c("(a)", "(b)", "(c)"),
          rel_heights = c(1.85,0.7,1))
dev.off()

# Classification scheme
## 3 affinity clusters
scheme <- tibble(distanceCluster = c("No", "Low", "Moderate", "High"),
                 Agonist = c("No", "No", "No", "No"),
                 None = c("No", "No", "Low", "Moderate"),
                 Antagonist = c("No", "Low", "Moderate", "High"))
## 2 affinity clusters
scheme <- tibble(distanceCluster = c("No", "Low", "Moderate", "High"),
                 Agonist = c("No", "No", "No", "No"),
                 Antagonist = c("No", "Low", "Moderate", "High")) 

## 4 affinity clusters
scheme <- tibble(distanceCluster = c("No", "Low", "Moderate", "High"),
                 No = c("No", "No", "No", "No"),
                 Low = c("No", "Low", "Low", "Moderate"),  
                 Moderate = c("No", "Low", "Moderate", "High"),
                 High = c("No", "Moderate", "High", "High")) %>%
  pivot_longer(-distanceCluster, names_to = "affinityCluster", values_to = "cholinergicCluster_base")

result.w.anno.full <- result.w.anno.full %>% 
  left_join(., scheme, c("distanceCluster", "affinityCluster")) %>%
  mutate(distanceCluster = factor(distanceCluster, levels = c("No", "Low", "Moderate", "High"),
                                  ordered = T),
         affinityCluster = factor(affinityCluster, levels = c("No", "Low", "Moderate", "High"), 
                                  ordered = T))

# Include Tanimoto chemical similarity, drug target signature, and MOA signature into consideration
library(proxy)
mols <- rcdk::parse.smiles(result.w.anno.full$SMILES)
fps <- lapply(mols, rcdk::get.fingerprint, type = "circular") # Extended connectivity (ECFP)
sim_matrix <- fingerprint::fp.sim.matrix(fps, method = "tanimoto")
rownames(sim_matrix) <- result.w.anno.full$Drug
colnames(sim_matrix) <- result.w.anno.full$Drug

Drug.target.df <- read.csv("1_PDI/DrugTarget.csv") %>%
  mutate(x = 1,
         Drug = gsub("6,7-di(sulfanyl-?S)",
                     "6,7-di(sulfanyl-κS)",
                     Drug, fixed = T),
         Drug = gsub("(2S) N-acetyl-L-alanyl-?L-phenylalanyl-chloroethylketone",
                     "(2S) N-acetyl-L-alanyl-αL-phenylalanyl-chloroethylketone",
                     Drug, fixed = T),
         Drug = gsub("5,7,2?-trihydroxy-6,8-dimethoxyflavone",
                     "5,7,2′-trihydroxy-6,8-dimethoxyflavone",
                     Drug, fixed = T)
  ) 

Target.list <- Drug.target.df %>%
  dplyr::filter(Target %in% rownames(CHRM.exp)) %>%
  distinct(Drug) %>% pull(Drug)

Drug.target.mat <- Drug.target.df %>%
  pivot_wider(names_from = "Target", values_from = "x", values_fill = 0) %>%
  dplyr::filter(Drug %in% result.w.anno.full$Drug) %>%
  column_to_rownames("Drug")
jaccard_target <- simil(as.matrix(Drug.target.mat),   method = "Jaccard") %>% as.matrix() 
jaccard_target <- jaccard_target[rownames(sim_matrix), colnames(sim_matrix)]

MOA.df <- read.csv("1_PDI/MOA.csv") %>%
  dplyr::select(pert_iname, moa) %>%
  dplyr::filter(moa != "") %>%
  separate_longer_delim(moa, delim = " | ") %>%
  mutate(pert_iname = str_to_sentence(pert_iname)) %>%
  dplyr::filter(pert_iname %in% result.w.anno.full$Drug) 

MOA.drug.list <- MOA.df %>%
  dplyr::filter(str_detect(moa, "acetylchol")) %>%
  distinct(pert_iname) %>% pull(pert_iname)

MOA.mat <- MOA.df %>%
  mutate(x = 1) %>%
  pivot_wider(names_from = "moa", values_from = "x", values_fill = 0) %>%
  column_to_rownames("pert_iname")
jaccard_moa <- simil(as.matrix(MOA.mat),   method = "Jaccard") %>% as.matrix() 

moa_full <- matrix(NA, nrow = length(result.w.anno.full$Drug), ncol = length(result.w.anno.full$Drug),
                   dimnames = list(result.w.anno.full$Drug, result.w.anno.full$Drug))
moa_full[rownames(jaccard_moa), colnames(jaccard_moa)] <- jaccard_moa

identical(rownames(sim_matrix), rownames(jaccard_target))
identical(colnames(sim_matrix), colnames(jaccard_target))
identical(rownames(sim_matrix), rownames(moa_full))

# Weighted average, ignoring NA
# Where MOA available:   tanimoto=0.33, target=0.33, moa=0.33
# Where MOA missing:     tanimoto=0.50, target=0.50, moa=0.00
moa_weight  <- ifelse(is.na(moa_full), 0, 1/3)
tani_weight <- ifelse(is.na(moa_full), 0.50, 1/3)
targ_weight <- ifelse(is.na(moa_full), 0.50, 1/3)
moa_filled  <- replace(moa_full, is.na(moa_full), 0)

unified <- (sim_matrix * tani_weight) + (jaccard_target * targ_weight) + (moa_filled * moa_weight)
dist_matrix <- as.dist(1 - unified)

View(unified[,c("Atropine", "Amitriptyline", "Nortriptyline", "Lorlatinib", "Candesartan")])

moa_idx <- which(rownames(dist_matrix) %in% unique(c(Target.list, MOA.drug.list)))

dist_profile <- dist_matrix[, moa_idx]
moa_min_distance <- apply(dist_profile, 1, min)

dist_profile_summary <- data.frame(
  Drug              = rownames(dist_profile),
  min_dist          = apply(dist_profile, 1, min),
  mean_dist         = apply(dist_profile, 1, mean),
  median_dist       = apply(dist_profile, 1, median),
  sd_dist           = apply(dist_profile, 1, sd),
  q10_dist          = apply(dist_profile, 1, quantile, probs = 0.10),
  n_close_03        = apply(dist_profile, 1, function(x) sum(x < 0.3)),
  prop_close_03     = apply(dist_profile, 1, function(x) mean(x < 0.3))
)

Elbow3 <- ClusterR::Optimal_Clusters_KMeans(as.data.frame(moa_min_distance),
                                            max_clusters = 10, num_init = 10, seed = 123,
                                            criterion = "WCSSE",
                                            initializer = "kmeans++")
fit <- ClusterR::KMeans_rcpp(as.data.frame(moa_min_distance), 
                             clusters = 4, num_init = 10, max_iters = 100, seed = 123,
                             initializer = "kmeans++")
ord <- order(rowMeans(fit$centroids)) %>% rev()
drug_class <- factor(match(fit$cluster, ord),
                     levels = 1:4,
                     labels = c("Unlikely", "Possible", "Likely" , "Very likely"),
                     ordered = TRUE)

BR.drug.sim <- ggdraw() + 
  draw_image("Biorender/DrugSimiliarity.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

tiff("Results/SFigure9.tiff", units="in", height=10.5, width=8.5, res=300)
plot_grid(BR.drug.sim,
          plot_grid(data.frame(WCSSE = Elbow3,
                               K = seq(1,10,1)) %>%
                      ggplot(aes(x = K, y = WCSSE)) +
                      geom_line() + geom_point(size = 2, color = "grey20") +
                      geom_vline(aes(xintercept=4), colour="red", linetype = "dashed") + 
                      labs(x = "Cluster (k)") +
                      scale_x_continuous(breaks = seq(1,10,1)) +
                      mytheme + theme(legend.position = "none"),
                    as.data.frame(moa_min_distance) %>%
                      mutate(Confidence = drug_class) %>%
                      ggplot(aes(x = 1-moa_min_distance, fill = Confidence)) +
                      geom_histogram(color =  "black") +
                      labs(x = "Similarity",
                           y = "Number of drugs",
                           fill = "", title = "") +
                      scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
                      scale_fill_bmj() +
                      mytheme +
                      theme(legend.position = c(0.99, 0.98),
                            legend.background = element_blank(),
                            legend.justification = c(1, 0.85)),
                    nrow = 1,
                    labels = c("(b)", "(c)"),
                    align = "h", axis = "tb"),
          ncol = 1, labels = c("(a)", ""),
          rel_heights = c(2, 0.85))
dev.off()

confidence.scheme <- tibble(cholinergicCluster_base = c("No", "Low", "Moderate", "High"),
                            Unlikely = c("No", "No", "No", "No"),
                            Possible = c("No", "Low", "Moderate", "High"),
                            Likely = c("No", "Low", "Moderate", "High"),
                            `Very likely` = c("No", "Low", "Moderate", "High")) %>%
  pivot_longer(-cholinergicCluster_base, names_to = "Confidence", values_to = "cholinergicCluster")

result.w.anno.full <- result.w.anno.full %>% 
  # dplyr::select(-cholinergicCluster, -moa_min_distance, -Confidence) %>%
  left_join(., as.data.frame(moa_min_distance) %>%
              mutate(Confidence = drug_class) %>% rownames_to_column("Drug"),
            "Drug") %>%
  left_join(., confidence.scheme, c("cholinergicCluster_base", "Confidence")) %>%
  mutate(cholinergicCluster_base = factor(cholinergicCluster_base, ordered = T,
                                          levels = c("No", "Low", "Moderate", "High")),
         cholinergicCluster = factor(cholinergicCluster, ordered = T,
                                     levels = c("No", "Low", "Moderate", "High")),
         Confidence = factor(Confidence, ordered = T, 
                             levels = c(c("Unlikely", "Possible", "Likely" , "Very likely")))
  )
  
# Descriptive statistics
write.csv(result.w.anno.full, "ISACS.csv", row.names = F)

table(result.w.anno.full$ACSBC.Anticholinergic.Activity)
table(result.w.anno.full$ORCA)
table(result.w.anno.full$cholinergicCluster_base)
round(table(result.w.anno.full$cholinergicCluster_base)/nrow(result.w.anno.full)*100, 2) # proportion

table(result.w.anno.full$cholinergicCluster)
round(table(result.w.anno.full$cholinergicCluster)/nrow(result.w.anno.full)*100, 2) # proportion

table(result.w.anno.full$cholinergicCluster, result.w.anno.full$Ground_truth)
table(result.w.anno.full$cholinergicCluster, fct_rev(result.w.anno.full$ACSBC.Anticholinergic.Activity))
table(result.w.anno.full$cholinergicCluster, fct_rev(result.w.anno.full$ORCA))
table(fct_rev(result.w.anno.full$ORCA), fct_rev(result.w.anno.full$ACSBC.Anticholinergic.Activity))

result.w.anno.full %>%
  ggplot(aes(x = QSAR.PCA, y = distance.value*zscore, 
             fill = cholinergicCluster)) +
  stat_density_2d(geom = "polygon", bins = 150,
                  aes(alpha = ..level.., fill = cholinergicCluster)) +
  geom_point(alpha = 0.5, size = 2, color = "white", pch = 21) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_alpha_continuous(range = c(0.3, 0.9)) +
  scale_fill_jama() +
  geom_rug(aes(color = cholinergicCluster)) +
  labs(x = "Mean averaged PSICHIC score", y = "Network distance × z-score") +
  scale_color_jama() +
  mytheme

lab.drugs <- ATC.df.enrich %>% 
  dplyr::filter(ATC_3rd %in% c("ANTIHISTAMINES FOR SYSTEMIC USE (R06A)", 
                               "ANTIDEPRESSANTS (N06A)", "UROLOGICALS (G04B)", 
                               "ANTIPSYCHOTICS (N05A)", "ANTICHOLINERGIC AGENTS (N04A)")) %>%
  pull(name)

axis.adj <- c(0.1, 0.1)
clu.p1 <- result.w.anno.full %>%
  mutate(lab = case_when(Drug %in% lab.drugs & cholinergicCluster != "No"~Drug,
                         T~"")) %>%
  arrange(cholinergicCluster) %>%
  ggplot(aes(x = QSAR_score, y = zscore*distance.value, 
             fill = cholinergicCluster)) +
  stat_density_2d(geom = "polygon", bins = 150,
                  aes(alpha = ..level.., fill = cholinergicCluster)) +
  geom_point(alpha = 0.5, size = 2, color = "white", pch = 21) +
  scale_x_continuous(expand = expansion(mult = axis.adj)) +
  scale_y_continuous(expand = expansion(mult = axis.adj)) +
  scale_alpha_continuous(range = c(0.3, 0.9)) +
  scale_fill_jama() +
  geom_rug(aes(color = cholinergicCluster), data = ~ filter(.x, cholinergicCluster == "No"),  alpha = 0.3) +
  geom_rug(aes(color = cholinergicCluster), data = ~ filter(.x, cholinergicCluster != "No"),  alpha = 0.8) +
  labs(x = "Mean mAChR PSICHIC score", y = "Network distance × z-score") +
  scale_color_jama() +
  mytheme +
  theme(legend.position = "none")
clu.p2 <- result.w.anno.full %>%
  ggplot(aes(x = QSAR_score, color = cholinergicCluster)) + 
  geom_density(linewidth = 1) +
  scale_x_continuous(expand = expansion(mult = axis.adj)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_color_jama() +
  theme_void() +
  theme(legend.position = "none")
clu.p3 <- result.w.anno.full %>%
  ggplot(aes(x = zscore*distance.value, color = cholinergicCluster)) + 
  geom_density(linewidth = 1) +
  scale_x_continuous(expand = expansion(mult = axis.adj)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_color_jama() +
  theme_void() +
  theme(legend.position = "none") + 
  coord_flip()

clu.p <- clu.p2 + patchwork::plot_spacer() + clu.p1 + clu.p3 + 
  patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

clu.p <- insert_xaxis_grob(clu.p1, clu.p2, grid::unit(.2, "null"), position = "top")
clu.p <- insert_yaxis_grob(clu.p, clu.p3, grid::unit(.2, "null"), position = "right")
ggdraw(clu.p)

KW.test <- result.w.anno.full %>% kruskal_test(distance.value ~ cholinergicCluster)    
dunn.test <- result.w.anno.full %>%
  dunn_test(distance.value ~ cholinergicCluster, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "cholinergicCluster") %>%
  mutate(y.position = y.position + 0.05) %>% dplyr::select(-groups)
clu.p4 <- result.w.anno.full %>%
  ggplot(aes(x = cholinergicCluster, y = distance.value, 
             fill = cholinergicCluster)) + 
  geom_violin(alpha = 0.25) + geom_boxplot(width = 0.15) +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.05,
                     tip.length = 0.02, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_fill_jama() + labs(x = "", y = "Network distance", subtitle = get_test_label(KW.test, detailed = TRUE)) + 
  mytheme + theme(legend.position = "none",
                  axis.text.x = element_text(angle = 25, hjust = 1))
  
KW.test <- result.w.anno.full %>% kruskal_test(zscore ~ cholinergicCluster)    
dunn.test <- result.w.anno.full %>%
  dunn_test(zscore ~ cholinergicCluster, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "cholinergicCluster") %>%
  mutate(y.position = y.position + 0.05) %>% dplyr::select(-groups)
clu.p5 <- result.w.anno.full %>%
  ggplot(aes(x = cholinergicCluster, y = zscore, 
             fill = cholinergicCluster)) + 
  geom_violin(alpha = 0.25) + geom_boxplot(width = 0.15) +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.05,
                     tip.length = 0.02, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_fill_jama() + labs(x = "", y = "Network z-score", subtitle = get_test_label(KW.test, detailed = TRUE)) + 
  mytheme + theme(legend.position = "none",
                  axis.text.x = element_text(angle = 25, hjust = 1))

KW.test <- result.w.anno.full %>% kruskal_test(QSAR_score ~ cholinergicCluster)    
dunn.test <- result.w.anno.full %>%
  dunn_test(QSAR_score ~ cholinergicCluster, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "cholinergicCluster") %>%
  mutate(y.position = y.position + 0.05) %>% dplyr::select(-groups)
clu.p6 <- result.w.anno.full %>%
  ggplot(aes(x = cholinergicCluster, y = QSAR_score, 
             fill = cholinergicCluster)) + 
  geom_violin(alpha = 0.25) + geom_boxplot(width = 0.15) +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.05,
                     tip.length = 0.02, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_fill_jama() + labs(x = "", y = "Mean mAChR PSICHIC score", subtitle = get_test_label(KW.test, detailed = TRUE)) + 
  mytheme + theme(legend.position = "none",
                  axis.text.x = element_text(angle = 25, hjust = 1))

# ATC Enrichment
ORA.result1 <- list()
for(i in levels(result.w.anno.full$cholinergicCluster)) {
  ORA.result1[[i]] <- enricher(result.w.anno.full %>% dplyr::filter(cholinergicCluster %in% i) %>% pull(Drug), 
                               TERM2GENE = ATC.df.enrich)
}
ORA.result1 <- merge_result(ORA.result1)

ORA.result2 <- list()
for(i in levels(result.w.anno.full$ORCA)) {
  ORA.result2[[i]] <- enricher(result.w.anno.full %>% drop_na(ORCA) %>% dplyr::filter(ORCA %in% i) %>% pull(Drug), 
                               TERM2GENE = ATC.df.enrich)
}
ORA.result2 <- merge_result(ORA.result2)

ORA.result3 <- list()
for(i in levels(result.w.anno.full$ACSBC.Anticholinergic.Activity)) {
  ORA.result3[[i]] <- enricher(result.w.anno.full %>% drop_na(ACSBC.Anticholinergic.Activity) %>% dplyr::filter(ACSBC.Anticholinergic.Activity %in% i) %>% pull(Drug), 
                               TERM2GENE = ATC.df.enrich)
}
ORA.result3 <- merge_result(ORA.result3)

x_numeric <- c("No" = 4, "Low" = 3, "Moderate" = 2, "High" = 1)
category_order <- ORA.result1@compareClusterResult %>%
  dplyr::filter(p.adjust < 0.05) %>%
  group_by(Description) %>%
  dplyr::summarise(x_levels = list(unique(Cluster)),
                   n_levels  = n_distinct(Cluster),
                   mean_x    = mean(x_numeric[Cluster])) %>%
  arrange(mean_x) %>%
  pull(Description)

Fig2f <- ORA.result1@compareClusterResult %>% mutate(Scale = "ISACS") %>%
  # rbind(ORA.result2@compareClusterResult %>% mutate(Scale = "ORCA"),
  # ORA.result3@compareClusterResult %>% mutate(Scale = "ACSBC")) %>%
  dplyr::filter(Cluster != "No consensus") %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Scale = factor(Scale, 
                        levels = c("ACSBC", "ORCA", "ISACS")),
         Cluster = factor(Cluster, 
                          levels = levels(result.w.anno.full$cholinergicCluster)),
         Description = factor(Description, levels = category_order)) %>%
  ggplot(aes(x = Cluster, y = fct_rev(Description), 
             size = Count, fill = Cluster, alpha = -log10(p.adjust))) +
  geom_point(pch = 21, color = "white") + 
  labs(x = "", y = "", alpha = "-log10(FDR-adjusted p-value)",
       fill = "Level",
       title = "ATC 3rd Level Enrichment Analysis") + 
  scale_size(range = c(5,10)) +
  scale_alpha(range = c(0.5, 0.9)) +
  scale_fill_jama() +
  # facet_grid(~Scale) +
  guides(fill = guide_legend(override.aes = list(size = 5, color = "black")),
         size = guide_legend(override.aes = list(color = "black", pch = 16)),
         alpha = guide_legend(override.aes = list(color = "black", fill = "grey50", 
                                                  size = 5))) +
  mytheme +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1), 
        panel.grid.major = element_line(color = "white", linewidth = 1))
Fig2f

# Supplementary Table 7
write_clip(result.w.anno.full %>% 
             mutate(moa_min_distance = 1-moa_min_distance,
                    `Chew's SAA list` = SAA_class) %>%
             dplyr::select(DrugBank.ID, Drug, 
                           distance.value, pval, zscore, zscore_sig,
                           ends_with("_QSAR"),
                           distanceCluster, affinityCluster, cholinergicCluster_base,
                           moa_min_distance, Confidence,
                           cholinergicCluster,
                           ORCA, ACSBC.Anticholinergic.Activity, 
                           ADS, ABS, ALS, 
                           AEC,`Chew's SAA list`, ACB, AAS, 
                           ARS, ABC, ATS,
                           pgp_sub, BBB, Neurotoxicity.DI
             ))

result.w.anno.full %>% 
  pivot_longer(c(ACSBC.Anticholinergic.Activity, ORCA)) %>%
  drop_na(value) %>%
  ggplot(aes(x = value, y = 1-moa_min_distance,
             fill = value)) +
  geom_violin(alpha = 0.25) + 
  geom_boxplot(width = 0.1) +
  labs(y = "Similarity", x = "") + mytheme +
  facet_grid(~name) +
  scale_fill_manual(values = rev(pal_jama()(4))) +
  theme(legend.position = "none")

# Supplementary Table 8
write_clip(ORA.result1@compareClusterResult %>% dplyr::select(-ID))

BR.ISACS <- ggdraw() + 
  draw_image("Biorender/ISACS_development.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

Fig2ab <- plot_grid(BR.ISACS, clu.p, ncol = 1,
                    rel_heights = c(0.4, 1),
                    labels = c("(a)", "(b)"))
Fig2cde <- plot_grid(clu.p4, clu.p5, clu.p6, ncol = 1,
                     labels = c("(c)", "(d)", "(e)"),
                     label_x = -0.1)

tiff("Results/Figure2.tiff", units="in", height=21, width=14.5, res=300)
plot_grid(plot_grid(Fig2ab, 
                    Fig2cde, 
                    rel_widths = c(1, 0.55)),
          Fig2f,
          rel_heights = c(1.35, 1),
          labels = c("", "(f)"),
          ncol = 1)
dev.off()

# SIDER/ONSIDES analysis with scale comparisons

## Distribution across scales
antichol.scales <- c("ORCA", "ACSBC.Anticholinergic.Activity", # "ADS", 
                     "cholinergicCluster")

scale.comp.df <- result.w.anno.full %>% dplyr::select(Drug, all_of(antichol.scales)) %>%
  mutate(ORCA = factor(ORCA, levels = c("No", "Low", "Moderate", "High"), ordered = T),
         ACSBC.Anticholinergic.Activity = factor(ACSBC.Anticholinergic.Activity, levels = c("No", "Low", "Moderate", "High"), 
                                                 ordered = T)) %>%
  drop_na(ORCA, ACSBC.Anticholinergic.Activity)

### Apple-to-apple
nk <- 4
scale.comp.bar1 <- scale.comp.df %>%
  pivot_longer(cols = antichol.scales) %>% drop_na() %>%
  mutate(name = case_when(name == "ACSBC.Anticholinergic.Activity"~"ACSBC",
                          name == "cholinergicCluster"~"ISACS",
                          T~name),
         name = factor(name, 
                       levels = c("ACSBC","ADS", "ORCA", "ISACS"))) %>%
  group_by(name, value) %>%
  count(name, value) %>%
  group_by(name) %>%
  mutate(total = sum(n),
         percent = n/total) %>%
  ggplot(aes(x = name, y = n, fill = fct_rev(value))) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
            position = position_fill(vjust = 0.5),
            size = 4,
            color = "white") +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Proportion (%)", 
       title = "Complete case across the 3 scales",
       subtitle = paste0("n = ", dim(scale.comp.df)[1]),
       fill = "Levels") +
  scale_fill_manual(values = rev(pal_jama()(nk))) + # Reversed
  mytheme +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))

### Apple-to-banana
scale.comp.bar2 <- result.w.anno.full %>% dplyr::select(Drug, all_of(antichol.scales)) %>%
  mutate(ORCA = factor(ORCA, levels = c("No", "Low", "Moderate", "High"), ordered = T),
         ACSBC.Anticholinergic.Activity = factor(ACSBC.Anticholinergic.Activity, levels = c("No", "Low", "Moderate", "High"), 
                                                 ordered = T)) %>%
  pivot_longer(cols = antichol.scales) %>% drop_na() %>%
  mutate(name = case_when(name == "ACSBC.Anticholinergic.Activity"~"ACSBC",
                          name == "cholinergicCluster"~"ISACS",
                          T~name),
         name = factor(name, 
                       levels = c("ACSBC", "ADS", "ORCA", "ISACS"))) %>%
  group_by(name, value) %>%
  count(name, value) %>%
  group_by(name) %>%
  mutate(total = sum(n),
         percent = n/total) %>%
  ggplot(aes(x = name, y = n, fill = fct_rev(value))) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
            position = position_fill(vjust = 0.5),
            size = 4,
            color = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  geom_text(data = data.frame(name = c("ACSBC", "ORCA", # "ADS",
                                       "ISACS"),
                              n = c(sum(!is.na(result.w.anno.full$ACSBC.Anticholinergic.Activity)),
                                    sum(!is.na(result.w.anno.full$ORCA)),
                                    # sum(!is.na(result.w.anno.full$ADS)),
                                    sum(!is.na(result.w.anno.full$cholinergicCluster))
                              )), y = 1.03, size = 4.75,
            aes(x = name, label = paste0("n = ", n)), inherit.aes = F) +
  labs(x = "", y = "Proportion (%)", 
       fill = "Levels") +
  scale_fill_manual(values = rev(pal_jama()(nk))) + # Reversed
  mytheme +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        plot.margin = margin(t = 25, r = 0, b = 0, l = 0)) + 
  labs(title = "All drugs available per scale") +
  theme(plot.title = element_text(vjust = 7.25))

Fig2e <- plot_grid(scale.comp.bar1, scale.comp.bar2, labels = c("(e)"))

# tiff(paste0("Results/Props.tiff"), units="in", height=5, width=9.5, res=300)
Fig2e
# dev.off()

## SIDER processing
library(rms)
SIDER.list <- list.files(path = "10_SIDER_4.1", full.names = T,
                         pattern = "\\.tsv$")
SIDER.files <- list()
for (i in SIDER.list) {
  df <- read.delim(i, header = F)
  names(df) <- paste0("col", seq_along(df))
  SIDER.files[[i]] <- df
}
SIDER.df <- right_join(SIDER.files[["10_SIDER_4.1/drug_names.tsv"]],
                       SIDER.files[["10_SIDER_4.1/meddra_all_label_se.tsv"]][,c(2, 5,7)],
                       c("col1"="col2")) %>% distinct
colnames(SIDER.df) <- c("CompoundID", "SIDER.drug", "MedDRA_type", "SE_name")
dim(SIDER.df %>% dplyr::filter(SE_name %in% "Delirium"))
SIDER.df <- SIDER.df %>%
  mutate(PubChem.Compound.ID = as.numeric(sub("^CID10*", "", CompoundID))) %>%
  left_join(., read.csv("1_PDI/1_DrugBank_v5.1.13/drug links.csv") %>% 
              dplyr::select(Name, PubChem.Compound.ID), 
            c("PubChem.Compound.ID")) %>%
  mutate(SIDER.drug = str_to_sentence(SIDER.drug)) %>%
  mutate(Drug = case_when(SIDER.drug == Name~SIDER.drug,
                          is.na(Name)~SIDER.drug,
                          SIDER.drug != Name~Name
  )) %>% 
  dplyr::select(-SIDER.drug, -Name)

rm(SIDER.files)
head(SIDER.df %>% dplyr::filter(MedDRA_type %in% "PT") %>% group_by(SE_name) %>% count)

## ONSIDES processing
ONSIDES.df <- read.csv("10_ONSIDES/ONSIDES.processed.csv") %>%
  dplyr::select(meddra_name, meddra_term_type, rxnorm_name) %>% 
  mutate(Resource = "ONSIDES",
         rxnorm_name = str_to_sentence(rxnorm_name))
colnames(ONSIDES.df) <- c("SE_name", "MedDRA_type", "Drug", "Resource")

head(ONSIDES.df %>% dplyr::filter(MedDRA_type %in% "PT") %>% group_by(SE_name) %>% count)

## Anticholinergic effects
central.SEs <- c("Agitation", "Ataxia", 
                 "Confusional state", "Coordination abnormal", 
                 "Delirium", "Delusion", "Disturbance in attention", 
                 "Dizziness", "Hallucination", 
                 "Memory impairment", "Speech disorder") # n = 11
peripheral.SEs <- c("Body temperature increased",
                    "Constipation", "Dehydration", 
                    "Dry eye", "Dry mouth", "Dry skin",
                    "Flushing", "Hot flush", 
                    "Mydriasis", "Tachycardia", 
                    "Thirst", "Urinary retention", "Vision blurred") # n = 13
# Slow speech (n = 1)
# Fall (n = 14)
# Dry throat (n = 38)
# Cognitive disorder (n = 46)

merged.SE.df <- SIDER.df %>% 
  dplyr::select(SE_name, MedDRA_type, Drug) %>%
  mutate(Resource = "SIDER") %>%
  rbind(ONSIDES.df) %>%
  dplyr::filter(MedDRA_type %in% "PT") %>%
  group_by(SE_name, MedDRA_type, Drug) %>%
  summarise(Resource = paste(unique(Resource), collapse = "; "),
            .groups = "drop") 
length(unique(merged.SE.df$Drug)) # Total drugs listed in AE dbs

# Only get those that are overlapping with our data
merged.SE.df.final <- merged.SE.df %>% 
  dplyr::select(-Resource, -MedDRA_type) %>%
  dplyr::filter(Drug %in% unique(result.w.anno.full$Drug)) %>% 
  mutate(x = 1) %>%
  pivot_wider(names_from = "SE_name", 
              values_from = "x", 
              values_fill = 0) %>% 
  dplyr::select(Drug, c(peripheral.SEs, central.SEs)) %>%
  pivot_longer(cols = c(peripheral.SEs, central.SEs),
               names_to = "SE_name", values_to = "x") 

length(unique(merged.SE.df.final$Drug)) # Total drugs overlap between AE dbs & our data

SE.res <- merged.SE.df.final %>% inner_join(., result.w.anno.full, "Drug") %>%
  mutate(val = case_when(SE_name %in% central.SEs~"Central",
                         SE_name %in% peripheral.SEs~"Peripheral"),
         zscore_sig = factor(zscore_sig, 
                             levels = c("Significant", "NS"))) 

as.data.frame(table(SE.res$SE_name, SE.res$x)) %>% 
  dplyr::filter(Var2 == 1) %>% dplyr::select(-Var2) %>%
  arrange(Freq)

SE.res <- SE.res %>% left_join(., as.data.frame(table(SE.res$SE_name, SE.res$x)) %>% dplyr::filter(Var2 == 1) %>% dplyr::select(-Var2),
                               c("SE_name" = "Var1")) %>%
  dplyr::filter(Freq > 50) %>%
  mutate(facet_lab = paste0(SE_name, "\n(n = ", Freq, ")"))

# HT for viz
SIDER.ht <- SE.res %>%
  dplyr::select(Drug, SE_name, val,x, drugbank_id,
                zscore, cholinergicCluster, QSAR_score,
                ACSBC.Anticholinergic.Activity, ORCA
  ) %>% distinct() %>%
  mutate(val = case_when(x == 1~val,
                         x == 0~NA)) %>%
  pivot_wider(names_from = "SE_name", values_from = "val") %>%
  dplyr::filter(x != 0) %>% distinct %>% dplyr::select(-x) %>%
  column_to_rownames("Drug")

SIDER.ht.anno <- SIDER.ht %>% dplyr::select(zscore, QSAR_score, cholinergicCluster, # drugbank_id, 
                                            ACSBC.Anticholinergic.Activity, ORCA) %>%
  dplyr::rename("ISACS" = "cholinergicCluster") %>%
  dplyr::rename("ACSBC" = "ACSBC.Anticholinergic.Activity") %>%
  dplyr::rename("Mean pKi" = "QSAR_score") %>%
  mutate(ISACS = factor(as.character(ISACS),
                        levels = c("High", "Moderate", "Low", "No")))

SIDER.ht.mat <- SIDER.ht %>% dplyr::select(central.SEs, peripheral.SEs)

SIDER.ht.mat.clu <- SIDER.ht.mat %>%
  mutate(across(everything(), ~ if_else(!is.na(.) & . != "", 1, 0))) %>%
  as.matrix()

row_hclust = hclust(vegan::vegdist(SIDER.ht.mat.clu, method = "jaccard"),
                    method = "complete")
col_hclust = hclust(vegan::vegdist(t(SIDER.ht.mat.clu), method = "jaccard"),
                    method = "complete")

col_fun = circlize::colorRamp2(c(5, 0, -1.302, -1.96, -10), 
                               c("#f0dddd", "#d6a9ab",  "#824446",  #"#b57d7f", 
                                 "#6e2b2d", "#4a1712")) 
col_fun2 = circlize::colorRamp2(c(7, 8, 9), 
                                c("#002c58", "white","#580000"))
ha <- rowAnnotation(df = SIDER.ht.anno %>%
                      dplyr::rename("Network z-score" = zscore) %>%
                      dplyr::rename("Mean mAChR PSICHIC score" = `Mean pKi`),
                    annotation_name_rot = 45, 
                    annotation_name_gp = gpar(fontsize = 10),
                    show_legend = F,
                    col = list("ISACS" = setNames(pal_jama("default")(4),
                                                  c("No", "Low", "Moderate", "High")),
                               "ACSBC" = setNames(pal_jama("default")(4),
                                                  c("No", "Low", "Moderate", "High")),
                               "ORCA" = setNames(pal_jama("default")(4),
                                                 c("No", "Low", "Moderate", "High")),
                               "Network z-score" = col_fun,
                               "Mean mAChR PSICHIC score" = col_fun2 
                    ))
drugs_to_show <- c(true_pos)

ha2 = rowAnnotation(drugnames = anno_mark(at = which(rownames(SIDER.ht.mat) %in% drugs_to_show),
                                          labels = rownames(SIDER.ht.mat)[rownames(SIDER.ht.mat) %in% drugs_to_show],
                                          labels_gp = gpar(fontsize = 10))
                    )

SIDER.ht <- oncoPrint(SIDER.ht.mat,
                      show_heatmap_legend = F,
                      show_column_names = T, 
                      show_row_names = F,
                      show_pct = F,
                      left_annotation = ha,
                      alter_fun = list(
                        Central = function(x, y, w, h)
                          grid.rect(x, y, w, h, gp = gpar(fill = "#9c1806", col = "#9c1806")),
                        Peripheral = function(x, y, w, h)
                          grid.rect(x, y, w, h, gp = gpar(fill = "#415494", col = "#415494"))),
                      col = c(Central = "#9c1806", Peripheral = "#415494"),
                      cluster_rows = row_hclust, cluster_columns = col_hclust) + ha2
SIDER.ht

lgd_zscore <- Legend(col_fun = col_fun, 
                     at = c(5, 0, -1.302, -1.96, -5, -10),  
                     title = "Network z-score", 
                     # direction = "vertical", title_position = "topleft",
                     direction = "horizontal", title_position = "topcenter",
                     legend_height = unit(3, "cm"), legend_width = unit(6, "cm"))
lgd_pKi <- Legend(col_fun = col_fun2,
                  title = "Mean mAChR PSICHIC score", 
                  # direction = "vertical", title_position = "topleft",
                  direction = "horizontal", title_position = "topcenter",
                  legend_height = unit(3, "cm"), legend_width = unit(6, "cm"))
lgd_txt = Legend(labels = c("High", "Moderate", "Low", "No"), 
                 title = "Anticholinergic level", 
                 row_gap = unit(0.7, "mm"),
                 column_gap = unit(6, "mm"),
                 nc = 1, legend_gp = gpar(fill = rev(pal_jama("default")(4))))
lgd_txt2 = Legend(labels = c("CNS", "PNS"), 
                  title = "Anticholinergic effect", 
                  row_gap = unit(0.7, "mm"),
                  column_gap = unit(6, "mm"),
                  nc = 1, legend_gp = gpar(fill = c("#9c1806", "#415494")))

lgd_continuous <- packLegend(lgd_zscore, lgd_pKi,
                             direction = "vertical", 
                             gap = unit(4, "mm"))

BR.AEs <- ggdraw() + 
  draw_image("Biorender/AnticholinergicEffects.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

tiff("Results/SFigure10.tiff", units="in", height=9, width=15, res=300)
plot_grid(BR.AEs,
          grid.grabExpr(draw(SIDER.ht, 
                             annotation_legend_list = packLegend(lgd_continuous, lgd_txt, lgd_txt2,
                                                                 column_gap = unit(1, "cm"),
                                                                 max_width = unit(11, "in"), direction = "horizontal"),
                             padding = unit(c(1, 20, 1, 1), "mm"),  # ← 20mm on the left
                             annotation_legend_side = "bottom")),
          rel_widths = c(0.8, 1),
          nrow = 1, labels = c("(a)", "(b)"))
dev.off()

# OR
SE.glm.df <- merged.SE.df.final %>%
  pivot_wider(names_from = "SE_name", 
              values_from = "x", values_fill = 0) %>%
  dplyr::select(Drug, central.SEs, peripheral.SEs) %>%
  inner_join(., result.w.anno.full, "Drug") %>%
  mutate(ISACS = relevel(factor(as.character(cholinergicCluster),
                                levels = c("No", "Low", "Moderate", "High")),
                         ref = "No"),
         ORCA = relevel(factor(ORCA, 
                               levels = c("No", "Low", "Moderate", "High")), 
                        ref = "No"),
         ACSBC = relevel(factor(ACSBC.Anticholinergic.Activity,
                                levels = c("No", "Low", "Moderate", "High")),
                         ref = "No"),
         PSICHIC = relevel(factor(as.character(affinityCluster),
                                levels = c("No", "Low", "Moderate", "High")),
                         ref = "No"),
         Network = relevel(factor(as.character(distanceCluster),
                                levels = c("No", "Low", "Moderate", "High")),
                         ref = "No"),
         Confidence = relevel(factor(as.character(Confidence),
                                     levels = c("Unlikely", "Possible", "Likely", "Very likely")),
                              ref = "Unlikely")
         ) 

# PK parameters
tiff("Results/SFigure11.tiff", units="in", height=3.5, width=8.5, res=300)
n.bins <- 20 
additional.cutoffs <- c(0.3, 0.7) 

bins1 <- seq(min(SE.glm.df$BBB), max(SE.glm.df$BBB), length.out = n.bins)
bins1 <- c(bins1, additional.cutoffs) %>% sort()
bins2 <- seq(min(SE.glm.df$pgp_sub), max(SE.glm.df$pgp_sub), length.out = n.bins)
bins2 <- c(bins2, additional.cutoffs) %>% sort()
plot_grid(SE.glm.df %>% 
            mutate(BBB.class = factor(case_when(BBB > 0.7~"Strong",
                                                BBB < 0.3~"Weak/no",
                                                T~"Moderate"),
                                      levels = c("Strong", "Moderate", "Weak/no"))) %>%
            ggplot(aes(x = BBB, fill = BBB.class)) +
            geom_histogram(color = "black", alpha = 0.75, breaks = bins1) +
            geom_vline(xintercept = 0.3, linetype = "dashed") +
            geom_vline(xintercept = 0.7, linetype = "dashed") +
            labs(x = "Probability of BBB permeability",
                 y = "Number of drugs",
                 fill = "", title = "") +
            scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
            scale_fill_bmj() +
            mytheme +
            theme(legend.position = c(0.99, 0.98),
                  legend.background = element_blank(),
                  legend.justification = c(1, 0.75)),
          SE.glm.df %>% 
            mutate(pgp_sub.class = factor(case_when(pgp_sub > 0.7~"Strong",
                                                    pgp_sub < 0.3~"Weak/no",
                                                    T~"Moderate"),
                                          levels = c("Strong", "Moderate", "Weak/no"))) %>%
            ggplot(aes(x = pgp_sub, fill = pgp_sub.class)) +
            geom_histogram(color = "black", alpha = 0.75, breaks = bins2) +
            geom_vline(xintercept = 0.3, linetype = "dashed") +
            geom_vline(xintercept = 0.7, linetype = "dashed") +
            labs(x = "Probability of p-glycoprotein substrate",
                 y = "Number of drugs",
                 fill = "", title = "") +
            scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
            scale_fill_bmj() +
            mytheme +
            theme(legend.position = c(0.99, 0.98),
                  legend.background = element_blank(),
                  legend.justification = c(1, 0.75)),
          nrow = 1, labels = c("(a)", "(b)"),
          label_y = 1)
dev.off()

# No of drugs assessed per scale
table(!is.na(SE.glm.df$ACSBC)) # 451
table(!is.na(SE.glm.df$ORCA)) # 176
table(!is.na(SE.glm.df$ISACS)) # all 937

# Single SE
pred.var <- c("ISACS", 
              "ISACS+pgp_sub",
              "ISACS+BBB", 
              "ISACS+BBB+pgp_sub")

SE.glm <- pmap_dfr(tidyr::crossing(scale = pred.var,
                                   SE_var = c(peripheral.SEs, central.SEs)),
                   function(scale, SE_var) {
                     fit <- glm(reformulate(termlabels = scale, response = SE_var), 
                                data = SE.glm.df,
                                family = binomial)
                     pred_prob <- fitted(fit)
                     actual <- model.frame(fit)[[SE_var]]
                     
                     roc_obj <- roc(response= model.frame(fit)[[SE_var]],
                                    predictor = fitted(fit),
                                    quiet = TRUE)
                     brier_score <- mean((pred_prob - actual)^2)
                     calib_stats <- val.prob(pred_prob, actual, pl = FALSE)
                     
                     pr_obj <- pr.curve(scores.class0 = pred_prob[actual == 1],
                                        scores.class1 = pred_prob[actual == 0],
                                        curve = TRUE)
                     AUPR <- pr_obj$auc.integral
                     
                     lvls <- levels(factor(actual))
                     pred_df <- tibble(
                       !!paste0(".pred_", lvls[1]) := 1 - pred_prob,
                       !!paste0(".pred_", lvls[2]) := pred_prob,
                       .actual = factor(actual, levels = lvls)
                     )
                     
                     broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                       dplyr::filter(term != "(Intercept)") %>%
                       mutate(Model = scale,
                              SE_name = SE_var,
                              AUROC = as.numeric(roc_obj$auc),
                              roc = list(roc_obj),
                              Baseline_AUPRC = mean(actual),
                              AUPRC = pr_obj$auc.integral,
                              pr_curve = list(pr_obj),
                              
                              Brier = brier_score,
                              sig = case_when(p.value < 0.05 ~ "Significant",
                                              T ~ "NS"))
                   }) 

SE.glm.ptrend <- pmap_dfr(tidyr::crossing(scale = pred.var,
                                          SE_var = c(peripheral.SEs, central.SEs)),
                          function(scale, SE_var) {
                            fit <- glm(reformulate(termlabels = scale, response = SE_var),
                                       data = SE.glm.df %>% 
                                         mutate(ISACS = factor(as.character(ISACS),
                                                               levels = c("No", "Low", "Moderate", "High"),
                                                               ordered = T)),
                                       family = binomial)
                            pred_prob <- fitted(fit)
                            actual <- model.frame(fit)[[SE_var]]
                            
                            roc_obj <- roc(response= model.frame(fit)[[SE_var]],
                                           predictor = fitted(fit),
                                           quiet = TRUE)
                            brier_score <- mean((pred_prob - actual)^2)
                            calib_plot <- val.prob(pred_prob, actual, pl = FALSE)
                            
                            broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                              mutate(Model = scale,
                                     SE_name = SE_var,
                                     AUROC = as.numeric(roc_obj$auc),
                                     roc = list(roc_obj),
                                     Brier = brier_score,
                                     calib = list(calib_plot),
                                     sig = case_when(p.value < 0.05~"Significant",
                                                     T~"NS"))
                          }) 
format_p <- function(p) {
  ifelse(
    p < 1e-4,
    formatC(p, format = "e", digits = 2),
    formatC(p, format = "f", digits = 3)
  )
}

SE.glm %>%
  distinct(Model, SE_name, AUROC, AUPRC, Baseline_AUPRC) %>%
  mutate(Model = factor(Model, levels = pred.var)) %>%
  pivot_longer(cols = c(AUROC, AUPRC), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = Model, y = value, fill = Model)) +
  geom_col(position = "dodge") +
  geom_hline(data = . %>% distinct(SE_name, Baseline_AUPRC) %>% 
               mutate(metric = "AUPRC"),
             aes(yintercept = Baseline_AUPRC), 
             linetype = "dashed", color = "black") +
  facet_grid(metric~SE_name) +
  scale_fill_bmj() +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(y = "Score", title = "AUC-ROC vs AUC-PR by model and side effect",
       x = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

SE.glm.ptrend <- SE.glm.ptrend %>% dplyr::filter(term == "ISACS.L") %>%
  dplyr::select(Model, SE_name, p.value) %>%
  add_significance("p.value") %>%
  mutate(p.value.signif = gsub("ns", "",p.value.signif),
         p.trend = paste0(format_p(p.value), p.value.signif))

# Supplementary Table 10
write_clip(SE.glm[,c(1:10, 13)] %>%
             mutate(term = gsub("ISACS", "", term)) %>%
             mutate(SE.class = case_when(SE_name %in% central.SEs~"Central effect",
                                         SE_name %in% peripheral.SEs~"Peripheral effect"),
                    nDrugs = nrow(SE.glm.df)) %>%
             left_join(., SE.glm.ptrend, c("SE_name", "Model")) %>%
             dplyr::select(Model, SE_name, SE.class, nDrugs, everything()))

SIDER.forest <- SE.glm %>% 
  dplyr::filter(Model == "ISACS") %>%
  mutate(term = gsub("ISACS", "", term)) %>% 
  dplyr::filter(term %in% levels(SE.glm.df$ISACS)) %>% 
  mutate(term = factor(term,levels = c("High", "Moderate", "Low")),
         sig = factor(sig, levels = c("Significant", "NS")),
         SE.class = case_when(SE_name %in% central.SEs~"Central effect",
                              SE_name %in% peripheral.SEs~"Peripheral effect")) %>%
  left_join(., SE.glm.ptrend, c("SE_name", "Model")) %>%
  ggplot(aes(x = estimate, y = SE_name,
             color = term, shape = sig,
             xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(size = 0.65, position = position_dodge(width= 0.75)) +
  geom_text(aes(label = p.trend, x = 70), 
            size = 3.25,
            color = "black") +
  geom_stripped_rows(color = NA) +
  scale_shape_manual(values = c(17,19)) +
  scale_y_discrete(limits=rev) +
  scale_x_log10(guide = "axis_logticks",
                expand = expansion(mult = c(0.05, 0.2))) +
  scale_color_manual(values = rev(pal_jama("default")(4)[-1])) +
  labs(x = "Odds ratio (log scale)",y = "",
       title = "Single anticholinergic effect analysis", 
       subtitle = paste0("Analyzed ", nrow(SE.glm.df), " drugs"),
       color = "Anticholinergic level", shape = "Significance") +
  facet_wrap(SE.class~., scales = "free", ncol = 2) + mytheme
SIDER.forest

# tiff("Results/SIDER_forest.tiff", units="in", height=8, width=10.5, res=300)
SIDER.forest
# dev.off()

single.AE.auc <- SE.glm %>% 
  mutate(Model = factor(Model, levels = pred.var),
         SE.class = case_when(SE_name %in% central.SEs~"CNS",
                              SE_name %in% peripheral.SEs~"PNS")) %>%
  ggplot(aes(x = Model, y = AUROC, fill = SE.class)) +
  geom_line(aes(group = SE_name, color = SE.class), 
            linewidth = 1, alpha = 0.9) + 
  geom_point(size = 3, pch = 21) + 
  geom_boxplot(outliers = F, alpha = 0.5, width = 0.5) +
  labs(#title = "Single anticholinergic effect analysis (AUC)",
    #subtitle = "Incorporating PK parameters", 
    x = "", y = 'AUC-ROC') +
  facet_wrap(~SE.class) +
  scale_fill_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  scale_color_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  mytheme + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

deltaAUC.df <- SE.glm %>% distinct(Model, SE_name, roc, AUROC) %>%
  group_by(SE_name) %>%
  mutate(base_auc = AUROC[Model == "ISACS"]) %>%
  group_modify(~{
    roc_base <- .x$roc[.x$Model == "ISACS"][[1]]
    .x %>%
      dplyr::filter(Model != "ISACS") %>%
      rowwise() %>%
      mutate(delta_auc = AUROC - base_auc,
             p_value = pROC::roc.test(roc_base, roc, method = "delong", paired = T)$p.value) %>%
      ungroup()
  }) %>% dplyr::select(-roc) %>%
  add_significance("p_value") %>%
  mutate(SE.class = case_when(SE_name %in% central.SEs~"CNS",
                              SE_name %in% peripheral.SEs~"PNS"))

deltaAUC.df %>%
  distinct(Model, SE.class, SE_name, delta_auc)%>%
  group_by(Model, SE.class) %>%
  dplyr::summarize(AUC.diff = median(delta_auc, na.rm = T))

colFn <- colorRampPalette(c("#002c58", "#2166AC", "white", "#B2182B", "#580000"), interpolate = 'spline')
fill <- colFn(100)[findInterval(deltaAUC.df$delta_auc, seq(-0.15, 0.15, length = 100))]

# tiff("Results/Delta.AUC.tiff", units="in", height=5.75, width=8, res=300)
deltaAUC.p <- deltaAUC.df %>% 
  mutate(Model = factor(Model, levels = pred.var),
         p_value.signif = case_when(p_value.signif == "ns"~"",
                                    T~p_value.signif)) %>%
  ggplot(aes(x = Model, y = SE_name)) +
  geom_point(aes(size = delta_auc), fill = fill,
             pch = 21, alpha = 0.75) +
  geom_text(aes(label = p_value.signif),
            show.legend = F, size = 5,
            vjust =0.7) +
  labs(x = "", y = "",
       size = "ΔAUC") +
  scale_y_discrete(limits=rev) +
  scale_size_continuous(range = c(2, 12)) +
  # scale_fill_gradientn(colours = colorRampPalette(c('#709AE1', '#FFFFFF', '#FD7446'))(100)) +
  guides(size = guide_legend(override.aes = list(color = "black", pch = 16))) +
  ggh4x::facet_wrap2(SE.class~., scales = "free", ncol = 2,
                     strip = strip_themed(background_x = elem_list_rect(fill = c("#9c1806", "#415494")))) + 
  mytheme +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1), 
        strip.text = element_text(color = "white", face = "bold"),
        strip.background = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "white", linewidth = 1))
# dev.off()

# Aggregated
pred.SIDER <- c("ISACS", 
                "ISACS+pgp_sub",
                "ISACS+BBB", 
                "ISACS+BBB+pgp_sub")
outcome.SIDER <- "anticholinergic.SE" 
n.cSE <- 1:length(central.SEs)
n.pSE <- 1:length(peripheral.SEs)

cutoff_grid <- expand.grid(n.cSE = n.cSE, n.pSE = n.pSE)

compare.scale.res <- cutoff_grid %>%
  mutate(result = map2(n.cSE, n.pSE, ~{
    SE.glm.df %>% 
      mutate(ISACS = factor(as.character(ISACS),
                            levels = c("No", "Low", "Moderate", "High"),
                            ordered = T)
      ) %>%
      mutate(n.central.SEs = rowSums(across(c(central.SEs)), na.rm = T),
             n.peripheral.SEs = rowSums(across(c(peripheral.SEs)), na.rm = T),
             anticholinergic.SE = case_when(n.central.SEs >= .x & n.peripheral.SEs >= .y ~ 1, T~0)) %>%
      dplyr::select(c(Drug, ISACS, BBB, pgp_sub, anticholinergic.SE)) %>% 
      group_modify(~{
        map_dfr(pred.SIDER, function(mod){
          formula_str <- as.formula(paste("anticholinergic.SE ~", mod))
          fit <- glm(formula_str,
                     data = .x,
                     family = binomial)
          
          roc_obj <- pROC::roc(response= .x$anticholinergic.SE,
                               predictor = fitted(fit),
                               quiet = T)
          pred_prob <- fitted(fit)
          actual <- model.frame(fit)[["anticholinergic.SE"]]
          
          brier_score <- mean((pred_prob - actual)^2)
          
          pr_obj <- pr.curve(scores.class0 = pred_prob[actual == 1],
                             scores.class1 = pred_prob[actual == 0],
                             curve = TRUE)
          AUPR <- pr_obj$auc.integral
          
          coords_df <- pROC::coords(roc_obj, "all", ret = "all", transpose = F)
          coords_df$AUROC <- pROC::auc(roc_obj)[1]
          coords_df$Baseline_AUPRC <- mean(actual)
          coords_df$AUPRC <- pr_obj$auc.integral
          coords_df$Brier <- brier_score
          coords_df$model <- mod
          as_tibble(coords_df)
        })
      }) %>%
      mutate(n.cSE.cutoff = .x, n.pSE.cutoff = .y)
  })) %>%
  unnest(result)

# Supplementary table 11
write_clip(compare.scale.res %>%
             distinct(model, n.cSE.cutoff, n.pSE.cutoff,
                      AUROC, AUPRC) %>%
             mutate(nDrugs = nrow(SE.glm.df)))

compare.scale.res.count <- cutoff_grid %>%
  mutate(counts = map2(n.cSE, n.pSE, ~{
    SE.glm.df %>%
      mutate(n.central.SEs = rowSums(across(c(central.SEs)), na.rm = T),
             n.peripheral.SEs = rowSums(across(c(peripheral.SEs)), na.rm = T),
             anticholinergic.SE = case_when(n.central.SEs >= .x & n.peripheral.SEs >= .y ~ 1, T ~ 0)) %>%
      count(anticholinergic.SE) %>%
      pivot_wider(names_from = anticholinergic.SE, 
                  values_from = n, 
                  names_prefix = "class_")}),
    n.cSE.cutoff = n.cSE,
    n.pSE.cutoff = n.pSE
  ) %>%
  unnest(counts) %>%
  dplyr::select(n.cSE.cutoff, n.pSE.cutoff, class_0, class_1) %>%
  mutate(perc = class_1*100/(class_0+class_1))

scales.p3 <- compare.scale.res.count %>% 
  ggplot(aes(x = n.cSE.cutoff, y = n.pSE.cutoff, fill = perc)) +
  geom_point(pch = 21, color = "black", alpha = 0.9, size = 10) + 
  geom_text(aes(label = class_1),
            color = ifelse(compare.scale.res.count$perc < 10, "black", "white")) +
  scale_fill_gradientn(colors = c("white", "#bef4fa", "#547478", "#023338"),
                       values = scales::rescale(c(0, 5, 10, 20), to = c(0, 1), from = c(0, 20)),
                       limits = c(0, 20),
                       oob = scales::squish, # clamps anything >20 to the max color
                       breaks = c(0, 5, 10, 20),
                       labels = c("0", "5", "10", ">20")) +
  scale_x_continuous(limits = c(min(n.cSE), max(n.cSE)),
                     breaks = seq(min(n.cSE), max(n.cSE), length.out = round(length(n.cSE)/2))) +
  scale_y_continuous(limits = c(min(n.pSE), max(n.pSE)),
                     breaks = seq(min(n.pSE), max(n.pSE), length.out = ceiling(length(n.pSE)/2))) +
  mytheme +
  labs(title = " ",
       subtitle = "Proportion and number of drugs", fill = "Percentage (%)",
       y = "Number of peripheral\nanticholinergic effects considered", 
       x = "Number of central\nanticholinergic effects considered") +
  theme(panel.grid.major = element_line(color = "white", linewidth = 1),
        panel.grid.minor = element_line(color = "white", linewidth = 1),
        panel.spacing.x = unit(1.1, "lines"))

col_fun = circlize::colorRamp2(c(0.4, 0.525, 0.65, 0.775, 0.9), 
                               c("#d6a9ab", "#b57d7f", "#824446", "#6e2b2d", "#4a1712")) 
col_palette <- col_fun(seq(0, 1, length.out = 256))
scales.p1 <- compare.scale.res %>% 
  distinct(n.cSE, n.pSE, AUROC, model) %>% 
  mutate(model = factor(model, levels = pred.var)) %>%
  ggplot(aes(x = n.cSE, y = n.pSE, size = AUROC, fill = AUROC)) +
  geom_point(pch = 21, color = "black", alpha = 0.9) + 
  facet_wrap(~model) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_gradientn(colors = col_palette) +
  scale_x_continuous(limits = c(min(n.cSE), max(n.cSE)),
                     breaks = seq(min(n.cSE), max(n.cSE), length.out = round(length(n.cSE)/2))) +
  scale_y_continuous(limits = c(min(n.pSE), max(n.pSE)),
                     breaks = seq(min(n.pSE), max(n.pSE), length.out = ceiling(length(n.pSE)/2))) +
  mytheme +
  labs(title = "Aggregated anticholinergic effect analysis",
       subtitle = "Performance", size = "AUC-ROC", fill = "AUC-ROC",
       y = "Number of peripheral\nanticholinergic effects considered", 
       x = "Number of central\nanticholinergic effects considered") +
  theme(panel.grid.major = element_line(color = "white", linewidth = 1),
        panel.grid.minor = element_line(color = "white", linewidth = 1),
        panel.spacing.x = unit(1.1, "lines"))

scales.p2 <- compare.scale.res %>% 
  mutate(model = factor(model, levels = pred.var)) %>%
  distinct(n.cSE, n.pSE, AUROC, model) %>% 
  ggplot(aes(x = model, y = AUROC)) +
  geom_violin() +
  geom_boxplot(outliers = T, width = 0.15, fill = "grey70") +
  labs(title = " ",
       subtitle = " ",
       x = "", y = 'AUC-ROC') +
  
  mytheme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

aggregated.result <- plot_grid(scales.p3, scales.p1, scales.p2, nrow = 1, 
                               rel_widths = c(1.05,1.25,0.65),
                               align = "h", axis = "tb")

# Cumulative
SE.glm.df.cumsum <- SE.glm.df %>%
  mutate(n.central.SEs = rowSums(across(c(central.SEs)), na.rm = T),
         n.peripheral.SEs = rowSums(across(c(peripheral.SEs)), na.rm = T),
         n.all.SEs = n.central.SEs+n.peripheral.SEs) 

SE.glm.cumsum <- pmap_dfr(tidyr::crossing(scale = pred.var,
                                          SE_var = c("n.central.SEs", "n.peripheral.SEs", "n.all.SEs")),
                          function(scale, SE_var) {
                            fit <- MASS::glm.nb(reformulate(termlabels = scale, response = SE_var), 
                                                data = SE.glm.df.cumsum)
                            broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                              dplyr::filter(term != "(Intercept)") %>%
                              mutate(SE_name = SE_var,
                                     model = scale,
                                     sig = case_when(p.value < 0.05~"Significant",
                                                     T~"NS"))
                          }) %>%
  mutate(term = gsub("ISACS", "", term)) %>%
  mutate(term = factor(term, levels = c("High", "Moderate", "Low", "BBB", "pgp_sub")),
         SE_name = case_when(SE_name == "n.all.SEs"~"Number of all AEs",
                             SE_name == "n.central.SEs"~"Number of CNS AEs",
                             SE_name == "n.peripheral.SEs"~"Number of PNS AEs"),
         model = factor(model, levels = c(pred.var)))

SE.glm.cumsum.ptrend <- pmap_dfr(tidyr::crossing(scale = pred.var,
                                                 SE_var = c("n.central.SEs", "n.peripheral.SEs", "n.all.SEs")),
                                 function(scale, SE_var) {
                                   fit <- MASS::glm.nb(reformulate(termlabels = scale, response = SE_var), 
                                                       data = SE.glm.df.cumsum %>% 
                                                         mutate(ISACS = factor(as.character(cholinergicCluster),
                                                                               levels = c("No", "Low", "Moderate", "High"),
                                                                               ordered = T))
                                   )
                                   
                                   broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                                     dplyr::filter(term != "(Intercept)") %>%
                                     mutate(SE_name = SE_var,
                                            model = scale,
                                            sig = case_when(p.value < 0.05~"Significant",
                                                            T~"NS"))
                                 }) %>%
  dplyr::filter(term %in% c("ISACS.L", "ORCA.L", "ACSBC.L")) %>% 
  mutate(SE_name = case_when(SE_name == "n.all.SEs"~"Number of all AEs",
                             SE_name == "n.central.SEs"~"Number of CNS AEs",
                             SE_name == "n.peripheral.SEs"~"Number of PNS AEs"),
         model = factor(model, levels = pred.var)) %>%
  dplyr::select(SE_name, model, p.value) %>%
  add_significance("p.value") %>%
  mutate(p.value.signif = gsub("ns", "",p.value.signif),
         p.trend = paste0(format_p(p.value), p.value.signif))

SE.cum.plot1 <- SE.glm.df.cumsum %>% 
  pivot_longer(cols = starts_with("n.")) %>%
  mutate(name = case_when(name == "n.all.SEs"~"Number of all AEs",
                          name == "n.central.SEs"~"Number of CNS AEs",
                          name == "n.peripheral.SEs"~"Number of PNS AEs")) %>%
  pivot_longer(cols = "ISACS",
               values_to = "Level", names_to = "Scale") %>% drop_na(Level) %>% 
  mutate(Level.num = case_when(Level == "No"~1,
                               Level == "Low"~2,
                               Level == "Moderate"~3,
                               Level == "High"~4,
  )) %>%
  ggplot(aes(x = Level, y = value, fill = Level)) +
  geom_point(size = 2, alpha = 0.75, 
             position = position_jitter(width = 0.15)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(outliers = F, alpha = 0.85, width = 0.2) +
  geom_smooth(method = "lm", se = F, inherit.aes = F,
              aes(x = Level.num, y = value), 
              color = "darkgreen") +
  scale_fill_jama() +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  facet_grid(.~name) +
  labs(x = "", y = "Count") +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1),
                  legend.position = "none")

SE.cum.plot2 <- SE.glm.cumsum %>%
  dplyr::filter(term %in% c("Low", "Moderate", "High")) %>%
  left_join(., SE.glm.cumsum.ptrend, c("model", "SE_name")) %>%
  mutate(sig = factor(sig, levels = c("Significant", "NS"))) %>%
  ggplot(aes(x = estimate, y = model, color = term, shape = sig,
             xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(size = 0.65, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = p.trend, x = 3.15), 
            size = 3.25,
            color = "black") +
  geom_stripped_rows(color = NA) +
  scale_color_manual(values = rev(pal_jama("default")(4)[-1])) +
  scale_y_discrete(limits=rev) +
  scale_shape_manual(values = c(17,19)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_grid(~SE_name) +
  labs(title = "Cumulative anticholinergic effect analysis",
       color = "Anticholinergic level", y = "", x = "Estimate", shape = "Significance") +
  mytheme

cum.result <- plot_grid(SE.cum.plot1, SE.cum.plot2,
                        nrow = 1, align = "h", 
                        rel_widths = c(0.8, 1.1))

# Supplementary Table 12
write_clip(SE.glm.cumsum[,c(1:9)] %>%
             mutate(term = gsub("ISACS", "", term)) %>%
             mutate(nDrugs = nrow(SE.glm.df)) %>%
             left_join(., SE.glm.cumsum.ptrend, c("SE_name", "model")) %>%
             dplyr::select(model, SE_name, nDrugs, everything()))

tiff("Results/Figure3.tiff", units="in", height=12.5, width=13, res=300)
plot_grid(plot_grid(NULL, SIDER.forest, NULL, rel_widths = c(0.05, 1, 0.05), nrow = 1),
          plot_grid(single.AE.auc, deltaAUC.p, nrow = 1, rel_widths = c(0.5, 1),
                    align = "h", axis = "tb"),
          ncol = 1, rel_heights = c(0.6, 0.4),
          labels = c("(a)", "(b)"))
dev.off()

tiff("Results/SFigure12.tiff", units="in", height=10, width=17, res=300)
plot_grid(aggregated.result,
          cum.result,
          rel_heights = c(1.3, 1),
          labels = c("(a)", "(b)"), 
          ncol = 1)
dev.off()

# Drug tissue score across ISACS
drug_tissue_score_plot.df <- drug_tissue_score %>% as.data.frame %>%
  rownames_to_column("Drug") %>% 
  pivot_longer(-Drug) %>% 
  inner_join(., result.w.anno.full %>% distinct(Drug, cholinergicCluster) %>% drop_na(),
             "Drug") %>%
  mutate(Truth_simple = case_when(cholinergicCluster == "No"~"Group1",
                                  T~"Group2"),
         Tissue = case_when(name %in% colnames(CHRM.exp.brain)~"Brain tissue",
                            T~"Other organs"))

region_order <- drug_tissue_score_plot.df %>%
  group_by(name, Truth_simple) %>%
  dplyr::summarise(mean_val = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Truth_simple, values_from = mean_val) %>%
  mutate(diff = abs(Group1 - Group2)) %>%
  arrange(diff)

drug_tissue_score_plot.2 <- drug_tissue_score_plot.df %>%
  mutate(name = factor(name, levels = region_order$name)) %>%
  ggplot(aes(x = value, y = fct_rev(name), fill = cholinergicCluster)) +
  geom_density_ridges(alpha = 0.75) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_jama() +
  labs(x = "Tissue-specific drug exposure score", y = "", fill = "") +
  # facet_grid(~Tissue, scales = "free_y") +
  theme_ridges(grid = F, center_axis_labels = T) +
  theme(legend.position = "right",
        plot.title = element_text(face="bold",size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",size=13, hjust = 0.5),
        strip.text = element_text(size = 11, color = "black"),
        panel.border = element_rect(fill = NA), 
        strip.background = element_rect(fill = NA))
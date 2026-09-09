# Decoding Anticholinergic Effects
## Analysis on only network results
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
  mutate(QSAR_score = mean(c_across(all_of(Outcomes)), na.rm = T)) 

Fig1C <- result.w.anno.full %>%
  left_join(., ATC.df %>% distinct(name_1st, drugbank_id) %>% mutate(name_1st = str_to_sentence(name_1st)),
            "drugbank_id") %>%
  mutate(lab = case_when(!is.na(Ground_truth)~Drug,
                         T~NA),
         name_1st = case_when(!is.na(name_1st)~name_1st,
                              T~"No ATC code"),
         name_1st = factor(name_1st,
                           levels = lvl_order),
         zscore_sig = factor(zscore_sig, levels = c("Significant", "NS"))) %>%
  ggplot(aes(x = distance.value, y = zscore, 
             shape = zscore_sig, color = name_1st)) +
  geom_point(data = ~ dplyr::filter(.x, name_1st == "No ATC code"),
             alpha = 0.8, size = 3) +
  geom_point(data = ~ dplyr::filter(.x, name_1st != "No ATC code"),
             alpha = 0.9, size = 3) +
  geom_hline(yintercept = -1.96, linetype = "dashed") +
  scale_shape_manual(values = c(17, 16)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.04))) +
  scale_alpha_continuous(range = c(0.3, 0.9)) +
  guides(fill = guide_legend(override.aes = aes(size = 3.5, pch = 21, label = "")),
         shape = guide_legend(override.aes = aes(size = 3.5, fill = "grey50"))) +
  ATC.col +
  labs(x = "Network distance", y = "Network z-score",
       color = "ATC level 1", shape = "Statistical significance") +
  mytheme

tiff("Results/SFigure2.tiff", units="in", height=13, width=12.5, res=300)
plot_grid(result.w.anno.full %>%
            left_join(., ATC.df %>% distinct(name_1st, drugbank_id) %>% mutate(name_1st = str_to_sentence(name_1st)),
                      "drugbank_id") %>% 
            dplyr::filter(zscore_sig == "Significant") %>% 
            mutate(name_1st = case_when(!is.na(name_1st)~name_1st,
                                        T~"No ATC code")) %>%
            ggplot(aes(x = forcats::fct_infreq(name_1st), fill = name_1st)) + 
            geom_bar(stat = "count", width = 0.5) + 
            geom_text(stat = "count",
                      aes(label = after_stat(count)),
                      vjust = -0.3) +
            scale_y_continuous(expand = expansion(mult = c(0.0, 0.15))) +
            labs(x = "", y = "Number of drugs",
                 title = "ATC level 1 of topologically-associated drugs") +
            ATC.fill + mytheme +
            theme(axis.text.x = element_text(angle = 25, hjust = 1),
                  legend.position = "none"),
          result.w.anno.full %>%
            left_join(., ATC.df %>% distinct(name_1st, drugbank_id) %>% mutate(name_1st = str_to_sentence(name_1st)),
                      "drugbank_id") %>%
            mutate(lab = case_when(!is.na(Ground_truth)~Drug,
                                   T~NA),
                   name_1st = case_when(!is.na(name_1st)~name_1st,
                                        T~"No ATC code"),
                   name_1st = factor(name_1st %>% gsub(" and ", " and\n" ,.) %>% gsub(", excl.", ",\nexcl." ,.) %>% 
                                       gsub("hormones and\n", "hormones and " ,.) %>% 
                                       gsub("insecticides and\n", "\ninsecticides and " ,.),
                                     levels = lvl_order %>% gsub(" and ", " and\n" ,.) %>% gsub(", excl.", ",\nexcl." ,.) %>% 
                                       gsub("hormones and\n", "hormones and " ,.) %>% 
                                       gsub("insecticides and\n", "\ninsecticides and " ,.)),
                   zscore_sig = factor(zscore_sig, levels = c("Significant", "NS"))) %>%
            ggplot(aes(x = distance.value, y = zscore, 
                       shape = zscore_sig, color = name_1st)) +  
            geom_point(alpha = 0.9, size = 2.5) +
            scale_shape_manual(values = c(17,16)) + #  c(24, 21)
            scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.04))) +
            guides(color = "none",
                   shape = guide_legend(override.aes = aes(size = 3.5, color = "grey50"))) +
            scale_color_manual(values = c("#e8ce4d","#8a7c74","#D2AF81FF",
                                         "#3b2929","#ed5565","#FD7446",
                                         "#ed87c1","#008280","#4dc0e8",
                                         "#ac92ec","#008B45FF","#5de4d2",
                                         "#a0cecb","#c5de5f","#99958b")) +
            facet_wrap(~name_1st, ncol = 5) + 
            gghighlight::gghighlight(use_direct_label = F) +
            geom_hline(yintercept = -1.96, linetype = "dashed") +
            labs(x = "Network distance", y = "Network z-score",
                 shape = "Statistical significance") +
            mytheme + theme(legend.position = "bottom"),
          ncol = 1, rel_heights = c(0.6,1),
          labels = c("(a)", "(b)"))
dev.off()

# ATC enrichment
ORA.result  <- enricher(result.w.anno.full %>% dplyr::filter(zscore_sig %in% "Significant") %>% pull(Drug), 
                        minGSSize = 5, TERM2GENE = ATC.df.enrich)

Fig1D <- ORA.result@result %>%
  left_join(., ATC.df %>% mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")"),
                                 name_1st = str_to_sentence(name_1st)) %>%
              distinct(name_1st, ATC_3rd),
            c("ID" = "ATC_3rd")) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(Count >= 5) %>%
  mutate(log.p = -log10(p.adjust)) %>%
  ggdotchart(., x = "Description", y = "Count",
             add = "segments", 
             add.params = list(color = "name_1st", size = 1.25),
             color = "name_1st",   
             sorting = "descending", rotate = T,
             dot.size = "log.p",
             ggtheme = theme_pubr()) +
  scale_size(range = c(3,6.5)) +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0., 0.05))) +
  guides(color = "none") +
  ATC.col +
  labs(x = "", y = "Number of drugs",
       title = "ATC Level 3 Enrichment Analysis of\nTopologically Significant Drugs",
       size = "-log10(FDR-adjusted p-value)") +
  theme(legend.position = "bottom",
        plot.title = element_text(face="bold",  size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        strip.text = element_text(size = 11, color = "black"),
        axis.ticks.length=unit(.15, "cm")) 

## Compare with Youden index threshold
## This corresponds to 0.09 p-value (one-sided)
message("When using Youden index z < -1.302, there were: ", 
        length(result.w.anno.full %>% dplyr::filter(zscore < -1.302) %>% 
                 pull(Drug)), " significant topologically-associated drugs")

ORA.result.Youden  <- enricher(result.w.anno.full %>% dplyr::filter(zscore < -1.302) %>% 
                                 pull(Drug), 
                               minGSSize = 5, TERM2GENE = ATC.df.enrich)
Youden.p <- rbind(ORA.result@result %>% mutate(Sig = "Standard cut-off\n(z < -1.96)"),
                  ORA.result.Youden@result %>% mutate(Sig = "Youden index\n(z < -1.302)")) %>%
  left_join(., ATC.df %>% mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")"),
                                 name_1st = str_to_sentence(name_1st)) %>%
              distinct(name_1st, ATC_3rd),
            c("ID" = "ATC_3rd")) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(Count >= 5) %>%
  mutate(log.p = -log10(p.adjust)) %>%
  ggdotchart(., x = "Description", y = "Count",
             add = "segments", 
             add.params = list(color = "name_1st", size = 1.25),
             color = "name_1st", group = "name_1st",  
             facet.by = "Sig",
             sort.by.groups = T,
             sorting = "descending", rotate = T,
             dot.size = "log.p",
             ggtheme = theme_pubr()) +
  scale_size(range = c(3,6.5)) +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0., 0.1))) +
  # guides(color = "none") +
  ATC.col + 
  labs(x = "", y = "Number of drugs", color = "ATC level 1",
       title = "ATC Level 3 Enrichment Analysis of\nTopologically Significant Drugs",
       size = "-log10(FDR-adjusted p-value)") +
  theme(legend.position = "right",
        panel.background = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face="bold",  size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        strip.text = element_text(size = 11, color = "black"),
        axis.ticks.length=unit(.15, "cm")) 
tiff("Results/SFigure3.tiff", units="in", height=6.5, width=14.25, res=300)
Youden.p
dev.off()

result.w.anno.full %>% drop_na(Concordance.ACSBC.ORCA) %>%
  ggplot(aes(x = Concordance.ACSBC.ORCA, y = zscore*distance.value, 
             fill = Concordance.ACSBC.ORCA)) + 
  geom_violin(alpha = 0.25) + geom_boxplot(width = 0.15) +
  ggpubr::stat_compare_means(comparisons = combn(levels(result.w.anno.full$Concordance.ACSBC.ORCA), 
                                                 2, simplify = F), method = "t.test") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_fill_jama() + labs(x = "", y = "Network distance × zscore") + 
  mytheme + theme(legend.position = "none",
                  axis.text.x = element_text(angle = 25, hjust = 1))

## GSEA just for fun
drugList <- result.w.anno.full$distance.value * result.w.anno.full$zscore
names(drugList) <- result.w.anno.full$Drug
drugList <- sort(drugList, decreasing = T)
DSEA.result <- GSEA(drugList, minGSSize = 5,
                    TERM2GENE = ATC.df.enrich)

DSEA.result@result %>%
  left_join(., ATC.df %>% mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")"),
                                 name_1st = str_to_sentence(name_1st)) %>%
              distinct(name_1st, ATC_3rd),
            c("ID" = "ATC_3rd")) %>%
  dplyr::filter(NES < 0) %>%
  dplyr::filter(p.adjust < 0.001) %>%
  mutate(log.p = -log10(p.adjust)) %>%
  ggdotchart(., x = "Description", y = "NES",
             add = "segments", 
             add.params = list(color = "name_1st", size = 1.25),
             color = "name_1st",   
             sorting = "descending", rotate = T,
             dot.size = "log.p",
             ggtheme = theme_pubr()) +
  scale_size(range = c(3,6.5)) +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0., 0.05))) +
  guides(color = "none") +
  ATC.col +
  labs(x = "", y = "NES",
       title = "ATC Level 3 Set Enrichment Analysis of\nTopologically Significant Drugs",
       size = "-log10(FDR-adjusted p-value)") +
  theme(legend.position = "bottom",
        plot.title = element_text(face="bold",  size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        strip.text = element_text(size = 11, color = "black"),
        axis.ticks.length=unit(.15, "cm")) 

write.csv(result.w.anno.full, "result.w.anno.full.csv", row.names = F)
colSums(!is.na(result.w.anno.full))

### ASSESSING GROUND TRUTH ###
# tiff(paste0("Results/1_GroundTruth_final.tiff"), units="in", height=4, width=5.75, res=300)
KW.test <- result.w.anno.full %>% 
  drop_na(Ground_truth) %>% 
  mutate(val.temp = distance.value*zscore) %>%
  kruskal_test(val.temp ~ Ground_truth)    

dunn.test <- result.w.anno.full %>% 
  drop_na(Ground_truth) %>% 
  mutate(val.temp = distance.value*zscore) %>%
  dunn_test(val.temp ~ Ground_truth, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Ground_truth") %>%
  mutate(y.position = y.position + 2) %>% dplyr::select(-groups)

Fig1B <- result.w.anno.full %>% 
  mutate(zscore_sig = factor(zscore_sig, levels = c("Significant", "NS"))) %>%
  drop_na(Ground_truth) %>% 
  ggplot(aes(x = Ground_truth, y = distance.value*zscore)) +
  geom_violin() +
  geom_boxplot(outliers = F, alpha = 0.3, width = 0.15, fill = "grey70") +
  geom_point(size = 3.5, alpha = 0.65, color = "white",
             position = position_jitterdodge(jitter.width = 0.4,
                                             dodge.width = 0.5),
             aes(fill = Drug_category, shape = zscore_sig)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_shape_manual(values = c(24,21)) +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.08) +
  labs(x = "", y = "Network distance × z-score", shape = "Statistical significance",
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
         shape = "none" # guide_legend(override.aes = list(color = "black", pch = c(19,17)))
         ) +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))
# dev.off()

result.w.anno.full %>% 
  drop_na(Ground_truth) %>% 
  group_by(Ground_truth, Drug_category) %>% count() 

### Relationship with SAA ###
v1 <- unique(SAA.df.average$Medication)
v2 <- result.w.anno.full %>% drop_na(SAA_class) %>% 
  dplyr::select(Drug, SAA_class) %>% distinct() %>% pull(Drug)
setdiff(v1, v2) # Primarily just antibiotics

SAA.plot.1 <- result.w.anno.full %>% drop_na(SAA_class) %>%
  mutate(SAA_class.num = as.numeric(SAA_class)-1) %>%
  ggplot(aes(x = SAA_class.num, y = distance.value*zscore)) +
  geom_violin(aes(group = as.factor(SAA_class.num))) +
  geom_point(size = 3.5, alpha = 0.75, color = "white",
             aes(fill = Class, shape = zscore_sig),
             position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(group = as.factor(SAA_class.num)), fill = "grey70",
               alpha = 0.35, outliers = F, width = 0.125) +
  stat_cor(method= "spearman", label.x = 2.25) +
  labs(x = "SAA categories",
       y = "Network distance × z-score") + 
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_continuous(breaks = c(0:4), labels = levels(result.w.anno.full$SAA_class)) +
  scale_fill_manual(values = c("Antidepressant" = "#E76F51", 
                               "Antihistamine" = "#D45D5D",
                               "Antipsychotic" = "#b27ef2", 
                               "Anxiolytic and sedative" = "#4d3f5e",
                               "Gastrointestinal and bowel agent" = "#415494",
                               "Urinary incontinence agent" = "#F4A261",
                               "Others" = "grey20")) + mytheme +
  scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(color = "black", pch = c(21))),
         shape = guide_legend(override.aes = list(color = "black", pch = c(19,17)))) +
  theme(legend.position = "none")

SAA.plot.2 <- result.w.anno.full %>% 
  dplyr::filter(SAA_class %in% c("+", "++", "+++")) %>%
  mutate(zscore_sig = factor(zscore_sig, levels = c("Significant", "NS"))) %>%
  # INTERESTING NOTE!!!!
  dplyr::filter(!Drug %in% c("Lithium chloride", "Lithium carbonate",
                             "Lithium cation", "Lithium succinate", "Lithium citrate")) %>%
  ggplot(aes(x = SAA.per.dose, y = zscore*distance.value)) + 
  scale_x_log10(guide = "axis_logticks") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.16))) +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  ggrepel::geom_text_repel(aes(label = Drug), size = 3.5,
                           max.overlaps = 5) +
  geom_point(size = 3.5, alpha = 0.75, color = "white",
             aes(fill = Class, shape = zscore_sig)) + 
  labs(x = "SAA (pmol/L) normalized by dose (mg)",
       y = "Network distance × z-score", 
       fill = "Drug class",
       shape = "Statistical significance") +
  stat_cor(method= "spearman", 
           label.x = log10(result.w.anno.full$SAA.per.dose[result.w.anno.full$Drug == "Paroxetine"]), 
           label.y = result.w.anno.full$distance.value[result.w.anno.full$Drug == "Temazepam"]*result.w.anno.full$zscore[result.w.anno.full$Drug == "Temazepam"]) +
  scale_shape_manual(values = c(24,21)) +
  guides(fill = guide_legend(override.aes = list(color = "black", pch = c(21))),
         shape = "none" # guide_legend(override.aes = list(color = "black", pch = c(17,19)))
         ) +
  scale_fill_manual(values = c("Antidepressant" = "#E76F51", 
                               "Antihistamine" = "#D45D5D",
                               "Antipsychotic" = "#b27ef2", 
                               "Anxiolytic and sedative" = "#4d3f5e",
                               "Gastrointestinal and bowel agent" = "#415494",
                               "Urinary incontinence agent" = "#F4A261",
                               "Others" = "grey20")) + mytheme

Fig1EF <- plot_grid(SAA.plot.1, SAA.plot.2, 
                    align = "h", axis = "tb",
                    ncol = 2, rel_widths = c(1,1.65),
                    labels = c("(e)", "(f)"))

# tiff("Results/3_SAA_final.tiff", units="in", height=3.75, width=11.5, res=300) # 3_SAA_v1_Li
Fig1EF
# dev.off()

BR.network <- ggdraw() + 
  draw_image("Biorender/NetworkMedicine.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

tiff("Results/Figure1.tiff", units="in", height=21, width=12, res=300)
plot_grid(plot_grid(BR.network, Fig1B, nrow = 1, rel_widths = c(1,0.85),
                    labels = c("(a)", "(b)")),
          Fig1C,
          Fig1D, 
          Fig1EF, ncol = 1, 
          rel_heights = c(0.9, 1.25, 1.25, 0.75),
          labels = c("", "(c)", "(d)", ""))
dev.off()

SAA.plot.2v2 <- result.w.anno.full %>% 
  dplyr::filter(SAA_class %in% c("+", "++", "+++")) %>%
  # INTERESTING NOTE!!!!
  dplyr::filter(!Drug %in% c(#"Thioridazine", "Temazepam",
    "Lithium chloride", "Lithium carbonate",
    "Lithium cation", "Lithium succinate", "Lithium citrate"
  )) %>%
  ggplot(aes(x = SAA.min.ther.dose, y = zscore*distance.value)) +
  scale_x_log10(guide = "axis_logticks") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  geom_smooth(method = "lm", se = F) +
  ggrepel::geom_text_repel(aes(label = Drug), size = 3.5,
                           max.overlaps = 5) +
  geom_point(size = 3.5, alpha = 0.75, color = "white",
             aes(fill = Class, shape = zscore_sig)) +  
  labs(x = "SAA (pmol/L) normalized by minimum therapeutic dose (mg)",
       y = "Network distance × z-score", shape = "Statistical significance") +
  stat_cor(method= "spearman", 
           label.x = log10(result.w.anno.full$SAA.min.ther.dose[result.w.anno.full$Drug == "Tolterodine"]), 
           label.y = result.w.anno.full$distance.value[result.w.anno.full$Drug == "Temazepam"]*result.w.anno.full$zscore[result.w.anno.full$Drug == "Temazepam"]) +
  scale_fill_manual(values = c("Antidepressant" = "#E76F51", 
                               "Antihistamine" = "#D45D5D",
                               "Antipsychotic" = "#b27ef2", 
                               "Anxiolytic and sedative" = "#4d3f5e",
                               "Gastrointestinal and bowel agent" = "#415494",
                               "Urinary incontinence agent" = "#F4A261",
                               "Others" = "grey20")) + 
  scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(color = "black", pch = c(21))),
         shape = guide_legend(override.aes = list(color = "black", pch = c(19,17)))) + 
  mytheme

### Relationship with binding affinity ###
BindingAffinity.plot.1 <- distance.df %>%
  inner_join(., BA.mAChR, c("Drug" = "DRUGNAME")) %>%
  mutate(GENENAME = factor(GENENAME, 
                           levels = c(paste0("CHRM", seq(5,1))))) %>%
  ggplot(aes(x = pKi, y = GENENAME, fill = GENENAME)) +
  geom_density_ridges(jittered_points = T, 
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, 
                      point_alpha = 1, alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = rev(wes_palette("FantasticFox1"))) +
  labs(x = "Binding affinity (pKi)", y = "") +
  theme_ridges(grid = F, center_axis_labels = T) +
  theme(legend.position = "none",
        plot.title = element_text(face="bold",  size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
        strip.text = element_text(size = 11, color = "black"),
        panel.border = element_rect(fill = NA), 
        strip.background = element_rect(fill = NA))

BindingAffinity.plot.2 <- distance.df %>%
  dplyr::select(-starts_with("Closest_ab")) %>%
  inner_join(., BA.mAChR, c("Drug" = "DRUGNAME")) %>%
  dplyr::filter(binding_type %in% "Ki") %>%
  ggplot(aes(x = pKi, y = Closest_b_zscore*Closest_b)) + 
  stat_density_2d(geom = "polygon", #bins = 10,  
                  aes(color = GENENAME), fill = NA) +
  facet_wrap(.~GENENAME, scales = "free", ncol = 1) +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = wes_palette("FantasticFox1")) +
  scale_fill_manual(values = wes_palette("FantasticFox1")) +
  geom_point(size = 3.5, alpha = 0.85, color = "white",
             pch = 21, aes(fill = GENENAME)) + 
  labs(x = "Binding affinity (pKi)",
       y = "Network distance × z-score", fill = "mAChR subtypes") +
  ggpubr::stat_cor(method = "spearman", label.sep = "\n",
                   label.x = 8.75) +
  mytheme + theme(legend.position = "none")

BindingAffinity.ht.df <- distance.df %>%
  inner_join(., BA.mAChR, c("Drug" = "DRUGNAME")) %>%
  mutate(binding_type = gsub("\\.\\.nM\\.", "", binding_type)) %>%
  dplyr::filter(binding_type %in% "Ki") %>%
  mutate(binding_value_uM = log10(binding_value/1000)) %>%
  dplyr::select(GENENAME, Drug, pKi) %>% #binding_value_uM
  pivot_wider(names_from = GENENAME, values_from = pKi) %>% #binding_value_uM
  rowwise() %>%
  mutate(row_avg = mean(c_across(paste0("CHRM", seq(1,5))), na.rm = T)) %>%
  ungroup() %>%
  arrange(desc(row_avg)) %>% 
  column_to_rownames("Drug") %>% dplyr::select(-row_avg) %>% 
  dplyr::select(paste0("CHRM", seq(1,5))) %>% as.matrix
BindingAffinity.ht.anno <- result.w.anno.full %>%
  dplyr::filter(Drug %in% rownames(BindingAffinity.ht.df)) %>% distinct %>%
  column_to_rownames("Drug")

BindingAffinity.ht.anno <- BindingAffinity.ht.anno[rownames(BindingAffinity.ht.df),]

ha <- rowAnnotation("Z-score" = anno_lines(BindingAffinity.ht.anno$zscore),
                    "Distance" = anno_lines(BindingAffinity.ht.anno$distance.value),
                    "Degree" = anno_lines(BindingAffinity.ht.anno$DrugTarget_degree),
                    
                    show_legend = T, show_annotation_name = T, 
                    annotation_name_rot = 45, 
                    annotation_name_gp = gpar(fontsize = 9),
                    na_col = "white", annotation_name_side = "top")
ha2 <- HeatmapAnnotation(df = data.frame(row.names =  colnames(BindingAffinity.ht.df),
                                         Gene = colnames(BindingAffinity.ht.df)),
                         show_legend = F, show_annotation_name = F, 
                         col = list(Gene = setNames(wes_palette("FantasticFox1"), colnames(BindingAffinity.ht.df))))
col_fun = circlize::colorRamp2(c(4, 5, 7, 9, 11, 13),  #log10(c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
                               c("#f0dddd", "#d6a9ab", "#b57d7f", "#824446", "#6e2b2d", "#4a1712")) 
BindingAffinity.plot.3 <- ComplexHeatmap::Heatmap(BindingAffinity.ht.df, col=col_fun,
                                                  border = T, 
                                                  border_gp = gpar(col = "gray20"),
                                                  rect_gp = gpar(col = "gray20"),
                                                  row_names_gp = gpar(fontsize = 9),
                                                  column_names_gp = gpar(fontsize = 9),
                                                  row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
                                                  na_col = "white",
                                                  show_row_names= T, cluster_rows = F,
                                                  show_column_names = T, cluster_columns = F, 
                                                  row_names_side = "left",
                                                  show_heatmap_legend = F,
                                                  right_annotation = ha,
                                                  top_annotation = ha2,
                                                  column_names_rot = 35)
lgd_scale <- Legend(col_fun = col_fun, title = "pKi", 
                    # at = c(log10(c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001))),
                    # labels =c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001),
                    direction = "vertical", title_position = "topleft",
                    legend_height = unit(3, "cm"), legend_width = unit(6, "cm"))

tiff("Results/SFigure4.tiff", units="in", height=16, width=9, res=300)
plot_grid(plot_grid(BindingAffinity.plot.1, BindingAffinity.plot.2,
                    ncol = 1, rel_heights = c(0.5,1.2),
                    labels = c("(a)", "(b)")), 
          grid.grabExpr(ComplexHeatmap::draw(BindingAffinity.plot.3, 
                                             annotation_legend_list = packLegend(lgd_scale, direction = "vertical"),
                                             annotation_legend_side = "right",
                                             padding = unit(c(2, 2, 2, 5), "mm"))),
          labels = c("", "(c)"),
          ncol = 2, rel_widths = c(0.65,1))
dev.off()

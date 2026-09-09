# Decoding Anticholinergic Effects
## Only analyse those drugs with complete network and pKi results
# Benchmarking
if (!exists("result.w.anno.full")) {
  source("6_3_RedefiningAnticholinergics_Classification.R")
}

### Relationship with Anticholinergic Scales ###
anticholinergic.scales <- c("ACSBC", "ADS", "ABS", "ALS","ORCA", 
                            "AEC","Chew's SAA list", "ACB", "AAS", 
                            "ARS", "ABC")
antichol.scale.df <- result.w.anno.full %>% 
  mutate(ISACS = as.numeric(cholinergicCluster)-1,
         ACSBC = as.numeric(fct_rev(ACSBC.Anticholinergic.Activity))-1,
         ORCA = as.numeric(fct_rev(ORCA))-1,
         `Chew's SAA list` = as.numeric(SAA_class)-1,
         AEC = case_when(AEC == "No consensus"~NA,
                         T~as.numeric(AEC))
  )
ndrug.scales <- antichol.scale.df %>% 
  dplyr::select(Drug, all_of(anticholinergic.scales)) %>%
  pivot_longer(cols = all_of(anticholinergic.scales)) %>%
  group_by(name) %>% drop_na() %>%
  count() %>%
  ungroup %>% mutate(label = paste0(name, "\n(n = ", n, ")"))

scale.plot.1 <- antichol.scale.df %>% 
  pivot_longer(cols = all_of(anticholinergic.scales)) %>%
  left_join(., ndrug.scales, "name") %>%
  mutate(name = factor(name, levels = anticholinergic.scales)) %>%
  ggplot(aes(x = value, y = distance.value*zscore)) + 
  geom_violin(alpha = 0.5,
              aes(fill = as.factor(value))) +
  geom_point(size = 2, alpha = 0.2, pch = 21,
             position = position_jitter(width = 0.15),
             aes(fill = as.factor(value))) +
  geom_boxplot(width = 0.2, outliers = F,
               aes(fill = as.factor(value))) +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  scale_fill_jama() +
  stat_cor(method= "spearman", label.x=2, label.y = 15, label.sep = "\n") +
  # stat_poly_eq(use_label(c("R2")), label.x = "right", label.y = "top") +
  labs(x = "Anticholinergic scale", y = "Network distance × z-score",
       title = "Relationship with network result") +
  facet_wrap(~name, nrow = 2,
             labeller = labeller(name = ndrug.scales %>%
                                   distinct(name, label) %>%
                                   tibble::deframe())) + 
  mytheme +
  theme(legend.position = "none")

scale.plot.2 <- antichol.scale.df %>% 
  pivot_longer(cols = all_of(anticholinergic.scales)) %>%
  left_join(., ndrug.scales, "name") %>%
  mutate(name = factor(name, levels = anticholinergic.scales)) %>%
  ggplot(aes(x = value, y = QSAR_score)) + 
  geom_violin(alpha = 0.5,
              aes(fill = as.factor(value))) +
  geom_point(size = 2, alpha = 0.2, pch = 21,
             position = position_jitter(width = 0.15),
             aes(fill = as.factor(value))) +
  geom_boxplot(width = 0.2, outliers = F,
               aes(fill = as.factor(value))) +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  scale_fill_jama() +
  stat_cor(method= "spearman", label.x = 2, label.y = 6, label.sep = "\n") +
  # stat_poly_eq(use_label(c("R2")), label.x = "right", label.y = "bottom") +
  labs(x = "Anticholinergic scale", y = "Mean mAChR PSICHIC score",
       title = "Relationship with PSICHIC") +
  facet_wrap(~name, nrow = 2,
             labeller = labeller(name = ndrug.scales %>%
                                   distinct(name, label) %>%
                                   tibble::deframe())) + 
  mytheme +
  theme(legend.position = "none")


ATS.plot.1 <- result.w.anno.full %>% drop_na(ATS) %>%
  ggplot(aes(x = ATS, y = zscore*distance.value)) +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  stat_cor(method= "spearman", label.x = 3, label.y.npc = "top") +
  # stat_poly_eq(use_label(c("R2")), label.x = "right") +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = Drug), max.overlaps = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(y = "Network distance × z-score",
       title = "Relationship with network result",
       subtitle = "ATS (n = 25)",
       x = "ATS") +
  mytheme +
  theme(plot.subtitle = element_text(face = "plain"))

ATS.plot.2 <- result.w.anno.full %>% drop_na(ATS) %>%
  ggplot(aes(x = ATS, y = QSAR_score)) +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  stat_cor(method= "spearman",  label.y.npc = "top") +
  # stat_poly_eq(use_label(c("R2")), label.x = "left", label.y = "top") +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = Drug), max.overlaps = 4) +
  labs(y = "Mean mAChR PSICHIC score",
       title = "Relationship with PSICHIC",
       subtitle = "ATS (n = 25)",
       x = "ATS") +
  mytheme +
  theme(plot.subtitle = element_text(face = "plain"))

# Compare against all scales (due to small sample size; do cumulative)
scale.dev <- read.csv("5_AnticholinergicScaleTable/ScaleDevelopment.csv") %>%
  column_to_rownames("Scale")

colnames(scale.dev) <- gsub("\\.", " ", colnames(scale.dev))
colnames(scale.dev)[4] <- "Activity assay (SAA)"
colnames(scale.dev)[5] <- "Activity assay (pKi)"
colnames(scale.dev)[6] <- "BBB & PK consideration"

scale.dev[scale.dev == 1] <- "Method used"

scale.dev <- t(scale.dev)
scale.dev.clu <- scale.dev %>% as.data.frame %>%
  mutate(across(everything(), ~ if_else(!is.na(.) & . != "", 1, 0))) %>%
  as.matrix()

row_hclust_scale = hclust(vegan::vegdist(scale.dev.clu, method = "jaccard"),
                          method = "complete")
col_hclust_scale = hclust(vegan::vegdist(t(scale.dev.clu), method = "jaccard"),
                          method = "complete")

scale.dev.ht <- oncoPrint(scale.dev,
                          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(height = unit(1, "cm"))),
                          show_heatmap_legend = T,
                          show_column_names = T, 
                          show_row_names = T,
                          show_pct = F,
                          name = " ",
                          alter_fun = list(
                            `Method used` = function(x, y, w, h)
                              grid.rect(x, y, w, h, gp = gpar(fill = "#355e3b", col = "#355e3b"))),
                          col = c(`Method used` = "#355e3b"),
                          cluster_rows = row_hclust_scale, cluster_columns = col_hclust_scale)

tiff(paste0("Results/SFigure13.tiff"), units="in", height=19, width=12, res=300)
plot_grid(plot_grid(NULL, grid.grabExpr(draw(scale.dev.ht)), NULL,
                    nrow = 1, rel_widths = c(0.25, 1, 0.25)),
          scale.plot.1,
          scale.plot.2,
          plot_grid(NULL, ATS.plot.1, NULL, ATS.plot.2, NULL,
                    rel_widths = c(0.2, 1,0.1, 1,0.2),
                    nrow = 1),
          ncol = 1,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          rel_heights = c(1.1, 1.2, 1.2, 0.8))
dev.off()

# Weighted Kappa concordance for ordinal ratings
## Fleiss-Cohen = quadratic penalisation
## Do concordance only to scales with >200 drugs
library(vcd)
library(forcats)
library(exact2x2)

anticholinergic.scales <- c("ISACS", "ACSBC", "ADS", "ABS", "ALS","ORCA", 
                            "AEC","Chew's SAA list", "AAS") 
# Scales with no 0s "ACB", "ARS", "ABC"
antichol.scale.df <- result.w.anno.full %>% 
  column_to_rownames("Drug") %>%
  mutate(ISACS = as.numeric(cholinergicCluster)-1,
         Confidence = as.numeric(Confidence)-1,
         distanceCluster = as.numeric(distanceCluster)-1,
         affinityCluster = as.numeric(affinityCluster)-1,
         
         ACSBC = as.numeric(fct_rev(ACSBC.Anticholinergic.Activity))-1,
         ORCA = as.numeric(fct_rev(ORCA))-1,
         `Chew's SAA list` = case_when(SAA_class == "0/+"~0,
                                       SAA_class == "0"~0,
                                       T~as.numeric(SAA_class)-2),
         AEC = case_when(AEC == "No consensus"~NA,
                         T~as.numeric(AEC))
  ) %>%
  dplyr::select(all_of(anticholinergic.scales)) %>%
  mutate(across(all_of(anticholinergic.scales),
                ~ factor(if_else(as.numeric(as.character(.)) == 0, "No", "Yes"), 
                         levels = c("No", "Yes"))))
binary_levels <- c("No", "Yes")

make_table_binary <- function(df, row_col, col_col) {
  row_var <- df[[row_col]]
  col_var <- df[[col_col]]
  
  tab <- table(row_var, col_var)
  dimnames(tab) <- list(binary_levels, binary_levels)
  names(dimnames(tab)) <- c(row_col, col_col)
  tab
}

extract_kappa_binary <- function(kappa_obj, label) {
  ci   <- confint(kappa_obj)
  k    <- kappa_obj$Unweighted   
  data.frame(
    comparison = label,
    kappa      = round(k["value"], 4),
    ASE        = round(k["ASE"],   5),
    lower      = round(ci["Unweighted", "lwr"], 4),
    upper      = round(ci["Unweighted", "upr"], 4)
  )
}
extract_mcnemar <- function(tab, label) {
  b <- tab["No", "Yes"] 
  c <- tab["Yes", "No"]
  
  mn <- mcnemar.exact(tab)
  
  data.frame(
    comparison   = label,
    b_only       = b,                 
    c_only       = c,              
    discordant_n = b + c,
    odds_ratio   = round(mn$estimate,  4), 
    p_value      = round(mn$p.value,   4),
    significant  = mn$p.value < 0.05
  )
}

pairs <- combn(seq_along(anticholinergic.scales), 2, simplify = FALSE)

kappa_summary_binary <- do.call(rbind, lapply(pairs, function(idx) {
  a_col <- anticholinergic.scales[idx[1]]
  b_col <- anticholinergic.scales[idx[2]]
  
  pair_df <- antichol.scale.df %>%
    dplyr::select(all_of(c(a_col, b_col))) %>%
    dplyr::filter(!is.na(.[[1]]), !is.na(.[[2]]))
  
  if (nrow(pair_df) < 10) {
    message("Skipping ", a_col, " - ", b_col, ": only ", nrow(pair_df), " complete cases")
    return(NULL)
  }
  
  tab <- make_table_binary(pair_df, row_col = a_col, col_col = b_col)
  
  if (any(rowSums(tab) == 0) | any(colSums(tab) == 0)) {
    message("Skipping ", a_col, " - ", b_col, ": degenerate table")
    return(NULL)
  }
  
  label <- paste(a_col, b_col, sep = " - ")
  extract_kappa_binary(Kappa(tab), label) %>%
    mutate(nDrugs = nrow(pair_df))
}))
rownames(kappa_summary_binary) <- NULL

mcnemar_summary <- do.call(rbind, lapply(pairs, function(idx) {
  a_col <- anticholinergic.scales[idx[1]]
  b_col <- anticholinergic.scales[idx[2]]
  
  pair_df <- antichol.scale.df %>%
    select(all_of(c(a_col, b_col))) %>%
    filter(!is.na(.[[1]]), !is.na(.[[2]]))
  
  if (nrow(pair_df) < 10) {
    message("Skipping ", a_col, " - ", b_col, ": only ", nrow(pair_df), " complete cases")
    return(NULL)
  }
  
  tab <- make_table_binary(pair_df, row_col = a_col, col_col = b_col)
  
  if (any(rowSums(tab) == 0) | any(colSums(tab) == 0)) {
    message("Skipping ", a_col, " - ", b_col, ": degenerate table")
    return(NULL)
  }
  
  label <- paste(a_col, b_col, sep = " - ")
  extract_mcnemar(tab, label)
}))
rownames(mcnemar_summary) <- NULL
concordance_summary <- kappa_summary_binary %>%
  left_join(mcnemar_summary, by = "comparison") %>%
  mutate(kappa_interp = case_when(kappa >= 0.80 ~ "Almost perfect",
                                  kappa >= 0.60 ~ "Substantial",
                                  kappa >= 0.40 ~ "Moderate",
                                  kappa >= 0.20 ~ "Fair",
                                  TRUE          ~ "Slight/Poor"
                                  ))


anticholinergic.scales <- c("ISACS",
                            # "Confidence", "distanceCluster", "affinityCluster",
                            "ACSBC", "ADS", "ABS", "ALS","ORCA", 
                            "AEC","Chew's SAA list", "ACB", "AAS", 
                            "ARS", "ABC", "ATS")
antichol.scale.df <- result.w.anno.full %>% 
  column_to_rownames("Drug") %>%
  mutate(ISACS = as.numeric(cholinergicCluster)-1,
         Confidence = as.numeric(Confidence)-1,
         distanceCluster = as.numeric(distanceCluster)-1,
         affinityCluster = as.numeric(affinityCluster)-1,
         
         ACSBC = as.numeric(fct_rev(ACSBC.Anticholinergic.Activity))-1,
         ORCA = as.numeric(fct_rev(ORCA))-1,
         `Chew's SAA list` = case_when(SAA_class == "0/+"~0,
                                       SAA_class == "0"~0,
                                       T~as.numeric(SAA_class)-2),
         AEC = case_when(AEC == "No consensus"~NA,
                         T~as.numeric(AEC))
  ) %>%
  dplyr::select(all_of(anticholinergic.scales))

cor_mat <- cor(antichol.scale.df, method = "spearman", use = "pairwise.complete.obs")
n_mat <- crossprod(!is.na(antichol.scale.df))

col_cor <- circlize::colorRamp2(c(-1,0,1), 
                                c("#002c58", "white","#580000"))

col_n <- colorRamp2(c(0, 20, 50, 100, 200, 550),
                    c("#f7fcf5", "#c7e9c0", "#74c476", "#31a354", "#006d2c", "#00441b")
)

colour_mat <- matrix(NA, nrow = nrow(cor_mat), ncol = ncol(cor_mat))
colour_mat[upper.tri(colour_mat, diag = TRUE)] <- cor_mat[upper.tri(cor_mat, diag = TRUE)]
colour_mat[lower.tri(colour_mat)]<- n_mat[lower.tri(n_mat)]

cell_fun <- function(j, i, x, y, width, height, fill) {
  if (i == j) {
    grid.rect(x, y, width, height, gp = gpar(fill = "grey20", col = "white"))
  } else if (i < j) {
    val <- cor_mat[i, j]
    grid.rect(x, y, width, height, gp = gpar(fill = col_cor(val), col = "white"))
    grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 9))
  } else {
    val <- n_mat[i, j]
    grid.rect(x, y, width, height, gp = gpar(fill = col_n(val), col = "white"))
    grid.text(sprintf("%d", val), x, y, gp = gpar(fontsize = 9))
  }
}

scale.compare.ht <- Heatmap(matrix= cor_mat,
                            col = col_cor,
                            cell_fun= cell_fun,
                            rect_gp = gpar(type = "none"), 
                            cluster_rows= F,
                            cluster_columns = F,
                            show_heatmap_legend = TRUE,
                            name= "Correlation",
                            column_title= "Spearman correlation (upper) & N (lower)",
                            heatmap_legend_param = list(title = "Correlation")) %>% 
  draw(annotation_legend_list = list(Legend(col_fun = col_n, title = "N (pairwise)")))

## Comparing scales for single SEs
## Compare apple-to-apple or apple-to-banana  
SE.glm.df <- SE.glm.df %>% drop_na(all_of(antichol.scales)) # 164
dim(SE.glm.df)

antichol.scales <- c("ISACS", 
                     "ISACS+pgp_sub",
                     "ISACS+BBB", 
                     "ISACS+BBB+pgp_sub",
                     
                     "ACSBC", 
                     "ACSBC+pgp_sub",
                     "ACSBC+BBB", 
                     "ACSBC+BBB+pgp_sub",
                     
                     "ORCA", 
                     "ORCA+pgp_sub",
                     "ORCA+BBB", 
                     "ORCA+BBB+pgp_sub"
                     )


SE.glm <- pmap_dfr(tidyr::crossing(scale = antichol.scales,
                                   SE_var = c(peripheral.SEs, central.SEs)),
                   function(scale, SE_var) {
                     fit <- glm(reformulate(termlabels = scale, response = SE_var), 
                                data = SE.glm.df,
                                family = binomial)
                     pred_prob <- fitted(fit)
                     actual <- model.frame(fit)[[SE_var]]
                     
                     roc_obj <- roc(response = model.frame(fit)[[SE_var]],
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
                              # Emax = calib_stats[["Emax"]],
                              # calib = list(pred_df),
                              sig = case_when(p.value < 0.05 ~ "Significant",
                                                  T ~ "NS"))
                   })

SE.glm.ptrend <- pmap_dfr(tidyr::crossing(scale = antichol.scales,
                                   SE_var = c(peripheral.SEs, central.SEs)),
                   function(scale, SE_var) {
                     fit <- glm(reformulate(termlabels = scale, response = SE_var), 
                                data = SE.glm.df %>% 
                                  mutate(ISACS = factor(as.character(cholinergicCluster),
                                                          levels = c("No", "Low", "Moderate", "High"),
                                                          ordered = T),
                                         ACSBC = factor(as.character(ACSBC),
                                                        levels = c("No", "Low", "Moderate", "High"),
                                                        ordered = T),
                                         ORCA = factor(as.character(ORCA),
                                                       levels = c("No", "Low", "Moderate", "High"),
                                                       ordered = T)
                                  ),
                                family = binomial)
                     pred_prob <- fitted(fit)
                     actual <- model.frame(fit)[[SE_var]]
                     
                     roc_obj <- roc(response = model.frame(fit)[[SE_var]],
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
                              # Emax = calib_stats[["Emax"]],
                              # calib = list(pred_df),
                              sig = case_when(p.value < 0.05 ~ "Significant",
                                                  T ~ "NS"))
                   }) 

SE.glm %>%
  mutate(Model = factor(Model, levels = antichol.scales),
         SE.class = case_when(SE_name %in% central.SEs~"CNS",
                              SE_name %in% peripheral.SEs~"PNS")) %>% 
  distinct(Model, SE_name, SE.class, AUROC, Baseline_AUPRC, AUPRC, Brier) %>%
  pivot_longer(c(AUROC, AUPRC, Brier)) %>%
  ggplot(aes(x = Model, y = value, fill = SE.class)) +
  geom_line(aes(group = SE_name, color = SE.class), 
            linewidth = 1, alpha = 0.9) + 
  geom_point(size = 3, pch = 21) + 
  geom_boxplot(outliers = F, alpha = 0.5, width = 0.5) +
  labs(title = "Single adverse effect analysis (AUC-ROC)",
       subtitle = "Incorporating PK parameters", x = "") +
  facet_grid(name~SE.class, scales = "free_y") +
  scale_fill_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  scale_color_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  mytheme + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

SE.glm.ptrend %>%
  mutate(Model = factor(Model, levels = antichol.scales),
         SE.class = case_when(SE_name %in% central.SEs~"CNS",
                              SE_name %in% peripheral.SEs~"PNS")) %>% 
  distinct(Model, SE_name, SE.class, AUROC, Baseline_AUPRC, AUPRC, Brier) %>%
  pivot_longer(c(AUROC, AUPRC, Brier)) %>%
  ggplot(aes(x = Model, y = value, fill = SE.class)) +
  geom_line(aes(group = SE_name, color = SE.class), 
            linewidth = 1, alpha = 0.9) + 
  geom_point(size = 3, pch = 21) + 
  geom_boxplot(outliers = F, alpha = 0.5, width = 0.5) +
  labs(title = "Single adverse effect analysis (AUROC)",
       subtitle = "Incorporating PK parameters", x = "") +
  facet_grid(name~SE.class, scales = "free_y") +
  scale_fill_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  scale_color_manual(values = c("CNS" = "#9c1806", "PNS" = "#415494")) +
  mytheme + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

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

# Aggregated
pred.SIDER <- c("ISACS", "ACSBC", "ORCA")
n.cSE <- 1:length(central.SEs)
n.pSE <- 1:length(peripheral.SEs)

cutoff_grid <- expand.grid(n.cSE = n.cSE, n.pSE = n.pSE)

compare.scale.res <- cutoff_grid %>%
  mutate(result = map2(n.cSE, n.pSE, ~{
    SE.glm.df %>% 
      mutate(ISACS = factor(as.character(ISACS),
                              levels = c("No", "Low", "Moderate", "High"),
                              ordered = T),
             ORCA = factor(as.character(ORCA),
                           levels = c("No", "Low", "Moderate", "High"),
                           ordered = T),
             ACSBC = factor(as.character(ACSBC),
                            levels = c("No", "Low", "Moderate", "High"),
                            ordered = T)
      ) %>%
      mutate(n.central.SEs = rowSums(across(c(central.SEs)), na.rm = T),
             n.peripheral.SEs = rowSums(across(c(peripheral.SEs)), na.rm = T),
             anticholinergic.SE = case_when(n.central.SEs >= .x & n.peripheral.SEs >= .y ~ 1, T~0)) %>%
      dplyr::select(c(Drug, ISACS, ACSBC, ORCA, BBB, pgp_sub, anticholinergic.SE)) %>% 
      group_modify(~{
        map_dfr(antichol.scales, function(mod){
          formula_str <- as.formula(paste("anticholinergic.SE ~", mod))
          fit <- glm(formula_str,
                     data = .x,
                     family = binomial)
          
          roc_obj <- pROC::roc(response  = .x$anticholinergic.SE,
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

col_fun = circlize::colorRamp2(c(0.4, 0.525, 0.65, 0.775, 0.9), 
                               c("#d6a9ab", "#b57d7f", "#824446", "#6e2b2d", "#4a1712")) 
col_palette <- col_fun(seq(0, 1, length.out = 256))
scales.p1 <- compare.scale.res %>% 
  distinct(n.cSE, n.pSE, AUROC, model) %>% 
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
  labs(title = "Aggregated adverse effect analysis",
       y = "Number of peripheral adverse effects considered", 
       x = "Number of central adverse effects considered") +
  theme(panel.grid.major = element_line(color = "white", linewidth = 1),
        panel.grid.minor = element_line(color = "white", linewidth = 1),
        panel.spacing.x = unit(1.1, "lines"))

scales.p2 <- compare.scale.res %>% 
  distinct(n.cSE, n.pSE, AUROC, model) %>% 
  ggplot(aes(x = model, y = AUROC, fill = model)) +
  geom_violin(alpha = 0.25) + geom_boxplot(width = 0.2) +
  labs(x = "") +
  mytheme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")


# Summary overall:
library(ggupset)
SE.summary.df1 <- SE.glm %>%
  mutate(Model = factor(Model, levels = antichol.scales),
         SE.class = case_when(SE_name %in% central.SEs~"CNS single AE",
                              SE_name %in% peripheral.SEs~"PNS single AE")) %>% 
  distinct(Model, SE_name, SE.class, AUROC, Baseline_AUPRC, AUPRC, Brier) %>%
  dplyr::select(-SE_name)
  
SE.summary.df2 <- compare.scale.res %>%
  mutate(SE.class = "Aggregated AE") %>% 
  distinct(model, SE.class, n.cSE, n.pSE, AUROC, Baseline_AUPRC, AUPRC, Brier) %>%
  dplyr::select(-n.cSE, -n.pSE) %>%
  dplyr::rename(Model = model)


combos <- list(c("Base"),
               c("Base", "pgp_sub"),
               c("Base", "BBB"),
               c("Base", "BBB", "pgp_sub"))

model_df <- tidyr::crossing(base_model = pred.SIDER,
                            covariate  = combos) %>%
  mutate(combo_set = map2(base_model, covariate, ~ c(.x, .y[.y != "Base"])),
         label = map_chr(combo_set, ~ paste(.x, collapse = "+"))) %>%
  right_join(., rbind(SE.summary.df1, SE.summary.df2), c("label" = "Model")) %>%
  mutate(label = factor(label, levels = antichol.scales),
         SE.class = factor(SE.class, 
                           levels = c("CNS single AE",
                                      "PNS single AE",
                                      "Aggregated AE")))

burner <- data.frame(val = c(1:length(antichol.scales)),
                     label = antichol.scales) %>% 
  right_join(., model_df, "label") %>%
  rowwise() %>% 
  mutate(dup_times = length(antichol.scales)+1-val) %>%
  ungroup() %>%
  uncount(weights = dup_times) %>% 
  mutate(x = 0,
         AUROC = NA,
         Brier = NA) %>% 
  dplyr::select(-val)

summaryplot <- model_df %>% 
  mutate(x = 1) %>%
  bind_rows(., burner) %>%
  ggplot(aes(x = combo_set, y = AUROC, fill = SE.class)) +
  scale_x_upset(sets = c(pred.SIDER, "BBB", "pgp_sub"),
                order_by = "freq") +
  geom_rect(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_point(size = 3, pch = 21, alpha = 0.5) + 
  geom_boxplot(aes(group = label), outliers = F, alpha = 0.75, width = 0.5) +
  labs(title = "", x = "", y = "AUC-ROC") +
  facet_grid(.~SE.class, scales = "free_y") +
  scale_fill_manual(values = c("CNS single AE" = "#9c1806", "PNS single AE" = "#415494")) +
  scale_color_manual(values = c("CNS single AE" = "#9c1806", "PNS single AE" = "#415494")) +
  mytheme + 
  theme(legend.position = "none",
        panel.spacing = unit(3.5, "lines")) +
  theme_combmatrix(combmatrix.label.text = element_text(size=12),
                   combmatrix.label.extra_spacing = 5) 
 

rbind(SE.summary.df1, SE.summary.df2) %>%
  mutate(Model = factor(Model, levels = antichol.scales),
         SE.class = factor(SE.class, 
                           levels = c("CNS single AE",
                                      "PNS single AE",
                                      "Aggregated AE"))) %>% 
  ggplot(aes(x = Model, y = AUROC, fill = SE.class)) +
  geom_rect(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
            fill = "grey85", inherit.aes = F) +
  geom_point(size = 3, pch = 21) + 
  geom_boxplot(outliers = F, alpha = 0.5, width = 0.5) +
  labs(title = "", x = "") +
  facet_grid(.~SE.class, scales = "free_y") +
  scale_fill_manual(values = c("CNS single AE" = "#9c1806", "PNS single AE" = "#415494")) +
  scale_color_manual(values = c("CNS single AE" = "#9c1806", "PNS single AE" = "#415494")) +
  mytheme + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none") 

# Cumulative
SE.glm.df.cumsum <- SE.glm.df %>%
  mutate(n.central.SEs = rowSums(across(c(central.SEs)), na.rm = T),
         n.peripheral.SEs = rowSums(across(c(peripheral.SEs)), na.rm = T),
         n.all.SEs = n.central.SEs+n.peripheral.SEs) 

SE.cum.plot1 <- SE.glm.df.cumsum %>% 
  pivot_longer(cols = starts_with("n.")) %>%
  mutate(name = case_when(name == "n.all.SEs"~"Number of all AEs",
                          name == "n.central.SEs"~"Number of CNS AEs",
                          name == "n.peripheral.SEs"~"Number of PNS AEs")) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(alpha=..count.., y = ..density..), 
                 fill="#0AD1C1", color = "black",
                 binwidth = 2) +
  geom_density(size = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  facet_grid(~name) +
  labs(x = "Number of AEs", y = "Density") +
  mytheme

SE.glm.cumsum <- map_dfr(c(), function(SE_var) {
  fit <- glm(reformulate(c("ISACS"), response = SE_var), 
             data = SE.glm.df.cumsum,
             family = poisson)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(SE_name = SE_var,
           sig = case_when(p.value < 0.05~"Significant",
                           T~"NS"))
})

SE.glm.cumsum <- pmap_dfr(tidyr::crossing(scale = antichol.scales,
                                          SE_var = c("n.central.SEs", "n.peripheral.SEs", "n.all.SEs")),
                          function(scale, SE_var) {
                            fit <- MASS::glm.nb(reformulate(termlabels = scale, response = SE_var), 
                                                data = SE.glm.df.cumsum)
                            aic.val <- AIC(fit)
                            bic.val <- BIC(fit)
                            dispersion <- sum(residuals(fit, type = "pearson")^2)/df.residual(fit)
                            cv_dev <- -2 * as.numeric(logLik(fit)) / nobs(fit)
                            
                            broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                              dplyr::filter(term != "(Intercept)") %>%
                              mutate(SE_name = SE_var,
                                     Scale = scale,
                                     AIC = aic.val,
                                     BIC = bic.val,
                                     cv_dev = cv_dev, 
                                     dispersion = dispersion,
                                     sig = case_when(p.value < 0.05~"Significant",
                                                     T~"NS"))
                          }) %>% 
  mutate(term = gsub("ACSBC", "", term)) %>%
  mutate(term = gsub("ISACS", "", term)) %>%
  mutate(term = gsub("ORCA", "", term)) %>%
  mutate(term = factor(term, levels = c("High", "Moderate", "Low")),
         SE_name = case_when(SE_name == "n.all.SEs"~"Number of all AEs",
                             SE_name == "n.central.SEs"~"Number of CNS AEs",
                             SE_name == "n.peripheral.SEs"~"Number of PNS AEs"),
         Scale = factor(Scale, levels = c("ACSBC", "ISACS", "ORCA")))

SE.glm.cumsum.ptrend <- pmap_dfr(tidyr::crossing(scale = antichol.scales,
                                                 SE_var = c("n.central.SEs", "n.peripheral.SEs", "n.all.SEs")),
                                 function(scale, SE_var) {
                                   fit <- MASS::glm.nb(reformulate(termlabels = scale, response = SE_var), 
                                                       data = SE.glm.df.cumsum %>% 
                                                         mutate(ISACS = factor(as.character(cholinergicCluster),
                                                                                 levels = c("No", "Low", "Moderate", "High"),
                                                                                 ordered = T)) %>%
                                                         mutate(ORCA = factor(as.character(ORCA),
                                                                              levels = c("No", "Low", "Moderate", "High"),
                                                                              ordered = T)) %>%
                                                         mutate(ACSBC = factor(as.character(ACSBC),
                                                                               levels = c("No", "Low", "Moderate", "High"),
                                                                               ordered = T))
                                   )
                                   aic.val <- AIC(fit)
                                   bic.val <- BIC(fit)
                                   dispersion <- sum(residuals(fit, type = "pearson")^2)/df.residual(fit)
                                   cv_dev <- -2 * as.numeric(logLik(fit)) / nobs(fit)
                                   
                                   broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
                                     dplyr::filter(term != "(Intercept)") %>%
                                     mutate(SE_name = SE_var,
                                            Scale = scale,
                                            AIC = aic.val,
                                            BIC = bic.val,
                                            cv_dev = cv_dev, 
                                            dispersion = dispersion,
                                            sig = case_when(p.value < 0.05~"Significant",
                                                            T~"NS"))
                                 }) %>%  
  dplyr::filter(term %in% c("ISACS.L", "ORCA.L", "ACSBC.L")) %>% 
  mutate(SE_name = case_when(SE_name == "n.all.SEs"~"Number of all AEs",
                             SE_name == "n.central.SEs"~"Number of CNS AEs",
                             SE_name == "n.peripheral.SEs"~"Number of PNS AEs"),
         Scale = factor(Scale, levels = c("ACSBC", "ISACS", "ORCA"))) %>%
  dplyr::select(SE_name, Scale, p.value) %>%
  add_significance("p.value") %>%
  mutate(p.value.signif = gsub("ns", "",p.value.signif),
         p.trend = paste0(format_p(p.value), p.value.signif))

SE.cum.plot2 <- SE.glm.cumsum %>%
  left_join(., SE.glm.cumsum.ptrend, c("Scale", "SE_name")) %>%
  ggplot(aes(x = estimate, y = Scale, color = term, shape = sig,
             xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(size = 0.65, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = p.trend, x = 3.15), 
            size = 3.25,
            color = "black") +
  geom_stripped_rows(color = NA) +
  scale_color_manual(values = rev(pal_jama("default")(4)[-1])) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_grid(~SE_name) +
  labs(title = "Performance of the scales on cumulative\nanticholinergic-related adverse effects",
       color = "Level", y = "", x = "Estimate") +
  mytheme

SE.cum.plot3 <- SE.glm.cumsum %>% dplyr::select(SE_name, Scale, AIC, BIC) %>% distinct %>%
  pivot_longer(cols = c(AIC, BIC)) %>%
  ggplot(aes(x = Scale, y = value, group = SE_name, fill = SE_name)) +
  geom_line() +
  geom_point(pch = 21, size = 3) +
  scale_fill_futurama() +
  facet_grid(~name) +
  labs(y = "Performance metrics", fill = "", 
       x = "") +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Supplementary Table
write_clip(ORA.result1@compareClusterResult %>% mutate(Scale = "ISACS") %>%
             rbind(ORA.result2@compareClusterResult %>% mutate(Scale = "ORCA"),
                   ORA.result3@compareClusterResult %>% mutate(Scale = "ACSBC")) %>%
             dplyr::filter(p.adjust < 0.05) %>%
             dplyr::select(Scale, Cluster, Description, GeneRatio, BgRatio, 
                           pvalue, p.adjust, qvalue, geneID, Count))

CompScaleplot <- result.w.anno.full %>%
  mutate(ORCA = factor(ORCA, levels = c("No", "Low", "Moderate", "High"), ordered = T),
         ACSBC = factor(ACSBC.Anticholinergic.Activity,
                        levels = c("No", "Low", "Moderate", "High"),
                        ordered = T),
         ISACS = cholinergicCluster) %>%
  # dplyr::filter(if_any(c(ACSBC, ORCA), ~ !is.na(.) & . != "")) %>%
  pivot_longer(cols = c("ACSBC", "ORCA", "ISACS")) %>%
  mutate(name = factor(name, levels = pred.SIDER)) %>%
  left_join(., ATC.df, "drugbank_id") %>%
  mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")")) %>%
  drop_na(value, name_1st) %>%
  dplyr::filter(ATC_3rd %in% c(ORA.result1@compareClusterResult %>% mutate(Scale = "ISACS") %>%
                                 rbind(ORA.result2@compareClusterResult %>% mutate(Scale = "ORCA"),
                                       ORA.result3@compareClusterResult %>% mutate(Scale = "ACSBC")) %>%
                                 dplyr::filter(!Cluster %in% c("No consensus", "No")) %>%
                                 dplyr::filter(p.adjust < 0.01) %>% pull(Description))) %>%
  ggplot(aes(x = fct_rev(factor(ATC_3rd)), fill = fct_rev(factor(value)))) +
  geom_bar(stat = "count", position = "stack", color = "black", width = 0.8) +
  scale_fill_manual(values = rev(pal_jama()(4))) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(y = "Number of drugs", x = "", fill = "Anticholinergic level",
       title = "ATC 3rd Level Analysis") +
  facet_grid(~name) +
  mytheme

# ## Drug class
# tiff(paste0("Results/CompScale_DrugClass2.tiff"), units="in", height=10, width=15, res=300)
ORA.result1@compareClusterResult %>% mutate(Scale = "ISACS") %>%
  rbind(ORA.result2@compareClusterResult %>% mutate(Scale = "ORCA"),
        ORA.result3@compareClusterResult %>% mutate(Scale = "ACSBC")) %>%
  dplyr::filter(Cluster != "No consensus") %>%
  dplyr::filter(Cluster == "Low") %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Scale = factor(Scale, levels = pred.SIDER),
         Cluster = factor(Cluster,
                          levels = levels(result.w.anno.full$cholinergicCluster))) %>%
  ggplot(aes(x = Cluster, y = fct_rev(Description),
             size = Count, fill = Cluster, alpha = -log10(p.adjust))) +
  geom_point(pch = 21, color = "white") +
  labs(x = "", y = "", alpha = "-log10(adjusted p-value)",
       fill = "Level",
       title = "ATC 3rd Level Enrichment Analysis") +
  scale_size(range = c(5,10)) +
  scale_alpha(range = c(0.5, 0.9)) +
  scale_fill_jama() +
  facet_grid(~Scale) +
  guides(fill = guide_legend(override.aes = list(size = 5, color = "black")),
         size = guide_legend(override.aes = list(color = "black", pch = 16)),
         alpha = guide_legend(override.aes = list(color = "black", fill = "grey50",
                                                  size = 5))) +
  mytheme +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        panel.grid.major = element_line(color = "white", linewidth = 1))
# dev.off()

tiff("Results/Figure4.tiff", units="in", height=18, width=12, res=300)
plot_grid(plot_grid(NULL, grid.grabExpr(draw(scale.compare.ht)), NULL,
                    nrow = 1, rel_widths = c(0.2, 1, 0.2)),
          summaryplot,
          CompScaleplot,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          rel_heights = c(1, 0.9, 1.5),
          ncol = 1)
dev.off()

Label_lode <- SE.glm.df %>% 
  dplyr::select(Drug, ORCA, ACSBC, ISACS) %>%
  pivot_longer(pred.SIDER) %>%
  count(name, value) %>%
  group_by(name) %>%
  mutate(Perc = round(n/(sum(n)),4)*100) %>% ungroup %>%
  mutate(Perc = paste0(sprintf(" \n(%s", Perc), "%)")) %>%
  mutate(label = paste0(n, " ", Perc))

df_lode <- SE.glm.df %>% dplyr::select(Drug, ORCA, ACSBC, ISACS) %>%
  column_to_rownames("Drug") %>%
  to_lodes_form() %>%  
  left_join(., Label_lode,
            by = c("stratum" = "value",
                   "x" = "name")) %>%
  mutate(x = factor(x, levels = pred.SIDER),
         stratum = factor(stratum, levels = c("High", "Moderate", "Low", "No")),
         lab.col = case_when(stratum == "No"~"X",
                             T~"Y")
         ) %>% 
  drop_na()

tiff("Results/SFigure14.tiff", units="in", height=5.5, width=6.5, res=300)
ggplot(df_lode, aes(x = x, stratum = stratum, alluvium = alluvium, 
                    fill = stratum, label = label)) +
  labs(x = "", y = "Number of drugs", 
       title = "Drugs overlapping across the 3 scales") + 
  geom_flow(stat = "flow", 
            width = 0.6,  knot.pos = 1/3) +
  geom_stratum(width = 0.6, na.rm = F) + 
  geom_text(stat = "stratum", size = 3.5, show.legend = F,
            fontface = "bold",
            aes(color = lab.col)) +
  theme_sankey(base_size = 15) +
  scale_color_manual(values = c("white", "black")) +
  scale_fill_manual(values = rev(pal_jama("default")(4))) +
  theme(legend.position = "right", legend.title = element_blank(),
        plot.title = element_text(face="bold",  size=14, hjust = 0.5),
        plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        strip.text = element_text(size = 11, color = "black"))
dev.off()

# Margin of safety
## IC50
IC50 <- read.csv("1_PDI/BA.mAChR.IC50.csv") %>%
  dplyr::filter(DRUGNAME %in% result.w.anno.full$Drug) 

table(IC50$GENENAME)

## Cmax data
Cmax <- read.csv("11_Cmax/tx0c00294_si_002/supplementary_ChemResTox/DataFile_08_computed_median_total_and_unbound_plasma_conc_per_molregno.csv") %>%
  mutate(Drug = str_to_sentence(pref_name))
dim(Cmax)

Cmax.df <- Cmax %>%
  left_join(., result.w.anno.full, 
            by = c("Drug")) %>%
  left_join(., result.w.anno.full, 
            by = c("chembl_id" = "ChEMBL.ID"), suffix = c("", "_secondary")) %>%
  mutate(across(
    .cols = ends_with("_secondary"),
    .fns = ~ coalesce(get(str_remove(cur_column(), "_secondary")), .),
    .names = "{str_remove(.col, '_secondary')}"
  )) %>%
  select(-ends_with("_secondary"))


Cmax.df %>% 
  mutate(ActivityCall = median.pMolar.unbound.plasma.concentration/QSAR_score) %>%
  dplyr::filter(ActivityCall >= 1) %>%
  distinct(Drug, cholinergicCluster, ORCA, ACSBC.Anticholinergic.Activity, ActivityCall) %>%
  arrange(cholinergicCluster, ActivityCall)

Cmax.df %>% 
  pivot_longer(c(median.pMolar.total.plasma.concentration,
                 median.pMolar.unbound.plasma.concentration)) %>%
  drop_na(cholinergicCluster) %>%
  ggplot(aes(x = cholinergicCluster, y = value, 
             fill = cholinergicCluster)) +
  ggpubr::stat_compare_means(comparisons = combn(unique(result.w.anno.full$cholinergicCluster), 
                                                 2, simplify = F),
                             label = "p.signif", method = "wilcox.test") +
  geom_violin(alpha = 0.25) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(width = 0.2, alpha = 0.9) +
  scale_fill_jama() +
  facet_wrap(~name, scales = "free_y")

unbound.df <- Cmax.df %>% 
  mutate("Rat brain [Yamada et al. 2023 (ABS)]" = case_when(ABS > 0~Standardized_IC50_uM*1000,
                                                            ABS == 0~100000
  )) %>%
  left_join(., IC50, c("Drug"="DRUGNAME")) %>%
  pivot_wider(names_from = "GENENAME", values_from = "binding_value") %>%
  pivot_longer(c(unique(IC50$GENENAME), "Rat brain [Yamada et al. 2023 (ABS)]"),
               values_to = "binding_value",
               names_to = "GENENAME"
  ) %>%
  # Cmax in molar whereas binding value is in nM
  mutate(pIC50 = -log10(binding_value/10^9),
         # log_ratio = pIC50 - median.pMolar.unbound.plasma.concentration,
         log_ratio = log10(median.Molar.unbound.plasma.concentration/(binding_value/10^9))
  ) %>%
  drop_na(log_ratio) %>%
  mutate(Level.num = case_when(cholinergicCluster == "No"~1,
                               cholinergicCluster == "Low"~2,
                               cholinergicCluster == "Moderate"~3,
                               cholinergicCluster == "High"~4))
  
unbound.plot <- unbound.df %>%
  ggplot(aes(x = Level.num, y = log_ratio)) +
  geom_violin(alpha = 0.25, aes(fill = cholinergicCluster)) +
  geom_point(position = position_jitter(width = 0.2), aes(fill = cholinergicCluster)) +
  geom_boxplot(width = 0.2, alpha = 0.9, aes(fill = cholinergicCluster), outliers = F) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_jama() +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  scale_x_continuous(breaks = c(1:4), labels = c("No", "Low", "Moderate", "High")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  stat_cor(method = "spearman") +
  # stat_poly_eq(use_label(c("eq", "R2", "n")), label.x = "left", label.y = "top") +
  labs(x = "", y = "log10[(Cmax unbound)/IC50]") +
  facet_wrap(~GENENAME,
             labeller = as_labeller(function(x) paste0(x, "\n(n = ", table(unbound.df$GENENAME)[x], ")"))) +
  mytheme + theme(legend.position = "none")

total.df <- Cmax.df %>% 
  mutate("Rat brain [Yamada et al. 2023 (ABS)]" = case_when(ABS > 0~Standardized_IC50_uM*1000,
                                                            ABS == 0~100000
  )) %>%
  left_join(., IC50, c("Drug"="DRUGNAME")) %>%
  pivot_wider(names_from = "GENENAME", values_from = "binding_value") %>%
  pivot_longer(c(unique(IC50$GENENAME), "Rat brain [Yamada et al. 2023 (ABS)]"),
               values_to = "binding_value",
               names_to = "GENENAME"
  ) %>%
  # Cmax in molar whereas binding value is in nM
  mutate(pIC50 = -log10(binding_value/10^9),
         log_ratio = log10(median.Molar.total.plasma.concentration/(binding_value/10^9))
  ) %>%
  drop_na(log_ratio) %>%
  mutate(Level.num = case_when(cholinergicCluster == "No"~1,
                               cholinergicCluster == "Low"~2,
                               cholinergicCluster == "Moderate"~3,
                               cholinergicCluster == "High"~4)) 

total.plot <- total.df %>%
  ggplot(aes(x = Level.num, y = log_ratio)) +
  geom_violin(alpha = 0.25, aes(fill = cholinergicCluster)) +
  geom_point(position = position_jitter(width = 0.2), aes(fill = cholinergicCluster)) +
  geom_boxplot(width = 0.2, alpha = 0.9, aes(fill = cholinergicCluster), outliers = F) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_jama() +
  geom_smooth(method = "lm", se = F, color = "darkgreen") +
  scale_x_continuous(breaks = c(1:4), labels = c("No", "Low", "Moderate", "High")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  stat_cor(method = "spearman") +
  # stat_poly_eq(use_label(c("eq", "R2", "n")), label.x = "left", label.y = "top") +
  labs(x = "", y = "log10[(Cmax total)/IC50]") +
  facet_wrap(~GENENAME,
             labeller = as_labeller(function(x) paste0(x, "\n(n = ", table(total.df$GENENAME)[x], ")"))) +
  mytheme + theme(legend.position = "none")

KW.test <- result.w.anno.full %>% kruskal_test(Neurotoxicity.DI ~ cholinergicCluster)    

dunn.test <- result.w.anno.full %>%
  dunn_test(Neurotoxicity.DI ~ cholinergicCluster, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "cholinergicCluster") %>%
  mutate(y.position = y.position + 0.05) %>% dplyr::select(-groups)

neurotox <- result.w.anno.full %>%
  ggplot(aes(x = cholinergicCluster, y = Neurotoxicity.DI, fill = cholinergicCluster)) +
  geom_violin(alpha = 0.25) +
  geom_boxplot(width = 0.1) +
  scale_fill_jama() + mytheme +
  stat_pvalue_manual(dunn.test,
                     label = "p.adj.signif",  
                     step.increase = 0.05,
                     tip.length = 0.02, hide.ns = TRUE) +
  labs(x = "", y = "Neurotoxicity potential",
       subtitle = get_test_label(KW.test, detailed = TRUE)) +
  theme(legend.position = "none")

tiff(paste0("Results/SFigure15.tiff"), units="in", height=5, width=9, res=300)
total.plot + theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()

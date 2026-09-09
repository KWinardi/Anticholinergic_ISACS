# Cross-sectional pharmacoepidemiological analysis
## Data analysis and visualization
## Created by Kevin Winardi (26/08/2026)

if (!exists("phenoage.dat")) {
  source("8_1_NHANES_Setup.R")
}

#### Merge all into 1 dataset ####
NHANES.df <- demo_all %>%
  mutate(SEQN = as.character(SEQN)) %>%
  left_join(., phenoage.dat %>% select(-cycle), "SEQN") %>%
  left_join(., ISACS.NHANES.score, "SEQN") %>%
  left_join(., nMeds.3, "SEQN") %>%
  mutate(WTMECPOOL = as.numeric(WTMEC2YR) / 10) %>%
  left_join(., chronic_conditions, c("SEQN", "cycle_years", "RIDAGEYR", "RIAGENDR")) %>%
  left_join(., ALQ_df, c("SEQN", "cycle_years")) %>%
  left_join(., SMQ_df, c("SEQN", "cycle_years")) %>%
  left_join(., PAQ_df, c("SEQN", "cycle_years")) %>%
  left_join(., CFQ_df, c("SEQN", "cycle_years")) %>%
  # Covariate adjustments
  mutate(cycle_years = factor(cycle_years),
         gender = factor(case_when(RIAGENDR == 1~"Male",
                                   RIAGENDR == 2~"Female")),
         ethnicity = factor(case_when(ethnicity == 1~"White",
                                      ethnicity == 2~"Black",
                                      ethnicity == 3~"Hispanic",
                                      ethnicity == 4~"Others"
         ), levels = c("White", "Black", "Hispanic","Others")),
         bmi = ifelse(bmi > 60, 60, bmi), # Put a cap on humanly possible range
         educationlevel = factor(case_when(education >= 3~"≥12",
                                           education < 3~"<12")),
         maritalstatus = factor(case_when(DMDMARTL == 1~"Married/Cohabit",
                                          DMDMARTL == 6~"Married/Cohabit",
                                          DMDMARTL > 6~NA,
                                          T~"Unmarried")),
         annual_income_bin = factor(case_when(annual_income >= 5~"≥$20K",
                                              annual_income < 5~"<$20K")),
         smoking = factor(case_when(SMQ020 == 2 ~ "Never",
                                    SMQ020 == 1 & SMQ040 == 3 ~ "Former",
                                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
                                    TRUE ~ NA_character_
                                    ), levels = c("Never", "Former", "Current")),
         alcohol = factor(case_when(alcohol == 1 ~ "Current drinker",
                                    alcohol == 0 ~ "Not drinker",
                                    alcohol %in% c(7, 9)~NA
         ), levels = c("Not drinker", "Current drinker")),
         active = factor(case_when(active == 1~"Active",
                                   active == 0~"Inactive",
                                   T~NA)),
         n_medicines = RXDCOUNT
  ) %>%
  mutate(kappa_ckdepi = if_else(gender == "Female", 0.7, 0.9),
         alpha_ckdepi = if_else(gender == "Female", -0.241, -0.302),
         min_ratio = pmin(creat / kappa_ckdepi, 1),
         max_ratio = pmax(creat / kappa_ckdepi, 1),
         egfr_CKD_EPI = 142 * (min_ratio ^ alpha_ckdepi) * (max_ratio ^ -1.200) *
           (0.9938 ^ age) * if_else(gender == "Female", 1.012, 1)) %>%
  select(-kappa_ckdepi, -alpha_ckdepi, -min_ratio, -max_ratio) 

drug_sampleIDs <- ISACS.NHANES.score$SEQN

NHANES.df %>% 
  ggplot(aes(x = age, y = phenoage, color = nMeds.class)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_abline(slope = 1, color = "red") +
  scale_color_bmj() +
  facet_wrap(~cycle_years) +
  labs(x = "Chronological age (y)", y = "PhenoAge (y)", color = "") +
  mytheme 

# Determine variables
outcome_labels <- c(# Cognitive function
                    CFDDS = "DSST",
                    CFDAST = "Animal fluency test",
                    CFDCST_total = "CERAD-IR",
                    CFDCSR = "CERAD-DR",
                    # Biological age
                    kdm_advance0 = "KDM gap",
                    kdm_advance = "Modified KDM gap",
                    phenoage_advance0 = "Levine phenotypic age gap",
                    phenoage_advance  = "Modified Levine phenotypic age gap"
                    )

covars <- c("age", "gender", "ethnicity", "cycle_years",
            "educationlevel", "maritalstatus", "annual_income_bin",
            "bmi", "active", "smoking", "alcohol",
            "egfr_CKD_EPI", "n_conditions")

NHANES.df <- NHANES.df %>%
  # Cohort of interest for analysis
  mutate(inAnalysis = if_all(all_of(c(covars, "ISACS.score")), ~ !is.na(.) &
                               age >= 60 & n_conditions_missing != length(cond_cols)) & 
           if_any(names(outcome_labels),~!is.na(.)))

NHANES.sub <- NHANES.df %>% filter(inAnalysis == T)

# Supplementary Table 16
print(tableone::CreateTableOne(vars = c(covars, "n_medicines", "nMeds.class"), data = NHANES.sub),
      nonnormal = c("age", "n_conditions", "n_medicines")) %>% write_clip

print(tableone::CreateTableOne(vars = colnames(ISACS.NHANES.score)[-1],
                               strata = "nMeds.class",
                               data = NHANES.sub))

# Checking missingness
## Overall
NHANES.sub %>%
  summarise(across(all_of(covars),
                   list(missing = ~ sum(is.na(.)),
                        filled  = ~ sum(!is.na(.))
                        ))) %>%
  pivot_longer(everything(),
               names_to = c("variable", ".value"),
               names_pattern = "(.+)_(missing|filled)$") %>%
  mutate(pct_missing = round(missing/(missing+filled)*100, 2)) %>%
  arrange(pct_missing)

## Per variable stratified by cycle
NHANES.df %>% filter(age >= 20) %>%
  count(annual_income_bin, cycle_years)

tic("Association analysis...")
stat_result <- list()
for (outcome in names(outcome_labels)) {
  # Fixed outcome and covariates
  model_vars <- c(outcome, "ISACS.score", covars)
  
  sampleIDs <- NHANES.df %>%
    filter(if_all(all_of(model_vars), ~ !is.na(.) & age >= 60 & n_conditions_missing != length(cond_cols))) %>%
    pull(SEQN)
  
  NHANES.sub <- NHANES.df %>% mutate(inAnalysis = SEQN %in% sampleIDs)
  
  nhanes_design <- svydesign(id = ~SDMVPSU,
                             strata  = ~SDMVSTRA,
                             weights = ~WTMECPOOL,
                             nest = TRUE,
                             data = NHANES.sub)
  
  design_sub <- subset(nhanes_design, inAnalysis)
  NHANES.sub <- NHANES.sub %>% filter(inAnalysis == T)
  
  message("Analysing ", outcome, "...")
  # Loop
  exposure_vars <- colnames(ISACS.NHANES.score)[-1]
  vif_list <- list()
  
  results <- map_dfr(exposure_vars, function(var1) {
    # Add non-anticholinergic drug count here - varies with scales
    var1_covar <- paste0(gsub(".score", "", var1), ".nMeds_diff") 
    f <- reformulate(c(var1, covars, var1_covar), response = outcome) 
    model <- svyglm(f, design = design_sub, family = gaussian())
    
    tidy(model, conf.int = TRUE) %>% 
      mutate(exposure = var1, .before = 1)  %>%
      mutate(sample.size = n_distinct(NHANES.sub$SEQN),
             outcome = outcome_labels[[outcome]])
  })
  stat_result[[outcome]] <- results
  stat_vif[[outcome]] <- bind_rows(vif_list)
}

stat_result_all <- bind_rows(stat_result, .id = "outcome") %>%
  mutate(outcome_label = outcome_labels[outcome],
         exposure = gsub(".score", "", exposure),
         exposure = fct_rev(factor(exposure, 
                                   levels = c("ISACS","ACSBC", "ADS", "ABS", "ALS","ORCA", 
                                              "AEC","Chew", "ACB", "AAS", 
                                              "ARS", "ABC", "ATS"))),
         signif = factor(ifelse(p.value < 0.05, "Significant", "NS"),
                         levels = c("Significant", "NS")),
         estimate_color = if_else(p.value < 0.05, estimate, NA_real_),
         facet_lab = paste0(outcome_label, "\nn = ", sample.size, " participants"),
         forest.lab = paste0(round(estimate, 2), " [", round(conf.low, 2), ", ", round(conf.high, 2), "], p = ", round(p.value, 3)),
  )
toc()

# Supplementary table 17
stat_result_all %>%
  filter(term %in% exposure_vars) %>% 
  select(outcome, outcome_label, exposure, 
         estimate, sample.size, std.error, statistic, 
         p.value, signif, conf.low, conf.high) %>%
  dplyr::rename("Scale"=exposure) %>%
  filter(!str_detect(outcome, "kdm")) %>% write_clip

# Cognition
NHANES.cog.res <- stat_result_all %>%
  filter(term %in% exposure_vars) %>%
  filter(str_detect(outcome, "^CFD")) %>% 
  filter(outcome != "CFDRIGHT") %>%
  mutate(label_x = conf.high + 0.03, # if_else(estimate < 0, conf.low - 0.02, conf.high + 0.02),
         label_hjust = 0, # if_else(estimate < 0, 1, 0),
         outcome = factor(outcome, 
                          levels = c("CFDDS", "CFDAST", "CFDCST_total", "CFDCSR")))

## Plot
NHANES.forest.cog <- NHANES.cog.res %>%
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, 
             y = fct_reorder(exposure, estimate), 
             color = estimate_color,
             shape = signif, label = forest.lab)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(size = 1, linewidth = 1) +
  scale_shape_manual(values = c("Significant" = 17, "NS" = 19)) +
  geom_text(aes(x = label_x, hjust = label_hjust), size = 4) +
  facet_wrap(~outcome, scales = "free", 
             labeller = labeller(outcome = deframe(distinct(NHANES.cog.res, outcome, facet_lab)))) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 1.05))) +
  scale_color_gradient2(#low = "#6A4C93", mid = "#A8A878", high = "#8C9E5E",
    low = "#264653", mid = "#C9A66B", high = "#E76F51",
    # low = "#2A9D8F", mid = "#D4A5A5", high = "#E76F51",
    midpoint = 0, na.value = "grey65") +
  labs(x = "Estimate (95% CI)", y = "", shape = "Statistical significance", color = "Estimate",
       # title = paste0("Association with ", outcome_labels[[outcome]]),
       subtitle = paste0("n = ", n_distinct(NHANES.sub$SEQN), " participants"),
       title = "Survey-weighted association analysis: NHANES Cohort"
  ) +
  mytheme

# tiff("Results/NHANES_forest_cog.tiff", units="in", height=10, width=12.5, res=300)
NHANES.forest.cog
# dev.off()

# Biological age
NHANES.BA.res <- stat_result_all %>%
  filter(term %in% exposure_vars) %>%
  filter(str_detect(outcome, "_advance")) %>% 
  mutate(label_x = conf.high + 0.03, # if_else(estimate < 0, conf.low - 0.02, conf.high + 0.02),
         label_hjust = 0, # if_else(estimate < 0, 1, 0),
         outcome = factor(outcome, 
                          levels = c("phenoage_advance0", "phenoage_advance",
                                     "kdm_advance0", "kdm_advance")))

## Plot
NHANES.forest.BA <- NHANES.BA.res %>%
  filter(!str_detect(outcome, "kdm")) %>%
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, 
             y = fct_reorder(exposure, estimate), 
             color = estimate_color,
             shape = signif, label = forest.lab)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(size = 1, linewidth = 1) +
  scale_shape_manual(values = c("Significant" = 17, "NS" = 19)) +
  geom_text(aes(x = label_x, hjust = label_hjust), size = 4) +
  facet_wrap(~outcome, scales = "free", 
             labeller = labeller(outcome = deframe(distinct(NHANES.BA.res, outcome, facet_lab)))) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.8))) +
  scale_color_gradient2(low = "#264653", mid = "#C9A66B", high = "#E76F51",
                        midpoint = 0, na.value = "grey65") +
  labs(x = "Estimate (95% CI)", y = "", shape = "Statistical significance", color = "Estimate",
       # title = paste0("Association with ", outcome_labels[[outcome]]),
       # subtitle = paste0("n = ", n_distinct(NHANES.sub$SEQN), " participants"),
       # subtitle = "60 years or older",
       title = "Survey-weighted association analysis: NHANES Cohort"
  ) +
  mytheme

# tiff("Results/NHANES_forest_BA.tiff", units="in", height=6.25, width=12, res=300)
NHANES.forest.BA
# dev.off()

# Summary
stat_result_all %>%
  filter(term %in% exposure_vars) %>%
  filter(!str_detect(outcome, "kdm")) %>%
  filter(signif == "Significant") %>%
  mutate(signlogp = sign(estimate)*-log10(p.value)) %>%
  select(exposure, outcome_label, signlogp) %>%
  pivot_wider(names_from = outcome_label, values_from = signlogp) %>% write_clip

### Weighted prevalence ###
n_distinct(RXQ_RX_df$RXDDRUG)
n_distinct(ISACS.NHANES$Drug)

drug_indicators <- RXQ_RX_df %>%
  filter(RXDUSE %in% c(1,2)) %>%
  filter(RXDDRUG != "", !is.na(RXDDRUG)) %>%
  left_join(., ISACS.NHANES, "RXDDRUG") %>%
  filter(!is.na(Drug)) %>%
  distinct(SEQN, Drug) %>%              
  mutate(present = 1) %>%
  pivot_wider(names_from = Drug, values_from = present, values_fill = 0)

drug_name_lookup <- tibble(Drug = colnames(drug_indicators)[-1],
                           clean = make.names(colnames(drug_indicators)[-1]))
colnames(drug_indicators)[-1] <- drug_name_lookup$clean

NHANES_drugs <- NHANES.df %>%
  mutate(WTINT2YR = as.numeric(WTINT2YR),
         WTINTPOOL = WTINT2YR / 10) %>%
  left_join(drug_indicators, by = "SEQN") %>%
  mutate(across(all_of(colnames(drug_indicators)[-1]), ~ replace_na(., 0)))

drug_design <- svydesign(id = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTINTPOOL,   
                         nest = TRUE,
                         data = NHANES_drugs)

drug_names <- colnames(drug_indicators)[-1]
drug_design_sub <- subset(drug_design, inAnalysis)

tic("Drug prevalence calculation")
batches <- split(drug_names, ceiling(seq_along(drug_names)/50))
drug_prevalence <- map_dfr(batches, function(v) {
  message("Calculating for drug batch ", names(v))
  fmla <- as.formula(paste("~", paste(v, collapse = " + ")))
  est <- svymean(fmla, design = drug_design_sub, na.rm = TRUE)
  tibble(clean = v, weighted_prevalence = coef(est), se = SE(est))
})
toc()

top15_weighted <- drug_prevalence %>%
  left_join(., drug_name_lookup, "clean") %>% 
  left_join(., ISACS.NHANES, "Drug") %>%
  distinct(ISACS, Drug, weighted_prevalence, se) %>% 
  mutate(ISACS = ifelse(is.na(ISACS), 0, ISACS)) %>%
  group_by(ISACS) %>%
  slice_max(weighted_prevalence, n=15, with_ties = T) %>%
  arrange(ISACS, weighted_prevalence) %>%
  mutate(Drug = factor(Drug, Drug),
         ISACS = case_when(ISACS == 3~"High",
                           ISACS == 2~"Moderate",
                           ISACS == 1~"Low",
                           ISACS == 0~"No",
                           T~NA),
         ISACS = factor(ISACS, levels = c("High", "Moderate",
                                          "Low", "No"))) %>%
  ggplot(aes(x = weighted_prevalence, y = Drug, fill = ISACS)) +
  geom_bar(stat = 'identity', color = "black", width = 1) +
  scale_fill_manual(values = rev(pal_jama()(4))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_percent()) +
  scale_y_discrete(expand = expansion(mult = c(0.025, 0.025))) +
  labs(title = "Top 15 most commonly prescribed medicines in NHANES stratified by ISACS",
       # subtitle = "60 years or older",
       x = "Survey-weighted prevalence", y = "", fill = "") +
  mytheme +
  theme(legend.position = c(0.99, 0.96),
        legend.background = element_blank(),
        legend.justification = c(1, 0.75))

top15 <- RXQ_RX_df %>%
  left_join(., phenoage.dat %>% select(SEQN, age), "SEQN") %>% filter(age >= 60) %>%
  filter(RXDDRUG != "") %>%
  left_join(., ISACS.NHANES, "RXDDRUG") %>%
  group_by(ISACS, Drug) %>% count %>%
  drop_na(ISACS) %>%
  ungroup %>% group_by(ISACS) %>%
  slice_max(n , n=15, with_ties = T) %>%
  arrange(ISACS, n) %>%
  mutate(n = n/sum(NHANES_drugs$inAnalysis, na.rm = T),
         Drug = factor(Drug, Drug),
         ISACS = case_when(ISACS == 3~"High",
                           ISACS == 2~"Moderate",
                           ISACS == 1~"Low",
                           ISACS == 0~"No",
                           T~NA),
         ISACS = factor(ISACS, levels = c("High", "Moderate",
                                          "Low", "No"))) %>%
  ggplot(aes(x = n, y = Drug, fill = ISACS)) +
  geom_bar(stat = 'identity', color = "black", width = 1) +
  scale_fill_manual(values = rev(pal_jama()(4))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_percent()) +
  scale_y_discrete(expand = expansion(mult = c(0.025, 0.025))) +
  labs(title = "Top 15 most commonly prescribed\nmedicines in NHANES stratified by ISACS",
       subtitle = "60 years or older",
       x = "Prevalence (crude)", y = "", fill = "") +
  mytheme +
  theme(legend.position = c(0.99, 0.98),
        legend.background = element_blank(),
        legend.justification = c(1, 0.75))

# tiff("Results/NHANES_top15.tiff", units="in", height=9.5, width=11, res=300)
plot_grid(top15, top15_weighted, nrow = 1)
# dev.off()

# Figure for paper
top15_weighted_inverted <- top15_weighted + 
  coord_flip() +
  scale_y_discrete(limits = rev, expand = expansion(mult = c(0.015, 0.015))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tiff("Results/Figure6.tiff", units="in", height=20, width=13.5, res=300)
plot_grid(top15_weighted_inverted,
          NHANES.forest.BA + labs(title = "Survey-weighted association with biological age") + theme(plot.subtitle = element_blank()),
          NHANES.forest.cog + labs(title = "Survey-weighted association with cognitive function") + theme(plot.subtitle = element_blank()),
          ncol = 1, 
          labels = c("(a)", "(b)", "(c)"),
          rel_heights = c(0.75,1,2)
          )
dev.off()

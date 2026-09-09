# tAge visualization
# Kevin Winardi

# Figure 5 requires safety metrics plots
if (!exists("neurotox")) {
  source("6_4_RedefiningAnticholinergics_CompareScales.R")
}

tAge.df <- read.csv("12_DrugMatrix/tAge_result/ALL_tAge_results.csv") %>%
  mutate(treatment = str_to_sentence(treatment)) %>%
  left_join(., read.csv("12_DrugMatrix/standardized_drug_names.csv") %>% 
              drop_na, c("treatment")) %>%
  mutate(Drug = case_when(!is.na(standardized_drugbank_name)~standardized_drugbank_name,
                          T~str_to_sentence(treatment)),
         array_type = case_when(gse_id %in% c("GSE57800", "GSE57811", "GSE57815", "GSE57816")~"Affymetrix",
                                TRUE~"CodeLink"),
         tissue.ch1 = gsub("Thigh muscle", "Skeletal muscle", tissue.ch1)) %>%
  dplyr::select(-standardized_drugbank_name)

tAge.df %>% 
  distinct(Drug, treatment, vehicle.ch1) %>%
  left_join(., result.w.anno.full, "Drug") %>%
  dplyr::filter(is.na(cholinergicCluster)) %>% 
  pull(Drug) %>% unique

tAge.df %>% 
  mutate(shape.var = case_when(treatment == "Control"~"Control",
                               TRUE~"Drug")) %>%
  ggplot(aes(x = scaled_diff_BR_tAge, y = yugene_diff_BR_tAge)) +
  geom_point(aes(shape = shape.var, color = platform_id), 
             size = 2, alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~tissue.ch1, scales = "free")

tAge.col <- c("High" = "#B24745FF" ,
              "Moderate" = "#00A1D5FF", 
              "Low" = "#DF8F44FF", 
              "No" = "#374E55FF", 
              "Control" = "#79AF97FF")

tAge.clean <- tAge.df %>%
  left_join(., result.w.anno.full, "Drug") %>% 
  mutate(cholinergicCluster_control = case_when(treatment == "Control"~"Control",
                                                TRUE~cholinergicCluster),
         time.numeric = as.numeric(gsub("d", "", time.ch1)),
         treatment_full = relevel(factor(case_when(treatment == "Control"~"Control",
                                                   TRUE~paste0(Drug, " (", dose.ch1, ")"))),
                                  ref = "Control"),
         label = paste0("Organ: ", tissue.ch1, "; Array: ", array_type,
                        "\nRoute: ", route.ch1, "; Vehicle: ", vehicle.ch1)
         
  ) %>%
  drop_na(cholinergicCluster_control) %>%
  mutate(cholinergicCluster_control = factor(cholinergicCluster_control, 
                                             levels = c("High",
                                                        "Moderate", 
                                                        "Low", 
                                                        "No", 
                                                        "Control")))

# Filter before analysis
match_keys <- c("vehicle.ch1", "gse_id", "route.ch1", "tissue.ch1", "time.ch1", "array_type")

sample_count <- tAge.clean %>%
  count(treatment_full, vehicle.ch1, route.ch1, gse_id, tissue.ch1, time.ch1, array_type) %>%
  dplyr::filter(n >= 3)

drug_strata <- sample_count %>%
  dplyr::filter(treatment_full != "Control") %>%
  distinct(across(all_of(match_keys)))

control_strata <- sample_count %>%
  dplyr::filter(treatment_full == "Control") %>%
  distinct(across(all_of(match_keys)))

valid_strata <- inner_join(drug_strata, control_strata, by = match_keys)

tAge.clean <- tAge.clean %>%
  semi_join(sample_count, by = c("treatment_full", match_keys)) %>%
  semi_join(valid_strata, by = match_keys) 

length(unique(tAge.clean$Drug))-1 # Control
length(unique(tAge.clean$route.ch1))
length(unique(tAge.clean$tissue.ch1))
length(unique(tAge.clean$time.ch1))
tAge.clean %>%
  dplyr::filter(treatment_full != "Control") %>% 
  distinct(Drug, tissue.ch1) %>% group_by(tissue.ch1) %>% count() %>% arrange(n)

# Descriptive statistics (Alluvial)
library(ggalluvial)
desc.stat <- tAge.clean %>%
  dplyr::filter(treatment_full != "Control") %>%
  distinct(Drug, dose.ch1, array_type, tissue.ch1, vehicle.ch1, route.ch1, time.ch1) %>%
  group_by(Drug, array_type, tissue.ch1, vehicle.ch1, route.ch1) %>%
  mutate(nDose = n_distinct(dose.ch1),
         nTime = n_distinct(time.ch1)) %>% ungroup

upset.df <- desc.stat %>% 
  distinct(Drug, tissue.ch1) %>%
  mutate(x = 1) %>% pivot_wider(names_from = "tissue.ch1", values_from = "x", values_fill = 0)

desc.stat.plot1 <- ComplexUpset::upset(upset.df, names(upset.df)[-1], name = "",
                                       base_annotations=list('Number of\nintersections'=(
                                         intersection_size(width=0.75, bar_number_threshold=1) +
                                           scale_y_continuous(expand=expansion(mult=c(0, 0.15))) +
                                           theme(panel.grid.major=element_blank(),
                                                 panel.grid.minor=element_blank(),
                                                 axis.line=element_line(colour='black')))),
                                       matrix=intersection_matrix(geom=geom_point(shape='circle filled', size=3.5, stroke=0.75)),
                                       set_sizes=(upset_set_size(geom=geom_bar(width = 0.5)) +
                                                    theme(axis.line.x=element_line(colour='black'),
                                                          axis.ticks.x=element_line(),
                                                          axis.title.x = element_text(size = 10)) + 
                                                    ylab('Number of drugs')),
                                       queries=list(upset_query(set="Brain", fill="#7B5EA7"),
                                                    upset_query(set="Heart", fill="#A02C3A"),
                                                    upset_query(set="Intestine", fill="#5A9E6F"),
                                                    upset_query(set="Kidney", fill="#2E6EA6"),
                                                    upset_query(set="Bone marrow", fill="#D4A72C"),
                                                    upset_query(set="Spleen", fill="#3A8FA8"),
                                                    upset_query(set="Liver", fill="#C45C2E"),
                                                    upset_query(set="Skeletal muscle", fill="#C47A8A")
                                       ),
                                       width_ratio = 0.2, height_ratio = 1, 
                                       mode = 'distinct', sort_sets="descending",
                                       sort_intersections='descending')

desc_lode <- desc.stat %>%
  group_by(Drug) %>%
  dplyr::summarise(
    Array_type = if (n_distinct(array_type) > 1) "Both arrays" else unique(array_type),
    Organ = if (n_distinct(tissue.ch1)  > 1) "Multiple organs" else unique(tissue.ch1),
    Vehicle = if (n_distinct(vehicle.ch1) > 1) "Multiple vehicles" else unique(vehicle.ch1),
    Route = if (n_distinct(route.ch1)   > 1) "Multiple routes" else unique(route.ch1),
    Number_of_doses_tested = max(nDose, na.rm = T),
    Number_of_timepoints_assessed = max(nTime, na.rm = T),
    .groups = "drop"
  ) %>%
  column_to_rownames("Drug") %>% 
  to_lodes_form(.) %>%
  mutate(x = gsub("_", " ", x)) %>% 
  mutate(x = gsub("of ", "of\n", x)) %>%
  mutate(x = factor(x, 
                    levels = c("Array type", 
                               "Organ", 
                               "Vehicle", 
                               "Route", 
                               "Number of\ndoses tested",
                               "Number of\ntimepoints assessed")))

Label_lode <- data.frame(table(desc_lode$stratum, desc_lode$x)) %>%
  dplyr::filter(Freq > 0) %>%
  mutate(Perc = round(Freq/(sum(Freq)/nlevels(desc_lode$x)),4)*100) %>%
  mutate(Perc = paste0(sprintf(" (%s", Perc), "%)")) %>%
  mutate(label = case_when(Freq > 25~paste0("**", Var1, "**", "<br>",
                                            Freq, " ", Perc),
                           TRUE~paste0("**", Var1, "**"))
         ) 
desc_lode <- left_join(desc_lode, Label_lode, c(c("stratum" = "Var1"), c("x" = "Var2"))) %>%
  mutate(label = ifelse(x == "Protein", paste0("**",as.character(stratum), "**"), label))

array_colors <- c("CodeLink"    = "#C8CDD0",
                  "Affymetrix"  = "#343E45",
                  "Both arrays" = "#7E8C94")
organ_colors <- c("Brain" = "#7B5EA7",
                  "Heart" = "#A02C3A",
                  "Liver" = "#C45C2E",
                  "Intestine" = "#5A9E6F",
                  "Kidney" = "#2E6EA6",
                  "Bone marrow" = "#D4A72C",
                  "Spleen" = "#3A8FA8",
                  "Multiple organs" = "#8C9E98")
vehicle_colors <- c("CMC" = "#C8894A", "Corn Oil" = "#E8C97A",
                    "Saline" = "#9BC4DD", "Water" = "#B8D4E8",
                    "Multiple vehicles" = "#A89880")
route_colors <- c("Oral Gavage" = "#E8A598", "Subcutaneous" = "#C75B6A",
                  "Intravenous" = "#8B2035", "Multiple routes" = "#9E8A8E")
dose_colors <- c("1" = "#EDE0F5", "2" = "#C9A8E0", "3" = "#9E6DC4",
                 "4" = "#7040A8", "5" = "#421D72")
all_colors <- c(array_colors,organ_colors,vehicle_colors,route_colors,dose_colors)

desc.stat.plot2 <- desc_lode %>% 
  mutate(stratum = factor(stratum, levels = names(all_colors))) %>% 
  ggplot(aes(x = x, stratum = stratum, alluvium = alluvium, 
             fill = stratum)) +
  labs(x = "") + 
  geom_stratum(width = 0.6) + 
  geom_flow(stat = "alluvium", 
            width = 0.6,  knot.pos = 1/3, cement.alluvia = T) +
  scale_y_discrete() + scale_x_discrete(expand = expansion(mult = c(0.15, 0.05)), position = "top") +
  theme_bw() +
  scale_fill_manual(values = all_colors) +
  ggtext::geom_richtext(aes(label = ifelse(Freq > 5, label, NA)),
                        stat = "stratum", size = 3, show.legend = F, 
                        fill = NA, label.color = NA,
                        label.padding = grid::unit(rep(0, 4), "pt")) +
  ggtext::geom_richtext(aes(label = ifelse(Freq < 5, label, NA)),
                        fill = NA, label.color = NA, show.legend = F, 
                        stat = "stratum", size = 3, nudge_x = 0.5) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 25, b = 10, l = 1, unit = "pt"),
        legend.position = "none",
        plot.title = element_text(face = "bold", color = "black", size = 12, hjust = 0.5),
        panel.border = element_blank(),
        axis.text = element_text(face = "bold", color = "black", size = 11),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

tiff(paste0("Results/SFigure16.tiff"), units="in", height=13, width=12, res=300)
plot_grid(desc.stat.plot1,
          desc.stat.plot2,
          ncol = 1, rel_heights = c(0.75, 1.25),
          labels = c("(a)", "(b)"))
dev.off()

analyze_drug_vs_vehicle <- function(drug_name, organ_name, array_name,
                                    vehicle_name, route_name,
                                    full_df, vehicle_label = "Control") {
  
  drug_rows <- full_df %>%
    filter(Drug == drug_name, tissue.ch1 == organ_name, gse_id == array_name,
           vehicle.ch1 == vehicle_name, route.ch1 == route_name
           )
  
  if (nrow(drug_rows) == 0) {
    return(tibble(Drug = drug_name, tissue.ch1 = organ_name, gse_id == array_name, 
                  vehicle.ch1 = vehicle_name, route.ch1 = route_name, 
                  term = NA, estimate = NA, p.value = NA, note = "no drug rows"))
  }
  
  strata <- drug_rows %>% distinct(route.ch1, vehicle.ch1, gse_id, time.ch1)
  
  control_rows <- full_df %>%
    filter(Drug == vehicle_label, tissue.ch1 == organ_name, gse_id == array_name,
           vehicle.ch1 == vehicle_name, route.ch1 == route_name
           ) %>%
    semi_join(strata, by = c("route.ch1", "vehicle.ch1", "gse_id", "time.ch1"))
  
  if (nrow(control_rows) == 0) {
    return(tibble(Drug = drug_name, tissue.ch1 = organ_name, gse_id = array_name, 
                  vehicle.ch1 = vehicle_name, route.ch1 = route_name,
                  term = NA, estimate = NA, p.value = NA, note = "no matched vehicle"))
  }
  
  one_drug_df <- bind_rows(drug_rows, control_rows) %>%
    mutate(time.numeric = as.numeric(gsub("d", "", time.ch1)),
           dose.ch1 = relevel(factor(dose.ch1), ref = "0 mg/kg"))
  
  n_doses <- n_distinct(one_drug_df$dose.ch1)
  n_times <- n_distinct(as.character(one_drug_df$time.numeric))
  
  counts_df <- one_drug_df %>%
    count(dose.ch1, time.numeric) %>%
    mutate(dose.ch1 = paste0("dose.ch1", dose.ch1)) %>%
    group_by(dose.ch1) %>% dplyr::summarize(n_min_treated = min(n), 
                                            n_average_treated = mean(n)) %>%
    ungroup %>%
    dplyr::rename("term"="dose.ch1")
  counts_df$n_min_control = counts_df %>% dplyr::filter(term == "dose.ch10 mg/kg") %>% pull(n_min_treated)
  counts_df$n_average_control = counts_df %>% dplyr::filter(term == "dose.ch10 mg/kg") %>% pull(n_average_treated)
  
  if (n_doses > 1 & n_times > 1) {
    formula_used <- yugene_diff_BR_tAge ~ dose.ch1 + time.numeric
  } else if (n_doses > 1) {
    formula_used <- yugene_diff_BR_tAge ~ dose.ch1
  } else if (n_times > 1) {
    formula_used <- yugene_diff_BR_tAge ~ time.numeric
  } else {
    formula_used <- yugene_diff_BR_tAge ~ 1
  }
  
  model <- tryCatch(lm(formula_used, data = one_drug_df), error = function(e) NULL)
  
  if (is.null(model)) {
    return(tibble(Drug = drug_name, tissue.ch1 = organ_name, gse_id = array_name,
                  vehicle.ch1 = vehicle_name, route.ch1 = route_name,
                  term = NA, estimate = NA, p.value = NA, note = "model failed"))
  }
  
  tidy(model) %>%
    mutate(Drug = drug_name,
           tissue.ch1 = organ_name,
           gse_id = array_name,
           vehicle.ch1 = vehicle_name, 
           route.ch1 = route_name,
           n_doses = n_doses,
           n_times = n_times,
           n_total_treated = nrow(drug_rows),
           n_total_control = nrow(control_rows),
           note = NA_character_) %>%
    dplyr::filter(str_detect(term, "dose.ch1")) %>%
    relocate(Drug, tissue.ch1, gse_id, vehicle.ch1, route.ch1) %>%
    left_join(., counts_df, "term")
}

combos <- tAge.clean %>%
  dplyr::filter(Drug != "Control") %>%
  distinct(Drug, tissue.ch1, gse_id, vehicle.ch1, route.ch1)

stat.results <- pmap_dfr(combos,
                         ~ analyze_drug_vs_vehicle(..1, ..2, ..3, ..4, ..5, tAge.clean)
                         ) %>%
  group_by(tissue.ch1) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         signlogp = sign(estimate)*-log10(p.adj)) %>%
  ungroup() %>%
  left_join(., result.w.anno.full %>% dplyr::select(Drug, cholinergicCluster), "Drug") %>%
  left_join(., tAge.clean %>% distinct(gse_id, array_type), "gse_id") %>%
  mutate(text.label = case_when(p.adj < 0.05 & cholinergicCluster != "No" & estimate > 0 ~ paste0(Drug, " (", gsub("dose.ch1", "", term), ")"),
                                TRUE ~ ""))

# Statistics
## Dose & platform level
stat.results %>% 
  dplyr::filter(p.adj < 0.05) %>%
  count(tissue.ch1, cholinergicCluster) %>% 
  pivot_wider(names_from = "cholinergicCluster", values_from = "n")

## Drug-level: at least one dose-platform combo tested is significant
stat.results %>% 
  dplyr::filter(p.adj < 0.05) %>%
  distinct(Drug, tissue.ch1, cholinergicCluster) %>%
  count(tissue.ch1, cholinergicCluster) %>% 
  pivot_wider(names_from = "cholinergicCluster", values_from = "n")

# 353 datapoints drugs-dose-route-vehicle
tiff(paste0("Results/SFigure17.tiff"), units="in", height=4.75, width=5.5, res=300)
stat.results %>%
  distinct(Drug, term, tissue.ch1, vehicle.ch1, array_type, estimate) %>%
  pivot_wider(names_from = array_type, values_from = estimate) %>%
  drop_na() %>%
  ggplot(aes(x = Affymetrix, y = CodeLink, color = tissue.ch1)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Concordance between platforms", color = "Organ",
       x = "Beta coefficient from Affymetrix",
       y = "Beta coefficient from CodeLink") +
  stat_cor(show.legend = FALSE, method = "spearman") + 
  scale_color_manual(values = organ_colors) +
  mytheme 
dev.off()

# Have to put strips manually
x_levels <- levels(factor(stat.results$tissue.ch1))
n <- length(x_levels)

stripe_df <- data.frame(
  xmin = seq(1, n, by = 2) - 0.5,
  xmax = seq(1, n, by = 2) + 0.5
)

plot.seed <- 88
tAge.summary.plot <- stat.results %>%
  ggplot(aes(x = tissue.ch1, y = signlogp)) +
  geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.5) +
  geom_point(pch = 21, aes(fill = cholinergicCluster),
             position = position_dodge(width = 1, preserve = "single")) +
  geom_boxplot(width = 0.5, aes(fill = cholinergicCluster), 
               position = position_dodge(width = 1, preserve = "single"),
               outliers = FALSE, alpha = 0.5) +
  # geom_stripped_cols() +
  facet_grid(array_type~., scales = "fixed") +
  geom_hline(linetype = "dashed", yintercept = log10(0.05)) +
  geom_hline(linetype = "dashed", yintercept = -log10(0.05)) +
  labs(x = "", y = "Signed log10(FDR-adjusted p-value)\n← Anti-aging     Pro-aging →",
       fill = "ISACS") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_jama() + 
  scale_color_jama() +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Stacked bar
## Descriptive statistics
stat.results %>% 
  dplyr::filter(p.adj < 0.05) %>%
  mutate(age_effect = case_when(estimate > 0~"Accelerated age",
                                T~"Decelerated age")) %>% 
  group_by(cholinergicCluster, age_effect, tissue.ch1, array_type) %>%
  count() %>%
  pivot_wider(names_from = "cholinergicCluster", values_from = "n")

stat.results %>% ggplot(aes(x = cholinergicCluster, y = estimate, fill = cholinergicCluster)) +
  geom_violin(alpha = 0.5) + geom_boxplot(width = 0.25) +
  scale_fill_jama() +
  facet_wrap(tissue.ch1~array_type) +
  mytheme

## Version 1
tAge.summary.plot.sig <- stat.results %>% 
  dplyr::filter(p.adj < 0.05) %>%
  mutate(age_effect = case_when(estimate > 0~"Accelerated age",
                                T~"Decelerated age")) %>%
  ggplot(aes(x = age_effect, fill = fct_rev(cholinergicCluster))) +
  geom_bar(position = "stack", width= 0.5, color = "black") +
  scale_fill_manual(values = rev(pal_jama()(4))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  facet_grid(array_type~tissue.ch1#, scales = "free_x"
  ) +
  labs(x = "", y = "Number of significant associations", fill = "ISACS") +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))

## Version 2
stat.results %>% 
  dplyr::filter(p.adj < 0.05) %>%
  mutate(age_effect = case_when(estimate > 0~"Accelerated age",
                                T~"Decelerated age")) %>%
  ggplot(aes(x = tissue.ch1, fill = fct_rev(cholinergicCluster))) +
  geom_bar(position = "stack", width= 0.5) +
  scale_fill_manual(values = rev(pal_jama()(4))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  facet_grid(array_type~age_effect, scales = "free_x"
  ) +
  labs(x = "", y = "Number of significant associations", fill = "ISACS") +
  mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Volcano 
vol.plot <- stat.results %>%
  arrange(cholinergicCluster) %>%
  # dplyr::filter(tissue.ch1 == "Liver") %>%
  ggplot(aes(x = estimate, y = -log10(p.adj), color = cholinergicCluster)) +
  geom_point(size = 2, alpha = 0.75) +
  # geom_text_repel(aes(label = text.label), show.legend = F) +
  facet_wrap(tissue.ch1~array_type, scales = "free") +
  geom_hline(linetype = "dashed", yintercept = -log10(0.05)) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_color_jama() +
  labs(color = "ISACS", y = "-log10(FDR-adjusted p-value)", x = "β estimate") +
  mytheme 

# tiff(paste0("Results/tAge1.1.tiff"), units="in", height=5, width=12, res=300)
tAge.summary.plot
# dev.off()

stat.results.mat <- stat.results %>%
  dplyr::filter(p.adj < 0.05) %>%
  dplyr::filter(cholinergicCluster != "No") %>%
  dplyr::select(Drug, term, tissue.ch1, array_type, cholinergicCluster, signlogp) %>%
  mutate(text.label = paste0(Drug, " (", gsub("dose.ch1", "", term), ")"),
         organ.array = paste0(tissue.ch1, " (", array_type, ")")) %>%
  dplyr::select(-Drug, -term, -tissue.ch1, -array_type) %>%
  arrange(organ.array) %>%
  pivot_wider(names_from = "organ.array", values_from = "signlogp") %>%
  mutate(cholinergicCluster = fct_rev(cholinergicCluster)) %>%
  arrange(cholinergicCluster, text.label) %>% 
  column_to_rownames("text.label")

stat.results.mat.anno <- stat.results.mat %>%
  dplyr::select(cholinergicCluster) %>%
  dplyr::rename(ISACS = cholinergicCluster)

stat.results.mat <- stat.results.mat %>%
  dplyr::select(-cholinergicCluster)

ha <- rowAnnotation(df = stat.results.mat.anno,
                    annotation_name_rot = 90, 
                    annotation_name_gp = gpar(fontsize = 10),
                    show_legend = T,
                    col = list("ISACS" = setNames(pal_jama("default")(4),
                                                  c("No", "Low", "Moderate", "High"))))
col_fun = circlize::colorRamp2(c(-15, -5, -1.3, 0, 1.3, 5, 15),
                               c("#002c58", "#3a6488", "#a8c4dc", "white", "#dcb3b3", "#a83838", "#580000"))

tAge.ht <- ComplexHeatmap::Heatmap(as.matrix(stat.results.mat), col = col_fun,
                                   cluster_rows = F, cluster_columns = F,
                                   heatmap_legend_param = list(title = "Signed log10(FDR-adjusted p-value)"),
                                   # name = "Signed log10(FDR-adjusted p-value)",
                                   right_annotation = ha)

# tiff(paste0("Results/tAge1.1.tiff"), units="in", height=9.5, width=7.5, res=300)
tAge.ht
# dev.off()

# Plot drug of interest
drug.of.interest <- "Chlorpromazine"  

drug_rows <- tAge.clean %>%
  dplyr::filter(Drug == drug.of.interest)

strata <- drug_rows %>%
  distinct(route.ch1, tissue.ch1, time.ch1, vehicle.ch1, array_type)

control_rows <- tAge.clean %>%
  dplyr::filter(treatment %in% "Control") %>%  
  semi_join(strata, by = c("route.ch1", "tissue.ch1", "time.ch1", "vehicle.ch1", "array_type"))

tAge.DOI.1 <- bind_rows(drug_rows, control_rows) %>%
  mutate(dose.label = case_when(treatment_full == "Control"~"Control",
                                TRUE~dose.ch1)) %>%
  mutate(dose.label = relevel(factor(dose.label), ref = "Control")) %>%
  ggplot(aes(x = time.numeric, y = yugene_diff_BR_tAge, 
             group = dose.label, color = dose.label)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3.75) +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~label) + 
  scale_x_continuous(breaks = unique(control_rows$time.numeric), labels = unique(control_rows$time.numeric)) +
  scale_color_igv() +
  labs(title = drug.of.interest, 
       color = "", x = "Duration of exposure (days)",
       y = "Relative tAge") +
  mytheme +
  theme(legend.position = "bottom")

# tiff(paste0("Results/tAge2.tiff"), units="in", height=5.5, width=n_distinct(control_rows$label)*6.5, res=300)
tAge.DOI.1
# dev.off()

drug.of.interest <- "Promethazine"  

drug_rows <- tAge.clean %>%
  dplyr::filter(Drug == drug.of.interest, array_type == "Affymetrix") 

strata <- drug_rows %>%
  distinct(route.ch1, tissue.ch1, time.ch1, vehicle.ch1, array_type)

control_rows <- tAge.clean %>%
  dplyr::filter(treatment %in% "Control") %>%  
  semi_join(strata, by = c("route.ch1", "tissue.ch1", "time.ch1", "vehicle.ch1", "array_type"))

tAge.DOI.2 <- bind_rows(drug_rows, control_rows) %>%
  mutate(dose.label = case_when(treatment_full == "Control"~"Control",
                                TRUE~dose.ch1)) %>%
  mutate(dose.label = relevel(factor(dose.label, 
                                     levels = c("2.3 mg/kg", "113 mg/kg", "Control")), ref = "Control")) %>%
  ggplot(aes(x = time.numeric, y = yugene_diff_BR_tAge, 
             group = dose.label, color = dose.label)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3.75) +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~label) + 
  scale_x_continuous(breaks = unique(control_rows$time.numeric), labels = unique(control_rows$time.numeric)) +
  scale_color_igv() +
  labs(title = drug.of.interest, 
       color = "", x = "Duration of exposure (days)",
       y = "Relative tAge") +
  mytheme +
  theme(legend.position = "bottom")
tAge.DOI.2

drug.of.interest <- "Amoxapine"  

drug_rows <- tAge.clean %>%
  dplyr::filter(Drug == drug.of.interest, tissue.ch1 == "Liver") 

strata <- drug_rows %>%
  distinct(route.ch1, tissue.ch1, time.ch1, vehicle.ch1, array_type)

control_rows <- tAge.clean %>%
  dplyr::filter(treatment %in% "Control") %>%  
  semi_join(strata, by = c("route.ch1", "tissue.ch1", "time.ch1", "vehicle.ch1", "array_type"))

tAge.DOI.3 <- bind_rows(drug_rows, control_rows) %>%
  mutate(dose.label = case_when(treatment_full == "Control"~"Control",
                                TRUE~dose.ch1)) %>%
  mutate(dose.label = relevel(factor(dose.label), ref = "Control")) %>%
  ggplot(aes(x = time.numeric, y = yugene_diff_BR_tAge, 
             group = dose.label, color = dose.label)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3.75) +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~label) + 
  scale_x_continuous(breaks = unique(control_rows$time.numeric), labels = unique(control_rows$time.numeric)) +
  scale_color_igv() +
  labs(title = drug.of.interest, 
       color = "", x = "Duration of exposure (days)",
       y = "Relative tAge") +
  mytheme +
  theme(legend.position = "bottom")
tAge.DOI.3

drug.of.interest <- "Alprazolam"  

drug_rows <- tAge.clean %>%
  dplyr::filter(Drug == drug.of.interest) 

strata <- drug_rows %>%
  distinct(route.ch1, tissue.ch1, time.ch1, vehicle.ch1, array_type)

control_rows <- tAge.clean %>%
  dplyr::filter(treatment %in% "Control") %>%  
  semi_join(strata, by = c("route.ch1", "tissue.ch1", "time.ch1", "vehicle.ch1", "array_type"))

tAge.DOI.4 <- bind_rows(drug_rows, control_rows) %>%
  mutate(dose.label = case_when(treatment_full == "Control"~"Control",
                                TRUE~dose.ch1)) %>%
  mutate(dose.label = relevel(factor(dose.label), ref = "Control")) %>%
  ggplot(aes(x = time.numeric, y = yugene_diff_BR_tAge, 
             group = dose.label, color = dose.label)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3.75) +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~label) + 
  scale_x_continuous(breaks = unique(control_rows$time.numeric), labels = unique(control_rows$time.numeric)) +
  scale_color_igv() +
  labs(title = drug.of.interest, 
       color = "", x = "Duration of exposure (days)",
       y = "Relative tAge") +
  mytheme +
  theme(legend.position = "bottom")
tAge.DOI.4
     
BR.tAge <- ggdraw() + 
  draw_image("Biorender/tAgeAnalysis.png") + 
  theme(plot.margin = margin(t=0.15, l=0.15, r=0.15, b=0.15, unit = "cm"))

# Plot all as figure 5
tiff(paste0("Results/Figure5.tiff"), units="in", height=20, width=16, res=300)
plot_grid(
  # Safety parameters
  plot_grid(unbound.plot + theme(axis.text.x = element_text(angle = 25, hjust = 1)), 
            neurotox + theme(axis.text.x = element_text(angle = 25, hjust = 1)), 
            nrow = 1, 
            labels = c("(a)", "(b)"),
            rel_widths = c(1.35, 0.65),
            align = "h", axis = "tb"),
  # tAge analysis 1: big picture 
  plot_grid(plot_grid(BR.tAge, vol.plot,
                      labels = c("(c)", "(d)"),
                      ncol = 1, rel_heights = c(1,1.75)), NULL,
            grid.grabExpr(draw(tAge.ht, newpage = FALSE)),
            nrow = 1, labels = c("", "(e)", ""),
            rel_widths = c(1.4,0.05, 1)),
  # tAge analysis 2: specific examples 
  plot_grid(tAge.DOI.1, tAge.DOI.2, tAge.DOI.3, tAge.DOI.4,
            nrow = 1, 
            labels = c("(f)", ""),
            align = "h", axis = "tb"),
  ncol = 1, rel_heights = c(5,7.25,3)
  )
dev.off()

tiff(paste0("Results/SFigure18.tiff"), units="in", height=10, width=9, res=300)
plot_grid(tAge.summary.plot + 
            theme(plot.margin = margin(t=0.1, l=0.5, r=0.1, b=0.1, unit = "cm")),
          tAge.summary.plot.sig + 
            theme(plot.margin = margin(t=0.1, l=0.5, r=0.1, b=0.1, unit = "cm")),
          rel_heights = c(1, 1.25), ncol = 1, 
          labels = c("(a)", "(b)"))
dev.off()

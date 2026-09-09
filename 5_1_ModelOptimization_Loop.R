# Model optimization
source("4_NetworkSetup.R")
set.seed(123)

# Ground truth drugs
prototypical.anticholinergic <- c("Scopolamine", "Butylscopolamine",  
                                  "Methscopolamine", "Methscopolamine bromide", 
                                  "Homatropine", "Homatropine methylbromide",
                                  "Atropine", "Ipratropium", "Tiotropium", "Tropicamide",
                                  "Homatropine", "Homatropine methylbromide",
                                  "Pirenzepine")
OAB.anticholinergic <- c("Oxybutynin", "Tolterodine",
                         "Solifenacin", "Darifenacin", "Fesoterodine", "Trospium")
TCA <- c("Nortriptyline", "Amitriptyline", "Imipramine", "Doxepin", "Dosulepin",
         "Amoxapine", "Clomipramine", "Protriptyline")
Sedative.antihistamine <- c("Chlorpheniramine", "Brompheniramine", "Diphenhydramine",
                            "Doxylamine", "Promethazine", "Hydroxyzine", "Triprolidine")
Bisphosphonate <- c("Alendronic acid", "Risedronic acid", "Pamidronic acid", "Etidronic acid")
Statin <- c("Atorvastatin", "Fluvastatin", "Simvastatin", "Lovastatin",
            "Rosuvastatin", "Pravastatin")
ACEi <- c("Enalapril", "Benazepril", "Captopril", "Lisinopril", "Fosinopril",
          "Trandolapril", "Moexipril")
DOAC <- c("Dabigatran", "Apixaban", "Betrixaban", "Edoxaban", "Rivaroxaban")
AChEi <- c("Donepezil", "Galantamine", "Edrophonium",
           "Rivastigmine", "Pyridostigmine", "Neostigmine", "Physostigmine")
Organophosphate <- c("Malathion", "Isoflurophate", "Paraoxon")
Cholinergic.Agonist <- c("Aceclidine", "Acetylcholine", "Carbamoylcholine",
                         "Arecoline", "Bethanechol", "Cevimeline", "Carbamoylcholine",
                         "Pilocarpine", "Xanomeline", "Methacholine", "Bethanechol")
true_pos <- c(prototypical.anticholinergic, OAB.anticholinergic, TCA, Sedative.antihistamine)
true_neg <- c(Bisphosphonate, Statin, ACEi, DOAC)
false_pos <- c(AChEi, Organophosphate, Cholinergic.Agonist)

drug.List <- unique(c(true_pos,true_neg,false_pos))

# Different cholinergic module
## CAP
p1 <- c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", 
        "ACHE", "CHAT", "SLC5A7", "SLC18A3")

# CAP + downstream
pathway_info <- keggGet("hsa04725")
gene_entries <- pathway_info[[1]]$GENE
gene_symbols <- sapply(gene_entries[seq(2, length(gene_entries), by = 2)],
                       function(x) {strsplit(x, ";")[[1]][1]}) %>%
  unname() %>% sort()
gene_symbols <- gene_symbols[!grepl("CHRN|CACN|JAK2|FYN|PIK|
                                    P3R3URF|AKT|BCL|CAMK",
                                    gene_symbols)] # Remove nAChR & downstream
p2 <- c(gene_symbols,
        "PIK3R6", "PIK3R5", "PIK3CG",
        "CACNA1A", "CACNA1B")

# CAP + effector
p3 <- c(p1,
        # Delirium (CXCL8 = IL8; IL1RN = IL1RA; CCL2 = MCP1)
        # n = 11
        "IL1B", "IL1RN", "IL6", " CXCL8", "IL10", "S100B", "TNF", 
        "CCL2", "GFAP", "IGF1",  "NEFL") 

# CAP + downstream + effector
p4 <- c(p2,
        # Delirium (CXCL8 = IL8; IL1RN = IL1RA; CCL2 = MCP1)
        # n = 11
        "IL1B", "IL1RN", "IL6", " CXCL8", "IL10", "S100B", "TNF", 
        "CCL2", "GFAP", "IGF1",  "NEFL") 

p.upset <- list(p1 = p1,
                p2 = p2,
                p3 = p3,
                p4 = p4) %>%
  imap_dfr(~ tibble(GENE = .x, list = .y)) %>%
  distinct() %>%                
  mutate(value = 1) %>%
  pivot_wider(names_from  = list,
              values_from = value,
              values_fill = 0) %>%
  column_to_rownames("GENE")

# Different weights
## None
g1 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                     Protein_B_SYMBOL), 
                           directed = F)
g1 <- igraph::simplify(g1)

## Hyperbola
g2 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g2)$weight <- PPI.df2$Evidence.weight
g2 <- igraph::simplify(g2)
E(g2)$weight <- 1 + (1/E(g2)$weight)

## Linear
g3 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g3)$weight <- PPI.df2$Evidence.weight 
g3 <- igraph::simplify(g3)
E(g3)$weight <- 1 + (20 - E(g3)$weight) / 20

## Linear piecewise
g4 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g4)$weight <- PPI.df2$Evidence.weight 
g4 <- igraph::simplify(g4)
E(g4)$weight <- ifelse(E(g4)$weight > 5, 1, 1 + (5 - E(g4)$weight)/5)

# Add names
graph.list <- list("None" = g1,
                   "Hyperbola" = g2,
                   "Linear" = g3,
                   "Linear piecewise" = g4)
protein.list <- list("CAP" = p1,
                     "CAP + downstream" = p2,
                     "CAP + effector" = p3,
                     "CAP + downstream + effector" = p4)

# Run analysis 
distance.metric <- c("Closest_b", "Closest_ab", "Shortest", "Kernel_b", "Kernel_ab", "Centre")

for (gname in names(graph.list)) {
  for (pname in names(protein.list)) {
    
    cat("Running: Graph weight =", gname, "& Pathway =", pname, "\n")
    result.w.sig <- compute_distance_significance(g = graph.list[[gname]],
                                                  drug_targets = drug.List,
                                                  pathway_proteins = protein.list[[pname]],
                                                  randomize_drugs = T,
                                                  distance_param = distance.metric,
                                                  n_perm = 1000,
                                                  seed = 123,
                                                  n_cores = n_cores)
    result.w.anno <- result.w.sig %>%
      pivot_longer(cols = -c(Drug, DrugTarget_degree),
                   names_to = c("distance.metric", "stat", "X"),
                   names_pattern = "(.+)_(pval|zscore)|(.+)") %>%
      mutate(distance.metric = paste0(distance.metric, X),
             stat = ifelse(stat == "", "distance.value", stat)) %>% dplyr::select(-X) %>%
      pivot_wider(values_from = "value", names_from = "stat") %>% unique %>%
      mutate(zscore_sig = factor(case_when(zscore < -1.96~"Significant",
                                           T~"NS"),
                                 levels = c("Significant", "NS")),
             Ground_truth = case_when(Drug %in% true_pos~"True positive",
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
                                              "False positive")))
    
    plot.A <- result.w.anno %>% 
      ggplot(aes(x = Ground_truth, y = distance.value)) +
      geom_violin() +
      geom_boxplot(outliers = F, alpha = 0.3, width = 0.5) +
      geom_point(size = 3, alpha = 0.65, shape = 21, color = "white",
                 position = position_jitterdodge(jitter.width = 0.36,
                                                 dodge.width = 0.85),
                 aes(fill = Drug_category)) + # , shape = zscore_sig
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
      ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", #"p.format",
                                 comparison = list(c("True negative", "True positive"),
                                                   c("False positive", "True positive"),
                                                   c("False positive", "True negative")
                                 ),
                                 step.increase = 0.2) +
      facet_grid(~distance.metric, scales = "free_y") +
      labs(x = "", y = "Proximity", fill = "Drug class",
           title = "Comparing distance values") +
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
      mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    
    plot.B <- result.w.anno %>%
      ggplot(aes(x = Ground_truth, fill = zscore_sig)) +
      geom_bar(position="fill", stat = "count", color = "black") + 
      scale_y_continuous(labels = scales::percent,
                         expand = expansion(mult = c(0, 0))) +
      scale_fill_manual(values = c("NS" = "#4DAA57",
                                   "Significant" = "#F26157")) +
      facet_grid(~distance.metric) +
      labs(x = "", y = "Proportion", fill = "Permutation test (n = 1000)",
           title = "Comparing statistical significance") + 
      mytheme + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    
    # Save result CSV
    csv.name <- paste0("Results/Optimization/", gname, "_", pname, "_res.csv")
    write.csv(result.w.anno, file = csv.name, row.names = F)
    
    subtitle <- paste("Graph weight =", gname, "& Protein =", pname)
    
    # Save plot
    final.plot <- plot_grid(plot.A, plot.B, nrow = 2, 
                            axis = "lr", align = "v", rel_heights = c(1, 0.8))
    fig.name <- paste0("Results/Optimization/", gname, "_", pname, "_plot.png")
    png(fig.name, width = 13, height = 7, units = "in", res = 600)
    print(final.plot)
    grid::grid.text(subtitle, x = 0.5, y = 0.02, gp = grid::gpar(fontsize = 10))
    dev.off()
  }
}

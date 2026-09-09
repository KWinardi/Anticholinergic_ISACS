# Calculate distance and significance
## After determining the best model, run for all ~7,000 drugs
source("4_NetworkSetup.R")

set.seed(123)
detectCores()
n_cores <- max(1, parallel::detectCores() - 2)

# Graph setup
g_final <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                           Protein_B_SYMBOL), 
                                 directed = F)
E(g_final)$weight <- PPI.df2$Evidence.weight
g_final <- igraph::simplify(g_final)
E(g_final)$weight <- 1 + (1/E(g_final)$weight)

# Anticholinergic module
p_final <- c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", 
             "ACHE", "CHAT", "SLC5A7", "SLC18A3")
  
# All drugs
drug.List <- unique(merged.PDI.df2$Protein_A_SYMBOL)
length(drug.List)

# Distance
distance.metric <- c("Closest_b", "Closest_ab", "Shortest", "Kernel_b", "Kernel_ab", "Centre")

# Run analysis
## Split into chunks 
chunk_size <- 2000

chunks <- split(drug.List, ceiling(seq_along(drug.List)/chunk_size))

for (i in seq_along(chunks)) {
  result <- compute_distance_significance(g = g_final,
                                          drug_targets = chunks[[i]],
                                          pathway_proteins = p_final,
                                          randomize_drugs = TRUE,
                                          distance_param = distance.metric,
                                          n_perm = 1000,
                                          seed = 123,
                                          n_cores = n_cores)
  end_index <- max(seq_along(chunks[[i]])) + (i - 1) * chunk_size
  file_name <- paste0("3_NetworkResults/AllDrugs.", end_index, ".csv")
  write.csv(result, file_name, row.names = FALSE)
  message("Saved ", file_name)
}
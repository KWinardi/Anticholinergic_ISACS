# tAge analysis on DrugMatrix dataset
# Kevin Winardi
library(tAge)
library(ggplot2)
library(tidyverse)
library(GEOquery)
options(timeout = 1200)

# Adjust function for microarrays
Sys.setenv(RETICULATE_PYTHON = paste0(getwd(), "/.venv/Scripts/python.exe"))

align_to_gene_list <- function(eset, gene_list) {
  gene_list <- as.character(gene_list)
  current <- rownames(eset)
  
  missing_genes <- setdiff(gene_list, current)
  
  if (length(missing_genes) > 0) {
    n_samples <- ncol(eset)
    zero_matrix <- matrix(NA_real_, nrow = length(missing_genes), ncol = n_samples,
                          dimnames = list(missing_genes, colnames(eset)))
    
    zero_fdata <- data.frame(row.names = missing_genes)
    
    zero_eset <- Biobase::ExpressionSet(assayData = zero_matrix,
                                        phenoData = Biobase::phenoData(eset),
                                        featureData = methods::new("AnnotatedDataFrame", data = zero_fdata),
                                        experimentData = Biobase::experimentData(eset),
                                        annotation = Biobase::annotation(eset))
    
    eset <- Biobase::combine(eset[intersect(gene_list, current), ], zero_eset)
  }
  
  eset[match(gene_list, rownames(eset)), ]
}

tAge_preprocessing <- function(
    eset,
    species = "mouse",
    gene_mapping_type = "Gene.Symbol",
    verbose = TRUE,
    control_group_column = NULL,
    control_group_label = NULL,
    count_threshold = 4,
    percent_threshold = 20
) {
  gene_list <- load_gene_list()
  
  eset_filtered <- filter_genes(eset, count_threshold = count_threshold, percent_threshold = percent_threshold, verbose = verbose)
  eset_genes_converted <- map_genes(eset_filtered, species, gene_mapping_type, verbose)
  eset_RLE <- eset_genes_converted
  eset_log_transformed <- eset_genes_converted
  
  eset_scaled <- scale_eset(eset_log_transformed, verbose = verbose)
  eset_yugene <- YuGene(eset_scaled, verbose = verbose)
  
  eset_scaled_aligned <- align_to_gene_list(eset_scaled, gene_list)
  eset_yugene_aligned <- align_to_gene_list(eset_yugene, gene_list)
  
  eset_scaled_diff <- control_subtraction(eset_scaled_aligned, column_name = control_group_column, control_label = control_group_label, verbose = verbose)
  eset_yugene_diff <- control_subtraction(eset_yugene_aligned, column_name = control_group_column, control_label = control_group_label, verbose = verbose)
  
  return(list(
    RLE_normalized = eset_RLE,
    log_transformed = eset_log_transformed,
    scaled = eset_scaled,
    scaled_diff = eset_scaled_diff,
    yugene = eset_yugene,
    yugene_diff = eset_yugene_diff
  ))
}

# Analyze tAge from microarray data

# Configurate GSE to the appropriate ML model for each organ
gse.config <- list(
  # Affy
  GSE57800 = list(organ = "Heart",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  GSE57811 = list(organ = "Kidney",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Kidney_WT_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Kidney_WT_yugenediff.pkl"),
  GSE57815 = list(organ = "Liver",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Liver_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Liver_yugenediff.pkl"),
  GSE57816 = list(organ = "Skeletal_muscle",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  # Codelink
  GSE59894 = list(organ = "Bone_marrow",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  GSE59895 = list(organ = "Brain",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Brain_WT_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Brain_WT_yugenediff.pkl"),
  GSE59905 = list(organ = "Heart",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  GSE59907 = list(organ = "Intestine",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  GSE59913 = list(organ = "Kidney",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Kidney_WT_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Kidney_WT_yugenediff.pkl"),
  GSE59923 = list(organ = "Liver",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Liver_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Liver_yugenediff.pkl"),
  GSE59925 = list(organ = "Spleen",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl"),
  GSE59926 = list(organ = "Skeletal_muscle",
                  scaled_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
                  yugene_diff = "12_DrugMatrix/tAge_models/BR_Chronoage_Rodents_Multitissue_yugenediff.pkl")
)

# List all DrugMatrix files (expression data and metadata)
DrugMatrixFiles.list <- list.files("12_DrugMatrix/Processed.Data", pattern = "\\.csv$", full.names = F)

DrugMatrixFiles <- data.frame(Filename = DrugMatrixFiles.list) %>%
  mutate(Data_type = if_else(str_detect(Filename, ".metadata"), "Metadata", "Matrix"),
         has_glp = grepl("_", Filename),
         temp_col = gsub(".metadata.csv", "", Filename),
         temp_col = gsub(".processedMat.csv", "", temp_col),
         gse = ifelse(has_glp, sub("_.*", "", temp_col), temp_col),
         glp = ifelse(has_glp, sub(".*_", "", temp_col), NA),
         fullpath = paste0("12_DrugMatrix/Processed.Data/", Filename)) %>% 
  dplyr::select(-temp_col)

# Perform tAge analysis
results.list <- list()
for (GSE in names(gse.config)) {
  message("\n========== Processing: ", GSE, " ==========\n")
  tryCatch({
    # Load data here
    gse.files <- DrugMatrixFiles %>% dplyr::filter(gse == GSE)
    
    meta.file <- gse.files %>% dplyr::filter(Data_type == "Metadata")
    if (nrow(meta.file) != 1) {
      stop("Expected exactly 1 metadata file for ", GSE, ", found ", nrow(meta.file))
    }
    sample.metadata <- read.csv(meta.file$fullpath, row.names = 1) 
    
    mat.files <- gse.files %>% dplyr::filter(Data_type == "Matrix")
    if (nrow(mat.files) == 0) {
      stop("No matrix files found for ", GSE)
    }
    
    sample.metadata <- sample.metadata %>%
      dplyr::select(-hybridization.date.ch1, -rna.extraction.date.ch1) %>%
      dplyr::filter(Sex.ch1 == "Male") %>%
      mutate(loop_var = paste(tissue.ch1, vehicle.ch1, route.ch1, sep = "_")) # time.ch1, 
    
    
    if (nrow(mat.files) == 1) {
      # Single matrix, straightforward load
      exp.data.gene <- read.csv(mat.files$fullpath[1], row.names = 1)
    } else {
      # Multiple GLP matrices
      mat.list <- lapply(mat.files$fullpath, function(f) {
        read.csv(f, row.names = 1) %>% tibble::rownames_to_column(var = "gene_id")
      })
      
      exp.data.gene <- purrr::reduce(mat.list, dplyr::inner_join, by = "gene_id") %>%
        tibble::column_to_rownames(var = "gene_id")
      
      message("    Merged ", nrow(mat.files), " GLP matrices (",
              paste(mat.files$glp, collapse = ", "), ") -> ",
              nrow(exp.data.gene), " shared genes")
    }
    
    # kNN-imputation if required
    n.missing <- sum(is.na(exp.data.gene))
    if (n.missing > 0) {
      message("    Found ", n.missing, " missing values (",
              round(100 * n.missing / (nrow(exp.data.gene) * ncol(exp.data.gene)), 2),
              "%) — running KNN imputation")
      knn.result <- impute::impute.knn(as.matrix(exp.data.gene), k = 10, rowmax = 0.5, colmax = 0.8)
      exp.data.gene <- as.data.frame(knn.result$data)
      
      message("    Imputation complete")
    } else {
      message("    No missing values detected — skipping imputation")
    }
    
    # tAge estimate nested per organ, vehicle, and ROA
    loop_var.levels <- unique(sample.metadata$loop_var)
    valid.levels <- purrr::keep(loop_var.levels, function(lv) {
      meta.sub <- sample.metadata %>% dplyr::filter(loop_var == lv)
      n.control <- sum(meta.sub$treatment == "Control")
      n.treated <- sum(meta.sub$treatment != "Control")
      if (n.control == 0) {
        message("  Skipping ", lv, " — no control samples found")
        return(FALSE)
      }
      if (n.treated == 0) {
        message("  Skipping ", lv, " — no treated samples found")
        return(FALSE)
      }
      TRUE
    })
    
    gse.results.list <- list()
    
    for (lv in valid.levels) {
      
      message("  -- loop_var: ", lv)
      
      tryCatch({
        meta.sub <- sample.metadata %>% dplyr::filter(loop_var == lv)
        n.control <- sum(meta.sub$treatment == "Control")
        n.treated <- sum(meta.sub$treatment != "Control")
        
        message("    Controls: ", n.control, " | Treated: ", n.treated)
        shared.samples <- intersect(colnames(exp.data.gene), meta.sub$geo_accession)
        exp.sub <- exp.data.gene[, shared.samples]
        meta.sub <- meta.sub %>%
          dplyr::filter(geo_accession %in% shared.samples)
        
        # tAge analysis
        eset.sub <- make_ExpressionSet(exp.sub, meta.sub)
        tAge_eset <- tAge_preprocessing(eset.sub,
                                        species = "rat",
                                        gene_mapping_type = "Gene.Symbol",
                                        control_group_column = "treatment",
                                        count_threshold = 4,
                                        percent_threshold = 20,
                                        control_group_label = "Control")
        
        model_paths <- list(scaled_diff = gse.config[[GSE]]$scaled_diff,
                            yugene_diff = gse.config[[GSE]]$yugene_diff)
        
        res <- predict_tAge(tAge_eset, model_paths, species = "rat", mode = "BR")
        
        gse.results.list[[lv]] <- res
        message("Done: ", lv)
        
      }, error = function(e) {
        message(" ERROR in loop_var ", lv, ": ", e$message)
      })
    }
    
    # Combine all results for a given GSE
    gse.results <- bind_rows(gse.results.list)
    results.list[[GSE]] <- gse.results
    
    write.csv(gse.results,
              paste0("12_DrugMatrix/tAge_result/", GSE, "_tAge.csv"),
              row.names = FALSE)
    
    message("Completed GSE: ", GSE)
    
  }, error = function(e) {
    message("ERROR in GSE ", GSE, ": ", e$message)
  })
}

# Export all result as one .csv file
results.all <- bind_rows(results.list, .id = "gse_id")
write.csv(results.all, "12_DrugMatrix/tAge_result/ALL_tAge_results.csv", row.names = FALSE)

message("\nAll done! Total samples processed: ", nrow(results.all))

# Consolidating DrugMatrix datasets from GEO
# Kevin Winardi
library(tidyverse)
library(GEOquery)
library(limma)
library(ggplot2)
options(timeout = 1200)

# Helper functions
parse_filename <- function(filename) {
  name <- tools::file_path_sans_ext(filename)
  has_glp <- grepl("-", name)
  gse <- ifelse(has_glp, sub("-.*", "", name), name)
  glp <- ifelse(has_glp, sub(".*-", "", name), NA)
  data.frame(gse = gse, glp = glp, stringsAsFactors = FALSE)
}

load_and_map_genes <- function(filepath, gse_metadata_cache) {
  filename <- basename(filepath)
  parsed <- parse_filename(filename)
  GSE <- parsed$gse[1]
  GLP <- parsed$glp[1]
  
  message("  Processing: ", filename)
  
  tryCatch({
    gse_obj <- gse_metadata_cache[[GSE]]
    
    if (length(gse_obj) == 1) {
      # Single platform: straightforward
      gse_element <- gse_obj[[1]]
      message("    Single GLP: ", names(gse_obj)[1])
      
    } else {
      # Multiple platforms
      cache_names <- names(gse_obj)
      expected_name <- gsub("\\.csv$", "_series_matrix.txt.gz", filename)
      match_idx <- which(cache_names == expected_name) 
      
      if (length(match_idx) == 0) {
        stop("Could not match GLP '", GLP, "' to any cache entry. Available: ",
             paste(cache_names, collapse = ", "))
      }
      if (length(match_idx) > 1) {
        warning("Multiple matches for GLP '", GLP, "' — using first match: ", cache_names[match_idx[1]])
        match_idx <- match_idx[1]
      }
      
      gse_element <- gse_obj[[match_idx]]
      message("    Matched: ", cache_names[match_idx])
    }
    
    # Feature data
    feature.data <- gse_element@featureData@data %>%
      mutate(ID = as.character(ID))
    
    # Standardized gene symbol column
    possible_names <- c("Gene Symbol", "GENE_SYMBOL", "Gene_Symbol", 
                        "gene_symbol", "Symbol", "SYMBOL", "GeneSymbol")
    
    match_col <- intersect(possible_names, colnames(feature.data))
    
    if (length(match_col) == 0) {
      stop("Could not find a gene symbol column. Available columns: ",
           paste(colnames(feature.data), collapse = ", "))
    }
    if (length(match_col) > 1) {
      warning("Multiple gene symbol-like columns found (", 
              paste(match_col, collapse = ", "), ") — using first: ", match_col[1])
      match_col <- match_col[1]
    }
    
    feature.data <- feature.data %>% rename(`Gene Symbol` = all_of(match_col))
    
    message("    Gene symbol column: '", match_col, "' → standardised")
    
    # Sample metadata
    sample.meta <- as.data.frame(gse_element@phenoData@data) %>%
      mutate(gse_id = GSE,
             treatment = case_when(is.na(`compound:ch1`) ~ "Control",
                                   TRUE~`compound:ch1`)) %>%
      dplyr::select(geo_accession, gse_id,
                    source_name_ch1, organism_ch1,
                    platform_id, treatment,
                    ends_with(":ch1"))
    
    # Expression data
    exp.data <- read.csv(filepath) %>%
      column_to_rownames("ID_REF") %>%
      as.data.frame() %>%
      rownames_to_column("ID") %>%
      left_join(., feature.data, by = "ID") %>%
      dplyr::filter(`Gene Symbol` != "") %>%
      dplyr::filter(!grepl(" /// ", `Gene Symbol`))
    
    exp.mat <- exp.data %>% dplyr::select(starts_with("GSM"))
    exp.anno <- exp.data %>% dplyr::select(ID, `Gene Symbol`)
    exp.data.gene <- limma::avereps(exp.mat, ID = exp.anno$`Gene Symbol`)
    write.csv(exp.data.gene, paste0(GSE, "_", GLP, ".processedMat.csv"))
    message("    ✓ ", nrow(exp.data.gene), " genes x ", ncol(exp.data.gene), " samples")
    
    list(data = exp.data.gene,
         sample_meta = sample.meta,
         filename = filename,
         gse_id = GSE,
         glp_id = GLP,
         success = TRUE)
    
  }, error = function(e) {
    message("    ✗ ERROR: ", e$message)
    list(filename = filename, gse_id = GSE, success = FALSE)
  })
}

# Identify all normalized matrix files
affy_files <- list.files("12_DrugMatrix/Affymetrix", pattern = "\\.csv$", full.names = TRUE)
codelink_files <- list.files("12_DrugMatrix/Codelink",   pattern = "\\.csv$", full.names = TRUE)

all_file_info <- bind_rows(data.frame(filepath = affy_files,     platform = "Affymetrix", stringsAsFactors = FALSE),
                           data.frame(filepath = codelink_files, platform = "Codelink",   stringsAsFactors = FALSE))

# Fetch GEO metadata for all unique GSEs (once, cached) 
unique_gses <- unique(sapply(all_file_info$filepath, function(f) {
  parse_filename(basename(f))$gse
}))

message("Fetching GEO metadata for ", length(unique_gses), " GSEs...")

destdir <- "geo_cache"
if (!dir.exists(destdir)) dir.create(destdir)

gse_metadata_cache <- list()
metadata_list <- list()

for (GSE in unique_gses) {
  tryCatch({
    gse_obj <- getGEO(GSE, GSEMatrix = TRUE)
    gse_metadata_cache[[GSE]] <- gse_obj
    
    meta <- lapply(gse_obj, function(gse_element) {
      as.data.frame(gse_element@phenoData@data) %>%
        mutate(gse_id = GSE,
               treatment = case_when(is.na(`compound:ch1`)~"Control",
                                     TRUE~`compound:ch1`)) %>%
        dplyr::select(geo_accession, gse_id,
                      source_name_ch1, organism_ch1,
                      platform_id, treatment,
                      ends_with(":ch1"))
      }) %>%
      bind_rows() %>%
      distinct(geo_accession, .keep_all = TRUE)  # deduplicate if sample appears in multiple GLPs
    
    metadata_list[[GSE]] <- meta
    write.csv(meta, paste0(GSE, ".metadata.csv"))
    message("  ✓ ", GSE, " (", length(gse_obj), " GLP(s), ", nrow(meta), " samples)")
    
  }, error = function(e) {
    message("  ✗ ERROR fetching ", GSE, ": ", e$message)
  })
}

sample_metadata_all <- bind_rows(metadata_list)

all_data <- mapply(load_and_map_genes,
                   filepath = all_file_info$filepath,
                   MoreArgs = list(gse_metadata_cache = gse_metadata_cache),
                   SIMPLIFY = FALSE)

# Remove samples that are not going to be analyze further
sample_metadata_filtered <- sample_metadata_all %>%
  dplyr::filter(`Sex:ch1`=="Male") %>%
  group_by(`vehicle:ch1`) %>%
  mutate(has_control   = any(treatment == "Control"),
         has_treatment = any(treatment != "Control")) %>%
  ungroup() %>%
  dplyr::filter(has_control & has_treatment) %>%
  dplyr::select(-has_control, -has_treatment)

write.csv(sample_metadata_filtered, "DrugMatrix_metadata.csv", row.names = FALSE) 

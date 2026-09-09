# DrugBank Integration
# 25/09/2025 
# Kevin Winardi

# Load library
library(tidyverse)

# Load each file and merge per folder
folder.path <- "1_PDI/1_DrugBank_v5.1.13"

subfolder.list <- c("Carrier", "Enzyme", "Target", "Transporter")

rename.cols <- list("all.csv" = c("Name" = "UniProt.Name",
                                  "Drug.IDs" = "Drug.Interactors"), # Format is Old = New
                    "uniprot link.csv" = c("Name" = "Drug.Name"))

for (protein in subfolder.list) {
  message("Processing: ", protein, " files...")
  DrugBank.files <- list.files(path = paste0(folder.path, "/", protein), 
                               full.names = T)
  DrugBank.files <- DrugBank.files[basename(DrugBank.files) != "pharmacologically_active.csv"]
  
  DrugBank.list <- lapply(DrugBank.files, function(f) {
    fname <- basename(f)
    rename_map <- rename.cols[[fname]]
    
    df <- read.csv(f, stringsAsFactors = F) %>%
      mutate(across(everything(), as.character))
    
    if(!is.null(rename_map)){
      df <- df %>% dplyr::rename(!!!setNames(names(rename_map), rename_map))
    }
    return(df)
  })
  
  merged.df <- purrr::reduce(DrugBank.list, full_join) %>%
    mutate(Pharm.Protein = protein)
  assign(paste0(protein, ".df"), merged.df)
}

# Merge all dfs
DrugBank.df <- do.call(rbind, mget(paste0(subfolder.list, ".df")))
rownames(DrugBank.df) <- NULL

# Include more drug identifiers
DrugBank.df <- DrugBank.df %>%
  left_join(., read.csv(paste0(folder.path, "/drug links.csv")), "DrugBank.ID") %>%
  left_join(., read.csv(paste0(folder.path, "/drugbank vocabulary.csv")), "DrugBank.ID")

dim(DrugBank.df)

write.csv(DrugBank.df, "DrugBank.df.csv")
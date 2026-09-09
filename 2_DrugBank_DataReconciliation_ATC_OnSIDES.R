# Process ATC and OnSIDES
library(dbparser)
library(XML)
library(progress)
library(tidyverse)

# Copied helpers from dbparser Github
# https://github.com/ropensci/dbparser
AbstractParser <- R6::R6Class(
  "AbstractParser",
  public = list(
    initialize = function(tibble_name = NULL,
                          object_node = NULL,
                          main_node = NULL,
                          secondary_node = NULL,
                          id = NULL) {
      private$tibble_name    <- tibble_name
      private$object_node    <- object_node
      private$main_node      <- main_node
      private$secondary_node <- secondary_node
      private$id             <- id
      
    },
    parse = function() {
      parsed_tbl        <- private$parse_record()
      parsed_tbl        <- as_tibble(parsed_tbl)
      names(parsed_tbl) <- gsub(pattern     = "-",
                                replacement = "_",
                                x           =  names(parsed_tbl))
      parsed_tbl
    }
  ),
  private = list(
    tibble_name = NULL,
    object_node = NULL,
    main_node = NULL,
    secondary_node = NULL,
    id = NULL,
    parse_record = function() {
      message("I am the abstract parser, please use proper parser")
    }
  )
)

init_dvobject <- function() {
  dvobject        <- list()
  class(dvobject) <- "dvobject"
  
  attr(dvobject, "original_db_info") <- list()
  dvobject
}

ATCParser <- R6::R6Class(
  "ATCParser",
  inherit = AbstractParser,
  private = list(
    parse_record = function() {
      drugs <-  xmlChildren(pkg_env$root)
      return(as_tibble(do.call(
        rbind,
        xmlApply(pkg_env$root, private$atc_recs)
      )))
    },
    atc_recs = function(rec) {
      if (xmlSize(rec[["atc-codes"]]) < 1) {
        return()
      }
      atcs <- xmlApply(rec[["atc-codes"]], private$atc_rec)
      atcs_tibble <-
        as_tibble(do.call(rbind, atcs))
      atcs_tibble[["drugbank-id"]] = xmlValue(rec[["drugbank-id"]])
      return(atcs_tibble)
    },
    atc_rec = function(atc) {
      c(
        alpha_5th = xmlGetAttr(atc, name = "code"),
        name_4th = xmlValue(atc[[1]]),
        alpha_4th = xmlGetAttr(atc[[1]], name = "code"),
        name_3rd = xmlValue(atc[[2]]),
        alpha_3rd = xmlGetAttr(atc[[2]], name = "code"),
        name_2nd = xmlValue(atc[[3]]),
        alpha_2nd = xmlGetAttr(atc[[3]], name = "code"),
        name_1st = xmlValue(atc[[4]]),
        alpha_1st = xmlGetAttr(atc[[4]], name = "code")
      )
    }
  )
)

drug_atc_codes <- function() {
  ATCParser$new("drug_atc_codes")$parse()
}

add_database_info <- function(dvobject,
                              db_type          = "DrugBank",
                              db_version       = NULL,
                              db_exported_date = NULL) {
  db_info <- attr(dvobject, db_type)
  
  db_info[["db_type"]]               <- db_type
  db_info[["db_version"]]            <- db_version
  db_info[["db_exported_date"]]      <- db_exported_date
  attr(dvobject, "original_db_info") <- db_info
  class(dvobject)                    <- "dvobject"
  dvobject
}

parseOnSIDES <- function(dataDir,
                         include_high_confidence = TRUE,
                         db_version              = NULL,
                         db_exported_date        = NULL) {
  # The 7 canonical files
  core_files <- c(
    "product_label.csv",
    "product_adverse_effect.csv",
    "product_to_rxnorm.csv",
    "vocab_meddra_adverse_effect.csv",
    "vocab_rxnorm_ingredient.csv",
    "vocab_rxnorm_product.csv",
    "vocab_rxnorm_ingredient_to_product.csv"
  )
  
  file_paths        <- file.path(dataDir, core_files)
  names(file_paths) <- gsub("\\.csv$", "", core_files)
  
  if (!all(file.exists(file_paths))) {
    missing <- core_files[!file.exists(file_paths)]
    stop("Core files not found in '", dataDir, "':\n", paste(missing, collapse = "\n"))
  }
  
  message("Parsing the 7 core OnSIDES database tables...")
  db_tables <- lapply(names(file_paths), function(name) {
    message("Reading ", name, ".csv ...")
    data.table::fread(file_paths[[name]])
  })
  names(db_tables) <- names(file_paths)
  
  # --- Optional Handling of high_confidence.csv ---
  if (include_high_confidence) {
    hc_path <- file.path(dataDir, "high_confidence.csv")
    if (file.exists(hc_path)) {
      message("Reading high_confidence.csv ...")
      hc_table <- data.table::fread(hc_path)
      db_tables$high_confidence <- hc_table
    } else {
      warning("`include_high_confidence` was TRUE, but 'high_confidence.csv' was not found.")
    }
  }
  
  message("Successfully parsed OnSIDES database.")
  db_tables <- add_database_info(dvobject         = db_tables,
                                 db_type          = "OnSIDES",
                                 db_version       = db_version,
                                 db_exported_date = db_exported_date)
  db_tables
}
merge_drugbank_onsides <- function(db_object, onsides_db) {
  if ("drugbank" %in% names(db_object)) {
    drugbank_db   <- db_object$drugbank
    merged_object <- db_object 
  } else {
    drugbank_db            <- db_object
    merged_object          <- init_dvobject()
    merged_object$drugbank <- db_object
    attr(merged_object, "DrugBankDB") <- attr(drugbank_db, "original_db_info")
  }
  if (!inherits(drugbank_db, "dvobject") ||
      (!"drugs" %in% names(drugbank_db))) {
    stop("`db_object` must contain a valid DrugBank dvobject.")
  }
  
  if (!inherits(drugbank_db, "dvobject") ||
      (!"external_identifiers" %in% names(drugbank_db$drugs))) {
    stop("`drugbank_db` must contain external_identifiers data.")
  }
  if (!inherits(onsides_db, "dvobject") ||
      (!"vocab_rxnorm_ingredient" %in% names(onsides_db))) {
    stop("`onsides_db` must be a valid dvobject from parseOnSIDES().")
  }
  message("Creating DrugBank ID <-> RxCUI mapping table...")
  rxcui_mapping_df <- drugbank_db$drugs$external_identifiers %>%
    dplyr::filter(.data$resource == "RxCUI") %>%
    dplyr::select(all_of("drugbank_id"), rxcui = .data$identifier) %>%
    dplyr::distinct()
  message("Enriching OnSIDES tables with DrugBank IDs...")
  onsides_ingredient_enriched <- onsides_db$vocab_rxnorm_ingredient %>%
    dplyr::left_join(rxcui_mapping_df, by = c("rxnorm_id" = "rxcui"))
  if ("high_confidence" %in% names(onsides_db)) {
    onsides_hc_enriched <- onsides_db$high_confidence %>%
      dplyr::mutate(ingredient_id = as.character(.data$ingredient_id)) %>%
      dplyr::left_join(rxcui_mapping_df, by = c("ingredient_id" = "rxcui"))
  }
  message("Assembling final merged object...")
  if (is.null(merged_object$onsides)) {
    merged_object$onsides <- list()
  }
  for (name in names(onsides_db)) {
    merged_object$onsides[[name]] <- onsides_db[[name]]
  }
  
  # Ensure integrated_data list exists
  if (is.null(merged_object$integrated_data)) {
    merged_object$integrated_data <- list()
  }
  
  merged_object$integrated_data[["vocab_rxnorm_ingredient_enriched"]] <- onsides_ingredient_enriched
  
  if (exists("onsides_hc_enriched")) {
    merged_object$integrated_data[["high_confidence_enriched"]] <- onsides_hc_enriched
  }
  
  # Add the mapping table itself for user reference
  merged_object$integrated_data[["DrugBank_RxCUI_Mapping"]] <- rxcui_mapping_df
  
  # Update metadata
  attr(merged_object, "onSideDB") <- attr(onsides_db, "original_db_info")
  
  # Assign a new class (Prepend to keep existing classes like DrugBankTWOSIDESDb)
  class(merged_object) <- unique(c("DrugBankOnSIDESDb", class(merged_object)))
  
  message("Merge complete.")
  merged_object
}

#### To link and extract ATC with DrugBank ####
pkg_env <- new.env(parent = emptyenv())
pkg_env$root <- xmlRoot(xmlParse("1_PDI/1_DrugBank_v5.1.13/drugbank_all_full_database.xml/full database.xml"))
atc_df <- ATCParser$new("1_PDI/1_DrugBank_v5.1.13/drugbank_all_full_database.xml/full database.xml")$parse()
write.csv(atc_df, "8_ATC/dbparser_atc_df.csv", row.names = F)

#### To link and extract OnSIDES with DrugBank ####
# First, parse the individual databases
drugbank <- parseDrugBank("1_PDI/1_DrugBank_v5.1.13/drugbank_all_full_database.xml/full database.xml",
                          drug_options = c("external_identifiers"))
onsides <- parseOnSIDES("10_ONSIDES/onsides-v3.1.0/csv")
# Merge them into a single, powerful object
merged_db <- merge_drugbank_onsides(drugbank, onsides)
# Create a clean dataframe
onsides.df <- left_join(onsides[["product_adverse_effect"]], onsides[["vocab_meddra_adverse_effect"]],
                        c("effect_meddra_id" = "meddra_id")) %>%
  drop_na(meddra_name) %>% dplyr::filter(meddra_term_type == "PT") %>%
  full_join(., full_join(onsides[["product_label"]], 
                         onsides[["product_to_rxnorm"]], "label_id") %>%
              full_join(., merged_db[["integrated_data"]][["vocab_rxnorm_ingredient_enriched"]], 
                        c("rxnorm_product_id" = "rxnorm_id")) %>% drop_na(drugbank_id) %>%
              dplyr::select(label_id, drugbank_id, rxnorm_name),
            c("product_label_id" = "label_id")) %>% drop_na(drugbank_id)

write.csv(onsides.df, "10_ONSIDES/ONSIDES.processed.csv", row.names = F)
write.csv(drugbank[["drugs"]][["general_information"]] %>% dplyr::select(drugbank_id, name),
          "1_PDI/DrugBankID.csv", row.names = F)

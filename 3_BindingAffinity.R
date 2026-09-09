# Binding Affinity: Merging BDB, TDD, and ChemBL
# 25/09/2025 
# Kevin Winardi

# Load library
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Drug identifiers
drug.structure.df <- read.csv("1_PDI/1_DrugBank_v5.1.13/Structure/structure links.csv") %>%
  dplyr::select(DrugBank.ID, Name, Drug.Groups, SMILES, ChEMBL.ID) %>% distinct()
drug.structure.df[drug.structure.df == ""] <- NA

# Drug name harmonization
harmonise_drugnames <- function(df) {
  df %>% mutate(DRUGNAME = case_when(
    DRUGNAME %in% c("Aclidinium bromide", "Cetylpyridinium bromide",
                    "Homatropine hydrobromide", "Mepenzolate bromide",
                    "Methylbenactyzium bromide", "Tiotropium bromide",
                    "Umeclidinium bromide", "Valethamate bromide") ~ gsub(" bromide", "", DRUGNAME),
    DRUGNAME %in% c("Benzethonium chloride",
                    "Trospium chloride", "Berberine chloride")  ~ gsub(" chloride", "", DRUGNAME),
    DRUGNAME == "Benztropine"                    ~ "Benzatropine",
    DRUGNAME == "Betamethasone phosphoric acid"  ~ "Betamethasone phosphate",
    DRUGNAME == "Carbachol"                      ~ "Carbamoylcholine",
    DRUGNAME == "Carbetapentane"                 ~ "Pentoxyverine",
    DRUGNAME == "Diphenidol hydrochloride"        ~ "Diphenidol",
    DRUGNAME == "Homatropine hydrobromide"        ~ "Homatropine",
    T ~ DRUGNAME
  ))
}

# BDB
BDB.df <- read.csv("1_PDI/4_BindingDB_v2025-09-01/BindingDB_filtered3.csv") %>%
  dplyr::select(-c("kon..M.1.s.1.", "koff..s.1.", "pH", "Temp..C.",
                   "Number.of.Protein.Chains.in.Target...1.implies.a.multichain.complex."))

Gene.anno <- bitr(BDB.df$UniProt..SwissProt..Primary.ID.of.Target.Chain.1,
                  fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

char.to.num <- c("Ki..nM.", "IC50..nM.", "Kd..nM.", "EC50..nM.")

BDB.df <- left_join(BDB.df, Gene.anno,
                    c("UniProt..SwissProt..Primary.ID.of.Target.Chain.1" = "UNIPROT")) %>%
  dplyr::select(-starts_with("UniProt")) %>% distinct() %>%
  mutate(across(all_of(char.to.num), ~ ifelse(grepl("^>", .x),
                                              as.numeric(sub(">", "", .x)),
                                              as.numeric(.x))))

BDB.mAChR <- BDB.df %>%
  dplyr::filter(SYMBOL %in% paste0("CHRM", 1:5)) %>%
  inner_join(., drug.structure.df, c("DrugBank.ID.of.Ligand" = "DrugBank.ID")) %>%
  pivot_longer(cols = all_of(char.to.num), values_to = "binding_value", names_to = "binding_type") %>%
  mutate(binding_type = gsub("\\.\\.nM\\.", "", binding_type)) %>%
  drop_na(binding_value) %>%
  group_by(Name, SYMBOL, binding_type, SMILES) %>%
  dplyr::summarize(binding_value = median(binding_value, na.rm = T), .groups = "drop")

colnames(BDB.mAChR) <- c("DRUGNAME", "GENENAME", "binding_type", "SMILES", "binding_value")

# TTD
TTD.df <- read.csv("1_PDI/3_TTD/P1-07-Drug-TargetMapping.csv")
TTD.Activity <- read.delim("1_PDI/3_TTD/P1-09-Target_compound_activity.txt") %>%
  dplyr::filter(str_ends(Activity, "nM")) %>%
  mutate(binding_type  = str_extract(Activity, "^[A-Za-z0-9]+"),
         binding_value = as.numeric(str_extract(Activity, "[0-9.]+(?=\\s*nM)")))

TTD.Prot.Anno <- read.csv("1_PDI/3_TTD/P1-01-TTD_target_download.csv") %>%
  dplyr::select(-D, -E) %>%
  filter(!if_all(everything(), ~ .x == "")) %>%
  group_by(A, B) %>%
  summarise(C = paste(C, collapse = " /// "), .groups = "drop") %>%
  pivot_wider(names_from = "B", values_from = "C") %>% dplyr::select(-A)

ENTREZ.anno <- bitr(TTD.Prot.Anno$GENENAME, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

TTD.Drug.Anno <- full_join(
  read.csv("1_PDI/3_TTD/P1-04-Drug_synonyms.csv") %>%
    filter(!if_all(everything(), ~ .x == "")) %>%
    group_by(A, B) %>% summarise(C = paste(C, collapse = " /// "), .groups = "drop") %>%
    pivot_wider(names_from = "B", values_from = "C") %>% dplyr::select(-A),
  read.csv("1_PDI/3_TTD/P1-02-TTD_drug_download.csv") %>%
    filter(!if_all(everything(), ~ .x == "")) %>%
    pivot_wider(names_from = "B", values_from = "C") %>% dplyr::select(-A),
  c("TTDDRUID" = "DRUG__ID"))

TTD.mAChR <- TTD.df %>%
  full_join(., TTD.Activity, by = c("TargetID" = "TTD.Target.ID", "DrugID" = "TTD.Drug.Compound.ID")) %>%
  left_join(., TTD.Drug.Anno, c("DrugID" = "TTDDRUID")) %>%
  left_join(., TTD.Prot.Anno, c("TargetID" = "TARGETID")) %>%
  dplyr::rename("DRUG.SYNONYM" = "SYNONYMS.x", "PROTEIN.SYNONYM" = "SYNONYMS.y") %>%
  left_join(., ENTREZ.anno %>% filter(ENTREZID != 100124696), by = c("GENENAME" = "SYMBOL")) %>%
  distinct() %>%
  mutate(DRUGNAME = gsub("Oxybutynine", "Oxybutynin", DRUGNAME)) %>%
  drop_na(Activity, DRUGNAME) %>%
  dplyr::filter(GENENAME %in% paste0("CHRM", 1:5)) %>%
  dplyr::select(DRUGNAME, DRUGSMIL, GENENAME, binding_type, binding_value) %>%
  mutate(DRUGNAME = str_to_sentence(DRUGNAME))

colnames(TTD.mAChR) <- c("DRUGNAME", "SMILES", "GENENAME", "binding_type", "binding_value")

# Load each file and process independently
process_CHEMBL <- function(pattern, binding_label) {
  
  files <- list.files(path = "11_QSAR/CHEMBL", pattern = pattern, full.names = T)
  
  dat <- files %>% lapply(read.csv) %>% bind_rows() %>%
    dplyr::filter(Target.Name != "", !is.na(Standard.Value)) %>%
    mutate(Molecule.Name = str_to_sentence(Molecule.Name),
           Target.Name   = gsub("Muscarinic acetylcholine receptor M", "CHRM", Target.Name)) %>%
    dplyr::select(Molecule.ChEMBL.ID, Molecule.Name, Smiles, Target.Name,
                  Standard.Type, Standard.Relation, Standard.Value, Standard.Units) %>%
    left_join(., drug.structure.df[, c("ChEMBL.ID", "DrugBank.ID")],
              by = c("Molecule.ChEMBL.ID" = "ChEMBL.ID")) %>%
    left_join(., drug.structure.df[, c("SMILES", "DrugBank.ID")],
              by = c("Smiles" = "SMILES")) %>%
    left_join(., drug.structure.df[, c("Name", "DrugBank.ID")],
              by = c("Molecule.Name" = "Name")) %>%
    mutate(DrugBank.ID = coalesce(DrugBank.ID, DrugBank.ID.x, DrugBank.ID.y)) %>%
    dplyr::select(-DrugBank.ID.x, -DrugBank.ID.y) %>%
    left_join(., drug.structure.df[, c("DrugBank.ID", "Name", "SMILES")], by = "DrugBank.ID") %>%
    distinct() %>%
    mutate(
      Name   = coalesce(Name, Molecule.Name, Molecule.ChEMBL.ID),
      Smiles = ifelse(Smiles == "" | is.na(Smiles), SMILES, Smiles)
    ) %>%
    dplyr::select(-SMILES) %>% distinct()
  
  dat %>%
    group_by(Molecule.ChEMBL.ID, Name, DrugBank.ID, Smiles, Target.Name, Standard.Type) %>%
    dplyr::summarize(binding_value = median(Standard.Value), .groups = "drop") %>%
    dplyr::select(Name, Smiles, Target.Name, Standard.Type, binding_value)
}

# Merge databases
merge_and_summarise <- function(CHEMBL.mAChR, binding_label) {
  
  colnames(CHEMBL.mAChR) <- c("DRUGNAME", "SMILES", "GENENAME", "binding_type", "binding_value")
  
  BA.df <- rbind(BDB.mAChR, TTD.mAChR, CHEMBL.mAChR) %>%
    harmonise_drugnames()
  
  BA <- BA.df %>%
    group_by(DRUGNAME, GENENAME, binding_type) %>%
    dplyr::summarize(binding_value = median(binding_value, na.rm = T), .groups = "drop") %>%
    mutate(across(everything(), ~ ifelse(. == "NaN", NA, .))) %>%
    dplyr::filter(binding_type == binding_label) %>%
    mutate(
      BDB    = as.integer(DRUGNAME %in% (BDB.mAChR    %>% filter(binding_type == binding_label) %>% pull(DRUGNAME))),
      TTD    = as.integer(DRUGNAME %in% (TTD.mAChR    %>% filter(binding_type == binding_label) %>% pull(DRUGNAME))),
      CHEMBL = as.integer(DRUGNAME %in% (CHEMBL.mAChR %>% filter(binding_type == binding_label) %>% pull(DRUGNAME)))
    )
  
  BA
}

# pKi
CHEMBL.Ki <- process_CHEMBL("_Ki\\.csv$", "Ki")
BA.mAChR.Ki <- merge_and_summarise(CHEMBL.Ki, "Ki")

print(table(BA.mAChR.Ki$GENENAME, BA.mAChR.Ki$binding_type))

BA.mAChR.Ki %>%
  mutate(binding_value = -log10(binding_value / 1000)) %>%
  ggplot(aes(x = binding_value, color = GENENAME)) +
  labs(title = "pKi") +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = -log10(10), linetype = "dashed")

write.csv(BA.mAChR.Ki, "1_PDI/BA.mAChR.pKi.csv", row.names = F)

# IC50
CHEMBL.IC50 <- process_CHEMBL("_IC50\\.csv$", "IC50")
BA.mAChR.IC50 <- merge_and_summarise(CHEMBL.IC50, "IC50")

print(table(BA.mAChR.IC50$GENENAME, BA.mAChR.IC50$binding_type))

write.csv(BA.mAChR.IC50, "1_PDI/BA.mAChR.IC50.csv", row.names = F)
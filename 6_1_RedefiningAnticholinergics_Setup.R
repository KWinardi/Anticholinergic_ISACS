# Decoding Anticholinergic Effects

# Load library
library(ComplexHeatmap)
library(igraph)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Matrix)
library(ggplot2)
library(ggpmisc)
library(clipr)
library(utils)
library(parallel)
library(KEGGREST)
library(ggpmisc)
library(cowplot)
library(ggrepel)
library(pbapply)
library(ggridges)
library(ggsci)
library(circlize)
library(glmnet)
library(wesanderson)
library(ggpubr)
library(ggsankey)
library(ggalluvial)
library(ComplexUpset)
library(ggsci)
library(rstatix)
library(pROC)
library(rsample)
library(ggh4x)
library(ggstats)
library(PRROC)

mytheme <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.key = element_blank(),
                 legend.box.background = element_blank(),
                 legend.title=element_text(size=11.5), 
                 legend.text=element_text(size=10.5),
                 plot.title = element_text(face="bold",  size=14, hjust = 0.5),
                 plot.subtitle = element_text(face="italic",  size=13, hjust = 0.5),
                 axis.title = element_text(size = 12, color = "black"),
                 axis.text = element_text(size = 11, color = "black"),
                 strip.text = element_text(size = 11, color = "black"),
                 panel.border = element_rect(fill = NA), 
                 strip.background = element_rect(fill = NA), 
                 axis.ticks.length=unit(.15, "cm"))

# Load network results
file.list <- list.files(path = "3_NetworkResults/", pattern = "^AllDrugs\\..*\\.csv$", 
                        full.names = T)

distance.df <- read.csv("1_PDI/DrugBankID.csv") %>%
  right_join(., map_dfr(file.list, read.csv),
             c("name" = "Drug")) %>% 
  dplyr::rename(Drug = name)

distance.df.long <- distance.df %>%
  pivot_longer(cols = -c(Drug, drugbank_id, DrugTarget_degree),
               names_to = c("distance.metric", "stat", "X"),
               names_pattern = "(.+)_(pval|zscore)|(.+)") %>%
  mutate(distance.metric = paste0(distance.metric, X),
         stat = ifelse(stat == "", "distance.value", stat)) %>% dplyr::select(-X) %>%
  pivot_wider(values_from = "value", names_from = "stat") %>% distinct %>%
  mutate(zscore_sig = case_when(zscore < -1.96~"Significant",
                                T~"NS"))

# Descriptive statistics
distance.df.long %>% 
  group_by(distance.metric, zscore_sig) %>%
  count

### Structural data
drug.structure.df <- read.csv("1_PDI/1_DrugBank_v5.1.13/Structure/structure links.csv") %>%
  dplyr::filter(Name %in% unique(distance.df$Drug)) %>% 
  dplyr::select(DrugBank.ID, Name, Drug.Groups, SMILES, ChEMBL.ID) %>% distinct

### ADMET
ADMET.metadata <- read.csv("9_ADMETLAB_3.0/ADMETLAB_Param.csv")

list.ADMET <- list.files(path = "9_ADMETLAB_3.0", full.names = T,
                         pattern = "\\.csv$")
list.ADMET <- list.ADMET[list.ADMET!="9_ADMETLAB_3.0/ADMETLAB_Param.csv"]
list.ADMET <- list.ADMET[list.ADMET!="9_ADMETLAB_3.0/compiled.ADMET.csv"]
ADMET.files <- list()
for (i in list.ADMET) {
  df <- read.csv(i)
  ADMET.files[[i]] <- df
}

ADMET.df <- do.call(rbind, ADMET.files) %>% 
  dplyr::filter(!MW %in% "Invalid Molecule") %>%
  dplyr::select(-molstr, -smiles) %>%
  dplyr::select(-c(ADMET.metadata %>% dplyr::filter(str_detect(Unit, "List of")) %>% pull(Abbreviation))) %>%
  mutate(across(!c(raw_smiles), as.numeric))
rownames(ADMET.df) <- NULL

drug.ADMET <- drug.structure.df %>% left_join(., ADMET.df, c("SMILES"="raw_smiles"))

### Anticholinergic scales, SAA, & ATS
anticholinergics.df <- read.csv("5_AnticholinergicScaleTable/Lozano-Ortega2020_AnticholinergicScale.csv") %>%
  mutate(Generic.Name = str_to_sentence(Generic.Name)) %>%
  mutate(Generic.Name = gsub("\xa0", " ", Generic.Name))

ACSBC.df <- read.csv("5_AnticholinergicScaleTable/Rihani2021_ACSBC.csv") %>% distinct
colnames(ACSBC.df)[-1] <- paste0("ACSBC.", colnames(ACSBC.df)[-1])

ABS.df <- read.csv("5_AnticholinergicScaleTable/Yamada2023_ABS_IC50.csv") %>%
  mutate(Drug = str_to_sentence(Drug))

ORCA.df <- read.csv("5_AnticholinergicScaleTable/ORCA.csv", fileEncoding = "Windows-1252") %>%
  mutate(Drug = Drug %>% stringi::stri_enc_toutf8() %>%      
           stringr::str_replace_all("[\u00A0\u200B\uFEFF]", " ") %>% 
           stringr::str_squish()) %>% 
  dplyr::filter(!Drug %in% c("Antimicrobial", "No DB data", "Plant", 
                             "Belladonna", "Lithium")) 

SAA.df <- read.csv("5_AnticholinergicScaleTable/Chew2008_SAA.csv") %>%
  mutate(SAA.numeric = as.numeric(ifelse(SAA == ">250", 251, SAA)),
         SAA.per.dose = case_when(is.na(Daily_dose)~0,
                                  T~SAA.numeric/Daily_dose),
         Medication = gsub("\xa0", " ", Medication)) %>%
  group_by(Medication) %>%
  mutate(Daily.dose.min = min(Daily_dose)) %>%
  ungroup() %>%
  mutate(SAA.min.ther.dose = case_when(is.na(Daily_dose)~0,
                                       T~SAA.numeric/(Daily_dose/Daily.dose.min)))

SAA.df.average <- SAA.df %>%
  mutate(SAA_class = factor(case_when(SAA_class == "No"~"0",
                                      SAA_class == "No or minimal"~"0/+",
                                      T~SAA_class),
                            levels = c("0", "0/+",
                                       "+", "++", "+++"))) %>% 
  group_by(Medication, Class, SAA_class) %>% 
  dplyr::summarize(SAA.per.dose = mean(SAA.per.dose, na.rm = T),
                   SAA.min.ther.dose = mean(SAA.min.ther.dose, na.rm = T),
                   .groups = "drop") 

ATS.df <- read.csv("5_AnticholinergicScaleTable/Xu2017_ATS.csv")
colnames(ATS.df)[2] <- "ATS.class"

AEC.df <- read.csv("5_AnticholinergicScaleTable/Bishara2016_AEC.csv")

### ATC code
ATC.df <- read.csv("8_ATC/dbparser_atc_df.csv") %>%
  dplyr::filter(!str_detect(name_4th, regex("combination", ignore_case = TRUE)))

### Binding affinity
BA.mAChR <- read.csv("1_PDI/BA.mAChR.pKi.csv") %>% 
  mutate(pKi = -log10(binding_value/10^9))
BA.mAChR.formatted <- BA.mAChR %>% 
  mutate(mAChR_activity = paste0(GENENAME, ".", binding_type)) %>%
  dplyr::select(DRUGNAME, mAChR_activity, pKi) %>% 
  pivot_wider(names_from = mAChR_activity, values_from = pKi) %>% 
  dplyr::select(order(colnames(.)))

BA.list <- BA.mAChR %>% pull(DRUGNAME) %>% unique

## QSAR from PSICHIC
### Only mAChR
# PSICHIC.df <- rbind(read.csv("4_QSAR/screening_CAP1.csv")) %>%
#   mutate(predicted_antagonist1 = predicted_antagonist/(predicted_antagonist+predicted_agonist),
#          predicted_agonist1 = predicted_agonist/(predicted_antagonist+predicted_agonist),
#          predicted_binding_affinity = predicted_binding_affinity*(1+predicted_antagonist1)
#          ) %>%
#   dplyr::select(Drug, Gene, predicted_binding_affinity) %>%
#   pivot_wider(names_from = "Gene", values_from = "predicted_binding_affinity") 

### ALL CAPs
PSICHIC.df <- rbind(read.csv("4_QSAR/screening_CAP1.csv"),
                    read.csv("4_QSAR/screening_CAP2.csv")) %>%
  mutate(predicted_antagonist1 = predicted_antagonist/(predicted_antagonist+predicted_agonist),
         predicted_agonist1 = predicted_agonist/(predicted_antagonist+predicted_agonist),
         predicted_binding_affinity = case_when(Gene == "ACHE"~predicted_binding_affinity*(1+predicted_agonist1),
                                                T~predicted_binding_affinity*(1+predicted_antagonist1))
  ) %>%
  dplyr::select(Drug, Gene, predicted_binding_affinity) %>%
  pivot_wider(names_from = "Gene", values_from = "predicted_binding_affinity")

colnames(PSICHIC.df)[-1] <- paste0(colnames(PSICHIC.df)[-1], "_QSAR")
Outcomes <- colnames(PSICHIC.df)[-1] 
QSAR.PCA <- prcomp(PSICHIC.df[,-1], scale. = T, center = T)
summary(QSAR.PCA)

PSICHIC.df$QSAR.PCA <- QSAR.PCA$x[,1]

## For ground truth assessment
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
false_pos <- c(AChEi, Organophosphate, Cholinergic.Agonist)
true_neg <- c(Bisphosphonate, Statin, ACEi, DOAC)

drug_category_map <- c(setNames(rep("Prototypical anticholinergics", length(prototypical.anticholinergic)), prototypical.anticholinergic),
                       setNames(rep("Urinary incontinence agent", length(OAB.anticholinergic)), OAB.anticholinergic),
                       setNames(rep("TCA", length(TCA)), TCA),
                       setNames(rep("Sedating antihistamine", length(Sedative.antihistamine)), Sedative.antihistamine),
                       setNames(rep("Bisphosphonate", length(Bisphosphonate)), Bisphosphonate),
                       setNames(rep("Statin", length(Statin)), Statin),
                       setNames(rep("ACE inhibitor", length(ACEi)), ACEi),
                       setNames(rep("DOAC", length(DOAC)), DOAC),
                       setNames(rep("Acetylcholinesterase inhibitor", length(AChEi)), AChEi),
                       setNames(rep("Organophosphate", length(Organophosphate)), Organophosphate),
                       setNames(rep("Cholinergic agonist", length(Cholinergic.Agonist)), Cholinergic.Agonist))

# ATC enrichment
ATC.df.enrich <- ATC.df %>% 
  mutate(ATC_3rd = paste0(name_3rd, " (", alpha_3rd, ")")) %>%
  left_join(., read.csv("1_PDI/DrugBankID.csv"), "drugbank_id") %>%
  distinct(ATC_3rd, name) %>% drop_na 

lvl_order <- c(c(ATC.df %>% distinct(name_1st) %>% 
                   mutate(name_1st = str_to_sentence(name_1st)) %>% arrange(name_1st) %>%
                   pull(name_1st)), "No ATC code")
ATC.col.code <- c("Alimentary tract and metabolism" = "#e8ce4d",
                  "Blood and blood forming organs" = "#ed5565",
                  "Cardiovascular system" = "#FD7446",
                  "Dermatologicals" = "#ed87c1",
                  "Genito urinary system and sex hormones" = "#008280",
                  "Systemic hormonal preparations, excl. Sex hormones and insulins" = "#a0cecb",
                  "Antiinfectives for systemic use" = "#8a7c74",
                  "Antiparasitic products, insecticides and repellents" = "#3b2929",
                  "Antineoplastic and immunomodulating agents" = "#D2AF81FF",
                  "Musculo-skeletal system" = "#4dc0e8",
                  "Nervous system" = "#ac92ec",
                  "Respiratory system" = "#008B45FF",
                  "Sensory organs" = "#5de4d2",
                  "Various" = "#c5de5f",
                  "No ATC code" = "#99958b")
ATC.col <- scale_color_manual(values = ATC.col.code,
                              na.value = "#99958b", breaks = lvl_order)
ATC.fill <- scale_fill_manual(values = ATC.col.code,
                              na.value = "#99958b", breaks = lvl_order)

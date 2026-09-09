# Cross-sectional pharmacoepidemiological analysis
## Data cleaning
## Created by Kevin Winardi (26/08/2026)
library(nhanesA)
library(tidyverse)
library(survey)
library(broom)
library(tableone)
library(tictoc)

source("6_1_RedefiningAnticholinergics_Setup.R")
result.w.anno.full <- read.csv("ISACS.csv")

cycle_info <- tibble::tibble(
  suffix = c("", "_B", "_C", "_D", "_E", "_F", "_G", "_H", "_I", "_J"),
  years  = c("1999-2000","2001-2002","2003-2004","2005-2006","2007-2008",
             "2009-2010","2011-2012","2013-2014","2015-2016","2017-2018"))

recode_yn <- function(x) {
  case_when(x == 1 ~ 1,
            x == 2 ~ 0,
            x %in% c(7, 9) ~ NA_real_,
            TRUE ~ NA_real_)
}

anticholinergic.scales <- c("ACSBC", "ADS", "ABS", "ALS","ORCA", 
                            "AEC","Chew's SAA", "ACB", "AAS", 
                            "ARS", "ABC", "ISACS", "ATS")

antichol.scale.df <- result.w.anno.full %>% 
  select(Drug, cholinergicCluster) %>%
  # Anticholinergic scale
  full_join(., anticholinergics.df, c("Drug" = "Generic.Name")) %>%
  # ORCA
  full_join(., ORCA.df, "Drug") %>%
  mutate(ORCA = factor(ORCA, levels = c("High", "Moderate", "Low", "No"))) %>%
  # ACSBC
  full_join(., ACSBC.df, "Drug") %>%
  mutate(ACSBC.Anticholinergic.Activity = factor(ACSBC.Anticholinergic.Activity,
                                                 levels = c("High", "Moderate", "Low", "No")),
         ACSBC.Sedative.Activity = factor(ACSBC.Sedative.Activity,
                                          levels = c("High", "Moderate", "Low", "No")),
         cholinergicCluster = factor(cholinergicCluster,
                                     levels = c("No", "Low", "Moderate", "High"))
         ) %>%
  # SAA
  full_join(., SAA.df.average, 
            c("Drug" = "Medication")) %>% 
  # ATS
  full_join(., ATS.df, c("Drug"="Medication")) %>% 
  # ABS 
  full_join(., ABS.df, "Drug") %>%
  # ATS
  full_join(., AEC.df, "Drug") %>%
  mutate(ISACS = as.numeric(cholinergicCluster)-1,
         ACSBC = as.numeric(fct_rev(ACSBC.Anticholinergic.Activity))-1,
         ORCA = as.numeric(fct_rev(ORCA))-1,
         `Chew's SAA` = as.numeric(SAA_class)-1,
         AEC = case_when(AEC == "No consensus"~NA,
                         T~as.numeric(AEC))
  ) %>% select(Drug, all_of(anticholinergic.scales))
antichol.scale.df[is.na(antichol.scale.df)] <- 0

# Extract demographic information
get_demo <- function(suffix) {
  tbl_name <- paste0("DEMO", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    mutate(cycle = tbl_name) 
}

demo_list <- map(cycle_info$suffix, get_demo)
names(demo_list) <- cycle_info$years
demo_all <- bind_rows(demo_list, .id = "cycle_years")

## Gender and age
table(demo_all$RIAGENDR, demo_all$RIDAGEYR)

# Cognitive test
get_CFQ <- function(suffix) {
  tbl_name <- paste0("CFQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

CFQ_list <- map(cycle_info$suffix, get_CFQ)
names(CFQ_list) <- cycle_info$years
CFQ_df <- bind_rows(CFQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN)) %>%
  unite("CFDDS", any_of(c("CFDDS", "CFDRIGHT")), sep = ";", na.rm = TRUE) %>%
  mutate(CFDCST_total = CFDCST1 + CFDCST2 + CFDCST3,
         CFDDS = as.numeric(CFDDS))
CFQ_df %>% drop_na(CFDDS) %>% count(cycle_years)
# CFDRIGHT: DSST 1999-2002
# CFDCST_total: CERAD immediate recall (IR) 2011-2014
# CFDCSR: CERAD delay recall (DR) 2011-2014
# CFDAST: Animal fluency test 2011-2014
# CFDDS: DSST 2011-2014

# Physical fitness
get_PAQ <- function(suffix) {
  tbl_name <- paste0("PAQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

PAQ_list <- map(cycle_info$suffix, get_PAQ)
names(PAQ_list) <- cycle_info$years
PAQ_df <- bind_rows(PAQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN)) %>% 
  unite("moderate_pa", any_of(c("PAD200", "PAQ665")), sep = ";", na.rm = TRUE) %>%
  unite("vigorous_pa", any_of(c("PAD320", "PAQ650")), sep = ";", na.rm = TRUE) %>%
  mutate(moderate_pa = recode_yn(as.numeric(moderate_pa)),
         vigorous_pa = recode_yn(as.numeric(vigorous_pa))) %>%
  mutate(active = case_when(moderate_pa == 1 | vigorous_pa == 1 ~ 1,
                            moderate_pa == 0 & vigorous_pa == 0 ~ 0,
                            TRUE ~ NA_real_)) %>%
  select(SEQN, cycle_years, moderate_pa, vigorous_pa, active)
# PAD200: Vigorous activity over 30 days 1999-2006
# PAD320: Moderate activity over 30 days 1999-2006
# PAQ650: Vigorous recreational activity 2007-2018
# PAQ665: Moderate recreational activity 2007-2018

# Smoking
get_SMQ <- function(suffix) {
  tbl_name <- paste0("SMQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

SMQ_list <- map(cycle_info$suffix, get_SMQ)
names(SMQ_list) <- cycle_info$years
SMQ_df <- bind_rows(SMQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN)) %>%
  select(SEQN, cycle_years, SMQ020, SMQ040)
# SMQ040

# Alcohol
get_ALQ <- function(suffix) {
  tbl_name <- paste0("ALQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

ALQ_list <- map(cycle_info$suffix, get_ALQ)
names(ALQ_list) <- cycle_info$years
ALQ_df <- bind_rows(ALQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN)) %>% 
  mutate(ALQ121 = case_when(ALQ111 == 1 & ALQ121 %in% c(1,2,3,4,5,6,7)~1,
                            ALQ111 == 1 & ALQ121 %in% c(0,8,9,10)~2,
                            ALQ111 == 2~2,
                            ALQ121 %in% c(77, 99)~NA)) %>%
  unite("alcohol", any_of(c("ALQ100", "ALD100", "ALQ101", "ALQ121")),
        sep = ";", na.rm = TRUE) %>%
  mutate(alcohol = recode_yn(alcohol)) %>%
  select(SEQN, cycle_years, alcohol)
# ALQ100 1999-2000
# ALD100 2001-2002
# ALQ101 2003-2016
# ALQ121 2017-2018

# Multimorbidity
get_BPQ <- function(suffix) {
  tbl_name <- paste0("BPQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

get_DIQ <- function(suffix) {
  tbl_name <- paste0("DIQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

get_MCQ <- function(suffix) {
  tbl_name <- paste0("MCQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

get_KIQ <- function(suffix) {
  tbl_name <- paste0("KIQ", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

get_KIQ_U <- function(suffix) {
  tbl_name <- paste0("KIQ_U", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F, cleanse_numeric = T), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

BPQ_list <- map(cycle_info$suffix, get_BPQ)
names(BPQ_list) <- cycle_info$years
BPQ_df <- bind_rows(BPQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN))
# BPQ020 Ever told HBP 1999-2018
# BPQ080 Doctor told you - high cholesterol level 1999-2018

DIQ_list <- map(cycle_info$suffix, get_DIQ)
names(DIQ_list) <- cycle_info$years
DIQ_df <- bind_rows(DIQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN))
# DIQ010 Doctor told you have diabetes 1999-2018

MCQ_list <- map(cycle_info$suffix, get_MCQ)
names(MCQ_list) <- cycle_info$years
MCQ_df <- bind_rows(MCQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN)) %>%
  unite("thyroid", any_of(c("MCQ160I", "MCD160M", "MCQ160M")),
        sep = ";", na.rm = TRUE) %>%
  mutate(thyroid = as.numeric(thyroid))
# MCQ010 Ever been told you have asthma 1999-2018
# MCQ100 Told have high blood pressure 1999-2000
# MCQ160A Doctor ever said you had arthritis 1999-2018
# MCQ160B Ever told had congestive heart failure 1999-2018
# MCQ160C Ever told you had coronary heart disease 1999-2018
# MCQ160D Ever told you had angina/angina pectoris 1999-2018
# MCQ160E Ever told you had heart attack 1999-2018
# MCQ160F Ever told you had a stroke 1999-2018
# MCQ160G Ever told you had emphysema 1999-2018
# MCQ160K Ever told you had chronic bronchitis 1999-2018
# MCQ160L Ever told you had any liver condition 1999-2018 
# MCQ220 Ever told you had cancer or malignancy 1999-2018 

# MCQ160I Ever told you had thyroid disease 1999-2000
# MCD160M Ever told you had a thyroid problem 2001-2002
# MCQ160M Ever told you had thyroid problem 2003-2018

KIQ_list <- map(cycle_info$suffix, get_KIQ)
names(KIQ_list) <- cycle_info$years
KIQ_df <- bind_rows(KIQ_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN))
# KIQ020 Ever told you had weak/failing kidneys 1999-2000

KIQ_U_list <- map(cycle_info$suffix, get_KIQ_U)
names(KIQ_U_list) <- cycle_info$years
KIQ_U_df <- bind_rows(KIQ_U_list, .id = "cycle_years") %>% mutate(SEQN = as.character(SEQN))

KIQ_all <- full_join(KIQ_df   %>% select(cycle_years, SEQN, KIQ020),
                     KIQ_U_df %>% select(cycle_years, SEQN, KIQ022),
                     by = c("cycle_years", "SEQN")) %>%
  unite("weak_kidney", KIQ020, KIQ022, sep = ";", na.rm = TRUE) %>%
  mutate(weak_kidney = as.numeric(weak_kidney))
# KIQ022 Ever told you had weak/failing kidneys 2001-2018

chronic_conditions <- demo_all %>%
  mutate(SEQN = as.character(SEQN)) %>%
  distinct(SEQN, cycle_years, RIDAGEYR, RIAGENDR) %>%
  left_join(BPQ_df %>% select(SEQN, BPQ020, BPQ080), by = "SEQN") %>%
  left_join(DIQ_df %>% select(SEQN, DIQ010), by = "SEQN") %>%
  left_join(MCQ_df %>% select(SEQN, MCQ010, MCQ160A, MCQ160B, MCQ160C, MCQ160D,
                              MCQ160E, MCQ160F, MCQ160G, MCQ160K, MCQ160L,
                              MCQ220, thyroid), by = "SEQN") %>%
  left_join(KIQ_all %>% select(SEQN, weak_kidney), by = "SEQN") %>%
  mutate(across(c(BPQ020, BPQ080, DIQ010, MCQ010, MCQ160A, MCQ160B, MCQ160C, MCQ160D,
                  MCQ160E, MCQ160F, MCQ160G, MCQ160K, MCQ160L, MCQ220, thyroid, weak_kidney),
                recode_yn))

cond_cols <- c("BPQ020", "BPQ080", "DIQ010", "MCQ010", "MCQ160A", "MCQ160B",
               "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ160G", "MCQ160K",
               "MCQ160L", "MCQ220", "thyroid", "weak_kidney")
# hypertension, hypercholesterolemia, diabetes, asthma, arthritis, CHF, CHD,
# angina, MI, stroke, emphysema, chronic bronchitis, liver disease, cancer,
# thyroid disorder, weak/failing kidney

cond_matrix <- as.matrix(chronic_conditions %>% select(all_of(cond_cols)))

chronic_conditions <- chronic_conditions %>%
  mutate(n_conditions_missing = rowSums(is.na(cond_matrix)),
         n_conditions = rowSums(cond_matrix, na.rm = TRUE))

# Medication list
get_RXQ_RX <- function(suffix) {
  tbl_name <- paste0("RXQ_RX", suffix)
  message("Fetching ", tbl_name)
  d <- tryCatch(nhanes(tbl_name, translated = F), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d %>%
    filter(SEQN %in% demo_all$SEQN) %>%
    mutate(cycle = tbl_name)
}

RXQ_RX_list <- map(cycle_info$suffix, get_RXQ_RX)
names(RXQ_RX_list) <- cycle_info$years
RXQ_RX_df <- bind_rows(RXQ_RX_list, .id = "cycle_years") %>%
  mutate(SEQN = as.character(SEQN)) %>%
  # Number of medicines
  unite("RXDCOUNT", RXD295, RXDCOUNT, sep = ";", na.rm = T) %>%
  # Past month
  unite("RXDUSE", RXD030, RXDUSE, sep = ";", na.rm = T) %>%
  # Generic name
  unite("RXDDRUG", RXD240B, RXDDRUG, sep = ";", na.rm = T) %>%
  # Duration of use
  unite("RXDDAYS", RXD260, RXDDAYS, sep = ";", na.rm = T) %>%
  # Verified by interviewer
  unite("RXQSEEN", RXQ250, RXQSEEN, sep = ";", na.rm = T) %>%
  mutate(across(c(RXDCOUNT, RXDUSE, RXDDAYS, RXQSEEN), as.numeric)) %>%
  select(cycle_years, SEQN, RXDCOUNT, RXDUSE, RXDDRUG, RXDDAYS, RXQSEEN)

nMeds <- RXQ_RX_df %>% distinct(SEQN, RXDCOUNT) %>%
  mutate(RXDCOUNT = case_when(is.na(RXDCOUNT)~0, T~as.numeric(RXDCOUNT)),
         nMeds.class = case_when(
           RXDCOUNT == 0 ~ "No medication use",
           RXDCOUNT >= 1 & RXDCOUNT <= 4 ~ "1-4 medications",
           RXDCOUNT >= 5 & RXDCOUNT <= 9 ~ "Polypharmacy (5-9)",
           RXDCOUNT >= 10 ~ "Hyperpolypharmacy (≥10)"
         ) %>% factor(levels = c("No medication use",
                                 "1-4 medications", 
                                 "Polypharmacy (5-9)",
                                 "Hyperpolypharmacy (≥10)")))

# More cleaning?
RXQ_RX_df <- RXQ_RX_df %>%
  filter(RXDUSE %in% c(1,2)) %>% 
  filter(!RXDDRUG %in% c("55555", "77777", "99999")) %>%
  filter(!(RXDUSE == 1 & RXDDRUG == "")) %>%
  filter(RXDDAYS >= 30 & RXDDAYS != 99999 & RXDDAYS != 77777 | is.na(RXDDAYS)) %>% 
  # filter(RXQSEEN %in% c(1, 3) | is.na(RXQSEEN)) %>%
  # Separate combination therapy
  separate_rows(RXDDRUG, sep = "; ")

nMeds.2 <- RXQ_RX_df %>% 
  count(SEQN) %>% 
  dplyr::rename("RXDCOUNT"=n) %>%
  left_join(., RXQ_RX_df %>% distinct(SEQN, RXDUSE), "SEQN") %>%
  mutate(RXDCOUNT = case_when(is.na(RXDCOUNT)~0, T~as.numeric(RXDCOUNT)),
         RXDCOUNT = case_when(RXDUSE==2~0, T~as.numeric(RXDCOUNT)),
         nMeds.class = case_when(
           RXDCOUNT == 0 ~ "No medication use",
           RXDCOUNT >= 1 & RXDCOUNT <= 4 ~ "1-4 medications",
           RXDCOUNT >= 5 & RXDCOUNT <= 9 ~ "Polypharmacy (5-9)",
           RXDCOUNT >= 10 ~ "Hyperpolypharmacy (≥10)"
         ) %>% factor(levels = c("No medication use",
                                 "1-4 medications", 
                                 "Polypharmacy (5-9)",
                                 "Hyperpolypharmacy (≥10)")))

salt_suffixes <- c("hydrochloride", "hydrobromide", "mesylate", "besylate", "tartrate",
                   "bitartrate", "succinate", "fumarate", "tannate", "citrate", "sulfate",
                   "maleate", "calcium", "bromide", "sodium", "potassium", "dipropionate",
                   "monohydrate", "terephthalate", "napsylate", "saccharate", "benzoate",
                   "pamoate", "hyclate", "nitrate", "gluconate", "trihydrate", "polistirex", "bismuth",
                   "carbonate", "bicarbonate", "hippurate", "lactate", "tromethamine", "dihydrochloride")
salt_pattern <- paste0(" (", paste(salt_suffixes, collapse = "|"), ")\\b")

# Matching with drugbank and anticholinergic scales
ISACS.NHANES <- RXQ_RX_df %>% distinct(RXDDRUG) %>% 
  mutate("Drug" = str_to_sentence(RXDDRUG)) %>%
  # Salts standardized
  mutate(Drug = gsub(salt_pattern, "", Drug)) %>%
  # Drug specific renaming to match DrugBank entries
  mutate(Drug = gsub("Aspirin", "Acetylsalicylic acid", Drug)) %>%
  mutate(Drug = gsub("Ursodiol", "Ursodeoxycholic acid", Drug)) %>%
  mutate(Drug = gsub("Divalproex", "Valproic acid", Drug)) %>%
  mutate(Drug = gsub("Esomeprazole magnesium", "Esomeprazole", Drug)) %>%
  mutate(Drug = gsub("Methscopolamine", "Methscopolamine bromide", Drug)) %>%
  mutate(Drug = gsub("Alloxanthine", "Oxypurinol", Drug)) %>%
  mutate(Drug = gsub("Beclomethasone", "Beclomethasone dipropionate", Drug)) %>%
  mutate(Drug = gsub("Azilsartan", "Azilsartan medoxomil", Drug)) %>%
  mutate(Drug = gsub("Alendronate", "Alendronic acid", Drug)) %>%
  mutate(Drug = gsub("Etidronate", "Etidronic acid", Drug)) %>%
  mutate(Drug = gsub("Risedronate", "Risedronic acid", Drug)) %>%
  mutate(Drug = gsub("Benztropine", "Benzatropine", Drug)) %>%
  mutate(Drug = gsub("Propoxyphene", "Dextropropoxyphene", Drug)) %>%
  mutate(Drug = gsub("Salmeterol xinafoate", "Salmeterol", Drug)) %>%
  mutate(Drug = gsub("Bisoprolol fumarate", "Bisoprolol", Drug)) %>%
  mutate(Drug = gsub("Fludrocortisone acetate", "Fludrocortisone", Drug)) %>%
  mutate(Drug = gsub("Pirbuterol acetate", "Pirbuterol", Drug)) %>%
  mutate(Drug = gsub("Desmopressin acetate", "Desmopressin", Drug)) %>%
  mutate(Drug = gsub("Triamcinolone acetate", "Triamcinolone", Drug)) %>%
  mutate(Drug = gsub("Triamcinolone diacetate", "Triamcinolone", Drug)) %>%
  mutate(Drug = gsub("Triamcinolone acetonide", "Triamcinolone", Drug)) %>%
  mutate(Drug = gsub("Clopidogrel bisulfate", "Clopidogrel", Drug)) %>%
  mutate(Drug = gsub("Codeine phosphate", "Codeine", Drug)) %>%
  mutate(Drug = gsub("Oxybutynin chloride", "Oxybutynin", Drug)) %>%
  mutate(Drug = gsub("Estrogens, conjugated", "Conjugated estrogens", Drug)) %>%
  mutate(Drug = gsub("Ethinyl estradiol", "Ethinylestradiol", Drug)) %>%
  mutate(Drug = gsub("Amphetamine aspartate", "Amphetamine", Drug)) %>%
  mutate(Drug = gsub("Doxycycline", "Doxycycline anhydrous", Drug)) %>%
  mutate(Drug = gsub("Torsemide", "Torasemide", Drug)) %>%
  mutate(Drug = gsub("Mesalamine", "Mesalazine", Drug)) %>%
  mutate(Drug = gsub("Nitrofurantoin", "Nitrofurantoin, macrocrystals", Drug)) %>%
  mutate(Drug = gsub("Pyrilamine", "Mepyramine", Drug)) %>%
  mutate(Drug = gsub("Clorazepate dipotassium", "Clorazepic acid", Drug)) %>%
  mutate(Drug = gsub("Omega-3 polyunsaturated fatty acids", "Omega-3 fatty acids", Drug)) %>%
  mutate(Drug = gsub("Dexamethasone phosphate", "Dexamethasone", Drug)) %>%
  mutate(Drug = gsub("Ethambutol", "Ethambutol hydrochloride", Drug)) %>%
  mutate(Drug = gsub("Clavulanate", "Clavulanic acid", Drug)) %>%
  mutate(Drug = gsub("Betamethasone valerate", "Betamethasone", Drug)) %>%
  mutate(Drug = gsub("Guanabenz acetate", "Guanabenz", Drug)) %>%
  mutate(Drug = gsub("Methylprednisolone acetate", "Methylprednisolone", Drug)) %>%
  mutate(Drug = gsub("Dimethyl", "Dimethyl fumarate", Drug)) %>%
  mutate(Drug = gsub("Bethanechol chloride", "Bethanechol", Drug)) %>%
  mutate(Drug = gsub("Clorazepate", "Clorazepic acid", Drug)) %>%
  mutate(Drug = gsub("Carbetapentane", "Pentoxyverine", Drug)) %>%
  mutate(Drug = gsub("Chlorpromamide", "Chlorpropamide", Drug)) %>%
  mutate(Drug = gsub("Clomiphene", "Enclomiphene", Drug)) %>%
  mutate(Drug = gsub("Colestipol", "Colestipol hydrochloride", Drug)) %>%
  mutate(Drug = gsub("Cefaclor", "Cefaclor anhydrous", Drug)) %>%
  mutate(Drug = gsub("Cromolyn", "Cromoglicic acid", Drug)) %>%
  mutate(Drug = gsub("Calcium acetate", "Calcium", Drug)) %>%
  mutate(Drug = gsub("Lithium", "Lithium cation", Drug)) %>%
  mutate(Drug = gsub("Cinarizine", "Cinnarizine", Drug)) %>%
  mutate(Drug = gsub("Ethacrynic acid", "Etacrynic acid", Drug)) %>%
  mutate(Drug = gsub("Fluocinolone acetonide", "Fluocinolone", Drug)) %>%
  mutate(Drug = gsub("Glycopyrrolate", "Glycopyrronium", Drug)) %>%
  mutate(Drug = gsub("Dexchlorpheniramine", "Dexchlorpheniramine maleate", Drug)) %>%
  mutate(Drug = gsub("Calcipotriene", "Calcipotriol", Drug)) %>%
  mutate(Drug = gsub("Estropipate", "Estrone sulfate", Drug)) %>%
  mutate(Drug = gsub("Calcipotriene", "Calcipotriol", Drug)) %>%
  mutate(Drug = gsub("Medroxyprogesterone acetate", "Medroxyprogesterone", Drug)) %>%
  mutate(Drug = gsub("Medroxyprogesterone", "Medroxyprogesterone acetate", Drug)) %>%
  mutate(Drug = gsub("Hydroxyprogesterone", "Hydroxyprogesterone caproate", Drug)) %>%
  mutate(Drug = gsub("Isometheptene", "Isometheptene mucate", Drug)) %>%
  mutate(Drug = gsub("Ergoloid mesylates", "Ergoloid mesylate", Drug)) %>%
  mutate(Drug = gsub("Clindamycin phosphate", "Clindamycin", Drug)) %>%
  mutate(Drug = gsub("Meclofenamate", "Meclofenamic acid", Drug)) %>%
  mutate(Drug = gsub("Norethindrone acetate", "Norethisterone", Drug)) %>%
  mutate(Drug = gsub("Norethindrone", "Norethisterone", Drug)) %>%
  mutate(Drug = gsub("Leuprolide acetate", "Leuprolide", Drug)) %>%
  mutate(Drug = gsub("Metaproterenol", "Orciprenaline", Drug)) %>%
  mutate(Drug = gsub("Demeclocylcine", "Demeclocycline",  Drug)) %>%
  mutate(Drug = gsub("Levalbuterol", "Levosalbutamol", Drug)) %>%
  mutate(Drug = gsub("Haloperidol decanoate", "Haloperidol", Drug)) %>%
  mutate(Drug = gsub("Ethynodiol", "Ethynodiol diacetate", Drug)) %>%
  mutate(Drug = gsub("Hyoscynamine", "Hyoscyamine", Drug)) %>%
  mutate(Drug = gsub("Flecainide acetate", "Flecainide", Drug)) %>%
  mutate(Drug = gsub("Cyproterone", "Cyproterone acetate", Drug)) %>%
  mutate(Drug = gsub("Tartaric acid", "D-tartaric acid", Drug)) %>%
  mutate(Drug = gsub("Isometheptene mucate", "Isometheptene", Drug)) %>%
  mutate(Drug = gsub(" desiccated", "", Drug)) %>%
  mutate(Drug = gsub("Thyroid", "Thyroid, porcine", Drug)) %>%
  mutate(Drug = gsub("Belladonna", "Atropine", Drug)) %>%
  mutate(Drug = gsub("Cannabis", "Cannabidiol", Drug)) %>%
  mutate(Drug = gsub("Opium", "Morphine", Drug)) %>%
  mutate(Drug = if_else(str_detect(Drug, regex("Erythromycin", ignore_case = TRUE)),
                        "Erythromycin", Drug)) %>%
  mutate(Drug = if_else(str_detect(Drug, regex("insulin", ignore_case = TRUE)),
                        "Insulin human", Drug)) %>%
  mutate(Drug = str_replace(Drug, "(?i)^(Sodium).*", "\\1")) %>%
  # ROA and unspecified agents
  filter(!str_detect(Drug, "topical|nasal|otic|ophthalmic|unspecified")) %>%
  left_join(., antichol.scale.df, "Drug")

ISACS.NHANES %>%
  summarise(across(ACSBC:ISACS,
                   list(missing = ~ sum(is.na(.)),
                        filled  = ~ sum(!is.na(.))))) %>%
  pivot_longer(everything(),
               names_to = c("variable", ".value"),
               names_sep = "_(?=missing$|filled$)")

View(ISACS.NHANES %>% filter(is.na(ISACS)))

ISACS.NHANES.score <- RXQ_RX_df %>%
  left_join(., ISACS.NHANES, "RXDDRUG") %>%
  group_by(SEQN) %>%
  summarise(ISACS.score = sum(ISACS, na.rm = TRUE),
            ACSBC.score = sum(ACSBC, na.rm = TRUE),
            ADS.score = sum(ADS, na.rm = TRUE),
            ABS.score = sum(ABS, na.rm = TRUE),
            ALS.score = sum(ALS, na.rm = TRUE),
            ORCA.score = sum(ORCA, na.rm = TRUE),
            AEC.score = sum(AEC, na.rm = TRUE),
            Chew.score = sum(`Chew's SAA`, na.rm = TRUE),
            ACB.score = sum(ACB, na.rm = TRUE),
            AAS.score = sum(AAS, na.rm = TRUE),
            ARS.score = sum(ARS, na.rm = TRUE),
            ABC.score = sum(ABC, na.rm = TRUE),
            ATS.score = sum(ATS, na.rm = TRUE))

# Calculate new medicine count as well as non-anticholinergic medicine counts
scales <- colnames(ISACS.NHANES[3:ncol(ISACS.NHANES)])

ISACS.NHANES.anticholinergic.only <- ISACS.NHANES %>% 
  mutate(across(all_of(scales), ~ ifelse(. != 0, 1, NA))) %>%
  filter(!if_all(all_of(scales), is.na)) %>% 
  select(-RXDDRUG) %>%
  distinct

anticholinergic.drugs <- RXQ_RX_df %>% 
  left_join(., ISACS.NHANES, "RXDDRUG") %>% distinct(SEQN, Drug) %>%
  left_join(., ISACS.NHANES.anticholinergic.only, "Drug") %>%
  group_by(SEQN) %>%
  summarise(ISACS.nMeds = sum(ISACS, na.rm = TRUE),
            ACSBC.nMeds = sum(ACSBC, na.rm = TRUE),
            ADS.nMeds = sum(ADS, na.rm = TRUE),
            ABS.nMeds = sum(ABS, na.rm = TRUE),
            ALS.nMeds = sum(ALS, na.rm = TRUE),
            ORCA.nMeds = sum(ORCA, na.rm = TRUE),
            AEC.nMeds = sum(AEC, na.rm = TRUE),
            Chew.nMeds = sum(`Chew's SAA`, na.rm = TRUE),
            ACB.nMeds = sum(ACB, na.rm = TRUE),
            AAS.nMeds = sum(AAS, na.rm = TRUE),
            ARS.nMeds = sum(ARS, na.rm = TRUE),
            ABC.nMeds = sum(ABC, na.rm = TRUE),
            ATS.nMeds = sum(ATS, na.rm = TRUE)
  ) 

nMeds.3 <- RXQ_RX_df %>% 
  left_join(., ISACS.NHANES, "RXDDRUG") %>% distinct(SEQN, Drug) %>%
  count(SEQN) %>% 
  dplyr::rename("RXDCOUNT"=n) %>%
  left_join(., RXQ_RX_df %>% distinct(SEQN, RXDUSE), "SEQN") %>%
  mutate(RXDCOUNT = case_when(is.na(RXDCOUNT)~0, T~as.numeric(RXDCOUNT)),
         RXDCOUNT = case_when(RXDUSE==2~0, T~as.numeric(RXDCOUNT)),
         nMeds.class = case_when(
           RXDCOUNT == 0 ~ "No medication use",
           RXDCOUNT >= 1 & RXDCOUNT <= 4 ~ "1-4 medications",
           RXDCOUNT >= 5 & RXDCOUNT <= 9 ~ "Polypharmacy (5-9)",
           RXDCOUNT >= 10 ~ "Hyperpolypharmacy (≥10)"
         ) %>% factor(levels = c("No medication use",
                                 "1-4 medications", 
                                 "Polypharmacy (5-9)",
                                 "Hyperpolypharmacy (≥10)"))) %>%
  left_join(., anticholinergic.drugs, "SEQN") %>%
  mutate(across(all_of(colnames(anticholinergic.drugs)[-1]), ~ RXDCOUNT - ., .names = "{.col}_diff")) %>%
  select(-ends_with(".nMeds"))

# Compare medicine count data
right_join(data.frame(table(nMeds$RXDCOUNT)), 
           data.frame(table(nMeds.2$RXDCOUNT)), 
           "Var1") %>%
  right_join(., data.frame(table(nMeds.3$RXDCOUNT)), "Var1")

plot_grid(nMeds %>% ggplot(aes(x = RXDCOUNT, fill = nMeds.class)) +
            geom_histogram(color = "black", alpha = 0.75, binwidth = 1) +
            geom_vline(xintercept = c(4.5, 9.5), linetype = "dashed", color = "black") +
            labs(x = "Number of medications",
                 y = "Number of participants",
                 fill = "", title = "Raw: Distribution of medicine users",
                 subtitle = paste0("n = ", nrow(nMeds))) +
            scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
            scale_fill_bmj() +
            mytheme +
            theme(legend.position = c(0.99, 0.98),
                  legend.background = element_blank(),
                  legend.justification = c(1, 0.75)), 
          nMeds.3 %>%
            ggplot(aes(x = RXDCOUNT, fill = nMeds.class)) +
            geom_histogram(color = "black", alpha = 0.75, binwidth = 1) +
            geom_vline(xintercept = c(4.5, 9.5), linetype = "dashed", color = "black") +
            labs(x = "Number of medications",
                 y = "Number of participants",
                 fill = "", title = "After cleaning: Distribution of medicine users",
                 subtitle = paste0("n = ", nrow(nMeds.2))) +
            scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
            scale_fill_bmj() +
            mytheme +
            theme(legend.position = c(0.99, 0.98),
                  legend.background = element_blank(),
                  legend.justification = c(1, 0.75)), ncol = 1)

# Calculate PhenoAge & PhenoAgeAccel
library(BioAge)

kdm = kdm_nhanes(biomarkers=c("albumin","alp","lncrp","totchol","lncreat","hba1c","sbp","bun","uap","lymph","mcv","wbc"))
phenoage = phenoage_nhanes(biomarkers=c("albumin_gL","alp","lncrp","totchol","lncreat_umol","hba1c","sbp","bun","uap","lymph","mcv","wbc"))

phenoage.dat <- phenoage$data %>%
  inner_join(., kdm$data %>% select(sampleID, kdm_advance), "sampleID") %>%
  select(-year) %>%
  separate(sampleID, c("cycle", "SEQN")) %>%
  filter(SEQN %in% demo_all$SEQN)

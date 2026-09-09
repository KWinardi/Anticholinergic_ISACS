# Network setup

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
library(readr)
library(ggridges)

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

source("0_HelperFunctions.R")

detectCores()
n_cores <- max(1, parallel::detectCores() - 2)

# Load network
## PPI network
PPI.df <- read.csv("2_PPI/2025_PPI_network.csv")

# PDI network
## DrugBank
DrugBank.df <- read.csv("1_PDI/1_DrugBank_v5.1.13/DrugBank.df.csv")
merged.PDI.df <- DrugBank.df %>%
  dplyr::filter(Pharm.Protein %in% "Target") %>% 
  dplyr::select(Name, Gene.Name) %>% 
  drop_na %>% unique %>% 
  dplyr::rename("Protein_A_SYMBOL" = "Name") %>% # Have all drugs on that column
  dplyr::rename("Protein_B_SYMBOL" = "Gene.Name")

## HPA (expression data)
HPA <- read.delim("6_HPA/tissue_category_rna_Any_Tissue.tsv") %>%
  dplyr::select(-c(starts_with("Cancer"), starts_with("RNA.pig"),
                   starts_with("RNA.mouse"), starts_with("RNA.cancer"),
                   starts_with("RNA.blood"), starts_with("RNA.cell"))) %>%
  dplyr::filter(Brain.expression.cluster != "")

HPA.anno <- bitr(HPA$Ensembl, 
                 fromType = "ENSEMBL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

HPA <- left_join(HPA, HPA.anno, c("Ensembl" = "ENSEMBL"))

## Include brain expression data?
PPI.df2 <- PPI.df %>% 
  left_join(., HPA %>% dplyr::select(Gene) %>% unique %>%
              mutate(Node1 = 1),
            c("Protein_A_SYMBOL" = "Gene")) %>%
  left_join(., HPA %>% dplyr::select(Gene) %>% unique %>%
              mutate(Node2 = 1),
            c("Protein_B_SYMBOL" = "Gene")) %>%
  mutate(Brain.Exp = Node1+Node2) %>%
  dplyr::filter(Brain.Exp == 2)
merged.PDI.df2 <- merged.PDI.df %>% 
  dplyr::filter(Protein_B_SYMBOL %in% 
                  unique(c(PPI.df2$Protein_A_SYMBOL, PPI.df2$Protein_B_SYMBOL))) %>% 
  left_join(., HPA %>% dplyr::select(Gene) %>% unique %>%
              mutate(Node2 = 1),
            c("Protein_B_SYMBOL" = "Gene")) %>%
  drop_na(Node2) 

## Important check!!!
all(unique(merged.PDI.df2$Protein_B_SYMBOL) %in% 
      unique(c(PPI.df2$Protein_A_SYMBOL, PPI.df2$Protein_B_SYMBOL)))

# No of edges & proteins
nrow(PPI.df)
length(unique(c(PPI.df$Protein_A_SYMBOL, PPI.df$Protein_B_SYMBOL)))

nrow(PPI.df2 %>%
       mutate(pair = paste(pmin(Protein_A_SYMBOL, Protein_B_SYMBOL), pmax(Protein_A_SYMBOL, Protein_B_SYMBOL), sep = "_")) %>%
       distinct(pair, .keep_all = TRUE))
length(unique(c(PPI.df2$Protein_A_SYMBOL, PPI.df2$Protein_B_SYMBOL)))

merged.PDI.df2 %>% group_by(Protein_A_SYMBOL) %>% count() %>%
  mutate(n = ifelse(n > 25, 25, n)) %>%
  ggplot(aes(x = n)) + 
  geom_histogram(aes(alpha=..count..), fill="#0AD1C1", bins = 25, color = "black") + 
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(x = "Number of drug interactors") + mytheme

# Graph setup
g <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                     Protein_B_SYMBOL), 
                           directed = F)
g <- igraph::simplify(g)

write_clip(merged.PDI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                            Protein_B_SYMBOL))

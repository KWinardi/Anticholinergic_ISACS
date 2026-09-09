# Constructing PPI network
# 25/09/2025 
# Kevin Winardi

# Load library
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load network
## PPI network
### Gross et al., 2022
PPI.df.Gross <- read.csv("2_PPI/2022Gross_HumanInteractome_SYMBOL/PPI_2022.csv") %>%
  mutate(Source = gsub("InateDB","InnateDB",Source),
         Source = gsub("Insider","INSIDER",Source),
         Source = gsub("Instruct","INstruct",Source),
         Source = gsub("Interactome3d","Interactome3D",Source),
         Source = gsub("BioGrid","BioGRID",Source),
         Source = gsub("InWeb","InWeb_IM",Source))
colnames(PPI.df.Gross) <- c("Protein_A_SYMBOL", "Protein_B_SYMBOL", "Gross.2022.Source")

### Zhou et al., 2020
PPI.df.Zhou <- read.delim("2_PPI/2020Zhou_HumanInteractome.tsv") %>%
  dplyr::select(-Type)
colnames(PPI.df.Zhou) <- c("Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "Source")
ENTREZ.anno <- bitr(unique(c(PPI.df.Zhou$Protein_A_Entrez_ID, PPI.df.Zhou$Protein_B_Entrez_ID)), 
                    fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
PPI.df.Zhou <- PPI.df.Zhou %>%
  mutate(Source = gsub("HI-II-14_Net","HI-Union",Source)) %>%
  mutate(across(everything(), as.character)) %>%
  left_join(., ENTREZ.anno, c("Protein_A_Entrez_ID" = "ENTREZID")) %>%
  dplyr::rename("Protein_A_SYMBOL" = "SYMBOL") %>%
  left_join(., ENTREZ.anno, c("Protein_B_Entrez_ID" = "ENTREZID")) %>%
  dplyr::rename("Protein_B_SYMBOL" = "SYMBOL") %>%
  dplyr::rename("Zhou.2020.Source" = "Source") %>%
  dplyr::select(Protein_A_SYMBOL, Protein_B_SYMBOL, Zhou.2020.Source) %>%
  group_by(Protein_A_SYMBOL, Protein_B_SYMBOL) %>%
  summarise(across(Zhou.2020.Source, ~ paste(unique(.), collapse = ",")), 
            .groups = "drop") %>%
  # Rename to match PPI.df.Gross
  mutate(Zhou.2020.Source = gsub("BioPlex2.0","BioPlex",Zhou.2020.Source),
         Zhou.2020.Source = gsub("PhosphositePlus","PhosphoSitePlus",Zhou.2020.Source),
         Zhou.2020.Source = gsub("SignaLink2.0","Signalink",Zhou.2020.Source))

PPI.df.Gross %>%
  pull(Gross.2022.Source) %>%                 
  str_split("\\|") %>%           
  unlist() %>%                   
  discard(~ .x == "") %>%        
  unique() %>%                   
  sort()  
PPI.df.Zhou %>%
  pull(Zhou.2020.Source) %>%                 
  str_split(",") %>%           
  unlist() %>%                   
  discard(~ .x == "") %>%        
  unique() %>%                   
  sort()

PPI.df <- bind_rows(PPI.df.Gross %>% dplyr::select(Protein_A_SYMBOL, Protein_B_SYMBOL, Gross.2022.Source),
                    PPI.df.Zhou %>% dplyr::select(Protein_A_SYMBOL, Protein_B_SYMBOL, Zhou.2020.Source)) %>%
  # Remove duplicates (A-B is the same as B-A)
  mutate(node_min = pmin(Protein_A_SYMBOL, Protein_B_SYMBOL),
         node_max = pmax(Protein_A_SYMBOL, Protein_B_SYMBOL),
         Gross.2022.Source = gsub("\\|",",", Gross.2022.Source)) %>%
  group_by(node_min, node_max) %>%
  summarise(Gross.2022.Source = paste(unique(na.omit(Gross.2022.Source)), collapse = ","),
            Zhou.2020.Source  = paste(unique(na.omit(Zhou.2020.Source)), collapse = ","),
            sources = list(unique(trimws(unlist(strsplit(paste(Gross.2022.Source, Zhou.2020.Source, sep = ","),
                                                         ","))))),
            Evidence.weight = length(sources[[1]]),
            .groups = "drop"
  ) %>%
  dplyr::rename(Protein_A_SYMBOL = node_min, Protein_B_SYMBOL = node_max) %>%
  mutate(sources = paste(sources[[1]], collapse = ",")) %>%
  drop_na(Protein_A_SYMBOL) %>%
  drop_na(Protein_B_SYMBOL)

# Number of nodes and edges
length(unique(c(PPI.df$Protein_A_SYMBOL, PPI.df$Protein_B_SYMBOL)))
nrow(PPI.df)

# Check for dups
PPI.df %>%
  group_by(Protein_A_SYMBOL, Protein_B_SYMBOL) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()

PPI.df %>% 
  ggplot(aes(x = Evidence.weight)) + 
  geom_histogram(aes(alpha= after_stat(count)), fill="#0AD1C1", bins = 25, color = "black") + 
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(x = "Number of evidence")

# Need to invert scale
PPI.df <- PPI.df %>% mutate(inv.Evidence.weight = 1/Evidence.weight)

write.csv(PPI.df, "2_PPI/2025_PPI_network.csv", row.names = F)
# write_clip(PPI.df)

# Weighted
g <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                     Protein_B_SYMBOL), 
                           directed = F)
is_simple(g)

# Unweighted (none)
E(g)$weight <- PPI.df2$Evidence.weight
g <- igraph::simplify(g)
# Inversed weight (linear)
g1 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g1)$weight <- 1 + (20 - PPI.df2$Evidence.weight)/20
g1 <- igraph::simplify(g1)
# Inversed weight (linear piecewise)
g2 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g2)$weight <- ifelse(PPI.df2$Evidence.weight > 5, 1, 1 + (5 - PPI.df2$Evidence.weight)/5)
g2 <- igraph::simplify(g2)
# Inversed weight (hyperbola)
g3 <- graph_from_data_frame(PPI.df2 %>% dplyr::select(Protein_A_SYMBOL,
                                                      Protein_B_SYMBOL), 
                            directed = F)
E(g3)$weight <- PPI.df2$inv.Evidence.weight + 1
g3 <- igraph::simplify(g3)

Drug <- "Allopurinol"
# Unweighted
distances(g, v = merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == Drug) %>% pull(Protein_B_SYMBOL),
          to = AChR.List, weights = NA)
# Weighted (raw)
distances(g, v = merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == Drug) %>% pull(Protein_B_SYMBOL),
          to = AChR.List)
# Weighted (inverse: linear)
distances(g1, v = merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == Drug) %>% pull(Protein_B_SYMBOL),
          to = AChR.List) %>% round(digits = 2)
# Weighted (inverse: linear piecewise)
distances(g2, v = merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == Drug) %>% pull(Protein_B_SYMBOL),
          to = AChR.List, weights = NA)
# Weighted (inverse: hyperbola)
distances(g3, v = merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == Drug) %>% pull(Protein_B_SYMBOL),
          to = AChR.List) %>% round(digits = 2)

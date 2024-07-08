library(tidyverse)
library(phyloseq)
library(microbiome)
library(parallel)
library(foreach)
library(doParallel)
library(MASS)
library(broom)
library(mediation)
library(wesanderson)
library(PNWColors)
library(gridExtra)
library(ggthemes)
library(ggpubr)
library(here)
library(randomForest)
library(caret)
library(ggbeeswarm)
library(randomForestExplainer)
library(ggsignif)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggforce)
library(ggnewscale)
library(cowplot)

# load physeq for positive metabolite data
met_pos_physeq <- readRDS("data/met_pos_ps.rds")
met_pos_ps <- microbiome::transform(met_pos_physeq, "compositional")

# select only the data that is from day 32 (after worm exposure, has defined worm burden, and was exposed)
D32_physeq_data <- subset_samples(met_pos_ps, day==32 & (!(is.na(total_worm))) & pcap==1)

# prune non-present taxa
cleaned_physeq_data <- prune_taxa(taxa_sums(D32_physeq_data)!=0, D32_physeq_data)

# extract the otu table, the taxonomy table, and the sample metadata
ASV_table <- data.frame(otu_table(cleaned_physeq_data))
ASV_tax_table <- data.frame(tax_table(cleaned_physeq_data)) %>%
  mutate(asv = rownames(ASV_tax_table))
positive_metadata <- data.frame(sample_data(cleaned_physeq_data))

# read in metabolite data for processing and selection
met_positive <- readRDS("data/met.pos.rds") %>%
  mutate_all(as.numeric)
metabolite_id_and_names <- read.csv(header=TRUE, "data/POS_MSMS_names.csv", na.strings=c("","NA")) %>%
  mutate_all(list(~na_if(.,""))) %>%
  mutate_all(list(~na_if(.," ")))

# manually investigated all metabolites to identify only those for which there is some biological information, otherwise
# they are functionally unannotated as nothing is known about them, and experimental validation may not be possible
poorly_understood_metabolites <- readRDS("data/poorly_understood_metabolites.rds")

unique_named_metabolites <- metabolite_id_and_names %>%
  filter(!(is.na(Accepted.Description))) %>%
  group_by(Accepted.Compound.ID) %>%
  mutate(unique_metab_id = seq(1, n())) %>%
  ungroup() %>%
  filter(unique_metab_id==1) %>%
  filter(!(Accepted.Description %in% poorly_understood_metabolites)) # extra step to remove poorly understood metabolites

nonunique_metabolites <- metabolite_id_and_names %>%
  filter(!(is.na(Accepted.Description))) %>%
  group_by(Accepted.Compound.ID) %>%
  mutate(unique_metab_id = seq(1, n())) %>%
  ungroup() %>%
  filter(unique_metab_id!=1)

# pseudo log transformation of metabolite data
met_positive_plus_1 <- met_positive + 1
met_pos <- met_positive_plus_1 %>%
  mutate_all(log) %>%
  mutate(sample_names = rownames(met_positive_plus_1)) %>%
  filter(sample_names %in% positive_metadata$sample_id) %>%
  dplyr::select(-sample_names)

met_pos_for_days <- met_positive_plus_1 %>%
  mutate_all(log) %>%
  mutate(sample_names = rownames(met_positive_plus_1)) %>%
  dplyr::select(-sample_names) %>%
  dplyr::select(all_of(unique_named_metabolites$Compound)) %>%
  dplyr::select(-"30.67_622.6098m/z",-"4.07_284.0996m/z")

# select only the named metabolites, this yields a set of 303 metabolites that are annotated, but we need to remove duplicate metabs
# note that these two metabolites were removed because they did not show up across any samples in the pared down data

metabolites_with_names <- met_pos %>%
  dplyr::select(all_of(unique_named_metabolites$Compound)) %>%
  dplyr::select(-"30.67_622.6098m/z",-"4.07_284.0996m/z")

worm_burden <- positive_metadata$total_worm
pcap_treatment <- positive_metadata$pcap

sample_data_all_days <- data.frame(sample_data(met_pos_ps)) %>%
  dplyr::select(sample_id, indiv_id, day, pcap, drug, total_worm) %>%
  mutate(actively_infected = ifelse(total_worm>0, "Actively Infected", "No Detected Worm Burden"))
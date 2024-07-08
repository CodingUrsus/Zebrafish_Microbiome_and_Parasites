## D3 analysis of NAEs (Day 0 in manuscript)
# clean D3 data
D3_physeq_data <- subset_samples(met_pos_ps, day==3) # select D3 data

cleaned_D3_data <- prune_taxa(taxa_sums(D3_physeq_data)!=0, D3_physeq_data) # remove non-present taxa

## extract the otu table, the taxonomy table, and the sample metadata
D3_otu_table <- data.frame(otu_table(cleaned_D3_data))
D3_tax_table <- data.frame(tax_table(cleaned_D3_data)) 
D3_metadata <- data.frame(sample_data(cleaned_D3_data)) %>%
  mutate(pcap_EXP = ifelse(group %in% c("PC", "PD"), 1, 0))

## Understand if D3 starting state impacts successional dynamics
bc_dist <- vegdist(D3_otu_table, method = "bray") # Calculate Bray-Curtis dissimilarity matrix

D3_metabs_with_id <- met_pos_for_days %>%
  rownames_to_column(var="sample_id") %>%
  filter(str_detect(sample_id, "_D3"))

# here are the NAE names and ids
D3_NAE <- data.frame("OEA" = D3_metabs_with_id$`26.12_325.2976n`,
                     "Anandamide (22:6, n-3)" = D3_metabs_with_id$`26.19_371.2817n`,
                     "LEA" = D3_metabs_with_id$`25.47_323.2821n`,
                     "alpha-LEA" = D3_metabs_with_id$`24.87_321.2659n`,
                     "Anandamide" = D3_metabs_with_id$`29.80_348.2890m/z`,
                     "Glycerophospho-N-Oleoyl Ethanolamine" = D3_metabs_with_id$`22.33_480.3068m/z`,
                     "Glycerophospho-N-Palmitoyl Ethanolamine" = D3_metabs_with_id$`21.70_453.2847n`,
                     "2-Linoleoyl Glycerol" = D3_metabs_with_id$`23.69_354.2761n`,
                     "sample_id" = D3_metabs_with_id$sample_id) %>%
  left_join(D3_metadata %>%
              dplyr::select(c("sample_id", "drug")))

# Assuming your data frame of NAEs is named 'D3_NAE' and the distance matrix is 'bc_dist':  
relevant_cols <- colnames(D3_NAE)[!(colnames(D3_NAE) %in% c("sample_id", "drug", "bc_dist"))] 

# Initialize an empty data frame to store results BEFORE the loop
D3_NAE_results_df <- data.frame(column = character(),
                                p_value = numeric(),
                                interaction_p_value = numeric())
# Loop through the columns
for (col in relevant_cols) {
  # Perform Adonis 2 test
  adonis_result <- extract_adonis_results(data = D3_NAE, column_name = col, bc_dist)
  D3_NAE_results_df <- rbind(D3_NAE_results_df, adonis_result)
}

# resulting D3 NAE associations are available at D3_NAE_results_df

# now do the same with D32 data
D32_physeq_data <- subset_samples(met_pos_ps, day==32) # select D322 data

cleaned_D32_data <- prune_taxa(taxa_sums(D32_physeq_data)!=0, D32_physeq_data) # remove taxa that sum to 0

## extract the otu table, the taxonomy table, and the sample metadata
D32_otu_table <- data.frame(otu_table(cleaned_D32_data))
D32_tax_table <- data.frame(tax_table(cleaned_D32_data)) 
D32_metadata <- data.frame(sample_data(cleaned_D32_data)) %>%
  mutate(pcap_EXP = ifelse(group %in% c("PC", "PD"), 1, 0))

## Understand if D32 starting state impacts successional dynamics

bc_dist32 <- vegdist(D32_otu_table, method = "bray") # Calculate Bray-Curtis dissimilarity matrix
drug_EXP_adonis <- adonis2(bc_dist32~drug*pcap_EXP, permutations=10000, data = D32_metadata) # look for interaction

D32_metabs_with_id <- met_pos_for_days %>%
  rownames_to_column(var="sample_id") %>%
  filter(str_detect(sample_id, "_Nx"))

D32_NAE <- data.frame("OEA" = D32_metabs_with_id$`26.12_325.2976n`,
                      "Anandamide (22:6, n-3)" = D32_metabs_with_id$`26.19_371.2817n`,
                      "LEA" = D32_metabs_with_id$`25.47_323.2821n`,
                      "alpha-LEA" = D32_metabs_with_id$`24.87_321.2659n`,
                      "Anandamide" = D32_metabs_with_id$`29.80_348.2890m/z`,
                      "Glycerophospho-N-Oleoyl Ethanolamine" = D32_metabs_with_id$`22.33_480.3068m/z`,
                      "Glycerophospho-N-Palmitoyl Ethanolamine" = D32_metabs_with_id$`21.70_453.2847n`,
                      "2-Linoleoyl Glycerol" = D32_metabs_with_id$`23.69_354.2761n`,
                      "sample_id" = D32_metabs_with_id$sample_id) %>%
  left_join(D32_metadata %>%
              dplyr::select(c("sample_id", "drug")))


# Assuming your data frame is named 'D32_NAE' and your distance matrix column is 'bc_dist':  
# Select potentially relevant columns 
relevant_cols <- colnames(D32_NAE)[!(colnames(D32_NAE) %in% c("sample_id", "drug", "bc_dist"))] 

# Initialize an empty data frame to store results BEFORE the loop
D32_NAE_results_df <- data.frame(column = character(),
                                 p_value = numeric(),
                                 interaction_p_value = numeric())

# Loop through the columns
for (col in relevant_cols) {
  # Perform Adonis 2 test
  adonis_result <- extract_adonis_results(data = D32_NAE, column_name = col, bc_dist32)
  D32_NAE_results_df <- rbind(D32_NAE_results_df, adonis_result)
}

colnames(D32_NAE_results_df) <- c("Compound","dist_NAE_p_value", "NAE:drug_p_value")
colnames(D3_NAE_results_df) <- c("Compound","dist_NAE_p_value", "NAE:drug_p_value")
combined_NAE_results_table <- rbind((D32_NAE_results_df %>%
                                       mutate(Day="29")), 
                                    (D3_NAE_results_df %>%
                                       mutate(Day="0")))

NAEs <- c("23.69_354.2761n", "21.70_453.2847n", "24.87_321.2659n", "26.19_371.2817n", "26.12_325.2976n", "25.47_323.2821n", "29.80_348.2890m/z", "22.33_480.3068m/z", "actively_infected")
my_NAEs <- c("Glycerophospho-N-Oleoyl Ethanolamine","Glycerophospho-N-Palmitoyl Ethanolamine", "Arachidonoyl Ethanolamide", "Anandamide (22:6, n-3)", "Oleoyl Ethanolamide", "Linoleoyl Ethanolamide", "2-Linoleoyl Glycerol", "α-Linolenoyl Ethanolamide")


names_and_ids <- unique_named_metabolites %>%
  dplyr::select(Compound, Accepted.Description)


NAE_df <- met_pos_for_days %>%
  mutate(sample_id = rownames(met_pos_for_days)) %>%
  left_join(sample_data_all_days) %>%
  filter(!(is.na(total_worm))) %>%
  filter(day==32) %>%
  dplyr::select(-c(day, sample_id, indiv_id, pcap, drug, total_worm)) %>%
  dplyr::select(NAEs) %>%
  pivot_longer(-actively_infected,
               names_to = "Compound",
               values_to="log_metab") %>%
  left_join(names_and_ids) %>%
  mutate(Accepted_Description = gsub(" ", "\n", Accepted.Description))

NAE_distinction_plot<- NAE_df %>%
  ggplot(aes(x=actively_infected, y=log_metab, color=actively_infected)) +
  geom_jitter(width=0.2) +
  geom_boxplot(alpha=0.1, outlier.shape = NA) +
  facet_wrap(~Accepted_Description, nrow=1,
             scales="free_y") +
  theme_clean() +
  geom_signif(comparisons = list(c("Actively Infected", "No Detected Worm Burden")),
              map_signif_level = T,
              y_position=c(10.7),
              color="black") +
  coord_cartesian(ylim = c(0, 11)) +
  scale_color_manual(values=c("#D95F02","#1B9E77")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(color = FALSE) +
  xlab(NULL) +
  ylab("Log Metabolite Abundance") +
  ggtitle("N-acylethanolamine and NAE-Precursor Abundance Distinguishes Zebrafish With Active Parasite Infection")

placeholder_metab <- met_positive_plus_1 %>%
  mutate_all(log)

D32_metabs <- placeholder_metab[rownames(placeholder_metab)%in%D32_metadata$sample_id,] %>%
  rownames_to_column(var="SampleID") %>%
  dplyr::select(c(unique_named_metabolites$Compound, "SampleID")) %>%
  pivot_longer(!SampleID, names_to="Compound", values_to="abund") %>%
  left_join(unique_named_metabolites %>%
              dplyr::select(Compound, Accepted.Description)) %>%
  dplyr::select(-Compound) %>%
  as_tibble() %>%
  left_join(D32_metadata %>%
              dplyr::select(total_worm, pcap, drug) %>%
              mutate("SampleID" = D32_metadata$sample_id))

NAE_data <- D32_metabs %>%
  filter(Accepted.Description %in% my_NAEs) %>%
  mutate(active_infection = ifelse(total_worm>0, "Active Infection", "No Active Infection"))

NAE_sample_data <- D32_metadata %>%
  mutate(PCAP = ifelse(pcap==0, "Control", "Parasite Exposed"))

NAE_pre_matrix <- NAE_data %>%
  pivot_wider(names_from = Accepted.Description,
              values_from = abund) %>%
  column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c("active_infection", "pcap", "drug", "total_worm"))

parasite_group <- D32_metadata %>%
  filter(pcap==1)

no_parasite_group <- D32_metadata %>%
  filter(pcap==0)

parasite_hclust = hclust(dist(as.matrix(NAE_pre_matrix)[rownames(NAE_pre_matrix) %in% parasite_group$sample_id,], method="canberra"))
no_parasite_hclust = hclust(dist(as.matrix(NAE_pre_matrix)[rownames(NAE_pre_matrix) %in% no_parasite_group$sample_id,], method="canberra"))

NAE_matrix <- (scale(as.matrix(NAE_pre_matrix)))

NAE_scaled <- NAE_matrix


parasite_sample_order <- parasite_hclust$labels[parasite_hclust$order]
no_parasite_sample_order <- no_parasite_hclust$labels[no_parasite_hclust$order]

external_row_order <- readRDS("/raid1/home/micro/hammera/darwinR/Analytical_Workflows_Mediation_Project/sample_order.rds")

row_order <- match(c(parasite_sample_order, rev(no_parasite_sample_order)), rownames(NAE_pre_matrix))

row_order <- match(external_row_order, rownames(NAE_pre_matrix))

ht_list = Heatmap(NAE_scaled, col = colorRamp2(c(-1.5, 0, 1.5), c("indianred3", "white", "midnightblue")), 
                  name = "scaled_expr", #column_title = ("Relative Log Metabolite Abundance"),
                  show_column_names = TRUE, width = ncol(NAE_scaled)*unit(4, "mm"), 
                  height = nrow(NAE_scaled)*unit(2, "mm"),
                  heatmap_legend_param = list(title = "Z-Scaled\nMetabolite Abundance"), show_column_dend = F, show_row_dend = F,
                  show_row_names = T, column_names_gp = gpar(fontsize = 11),
                  #row_order = row_order,
                  cluster_rows = T, cluster_columns = T) +
  Heatmap(D32_metadata$total_worm, name = "Worm Count", width = unit(3.5, "mm"),
          heatmap_legend_param = list(title = "Worm Count"), col = colorRamp2(c(0, 13), c("whitesmoke", "royalblue")), column_names_gp = gpar(fontsize = 11))

D3_physeq_data <- subset_samples(met_pos_ps, day==3) # select D3 data
cleaned_D3_data <- prune_taxa(taxa_sums(D3_physeq_data)!=0, D3_physeq_data) # remove taxa that sum to 
D3_physeq_only_abx_exposed <- subset_samples(cleaned_D3_data, (!(is.na(total_worm)) & (group=="PD" | group=="PC"))) # get just the parasite exposed samples and those that weren't NA for total_worm

D3_physeq_abx_asv_cleaned <- prune_taxa(taxa_sums(D3_physeq_only_abx_exposed)!=0, D3_physeq_only_abx_exposed)

# Prepare data
otu_table <- as.data.frame(otu_table(D3_physeq_abx_asv_cleaned))  # Get ASV abundances
response <- as.numeric(sample_data(D3_physeq_abx_asv_cleaned)$total_worm)  # Get total_worm values
D3_pcap_exposed_metadata <- data.frame(sample_data(D3_physeq_abx_asv_cleaned))


# Build random forest model
rf <- randomForest(response ~ ., data = otu_table, ntree = 20000)

# Select 10 most important ASVs
importance_df <- data.frame(IncNodePurity = importance(rf)) %>%
  rownames_to_column(var = "ASV") %>%
  left_join(data.frame(tax_table(D3_physeq_abx_asv_cleaned)) %>%
              rownames_to_column(var="ASV")) %>%
  arrange(desc(IncNodePurity))

# Calculate variance explained
rsq <- rf$rsq[20000]  # 

# Select top 10 ASVs
head(importance_df, 10) %>%
  mutate(asv_genus = paste0(ASV, " (", Genus, ")")) %>%
  ggplot(aes(x = IncNodePurity, y = reorder(asv_genus, IncNodePurity))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    x = "Increase in Node Purity",
    y = "ASV ID",
    title = "Top 10 Most Predictive ASVs for Total Worm Count\nEXCLUSIVELY in No-ABX and PCAP-Exposed Hosts"
  ) +
  theme_bw() +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  ggtitle("Top 10 Most Predictive ASVs for Total Worm Count\nEXCLUSIVELY in No-ABX and PCAP-Exposed Hosts",
          subtitle=paste0("Percentage of variance explained by the model:", as.character(round(rsq * 100, 2)), "%\n"))
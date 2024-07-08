D32_physeq_data <- subset_samples(met_pos_ps, day==32 & (!(is.na(total_worm))))

# Extract OTU table and metadata
my_otu_table <- otu_table(D32_physeq_data)
metadata <- data.frame(sample_data(D32_physeq_data))

# Calculate Bray-Curtis dissimilarity matrix
D32_bray_dist <- vegdist((my_otu_table), method = "bray", na.rm=TRUE)

# Perform PCoA
pcoa <- (pcoa(D32_bray_dist))

# PERMANOVA analysis
D32_adonis_exp_and_drug <- adonis2(D32_bray_dist ~ pcap*drug, data = metadata, method = "bray", permutations=10000)
D32_adonis_burden_and_drug <- adonis2(D32_bray_dist ~ total_worm*drug, data = metadata, method = "bray", permutations=10000)

# Create a ggplot object with desired features
pcoa_df <- data.frame(pcoa$vectors) %>%
  rownames_to_column(var="sample_id") %>%
  full_join(metadata) %>%
  mutate(ABX = ifelse(drug==0, "No Antibiotic\nExposure", "Antibiotic Exposed"))

pcoa_cent <- aggregate(cbind(Axis.1, Axis.2)~drug, data=pcoa_df, FUN = mean, na.rm=TRUE)

pcoa_scores <- merge(pcoa_df, setNames(pcoa_cent, c("drug", 'oAxis.1', 'oAxis.2')), by = 'drug', sort=FALSE)

# Perform envfit analysis
envfit_results <- envfit(pcoa$vectors ~ total_worm:drug, data = metadata)

# Extract vector coordinates for visualization
envfit_arrows <- data.frame(envfit_results$vectors$arrows) %>%
  mutate(mlabel = c("Worm Burden", "Drug Exposure"))

myp <- ggplot(data = pcoa_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = total_worm), size = 14, alpha = 0.30) +
  scale_color_gradient(low = "white", high = "#625a94",
                       breaks = c(0, 4, 8, 12),
                       labels = c("0", "4", "8", "12"),
                       name = "Total Worm\nCount") +  # Set legend title
  new_scale_color() +
  geom_point(aes(color = ABX), size = 5) +
  scale_color_manual(values = c("royalblue4", "indianred"),
                     name="Antibiotics") +
  labs(x = paste0("PCoA 1 (", round(pcoa$values$Relative_eig[1] * 100, 2), "% variance)"),
       y = paste0("PCoA 2 (", round(pcoa$values$Relative_eig[2] * 100, 2), "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Add envfit vectors to the plot
comp1 <- myp +
  geom_segment(
    aes(x = 0, y = 0, xend = envfit_arrows[1, 1] * 0.7, yend = envfit_arrows[1, 2] * 0.7),
    arrow = arrow(length = unit(0.2, "cm")), color = "#625a94", linewidth = 0.5
  ) +
  geom_text(data=envfit_arrows %>% filter(mlabel=="Worm Burden"),
            aes(x = Axis.1 * 0.7, y = Axis.2* 0.7, label = mlabel),
            hjust = "left", vjust = "center", nudge_x = 0.05,
            size = 6,     # Adjust font size
            #    family = "sans",  # Change font family
            color = "#625a94"  # Choose a color that complements the plot
  ) +
  geom_point(aes(x=0,y=0), color="#625a94", size=1)
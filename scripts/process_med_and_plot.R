# neg bin modeling of worm burden ~ metabs, this can be thought of as the direct effect
glmnb_worm_metabolite_set <- data.frame(t(sapply(1:ncol(metabolites_with_names), function(x) glmnb_total_worm_on_feature(x, metabolites_with_names, colnames(metabolites_with_names), worm_burden))))

colnames(glmnb_worm_metabolite_set) <- c("glm.nb_coefficients", "p_values", "edge_1")

# get metabolites to merge
named_metabs <- unique_named_metabolites %>%
  mutate(edge_1 = Compound) %>%
  dplyr::select(-c("Compound", "Accepted.Compound.ID", "Compound.Link"))

glmnb_worm_metabolites_and_ids <- left_join(glmnb_worm_metabolite_set, named_metabs) %>%
  dplyr::select(-unique_metab_id) %>%
  filter(!(Accepted.Description %in% poorly_understood_metabolites))

metabolite_glmnb_with_adjusted_pvals <- glmnb_worm_metabolites_and_ids %>%
  mutate(adjusted_p_values = p.adjust(p_values, method="fdr")) %>%
  filter(adjusted_p_values < 0.1) %>%
  arrange(adjusted_p_values)

test_met_set <- (metabolites_with_names %>% dplyr::select(metabolite_glmnb_with_adjusted_pvals$edge_1))

colnames(test_met_set) <- make.names(colnames(test_met_set))


# apply a sparsity filter to select prevalent ASVs
sparsity_filter <- apply(ASV_table, 2, function(x) sum(x>0))

asvs_more_than_x_obs <- sparsity_filter[sparsity_filter>8]

sparse_selected_ASV_table <- ASV_table %>%
  dplyr::select(all_of(names(asvs_more_than_x_obs)))

test_asv_set <- sparse_selected_ASV_table

set.seed(1001)  # for reproducibility

# Set up parallel cluster
cl <- makeCluster(40)
registerDoParallel(cl)
clusterEvalQ(cl, library("mediation", lib.loc = "/raid1/home/micro/hammera/R/"))
clusterEvalQ(cl, library("nptest", lib.loc = "/raid1/home/micro/hammera/R/"))
clusterEvalQ(cl, library("dplyr", lib.loc = "/raid1/home/micro/hammera/R/"))
clusterExport(cl, c("parallel_calculate_med_and_cors", "test_met_set", "test_asv_set", "worm_burden"))

# Apply the function for each metabolite
results <- parLapply(cl, 1:ncol(test_met_set), function(i) {
  parallel_calculate_med_and_cors(test_met_set[, i], colnames(test_met_set)[i], test_asv_set, worm_burden)
})

# Close the parallel cluster
stopCluster(cl)

# Convert results to dataframe
res_df <- do.call(rbind, lapply(results, as.data.frame))

structured_results <- as.data.frame(lapply(data.frame(t(data.frame(results))), unlist), stringsAsFactors = FALSE)

mediating_relationships <- structured_results %>%
  filter(!(is.na(ADE_pval))) %>%
  mutate(simple_fdr = p.adjust(simple_p_value, method="fdr"),
         pcor_fdr = p.adjust(partial_p_value, method="fdr")) %>%
  mutate(ACME_fdr = p.adjust(ACME_pval, method="fdr")) %>%
  filter(ACME_fdr<0.3 & pcor_fdr<0.3) %>%
  left_join(named_metabs %>%
              mutate(metab_name = make.names(edge_1)) %>%
              dplyr::select(c("metab_name", "Accepted.Description"))) %>%
  left_join(ASV_tax_table %>%
              mutate(asv_name = asv))
left_join(glmnb_worm_genus_test %>%
            mutate(asv_name = asv) %>%
            dplyr::select(c("asv", "worm_on_asv_coefficient")))


sr1 <- mediating_relationships 

sr2 <- sr1 %>%
  left_join(mediating_relationships, by = c("metab_name", "asv_name"))

sr2 %>%
  ggplot(aes(x=ACME_pval.x, y=ACME_pval.y)) +
  geom_point()


l <- 1
r <- 5
range <- c(l,r)

forward_mediation_results_df <- mediating_relationships %>%
  dplyr::mutate(asv_designation = paste0(ifelse(is.na(Genus), ifelse(is.na(Family), Order, Family), Genus), " | ", asv)) %>%
  mutate(x0 = l,
         y0 = (match(Accepted.Description, unique(Accepted.Description))),
         x1 = r,
         y1 = (match(asv_designation, unique(asv_designation)))*((length(unique(Accepted.Description)))/(length(unique(asv_designation)))),
         mid = mean(range))

blue_to_red <- c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")
gold_to_purple <- c("coral1", "white", "midnightblue")


forward_mediation_results_plot <- ggplot(data=forward_mediation_results_df) +
  geom_tile(aes(x=x1, y=(y1)), height = 0.85, width = 0.7, fill="white", color="black") +
  geom_text(data = unique(forward_mediation_results_df %>%
                            dplyr::select(x1, y1, asv_designation)), aes(x = x1, y = y1, label = asv_designation),
            size = 3, fontface = "italic") +
  geom_tile(data = forward_mediation_results_df %>%
              dplyr::select(mid, y0), aes(x=mid, y=(y0)), width = 1.1, color="black", fill="white") +
  geom_text(data = unique(forward_mediation_results_df %>%
                            dplyr::select(y0, mid, Accepted.Description)), aes(x = mid, y = y0, label = Accepted.Description),
            size = 3, fontface = "italic") +
  geom_segment(data=forward_mediation_results_df, aes(x=x1-0.35, xend=mid+0.55, y=y1, yend=y0, color=as.numeric(simple_estimate))) +
  scale_color_gradientn(name="Correlation Estimate", colors=gold_to_purple) +
  scale_size_continuous(name="Proportion ASV-Worm\nBurden Mediated", range = c(0.05, 2)) +
  scale_alpha_continuous(guide="none") +
  theme_void() +
  theme(legend.position="bottom") +
  scale_x_continuous(
    breaks = c(unique(forward_mediation_results_df$x1), unique(forward_mediation_results_df$mid)),
    labels = c("Taxa", "Metabolites"),
    position = "top") +
  theme(axis.text.x=element_text(size = 20, face = "bold"))+
  theme(plot.title = element_text(size = 15, face = "bold")) +
  geom_tile(data = forward_mediation_results_df %>% 
              mutate(fill_color = ifelse(asv_designation == "Pelomonas | asv_4", "Taxon Inversely Associated\nwith Worm Burden", "white")) %>%
              filter(asv_designation == "Pelomonas | asv_4"),
            aes(x=x1, y=(y1), fill=fill_color),
            height = 0.85, width = 0.7, color="black", alpha=0.1) +
  scale_fill_manual(values = c("darkorange3")) +
  labs(fill="")
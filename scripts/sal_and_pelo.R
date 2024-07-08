pelomonas_salicylaldehyde_df <- data.frame(pelomonas = test_asv_set$asv_4,
                                           salicylaldehyde = test_met_set$X2.70_123.0433m.z,
                                           worm_burden = worm_burden)

pelomonas_and_worm_burden_plot <- pelomonas_salicylaldehyde_df %>%
  ggplot(aes(x=pelomonas, y=worm_burden)) +
  geom_point(size=1.5, alpha=0.8) +
  geom_smooth(method="lm", se=T, color="blue") +
  ylab("Total Worm Burden") +
  xlab("Pelomonas Relative Abundance") +
  theme_bw()+

  coord_cartesian(ylim = c(-0, max(pelomonas_salicylaldehyde_df$worm_burden)))+
  annotate("text", x=Inf, y=Inf, hjust = 1.05, vjust = 2, label=paste("p-value = 0.002, glm.nb coefficient = -13.53"))

salicylaldehyde_and_worm_burden_plot <- pelomonas_salicylaldehyde_df %>%
  ggplot(aes(y=worm_burden, x=(salicylaldehyde))) +
  geom_point(size=1.5, alpha=0.8) +
  geom_smooth(method="lm", se=T, color="blue") +
  ylab("Total Worm Burden") +
  xlab("log(Salicylaldehyde Abundance)") +
  theme_bw() +

  coord_cartesian(ylim = c(-0, max(pelomonas_salicylaldehyde_df$worm_burden)))+
  annotate("text", x=Inf, y=Inf, hjust = 1.05, vjust = 2, label=paste("p-value = 0.01, glm.nb coefficient = -0.35"))

pelomonas_and_salicylaldehyde_plot <- pelomonas_salicylaldehyde_df %>%
  ggplot(aes(x=pelomonas, y=salicylaldehyde)) +
  geom_point(size=1.5, alpha=0.8) +
  geom_smooth(method="lm", se=T, color="blue") +
  xlab("Pelomonas Relative Abundance") +
  ylab("log(Salicylaldehyde Abundance)") +
  theme_bw() +
  annotate("text", x=Inf, y=Inf, hjust = 1.05, vjust = 2, label=paste("p-value = 0.002, rho = 0.62"))



# Create empty plots with labels
plot_A <- salicylaldehyde_and_worm_burden_plot
plot_B <- pelomonas_and_worm_burden_plot
plot_C <- pelomonas_and_salicylaldehyde_plot

# Arrange the plots
tripanel <-  plot_grid(plot_A, plot_B, plot_C, ncol = 3, nrow=1, labels = c("A", "B", "C"), label_size=18)



larvation_data <- data.frame("Group" = c("Control\n(DMSO)","Treatment\n(2mg/L)","Control\n(DMSO)","Treatment\n(2mg/L)"),
                             "Larvate" = c(68, 0, 186, 0),
                             "Unlarvate" = c(45, 90, 68, 131),
                             "Pass" = c("Second", "Second", "First", "First")) %>%
  group_by(Group, Pass) %>%
  mutate(Larvated = (round(Larvate/(Larvate+Unlarvate), 3)*100),
         Unlarvated = (round(Unlarvate/(Larvate+Unlarvate), 3)*100)) %>%
  pivot_longer(cols = c(Larvated, Unlarvated), names_to = "HatchType", values_to = "Percent") %>%
  mutate(Hatch = ifelse(HatchType=="Larvated", "Larvated", "Unlarvated or\nDead"))

f4c <- larvation_data %>%
  filter(Pass=="First" & HatchType=="Unlarvated") %>%
  mutate(SampleGroup = ifelse(Group=="Control\n(DMSO)", "Control (DMSO)", "Salicylaldehyde\n(2mg/L)")) %>%
  ggplot(aes(x = SampleGroup, y = Percent, fill = SampleGroup)) +
  geom_bar(stat = "identity", width=0.4) +
  labs(x = NULL,
       y = NULL) +
  scale_fill_manual(values = c("#625a94","#f57946")) +
  theme_bw() +
  #  geom_segment(aes(x = 1, xend = 2, y = 103, yend = 103), color = "black") +
  # geom_text(aes(x = 1.5, y = 103.2, label = "***"), size = 4, vjust = 0) +
  labs(fill=NULL) +
  theme(axis.text.x = element_text(color = "black", size=11),
        axis.title.y = element_text(size = 11)) +
  ylab("Percent Eggs Dead\nor Unlarvated")


sal_count_data <- readRDS("data/sal_mature_count_data.rds")

sal_worm_percent_df <- sal_count_data %>%
  group_by(Group, FemaleWorms) %>%
  dplyr::summarize(FemaleWormCount_Number = n()) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(GroupTotal = sum(FemaleWormCount_Number)) %>%
  ungroup() %>%
  mutate(Female_Worm_Percent = (round(FemaleWormCount_Number/GroupTotal, 3))*100)


f4b <- sal_count_data %>%
  mutate(SampleGroup = ifelse(Group=="Control", "Control (DMSO)", "Salicylaldehyde\n(2mg/L)")) %>%
  ggplot(aes(x=FemaleWorms, fill=SampleGroup)) +
  geom_histogram(bins=4) +
  facet_wrap(~SampleGroup, labeller = labeller(NULL)) +
  scale_fill_manual(values = c("#625a94","#f57946")) +
  theme_bw() +
  ylab("Sample Count") +
  xlab("# Mature (Egg-Producing) Female Worms in Sample") +
  labs(fill=NULL) +
  theme(legend.position = "none")

fig_de <- plot_grid(f4c, f4b, ncol = 2, nrow=1, labels = c("D", "E"), label_size=18)

combined_plot <- plot_grid(
  plot_grid(f4b, f4c, ncol = 2, nrow=1, labels = c("D", "E"), label_size=18),
  ncol = 2,
  nrow=1,
  rel_heights = c(1, 1),
  rel_widths = c(1, 1),
  hjust = -0.2, 
  label_size=18  # Adjust horizontal justification to move the labels closer to the plots
)

sal_worm_abund_data <- readRDS("data/sal_worm_abund.rds")

administered_sal_worm_abund_plot <- sal_worm_abund_data %>%
  ggplot(aes(x=group, y=worm_count, color=group)) +
  geom_boxplot(alpha=0) +
  geom_jitter() +
  theme_bw() +
  xlab(NULL) +
  ylab("Worm Count") +
  labs(color=NULL) +
  scale_color_manual(values = c("#625a94","#f57946"))  +
  geom_signif(comparisons = list(c("Control (DMSO)", "Salicylaldehyde (2mg/L)")),
              map_signif_level = T,
              y_position=c(55),
              color="black",
              textsize = 5) +
  coord_cartesian(ylim = c(0, 58))


test_sample <- plot_grid(administered_sal_worm_abund_plot, f4b, ncol = 2, labels = c("4d", "4e"), label_size=18)

plot_4f <- plot_grid(f4c, ncol = 1, labels = c("4f"), label_size = 16)

larvation_data <- readRDS("data/larvation_data.rds")

# Perform Chi-Square Test
chi2_result <- chisq.test(larvation_data)

# Create the contingency table
gravid_contingency_table <- (readRDS("data/gravid_contingency_table.rds"))

# Perform Fisher's exact test
fishers_exact_result <- fisher.test(gravid_contingency_table)
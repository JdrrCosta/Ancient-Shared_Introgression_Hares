####Final RESULTS 
library(ggplot2)
library(dplyr)

# List of files and model names
files <- c("1000_rep1_MSC_out.mcmc.txt", 
           "1000_rep1_MSC-I_proto_out.mcmc.txt", 
           "1000_rep1_MSC-I_cas_proto_out.mcmc.txt",
           "1000_rep1_MSC-I_cor_proto_out.mcmc.txt",
           "1000_rep1_MSC-I_cascor_proto.mcmc.txt",
           "1000_rep1_MSC-I_cas_out.mcmc.txt",
           "1000_rep1_MSc-I_cor_out.mcmc.txt",
           "1000_rep1_MSC-I_cascor_out.mcmc.txt",
           "1000_rep1_MSC-M_out.mcmc.txt",
           "1000_rep1_MSC-M_proto_cas_out.mcmc.txt",
           "1000_rep1_MSC-I_cascor_proto_more_chains.mcmc.txt",
           "1000_rep1_MSC-I_cascor_proto_flipped_more_chains.mcmc.txt")

models <- c("No_geneflow", "Proto", "Proto_cas","Proto_cor", 
            "Only_cas","Only_cor","Only_cascor","MSC-M","MSC-M_proto_cas",
            "MSC-M_proto_cascor","Proto_cascor", "Proto_cascor_flipped")

# Read only lnL from each file and add model column
lnL_data <- lapply(seq_along(files), function(i) {
  read.table(files[i], header = TRUE) %>%
    select(lnL) %>%
    mutate(model = models[i])
}) %>% bind_rows()

# Reorder models by descending median lnL
model_order <- lnL_data %>%
  group_by(model) %>%
  summarise(median_lnL = median(lnL)) %>%
  arrange(desc(median_lnL)) %>%
  pull(model)

lnL_data <- lnL_data %>%
  mutate(model = factor(model, levels = model_order))

# Boxplot of lnL for the models 
plot <- ggplot(lnL_data, aes(x = model, y = lnL, fill = model)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_minimal(base_size = 24) +
  labs(
    title = "lnL Distributions for BPP Models",
    x = "Model",
    y = "Log-Likelihood (lnL)"
  ) +
  scale_fill_manual(values = c(
    "No_geneflow"   = "#ffba08",
    "Proto"         = "#90be6d",
    "Proto_cas"     = "#3d348b",
    "Proto_cor"     = "#01665e",
    "Proto_cascor"  = "#43aa8b",
    "Only_cas"      = "#f94144",
    "Only_cor"      = "#f3722c",
    "Only_cascor"   = "#577590",
    "MSC-M"         = "#277da1",
    "MSC-M_proto_cas" = '#7209b7',
    "Proto_cascor_more" = '#9cc1f9',
    "Proto_cascor_flipped" = '#b56576'
  )) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.text.y = element_text(color = "black"),      # show y-axis labels
    axis.ticks.y = element_line(color = "black"),     # show y-axis ticks
    panel.grid.major.y = element_line(color = "grey80"),  # add horizontal grid lines
    panel.grid.minor.y = element_blank(),            # optional: remove minor y grids
    panel.grid.major.x = element_blank()             # optional: remove vertical grids
  )

ggsave(plot = plot, filename = "All_models_likelihood_comparisons.svg", 
       dpi = 300, width = 35, height = 10)

plot

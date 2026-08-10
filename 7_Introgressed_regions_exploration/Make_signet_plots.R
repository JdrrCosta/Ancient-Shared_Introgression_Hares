library(ggplot2)
library(dplyr)
library(patchwork)

# ---------------------------
# Load data
# ---------------------------
signet_out_kegg <- read.delim("Cactus_full_genes_Signet_output_kegg.tsv", header = TRUE)
signet_out_reactome  <- read.delim("Cactus_full_genes_Signet_output_reactome.tsv", header = TRUE)

# ---------------------------
# Function to prepare data
# ---------------------------
prep_data <- function(df) {
  df %>%
    filter(p.val < 0.01) %>%
    arrange(desc(subnet.score)) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
}

significant_kegg <- prep_data(signet_out_kegg)
significant_reactome  <- prep_data(signet_out_reactome)

# ---------------------------
# Plotting function
# ---------------------------
plot_pathways <- function(data, title) {
  ggplot(data, aes(x = pathway, y = subnet.score)) +
    geom_col(fill = "#4A90E2", width = 0.7) +
    coord_flip() +
    labs(x = NULL, y = "Subnet Score", title = title) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ---------------------------
# Generate plots
# ---------------------------
p_kegg <- plot_pathways(significant_kegg, "Kegg: Significant Pathways")
p_reactome  <- plot_pathways(significant_reactome, "Reactome: Significant Pathways")

# ---------------------------
# Combine using patchwork
# ---------------------------
combined_plot <- p_kegg | p_reactome

combined_plot
ggsave("Kegg_reactome_both_cactus.pdf",combined_plot, dpi = 300, width = 20, height = 10)

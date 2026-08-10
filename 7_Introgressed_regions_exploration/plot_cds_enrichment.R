library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# --- Read data ---
boot <- readr::read_tsv("CDS_bootstrap_counts_bp_no_overlaps.tsv")
empirical <- readr::read_tsv("CDS_empirical_summary_bp_no_overlaps.tsv")

# --- Convert bootstrap to long format ---
boot_long <- tidyr::pivot_longer(
  boot,
  cols = c(CDS, Non_CDS),
  names_to = "Region",
  values_to = "Count"
)

# --- Extract empirical values ---
empirical_counts <- as.list(empirical[1, c("CDS", "Non_CDS")])

# --- Function to create density plot ---
plot_region_density <- function(region_name, fill_color) {
  
  emp_val <- empirical_counts[[region_name]]
  
  region_data <- dplyr::filter(boot_long, Region == region_name)
  
  dens <- density(region_data$Count)
  
  ggplot(region_data, aes(x = Count)) +
    
    geom_density(
      fill = fill_color,
      alpha = 0.6,
      color = "black",
      linewidth = 1
    ) +
    
    geom_vline(
      xintercept = emp_val,
      color = "red",
      linetype = "dashed",
      linewidth = 1.2
    ) +
    
    annotate(
      "text",
      x = emp_val,
      y = max(dens$y) * 0.9,
      label = paste0(round(emp_val / 1000, 1), "k"),
      color = "red",
      fontface = "bold",
      vjust = -0.2
    ) +
    
    scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    
    labs(
      title = region_name,
      x = "Basepair Count",
      y = "Density"
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = element_text(
        color = "#333333"
      ),
      axis.text = element_text(
        color = "#333333"
      )
    )
}

# --- Create plots ---
p_cds <- plot_region_density("CDS", "#94d6c4")
p_non_cds <- plot_region_density("Non_CDS", "#f4b183")

# --- Combine plots ---
combined_plot <- p_cds | p_non_cds

# --- Display ---
print(combined_plot)

# --- Save figure ---
ggplot2::ggsave(
  "bootstrap_vs_empirical_CDS_vs_nonCDS_density.svg",
  combined_plot,
  width = 8,
  height = 4
)

library(ggplot2)
library(patchwork)

# Function
plot_overlap_gg <- function(file, empirical, title,
                            col_below = "lightblue",
                            col_above = "orange") {
  
  boot <- scan(file)
  df <- data.frame(value = boot)
  
  ggplot(df, aes(x = value)) +
    
    geom_histogram(aes(fill = value < empirical),
                   bins = 30,
                   color = "black",
                   alpha = 0.9) +
    
    scale_fill_manual(
      values = c("TRUE" = col_below, "FALSE" = col_above),
      guide = "none"
    ) +
    
    geom_vline(xintercept = empirical,
               color = "red",
               linewidth = 1.2) +
    
    labs(title = title,
         x = "Overlap (%)",
         y = "Bootstrap Replicates") +
    
    # Fixed axis for ALL plots
    coord_cartesian(xlim = c(2, 7.5)) +
    
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

# ---- Cactus4 ----
p1 <- plot_overlap_gg(
  "Percentage_overlaps_cactus4_bootstrap_ngsParalog_Lcor.txt",
  2.44,
  "Cactus4 - Italian",
  "#ff7f0e", "#ff7f0e"
)

p2 <- plot_overlap_gg(
  "percentage_overlaps_cactus4_bootstrap_ngsParalog_Lcas.txt",
  3.58,
  "Cactus4 - Broom",
  "#f4a6b8ff", "#f4a6b8ff"
)

# ---- Cactus12 ----
p3 <- plot_overlap_gg(
  "percentage_overlaps_bootstrap_cactus12_ngsParalog_Lcor.txt",
  3.10,
  "Cactus12 - Italian",
  "#ff7f0e", "#ff7f0e"
)

p4 <- plot_overlap_gg(
  "percentage_overlaps_bootstrap_cactus12_ngsParalog_Lcas.txt",
  5.10,
  "Cactus12 - Broom",
  "#f4a6b8ff", "#f4a6b8ff"
)

# Combine into 2x2 layout
final_plot <- (p1 | p2) / (p3 | p4)

# Save as SVG
ggsave("overlap_comparison_4panel.svg",
       plot = final_plot,
       width = 12,
       height = 10)

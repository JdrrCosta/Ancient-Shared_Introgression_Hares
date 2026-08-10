library(ggplot2)

# ----------------------------
# Load data
# ----------------------------
bootstrap_high <- scan("cactus4_boostrap_sites_with_DI_10_to_20.txt")
bootstrap_all  <- scan("cactus4_boostrap_sites_NO_DI.txt")

empirical_high <- scan("cactus4_empirical_sites_with_DI_10_to_20.txt")
empirical_all  <- scan("cactus4_empirical_sites_NO_DI.txt")

# ----------------------------
# Calculate proportions
# ----------------------------
bootstrap_ratio <- bootstrap_high / bootstrap_all
empirical_ratio <- empirical_high / empirical_all

df <- data.frame(ratio = bootstrap_ratio)

# Density for annotation placement
dens <- density(df$ratio)
density_max <- max(dens$y)

# One-sided p-value
p_value <- mean(bootstrap_ratio >= empirical_ratio)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(df, aes(x = ratio)) +
  
  geom_density(
    fill = "#4C78A8",
    color = "#1F3552",
    alpha = 0.4,
    linewidth = 1.2,
    adjust = 1.2
  ) +
  
  geom_vline(
    xintercept = empirical_ratio,
    color = "#D62728",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  annotate(
    "label",
    x = empirical_ratio,
    y = density_max * 0.9,
    label = paste0(
      "Empirical ratio = ", round(empirical_ratio, 4), "\n",
      "(", empirical_high, "/", empirical_all, ")\n",
      "P-value = ", signif(p_value, 3)
    ),
    color = "#D62728",
    fill = "white",
    fontface = "bold",
    size = 4,
    label.size = 0.25
  ) +
  
  labs(
    x = "Proportion of DI sites between -10 and -20",
    y = "Density (bootstrap distribution)",
    title = "Cactus 4 Bootstrap vs empirical: proportion of DI sites between -10 and -20"
  ) +
  scale_x_continuous(
    limits = c(0, 0.60),
    breaks = seq(0, 0.60, by = 0.05),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey92"),
    plot.margin = margin(15, 20, 15, 15)
  )

p

# ----------------------------
# Save
# ----------------------------
ggsave(
  "Cactus4_bootstrap_vs_empirical_DI_10_to_20.svg",
  plot = p,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300
)

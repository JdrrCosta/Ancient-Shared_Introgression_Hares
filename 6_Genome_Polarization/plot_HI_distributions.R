library(ggplot2)

# Load data
bootstrap_values <- read.table("Cactus12_bootstrap_averages_no_auto.txt", header = FALSE)[,1]
empirical_value <- as.numeric(read.table("Cactus12_empirical_average_no_auto.txt", header = FALSE)[1,1])

# Transform HI values
bootstrap_values <- 1 - (bootstrap_values / 2)
empirical_value <- 1 - (empirical_value / 2)

df <- data.frame(value = bootstrap_values)

# Density max for annotation positioning
density_max <- max(density(bootstrap_values)$y)

# Plot
p <- ggplot(df, aes(x = value)) +
  
  # Filled density curve
  geom_density(
    fill = "#4C78A8",
    color = "#1F3552",
    alpha = 0.4,
    linewidth = 1.2,
    adjust = 1.2
  ) +
  
  # Empirical value line
  geom_vline(
    xintercept = empirical_value,
    color = "#D62728",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  # Label for empirical value
  annotate(
    "label",
    x = empirical_value,
    y = density_max * 0.9,
    label = paste0("Empirical HI = ", round(empirical_value, 3)),
    color = "#D62728",
    fill = "white",
    fontface = "bold",
    size = 4,
    label.size = 0.25
  ) +
  
  # Axis labels
  labs(
    x = "Average HI",
    y = "Density",
    title = "Distribution of Bootstrap HI Values"
  ) +
  
  # Nice limits and breaks
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  # Clean theme
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey92"),
    plot.margin = margin(15, 20, 15, 15)
  )

p

# Save as SVG
ggsave(
  "Cactus12_HI_distribution_no_auto.svg",
  plot = p,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300
)

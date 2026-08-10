# Load libraries
library(ggplot2)

# Read bootstrap values
candidate_cov <- scan("Bootstrap_Cactus4_bp_overlaps_seixas_genes.txt")
total_cov <- scan("Bootstrap_Cactus4_bp_overlaps_all_genes.txt")

# Check lengths
stopifnot(length(candidate_cov) == length(total_cov))

# Calculate bootstrap ratios
bootstrap <- candidate_cov / total_cov

# Replace with your empirical ratio
# e.g. empirical_candidate_cov / empirical_total_gene_cov
empirical <- 0.01

# Summary statistics
bootstrap_mean <- mean(bootstrap)

# Empirical two-tailed p-value
p_value <- mean(
  abs(bootstrap - bootstrap_mean) >=
    abs(empirical - bootstrap_mean)
)

# Print results
cat("Number of bootstrap replicates:", length(bootstrap), "\n")
cat("Bootstrap mean ratio:", bootstrap_mean, "\n")
cat("Empirical ratio:", empirical, "\n")
cat("P-value:", p_value, "\n")

# Create dataframe for plotting
df <- data.frame(ratio = bootstrap)

# Plot
p <- ggplot(df, aes(x = ratio)) +
  geom_histogram(
    bins = 40,
    fill = "skyblue",
    color = "black"
  ) +
  geom_vline(
    xintercept = empirical,
    color = "red",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = empirical,
    y = Inf,
    label = paste("Empirical =", signif(empirical, 4)),
    vjust = 2,
    color = "red"
  ) +
  labs(
    title = "Bootstrap Distribution of Coverage Ratios",
    subtitle = paste("P-value =", signif(p_value, 4)),
    x = "Candidate gene coverage / Total gene coverage",
    y = "Frequency"
  ) +
  theme_minimal()

p

ggsave(
  "Seixas_bootstrap_bp_ratio.svg",
  p,
  dpi = 300
)
ggsave(
  "Seixas_bootstrap_bp_ratio.png",
  p,
  dpi = 300
)
ggsave(
  "Seixas_bootstrap_bp_ratio.pdf",
  p
)

# Load libraries
library(ggplot2)

# Input file
file <- "No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI.txt"

# Read data (skip comment lines starting with ## or #)
df <- read.table(file,
                 header = FALSE,
                 comment.char = "#",
                 stringsAsFactors = FALSE)

# Keep only real data rows (remove header lines manually)
df <- df[!grepl("^##|^#", df$V1), ]

# Assign column names based on your format
colnames(df) <- c("Chrom", "Start", "End", "Genotype", "Polarity", "DI")

# Convert DI to numeric
df$DI <- as.numeric(df$DI)

# Basic histogram of DI values
p_hist <- ggplot(df, aes(x = DI)) +
  geom_histogram(bins = 50, color = "black", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Diagnostic Index (DI)",
       x = "Diagnostic Index (DI)",
       y = "Frequency")

p_den <- ggplot(df, aes(x = DI)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "DI distribution (density)",
       x = "Diagnostic Index (DI)",
       y = "Density")

ggsave("DI_distribution_hist.svg", plot = p_hist, dpi = 300)
ggsave("DI_distribution_den.svg", plot = p_den, dpi = 300)

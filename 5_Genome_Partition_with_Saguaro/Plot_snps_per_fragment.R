# Read the data
df <- read.table("Snp_counts_per_fragment_cactus4.bed", header = FALSE, sep = "\t")
colnames(df) <- c("chr", "start", "end", "snp_count")

# Summary statistics
snp_stats <- df %>%
  summarise(
    mean = mean(snp_count),
    median = median(snp_count),
    min = min(snp_count),
    max = max(snp_count)
  )

#---------------------------
# Full histogram
#---------------------------
p1 <- ggplot(df, aes(x = snp_count)) +
  geom_histogram(binwidth = 1, fill = "#69b3a2", color = "black", alpha = 0.8) +
  annotate(
    "label",
    x = max(df$snp_count) * 0.6,
    y = max(ggplot_build(ggplot(df, aes(x=snp_count)) + geom_histogram(binwidth=1))$data[[1]]$count) * 0.8,
    label = paste0(
      "Mean: ", round(snp_stats$mean, 2), "\n",
      "Median: ", snp_stats$median, "\n",
      "Min: ", snp_stats$min, "\n",
      "Max: ", snp_stats$max
    ),
    size = 5,
    fill = "white",
    color = "black",
    fontface = "bold"
  ) +
  labs(
    title = "Histogram of SNPs per Fragment (Full Range)",
    x = "Number of SNPs",
    y = "Number of Fragments"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

# Print plots
print(p1)

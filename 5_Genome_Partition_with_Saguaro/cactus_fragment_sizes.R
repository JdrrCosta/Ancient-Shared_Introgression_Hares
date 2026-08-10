library(tidyverse)

# Read BED file
cactus12 <- read_tsv("cactus_iter12_cactus4_windows.bed", col_names = FALSE)

# Calculate fragment sizes
cactus12 <- cactus12 %>%
  mutate(size = X3 - X2)

# Summary statistics
stats <- cactus12 %>%
  summarise(
    mean = mean(size),
    median = median(size),
    min = min(size),
    max = max(size)
  )

# ---- Top 0.01% largest fragments ----
threshold <- quantile(cactus12$size, probs = 0.9999)

top_fragments <- cactus12 %>%
  filter(size >= threshold) %>%
  arrange(desc(size))

cat("99.99th percentile threshold:", threshold, "bp\n")
cat("Number of fragments in the top 0.01%:", nrow(top_fragments), "\n")

# Optional: save the top fragments
write_tsv(top_fragments, "Cactus4_top_0.01_percent_fragments.tsv")

# Histogram with nicer aesthetics
ggplot(cactus12, aes(x = size)) +
  geom_histogram(binwidth = 50, fill = "#69b3a2", color = "black", alpha = 0.8) +
  
  # Add summary stats in a white box
  annotate(
    "label",
    x = max(cactus12$size) * 0.6,
    y = max(ggplot_build(
      ggplot(cactus12, aes(x = size)) +
        geom_histogram(binwidth = 50)
    )$data[[1]]$count) * 0.8,
    label = paste0(
      "Mean: ", round(stats$mean, 2), "\n",
      "Median: ", round(stats$median, 2), "\n",
      "Min: ", stats$min, "\n",
      "Max: ", stats$max, "\n",
      "Top 0.01% cutoff: ", round(threshold, 1), " bp\n",
      "Top 0.01% n = ", nrow(top_fragments)
    ),
    size = 5,
    fill = "white",
    color = "black",
    fontface = "bold"
  ) +
  
  labs(
    title = "Distribution of Fragment Sizes",
    x = "Fragment Size (bp)",
    y = "Frequency"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

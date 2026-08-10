library(tidyverse)

# Find all files
files <- list.files(
  pattern = "^HI_masked_(cacti|leur)_no_auto_(DI_.*|no_threshold)\\.txt$",
  full.names = TRUE
)

# Extract threshold
get_threshold <- function(x) {
  if (str_detect(x, "no_threshold")) {
    return("No threshold")
  }
  str_extract(basename(x), "(?<=DI_)-?\\d+")
}

# Dataset label (used for linetype)
get_dataset <- function(x) {
  case_when(
    str_detect(basename(x), "masked_cacti") ~ "Masked",
    str_detect(basename(x), "masked_leur")  ~ "Active"
  )
}

# Read + summarise
get_species_means <- function(file) {
  dat <- read.delim(file)
  
  dat %>%
    mutate(
      Mean_HI = Mean_HI / 2,
      Species = str_to_lower(str_extract(Individual, "^[A-Za-z]+"))
    ) %>%
    group_by(Species) %>%
    summarise(mean_HI = mean(Mean_HI, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      threshold = get_threshold(file),
      dataset   = get_dataset(file)
    )
}

results <- map_dfr(files, get_species_means)

# Order thresholds
threshold_order <- c(
  "No threshold",
  sort(
    unique(results$threshold[results$threshold != "No threshold"]) %>%
      as.numeric(),
    decreasing = TRUE
  )
)

results$threshold <- factor(results$threshold, levels = as.character(threshold_order))

# Species colors
species_colors <- c(
  lgra = "#c5cca0ff",
  ltim = "#a3bedcff",
  lcas = "#f4a6b8ff",
  lcor = "#fca64dff"
)

# Common plot function
make_plot <- function(df, y_limits, title) {
  ggplot(
    df,
    aes(
      x = threshold,
      y = mean_HI,
      colour = Species,
      group = interaction(Species, dataset),
      linetype = dataset
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~Species, scales = "fixed") +
    scale_colour_manual(values = species_colors) +
    scale_linetype_manual(values = c("Masked" = "dashed", "Active" = "solid")) +
    coord_cartesian(ylim = y_limits) +
    labs(
      x = "Threshold",
      y = "Mean HI",
      colour = "Species",
      linetype = "Dataset",
      title = title
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# Split datasets
plot_low  <- make_plot(
  filter(results, Species %in% c("lgra", "ltim")),
  c(0, 0.5),
  "Mean HI (0–0.5 range): lgra & ltim"
)

plot_high <- make_plot(
  filter(results, Species %in% c("lcas", "lcor")),
  c(0.5, 1),
  "Mean HI (0.5–1 range): lcas & lcor"
)

# Save
ggsave("HI_low_range.svg", plot_low, dpi = 300, width = 10, height = 5)
ggsave("HI_high_range.svg", plot_high, dpi = 300, width = 10, height = 5)

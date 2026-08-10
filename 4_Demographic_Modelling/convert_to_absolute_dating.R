# --- Load libraries ---
library(dplyr)
library(readr)

# --- Parameters ---
mu <- 2.8e-9  # mutation rate per site per generation
g  <- 2       # generation time in years

# --- Read BPP table ---
bpp_data <- read.table("1000_rep1_MSC-I_proto_cascor_flipped.txt", header = TRUE)

# Keep only theta, tau, and phi rows
bpp_filtered <- bpp_data %>%
  filter(grepl("^theta|^tau|^phi", param))

# --- Conversion functions ---

# Convert theta to Ne
convert_theta <- function(x) x / (4 * mu)

# Convert tau to generations and years
convert_tau <- function(x) {
  T_gen <- x / mu
  T_years <- T_gen * g
  return(list(T_gen = T_gen, T_years = T_years))
}

# --- Absolute values table ---
results <- bpp_filtered %>%
  rowwise() %>%
  mutate(
    abs_mean = ifelse(grepl("^theta", param), convert_theta(mean),
                      ifelse(grepl("^tau", param), convert_tau(mean)$T_years, mean)),
    abs_median = ifelse(grepl("^theta", param), convert_theta(median),
                        ifelse(grepl("^tau", param), convert_tau(median)$T_years, median)),
    abs_2.5 = ifelse(grepl("^theta", param), convert_theta(X2.5.HPD),
                     ifelse(grepl("^tau", param), convert_tau(X2.5.HPD)$T_years, X2.5.HPD)),
    abs_97.5 = ifelse(grepl("^theta", param), convert_theta(X97.5.HPD),
                      ifelse(grepl("^tau", param), convert_tau(X97.5.HPD)$T_years, X97.5.HPD))
  ) %>%
  ungroup() %>%
  select(param, mean, median, X2.5.HPD, X97.5.HPD, abs_mean, abs_median, abs_2.5, abs_97.5)

# --- Pretty transformed table ---

# Helper to convert theta Ne into "K" notation
format_k <- function(x) {
  ifelse(x >= 1000, paste0(round(x / 1000, 1), "k"), round(x, 1))
}

# Helper to convert years into mya/kya
format_time <- function(x) {
  ifelse(x >= 1e6,
         paste0(round(x / 1e6, 2), " mya"),
         paste0(round(x / 1e3, 0), " kya"))
}

# Helper to convert phi into %
format_phi <- function(x) paste0(round(x * 100, 2), "%")

pretty_table <- results %>%
  mutate(
    transformed_mean = as.character(case_when(
      grepl("^theta", param) ~ format_k(abs_mean),
      grepl("^tau", param) ~ format_time(abs_mean),
      grepl("^phi", param) ~ format_phi(mean),
      TRUE ~ as.character(abs_mean)
    )),
    transformed_median = as.character(case_when(
      grepl("^theta", param) ~ format_k(abs_median),
      grepl("^tau", param) ~ format_time(abs_median),
      grepl("^phi", param) ~ format_phi(median),
      TRUE ~ as.character(abs_median)
    )),
    transformed_2.5 = as.character(case_when(
      grepl("^theta", param) ~ format_k(abs_2.5),
      grepl("^tau", param) ~ format_time(abs_2.5),
      grepl("^phi", param) ~ format_phi(X2.5.HPD),
      TRUE ~ as.character(abs_2.5)
    )),
    transformed_97.5 = as.character(case_when(
      grepl("^theta", param) ~ format_k(abs_97.5),
      grepl("^tau", param) ~ format_time(abs_97.5),
      grepl("^phi", param) ~ format_phi(X97.5.HPD),
      TRUE ~ as.character(abs_97.5)
    ))
  ) %>%
  select(param, transformed_mean, transformed_median, transformed_2.5, transformed_97.5)

# --- Save both tables ---
write_csv(results, "MSC-I_proto_cascor_flipped_absolute_values.csv")
write_csv(pretty_table, "MSC-I_proto_cascor_flipped_Pretty_table.csv")

# --- Print results ---
cat("\nAbsolute values table:\n")
print(results)

cat("\nTransformed (pretty) table:\n")
print(pretty_table)

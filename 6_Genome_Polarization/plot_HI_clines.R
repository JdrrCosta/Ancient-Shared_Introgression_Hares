library(data.table)
library(ggplot2)

# Read files
pos <- fread("Cactus4_Positive_positions_cline_HI_10bp_no_auto_DI_20.txt",
             col.names = c("chr", "start", "end", "HI", "distance"))

neg <- fread("Cactus4_Negative_positions_cline_HI_10bp_no_auto_DI_20.txt",
             col.names = c("chr", "start", "end", "HI", "distance"))

# Merge datasets
dat <- rbind(pos, neg)

# Transform HI values:
# 1) divide by 2
# 2) flip values so low becomes high and vice versa
# Example: 0.2 -> 0.8
dat[, HI := 1 - (HI / 2)]

# Compute summary statistics per distance
dat_summary <- dat[, .(
  mean_HI = mean(HI, na.rm = TRUE),
  sd_HI   = sd(HI, na.rm = TRUE),
  n       = .N
), by = distance]

# Standard error and 95% confidence interval
dat_summary[, se := sd_HI / sqrt(n)]

dat_summary[, `:=`(
  lower = mean_HI - qt(0.975, df = n - 1) * se,
  upper = mean_HI + qt(0.975, df = n - 1) * se
)]

# Plot
p <- ggplot(dat_summary, aes(x = distance, y = mean_HI)) +
  
  # 95% CI of the observed means
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.25) +
  
  # Observed means
  geom_line(color = "black", linewidth = 1) +
  geom_point(color = "black", size = 2) +
  
  # LOESS trend with its own 95% CI
 # geom_smooth(method = "loess",
 #             span = 0.3,
#              se = TRUE,
#              color = "red",
#              fill = "pink",
#              linewidth = 1) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1)) +
  
  theme_minimal() +
  
  labs(
    x = "Distance from introgressed region",
    y = "Mean HI",
    title = "Cactus 4 HI cline around introgressed region (10 bp)"
  )

p

ggsave("Cactus4_HI_Cline_10bp_confidence_interval_no_smooth.svg", plot = p, dpi = 300, width = 10, height = 10)

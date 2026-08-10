library(ggplot2)
library(dplyr)
library(patchwork)

# ---------- FUNCTION ----------
plot_dxy <- function(boot_lcas, boot_lcor, emp_lcas, emp_lcor, title, xlim_range) {
  
  # Bootstrap data
  bootstrap_df <- bind_rows(
    data.frame(value = boot_lcas, comparison = "Lcas vs Ltim"),
    data.frame(value = boot_lcor, comparison = "Lcor vs Ltim")
  )
  
  # Empirical data
  empirical_df <- data.frame(
    comparison = c("Lcas vs Ltim", "Lcor vs Ltim"),
    empirical = c(emp_lcas, emp_lcor)
  )
  
  ggplot(bootstrap_df, aes(x = value, fill = comparison, colour = comparison)) +
    geom_density(alpha = 0.4, linewidth = 1, trim = TRUE) +
    
    geom_vline(
      data = empirical_df,
      aes(xintercept = empirical, colour = comparison),
      linetype = "dashed",
      linewidth = 1
    ) +
    
    coord_cartesian(xlim = xlim_range) +  # ✅ SAME AXIS
    
    labs(
      title = title,
      x = "Mean Dxy",
      y = "Density"
    ) +
    
    scale_fill_manual(values = c(
      "Lcas vs Ltim" = "#f4a6b8ff",  # blue
      "Lcor vs Ltim" = "#ff7f0e"   # orange
    )) +
    
    scale_colour_manual(values = c(
      "Lcas vs Ltim" = "#f4a6b8ff",
      "Lcor vs Ltim" = "#ff7f0e"
    )) +
    
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# ---------- LOAD DATA ----------

# Cactus4
boot_lcas_4 <- scan("Circular_bootstrap_cactus4_Lcas_Ltim.dxy")
emp_lcas_4  <- scan("Saguaro_iter12_cactus4_Lcas_Ltim.dxy")

boot_lcor_4 <- scan("Circular_bootstrap_cactus4_Lcor_Ltim.dxy")
emp_lcor_4  <- scan("Saguaro_iter12_cactus4_Lcor_Ltim.dxy")

# Cactus12
boot_lcas_12 <- scan("Circular_bootstrap_Lcas_Ltim.dxy")
emp_lcas_12  <- scan("Saguaro_iter12_cactus12_Lcas_Ltim.dxy")

boot_lcor_12 <- scan("Circular_bootstrap_Lcor_Ltim.dxy")
emp_lcor_12  <- scan("Saguaro_iter12_cactus12_Lcor_Ltim.dxy")

# ---------- FIXED X-AXIS RANGE ----------
# Compute once across ALL datasets
all_values <- c(boot_lcas_4, boot_lcor_4, boot_lcas_12, boot_lcor_12,
                emp_lcas_4, emp_lcor_4, emp_lcas_12, emp_lcor_12)

xlim_range <- range(all_values)

# ---------- PLOTS ----------
p_cactus4 <- plot_dxy(
  boot_lcas_4, boot_lcor_4,
  emp_lcas_4, emp_lcor_4,
  "Cactus4",
  xlim_range
)

p_cactus12 <- plot_dxy(
  boot_lcas_12, boot_lcor_12,
  emp_lcas_12, emp_lcor_12,
  "Cactus12",
  xlim_range
)

# Combine (side-by-side)
final_plot <- p_cactus4 / p_cactus12

# ---------- SAVE ----------
ggsave("Dxy_distribution_cactus4_vs_cactus12.svg",
       plot = final_plot,
       width = 12,
       height = 5)

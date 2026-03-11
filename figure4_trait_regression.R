##' Produce preliminary figures for thesis proposal.
##'
##' @author [Nathan Malamud]
##' @date [2025-01-]

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(factoextra)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Import data
traits <- read_csv("./data/traits.csv")

# Ensure species is a factor
traits$species <- factor(traits$species, 
                         levels = c("R. sativus",
                                    "B. officinalis",
                                    "H. vulgare"))

# Convert to proper units
# In dataset, LMA is in g/m^2, LDMC in mg/g, CHL in µg/cm^2, and GRT in g/day
# Want:
# LMA in g/m^2 -> kg/m^2
# LDMC in mg/g -> the same
# CHL in µg/cm^2 -> mg/m^2
# Filter only by traits of interest
traits <- traits %>% select(
  species, LDMC, LMA, CHL, dry_whole_g, treatment_mmol
)

# Define growth period
growth_period_days <- 6 * 7  # 6-week experiment (42 days)

# Convert to proper units
# In dataset, LMA is in g/m^2, LDMC in mg/g, CHL in µg/cm^2, and GRT in g/day
# Want:
# LMA in g/m^2 -> kg/m^2
# LDMC in mg/g -> the same
# CHL in µg/cm^2 -> mg/m^2
# Filter only by traits of interest
traits <- traits %>% select(
  species, LDMC, LMA, CHL, dry_whole_g, treatment_mmol
)

# Define growth period
growth_period_days <- 6 * 7  # 6-week experiment (42 days)

# Convert units and compute growth rate
traits <- traits %>%
  select(species, LDMC, LMA, CHL, treatment_mmol, dry_whole_g) %>%
  mutate(
    LMA = LMA / 1000,  # Convert g/m² to kg/m²
    CHL = CHL * 10,     # Convert µg/cm² to mg/m²
    GRT = dry_whole_g / growth_period_days  # Compute growth rate per day
  ) %>%
  select(-dry_whole_g)  # Remove intermediate variable if not needed

# Styles ----
# Display first rows of the dataset
# Define consistent font size
base_font_size <- 12

# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = base_font_size),
    axis.text = element_text(size = 9),  
    legend.position = "bottom",
    
    # Title - left justified, 12 point font, bold
    plot.title = element_text(size = base_font_size, hjust = 0.5),
    
    # Aspect ratio 1:1
    aspect.ratio = 1,
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define units for variables with LaTeX formatting
label_units <- list(
  "LDMC" = expression("LDMC"~"("*mg~g^-1*")"),  # No change
  "LMA" = expression("LMA"~"("*kg~m^-2*")"),  # Converted from g/m² (measured) to kg/m²
  "CHL" = expression("Chlorophyll content"~"("*mg~m^-2*")"),  # Converted from µg/cm² (estimated) to mg/m²
  "treatment_mmol" = "Nitrogen (mM)"  # No change
)

# Log-Log Trait Regressions ----
# 2) Show pairwise responses of LDMC, LMA, CHL, and area_cm2 using glms (log-log sma)

## LMA vs LDMC----
# Compare trait-trait scaling relationships using SMA (standard major axis)
# regression models.
# Milos: OLS vs SMA, sum of squared residuals is changing. R2 is the same, trendlines different.
#   Syntax: * does slope, + does slope and elevation
lma_ldmc_sma <- smatr::sma(LMA ~ LDMC * species, log = "XY",
                           method = "SMA",
                           data = traits)

lma_ldmc_plot <- ggplot(traits, aes(y = LDMC, x = LMA)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(aes(color = species), method = "SMA", se=F) +  # SMA regression line
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["LMA"]], y = label_units[["LDMC"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme


## LMA vs CHL----
lma_chl_sma <- smatr::sma(LMA ~ CHL * species, log = "XY",
                          method = "SMA",
                          data = traits)

lma_chl_plot <- ggplot(traits, aes(y = LMA, x = CHL)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(data = subset(traits, species == "R. sativus"), 
               method = "SMA", se = F, color = josef_colors["R. sativus"], linetype = "solid") +
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["CHL"]], y = label_units[["LMA"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

## LDMC vs CHL----
ldmc_chl_sma <- smatr::sma(LDMC ~ CHL * species, log = "XY",
                           method = "SMA",
                           data = traits)

ldmc_chl_plot <- ggplot(traits, aes(y = LDMC, x = CHL)) +
  geom_point(size = 2, alpha = 0.5, aes(color = species, shape = species)) +  # Adjusted point size and opacity
  stat_ma_line(data = subset(traits, species == "R. sativus"), 
               method = "SMA", se = F, color = josef_colors["R. sativus"], linetype = "dashed") +  # SMA regression line for R. sativus
  scale_color_manual(values = josef_colors, name = "Species") +
  scale_shape_manual(values = c(16, 17, 18), name = "Species") +
  scale_y_log10() +  # Log scale for Y-axis
  scale_x_log10() +  # Log scale for X-axis
  labs(x = label_units[["CHL"]], y = label_units[["LDMC"]]) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 1,
    legend.position = "bottom"  # Legend position at the bottom
  ) +
  guides(
    color = guide_legend("Species", override.aes = list(shape = c(16, 17, 18), alpha=1.0, size = 4)),  # Bigger shapes
    shape = guide_legend("Species", override.aes = list(size = 4))
  ) + custom_theme

# # Apply ggarrange and save as figure 3.
# # Combine all plots
figure3_regression_plot <- ggarrange(
  lma_ldmc_plot, lma_chl_plot, ldmc_chl_plot,
  #   leaf_area_lma_plot, leaf_area_ldmc_plot, leaf_area_chl_plot,
  ncol = 3, nrow = 1,  # 2 rows, 3 columns
  labels = c("a", "b", "c"), #"d", "e", "f"),  # Subplot labels
  common.legend = TRUE,  # Combine legends
  legend = "bottom"  # Legend at the bottom
)

# Save as png
ggsave("figure3_trait_correlations.png", figure3_regression_plot, 
       width = 8.5, height = 3.75, dpi = 300, bg="white")

# Save as svg
library(svglite)
ggsave("figure3_trait_correlations.svg", figure3_regression_plot,
       width = 8.5, height = 3.75, bg="white")

# # Print to console
# print(figure3_regression_plot)
# 
# # Save the figure
ggsave(
  filename = "./figures/prelim/figure3_regression_plots.png",
  plot = figure3_regression_plot,
  width = 10, height = 5  # Adjust width and height for layout
)

# Fit models with species factor and compare regression model stats
summary(lma_ldmc_sma)
summary(lma_chl_sma)
summary(ldmc_chl_sma)

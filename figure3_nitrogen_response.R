##' Plot reciprocal transformeed
##' responses with random intercepts
##' 
##' @author [Nathan D. Malamud]
##' @date [2025-01-29]
##' 

# Libraries ----
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(scales)
library(smatr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(cowplot)

# Model fitting ----
library(lme4)
library(MuMIn)

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
  "treatment_mmol" = "Nitrogen addition (mM)"  # No change
)

# Import Data ----
# REMINDER: Set Working Directory -> Source File Location
# Define traits of interest
traits_of_interest <- c("LMA", "LDMC", "CHL")

# Apply specific transformations for each trait
traits <- traits %>%
   mutate(
     LMA_trans = I(LMA),
     LDMC_trans = I(LDMC),
     CHL_trans = I(CHL)
)

# Load required libraries
library(lme4)
library(ggplot2)
library(ggpubr)

# Fit mixed-effects models with appropriate structure
mod_LMA <- lm(LMA ~ treatment_mmol, data = traits)
mod_LDMC <- lmer(LDMC ~ treatment_mmol + (1 | species), data = traits, REML = FALSE)
mod_CHL <- lm(CHL ~ treatment_mmol, data = traits)

# Get predictions with standard errors
pred_LMA <- predict(mod_LMA, se.fit = TRUE)
pred_LDMC <- predict(mod_LDMC, se.fit = TRUE)
pred_CHL <- predict(mod_CHL, se.fit = TRUE)

# Compute 95% Confidence Intervals
traits$predicted_LMA <- pred_LMA$fit
traits$LMA_lower_CI <- pred_LMA$fit - 1.96 * pred_LMA$se.fit
traits$LMA_upper_CI <- pred_LMA$fit + 1.96 * pred_LMA$se.fit

traits$predicted_LDMC <- pred_LDMC$fit
traits$LDMC_lower_CI <- pred_LDMC$fit - 1.96 * pred_LDMC$se.fit
traits$LDMC_upper_CI <- pred_LDMC$fit + 1.96 * pred_LDMC$se.fit

traits$predicted_CHL <- pred_CHL$fit
traits$CHL_lower_CI <- pred_CHL$fit - 1.96 * pred_CHL$se.fit
traits$CHL_upper_CI <- pred_CHL$fit + 1.96 * pred_CHL$se.fit

# Compute 95% Prediction Intervals (Optional)
sigma_LMA <- sigma(mod_LMA)  # Residual SD
traits$LMA_lower_PI <- pred_LMA$fit - 1.96 * sqrt(pred_LMA$se.fit^2 + sigma_LMA^2)
traits$LMA_upper_PI <- pred_LMA$fit + 1.96 * sqrt(pred_LMA$se.fit^2 + sigma_LMA^2)

sigma_LDMC <- sigma(mod_LDMC)
traits$LDMC_lower_PI <- pred_LDMC$fit - 1.96 * sqrt(pred_LDMC$se.fit^2 + sigma_LDMC^2)
traits$LDMC_upper_PI <- pred_LDMC$fit + 1.96 * sqrt(pred_LDMC$se.fit^2 + sigma_LDMC^2)

sigma_CHL <- sigma(mod_CHL)
traits$CHL_lower_PI <- pred_CHL$fit - 1.96 * sqrt(pred_CHL$se.fit^2 + sigma_CHL^2)
traits$CHL_upper_PI <- pred_CHL$fit + 1.96 * sqrt(pred_CHL$se.fit^2 + sigma_CHL^2)

# Define individual plots

# LMA: Single trend line with CI and PI
p_LMA <- ggplot(traits, aes(x = treatment_mmol, y = LMA)) +
  geom_point(aes(color = species, alpha = 0.5, shape = species), size = 2) +
  #geom_ribbon(aes(ymin = LMA_lower_CI, ymax = LMA_upper_CI), fill = "#999999", alpha = 0.2) +  # 95% CI
  #geom_ribbon(aes(ymin = LMA_lower_PI, ymax = LMA_upper_PI), fill = "#CCCCCC", alpha = 0.2) +  # 95% PI
  geom_line(aes(y = predicted_LMA), linewidth = 1, color = "#666666") +  # Trend line
  scale_y_log10() +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["LMA"]], color = "Species", shape = "Species") +
  scale_color_manual(values = josef_colors) +
  custom_theme +
  guides(alpha = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4)), 
         shape = guide_legend(override.aes = list(alpha = 1, size = 4)))

# LDMC: Separate intercepts per species with CI and PI
p_LDMC <- ggplot(traits, aes(x = treatment_mmol, y = LDMC, color = species)) +
  geom_point(aes(alpha = 0.5, shape = species), size = 2) +
  #geom_ribbon(aes(ymin = LDMC_lower_CI, ymax = LDMC_upper_CI, fill = "#999999"), alpha = 0.2) +  # 95% CI
  #geom_ribbon(aes(ymin = LDMC_lower_PI, ymax = LDMC_upper_PI, fill = species), alpha = 0.2) +  # 95% PI
  geom_line(aes(y = predicted_LDMC), linewidth = 1) +  # Trend lines per species
  scale_y_log10() +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["LDMC"]], color = "Species", shape = "Species") +
  scale_color_manual(values = josef_colors) +
  scale_fill_manual(values = josef_colors) +  # Fill matches species colors
  custom_theme +
  guides(alpha = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4)), 
         shape = guide_legend(override.aes = list(alpha = 1, size = 4)))

# CHL: Single trend line with CI and PI
p_CHL <- ggplot(traits, aes(x = treatment_mmol, y = CHL)) +
  geom_point(aes(color = species, alpha = 0.5, shape = species), size = 2) +
  #geom_ribbon(aes(ymin = CHL_lower_CI, ymax = CHL_upper_CI), fill = "#999999", alpha = 0.2) +  # 95% CI
  #geom_ribbon(aes(ymin = CHL_lower_PI, ymax = CHL_upper_PI), fill = "#CCCCCC", alpha = 0.2) +  # 95% PI
  geom_line(aes(y = predicted_CHL), linewidth = 1, color = "#666666") +  # Trend line
  scale_y_log10() +
  labs(x = label_units[["treatment_mmol"]], y = label_units[["CHL"]], color = "Species", shape = "Species") +
  scale_color_manual(values = josef_colors) +
  custom_theme +
  guides(alpha = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4)), 
         shape = guide_legend(override.aes = list(alpha = 1, size = 4)))

# Arrange all plots in a grid with a common legend
final_plot <- ggarrange(p_LMA, p_LDMC, p_CHL,
                        labels = c("a", "b", "c"),
                        nrow = 1, ncol = 3,
                        common.legend = TRUE, legend = "bottom")

# Display the final plot
print(final_plot)

# Save as png
ggsave("figure2_nitrogen_response.png", final_plot,
       width = 8.5, height = 3.75, dpi = 300, bg="white")

# Save as svg
library(svglite)
ggsave("figure2_nitrogen_response.svg", final_plot,
       width = 8.5, height = 3.75, bg="white")

##  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Create an empty list to store p-values
p_values <- list()
r2_values <- list()

# Iterate over traits and fit models accordingly
for (trait in traits_of_interest) {
  
  # Determine the transformed trait name and formula
  transformed_trait <- paste0(trait, "_trans")
  formula <- as.formula(paste(transformed_trait,
                              "~ treatment_mmol + (1 | species)"))
  
  # Fit the mixed-effects model using REML = FALSE for hypothesis testing
  mod <- lmer(formula, data = traits, REML = FALSE)
  
  # Extract p-value for treatment effect
  p_value <- anova(mod)$`Pr(>F)`[1]
  # Extract R-squared value
  r2_value <- r.squaredGLMM(mod)[1]
  
  # Store in the list with trait name
  p_values[[trait]] <- p_value
  r2_values[[trait]] <- r2_value
}

# Print p-values to R console
print(p_values)

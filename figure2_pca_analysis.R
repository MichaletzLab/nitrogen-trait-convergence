##' Run PCA analysis of traits.
##' View differences across crops and nitrogen treatment.
##' 
##' @author [Nathan D. Malamud]
##' @date [2025-01-27] 
##'
##' Revised March 2026 to add convergence analysis
##' Sean Michaletz (sean.michaletz@ubc.ca)
##'

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
# Define factor levels as species
traits <-  read_csv("./data/traits.csv")
traits$species <- factor(traits$species,
                         levels=c("R. sativus", "B. officinalis", "H. vulgare"))

# Calculate rate of growth
growth_period_days <- 6 * 7 # 6 week experiment
traits$GRT <- (traits$dry_whole_g / growth_period_days)

# Define traits of interest here
traits_of_interest <- c("LDMC", "LMA", "CHL", "GRT")

# Filter by metrics of interest only
traits <- traits %>%
  select(barcodeID, species, treatment_mmol,
         all_of(traits_of_interest))

# Styling and aesthetics ----
# Define a custom theme for all plots
custom_theme <- theme_classic() +  # Base theme
  theme(
    # Text and font styling
    text = element_text(family = "sans", size = 12),
    axis.text = element_text(size = 10),  
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    
    # Axis labels
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # Margin for y-axis title
    
    # Panel and grid styling
    panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1),  
    panel.grid.minor = element_blank(),  # No minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
    
    # Aspect ratio
    aspect.ratio = 1  # 1:1 ratio
  )

# Custom colors advised by J. Garen
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Define units for variables
# TODO: format with Latex expressions
label_units <- c(
  "CHL" = "CHL (ug / cm²)",
  "LDMC" = "LDMC (mg / g)",
  "LMA" = "LMA (g / m²)",
  "area_cm2" = "Leaf Area (cm²)",
  "treatment_mmol" = "N (mM)"
)

# PCA Analysis ----
# Keep treatment_mmol numeric here. Convert to factor only inside plotting calls.
# This keeps the PCA trait inputs unchanged while allowing numeric trend tests later.

pca_data <- traits %>%
  select(barcodeID, species, treatment_mmol, LMA, LDMC, CHL, GRT)

pca_input <- pca_data %>%
  select(LMA, LDMC, CHL, GRT)

# Save the exact scaled trait matrix that enters the PCA so the 4-D sensitivity
# analysis uses the same rows and standardization as the ordination.
pca_input_scaled <- scale(pca_input, center = TRUE, scale = TRUE)

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)

# Calculate % variance explained by each PC
explained_variation <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

## Ordination biplot, groups by species ----
pca_plot_species <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = pca_data$species,  # Color individuals by species
  palette = josef_colors,  # Custom species colors
  addEllipses = FALSE,  # Add confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  legend.title = "Species",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

# Print the PCA plot
print(pca_plot_species)

## Ordination groups by treatment ----
# Define custom colors for specific nitrogen levels
custom_colors1 <- c(
  "0"  = "#D2B48C",  # Darker beige
  "5"  = "#E6A96B",  # Light brownish-orange
  "10" = "#FDBE85",  # Soft orange
  "15" = "#FDAE61",  # Bright orange
  "20" = "#F46D43",  # Deep orange
  "25" = "#E31A1C",  # Red
  "30" = "#BD0026",  # Dark red
  "35" = "#800026"   # Deep maroon
)

# Highlight key treatments using custom colors for specific nitrogen levels
# custom_colors2 <- c(
#   "0"  = "#D2B48C",  # Darker beige
#   "5"  = "grey",  
#   "10" = "grey", 
#   "15" = "#FDAE61",  # Bright orange
#   "20" = "grey",
#   "25" = "grey",
#   "30" = "#BD0026",  # Dark red
#   "35" = "grey" 
# )


# Define distinct shapes (must match the number of nitrogen levels)
custom_shapes <- c(16, 17, 15, 3, 8, 10, 7, 4)  # 8 values

pca_plot_treatment_ellipse <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = factor(pca_data$treatment_mmol),  # Ensure discrete coloring
  shape.ind = factor(pca_data$treatment_mmol),  # Ensure discrete shapes
  addEllipses = TRUE,  # Add confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  ellipse.alpha = 0.0,  # Set fill transparency to 0 for ellipses
  legend.title = "Nitrogen addition (mM)",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  scale_color_manual(values = custom_colors1) +  # Assign custom colors
  scale_shape_manual(values = custom_shapes) +  # Assign distinct shapes
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

pca_plot_treatment_no_ellipse <- fviz_pca_biplot(
  pca,
  geom.ind = "point",  # Plot individuals as points
  col.ind = factor(pca_data$treatment_mmol),  # Ensure discrete coloring
  shape.ind = factor(pca_data$treatment_mmol),  # Ensure discrete shapes
  addEllipses = FALSE,  # No confidence ellipses for each group
  ellipse.level = 0.95,  # Set confidence level for ellipses
  ellipse.alpha = 0.0,  # Set fill transparency to 0 for ellipses
  legend.title = "Nitrogen addition (mM)",  # Legend title
  repel = TRUE  # Avoid overlap of labels
) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variation[1], "% variance)"),
    y = paste0("PC2 (", explained_variation[2], "% variance)")
  ) +
  custom_theme +
  scale_color_manual(values = custom_colors1) +  # Assign custom colors
  scale_shape_manual(values = custom_shapes) +  # Assign distinct shapes
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

# Full plot with ellipses
pca_plot_treatment <- ggarrange(
  pca_plot_treatment_ellipse,
  pca_plot_treatment_no_ellipse,
  nrow = 1, labels = c("a", "b")
)

# Print the PCA plot
print(pca_plot_treatment)

# Create a scree plot for eigenvalues
pca_eigenplot <- fviz_eig(pca, addlabels = TRUE, geom = c("line", "point")) +
  custom_theme +
  geom_segment(
    aes(x = seq_along(pca$sdev), xend = seq_along(pca$sdev), 
        y = 0, yend = pca$sdev^2 / sum(pca$sdev^2) * 100),
    linetype = "dashed", color = "gray"
  ) +
  geom_point(
    aes(x = seq_along(pca$sdev), 
        y = pca$sdev^2 / sum(pca$sdev^2) * 100), color = "black"
  ) +
  theme(
    aspect.ratio = 1
  ) +
  labs(
    y = "% Variance",
    title = NULL
  )

# Full plot with ellipses
pca_plot_treatment <- ggarrange(
  pca_plot_treatment_ellipse,
  pca_plot_treatment_no_ellipse,
  nrow = 1, labels = c("", "")
)

print(pca_plot_treatment)

print(pca_eigenplot)

# Save the PCA plot as Figure 4
# TODO: investigated ignored "override.aes" warnings
ggsave(
  filename = "./figures/treatment_pca_plot.png",
  plot = pca_plot_treatment,
  width = 10, height = 5,
  bg = "white"
)

# Save as svg
library(svglite)
ggsave("./figures/treatment_pca_plot.svg", pca_plot_treatment,
       width = 16, height = 6, bg="white")

# Save scree plot
ggsave(
  filename = "./figures/pca_eigenplot.png",
  plot = pca_eigenplot,
  width = 5, height = 5,
  bg = "white"
)


# Convergence metrics for Fig. 4 ----
# Analytical rationale:
# 1. All convergence metrics are computed in PC1-PC2 space to match Fig. 4.
# 2. All metrics are computed at the nitrogen-treatment level, pooled across species.
# 3. All computations use individual PC scores from the PCA, not ellipse parameters.

# Step 2. Extract individual PC scores with metadata.
# This is built from the exact rows retained in the PCA so there is a one-to-one
# correspondence between the plotted points and the convergence metrics.
pc_scores <- pca_data %>%
  select(barcodeID, species, treatment_mmol) %>%
  bind_cols(as_tibble(pca$x[, c("PC1", "PC2"), drop = FALSE])) %>%
  mutate(treatment_mmol = as.numeric(treatment_mmol))

# Step 3-4. Compute treatment sample sizes, centroids, and mean distance to centroid.
# mean_dist_to_centroid is analogous to functional dispersion (FDis), but is
# computed in ordination space rather than the original trait space.
treatment_summary <- pc_scores %>%
  group_by(treatment_mmol) %>%
  mutate(
    centroid_pc1 = mean(PC1),
    centroid_pc2 = mean(PC2),
    dist_to_centroid = sqrt(
      (PC1 - centroid_pc1)^2 +
        (PC2 - centroid_pc2)^2
    )
  ) %>%
  summarise(
    sample_size = n(),
    centroid_pc1 = first(centroid_pc1),
    centroid_pc2 = first(centroid_pc2),
    mean_dist_to_centroid = mean(dist_to_centroid),
    .groups = "drop"
  ) %>%
  arrange(treatment_mmol)

# Helper functions ----
get_reference_centroid <- function(score_df, ref_levels) {
  score_df %>%
    filter(treatment_mmol %in% ref_levels) %>%
    summarise(
      ref_pc1 = mean(PC1),
      ref_pc2 = mean(PC2)
    )
}

add_centroid_distance <- function(summary_df, score_df, ref_levels) {
  ref_centroid <- get_reference_centroid(score_df, ref_levels)
  
  summary_df %>%
    mutate(
      centroid_distance_to_convergence_region = if_else(
        treatment_mmol %in% ref_levels,
        NA_real_,
        sqrt(
          (centroid_pc1 - ref_centroid$ref_pc1)^2 +
            (centroid_pc2 - ref_centroid$ref_pc2)^2
        )
      )
    )
}

run_trend_tests <- function(df, response_var) {
  x <- df$treatment_mmol
  y <- df[[response_var]]
  
  spearman_fit <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )
  
  lm_formula <- as.formula(paste(response_var, "~ treatment_mmol"))
  lm_fit <- lm(lm_formula, data = df)
  lm_coef <- summary(lm_fit)$coefficients["treatment_mmol", ]
  
  tibble(
    response = response_var,
    spearman_rho = unname(spearman_fit$estimate),
    spearman_p = spearman_fit$p.value,
    lm_slope = unname(lm_coef[["Estimate"]]),
    lm_r2 = unname(summary(lm_fit)$r.squared),
    lm_p = unname(lm_coef[["Pr(>|t|)"]])
  )
}

run_reference_sensitivity <- function(summary_df, score_df, ref_levels, label) {
  add_centroid_distance(summary_df, score_df, ref_levels) %>%
    filter(!is.na(centroid_distance_to_convergence_region)) %>%
    run_trend_tests("centroid_distance_to_convergence_region") %>%
    transmute(
      analysis = label,
      spearman_rho,
      spearman_p
    )
}

# Step 5-6. Define the pooled 35 mM reference centroid and compute
# centroid distance for non-reference treatments.
treatment_summary_main <- add_centroid_distance(
  summary_df = treatment_summary,
  score_df = pc_scores,
  ref_levels = c(35)
)

# Step 7. Test mean distance to centroid ~ nitrogen (n = 8 treatments).
mean_dist_tests_pc <- treatment_summary_main %>%
  run_trend_tests("mean_dist_to_centroid")

# Step 8. Test centroid distance ~ nitrogen (n = 7 treatments; 35 mM excluded).
centroid_dist_tests_main <- treatment_summary_main %>%
  filter(!is.na(centroid_distance_to_convergence_region)) %>%
  run_trend_tests("centroid_distance_to_convergence_region")

# Step 9. Sensitivity analyses (diagnostic only).
# These are robustness checks, not alternative analyses to report in the manuscript.

# 9a. Main reference for comparison table: 35 mM alone
centroid_dist_tests_ref_35 <- run_reference_sensitivity(
  summary_df = treatment_summary,
  score_df = pc_scores,
  ref_levels = c(35),
  label = "centroid_distance_pc_ref_35"
)

# 9b. Alternative reference: pooled 25 + 30 + 35 mM
centroid_dist_tests_ref_25_30_35 <- run_reference_sensitivity(
  summary_df = treatment_summary,
  score_df = pc_scores,
  ref_levels = c(25, 30, 35),
  label = "centroid_distance_pc_ref_25_30_35"
)

# 9c. Alternative reference: pooled 30 + 35 mM
centroid_dist_tests_ref_30_35 <- centroid_dist_tests_main %>%
  transmute(
    analysis = "centroid_distance_pc_ref_30_35",
    spearman_rho,
    spearman_p
  )

# 9d. Mean distance to centroid in the original 4-D standardized trait space.
# This applies only to the within-treatment dispersion metric, not to the
# directional centroid-distance metric.
trait_scores_4d <- pca_data %>%
  select(barcodeID, species, treatment_mmol) %>%
  bind_cols(as_tibble(pca_input_scaled)) %>%
  mutate(treatment_mmol = as.numeric(treatment_mmol))

trait_summary_4d <- trait_scores_4d %>%
  group_by(treatment_mmol) %>%
  mutate(
    centroid_LMA = mean(LMA),
    centroid_LDMC = mean(LDMC),
    centroid_CHL = mean(CHL),
    centroid_GRT = mean(GRT),
    dist_to_centroid_4d = sqrt(
      (LMA - centroid_LMA)^2 +
        (LDMC - centroid_LDMC)^2 +
        (CHL - centroid_CHL)^2 +
        (GRT - centroid_GRT)^2
    )
  ) %>%
  summarise(
    mean_dist_to_centroid_4d = mean(dist_to_centroid_4d),
    .groups = "drop"
  ) %>%
  arrange(treatment_mmol)

mean_dist_tests_4d <- trait_summary_4d %>%
  run_trend_tests("mean_dist_to_centroid_4d")

# Comparison table for sensitivity checks: Spearman rho and p-value only
sensitivity_table <- bind_rows(
  mean_dist_tests_pc %>%
    transmute(
      analysis = "mean_dist_to_centroid_pc",
      spearman_rho,
      spearman_p
    ),
  centroid_dist_tests_ref_30_35,
  centroid_dist_tests_ref_35,
  centroid_dist_tests_ref_25_30_35,
  mean_dist_tests_4d %>%
    transmute(
      analysis = "mean_dist_to_centroid_4d",
      spearman_rho,
      spearman_p
    )
)

# Step 10. Print summary table and main test results to console.
cat("\nTreatment-level summary table:\n")
print(treatment_summary_main, n = Inf)

cat("\nVariance explained by PC1 and PC2:\n")
print(explained_variation[1:2])

cat("\nMain trend tests: mean distance to centroid in PC1-PC2\n")
print(mean_dist_tests_pc, n = Inf)

cat("\nMain trend tests: centroid distance to 35 mM reference\n")
print(centroid_dist_tests_main, n = Inf)

cat("\nSensitivity checks (Spearman only):\n")
print(sensitivity_table, n = Inf)

# Step 11. PERMANOVA on PC1-PC2 scores across treatments.
library(vegan)

pc_matrix <- as.matrix(pc_scores[, c("PC1", "PC2")])
treatment_factor <- factor(pc_scores$treatment_mmol)

set.seed(42)
permanova_overall <- adonis2(
  pc_matrix ~ treatment_factor,
  method = "euclidean",
  permutations = 999
)

cat("\nOverall PERMANOVA (PC1 + PC2 ~ treatment):\n")
print(permanova_overall)

# Pairwise PERMANOVA with Bonferroni correction
treatment_levels <- sort(unique(pc_scores$treatment_mmol))
pairs <- combn(treatment_levels, 2, simplify = FALSE)

pairwise_permanova <- map_dfr(pairs, function(pair) {
  sub_df <- pc_scores %>% filter(treatment_mmol %in% pair)
  sub_matrix <- as.matrix(sub_df[, c("PC1", "PC2")])
  sub_treatment <- factor(sub_df$treatment_mmol)
  
  set.seed(42)
  fit <- adonis2(
    sub_matrix ~ sub_treatment,
    method = "euclidean",
    permutations = 999
  )
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    F_value = round(fit$F[1], 3),
    R2    = round(fit$R2[1], 3),
    p_value = fit$`Pr(>F)`[1]
  )
})

pairwise_permanova <- pairwise_permanova %>%
  mutate(p_adj_bonferroni = p.adjust(p_value, method = "bonferroni"))

cat("\nPairwise PERMANOVA (Bonferroni-corrected):\n")
print(pairwise_permanova, n = Inf)

# Step 12. Plots.
plot_mean_dist <- ggplot(
  treatment_summary_main,
  aes(x = treatment_mmol, y = mean_dist_to_centroid)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE) +
  labs(
    x = "Nitrogen addition (mM)",
    y = "Mean distance to centroid"
  ) +
  custom_theme +
  theme(aspect.ratio = 1)

plot_centroid_dist <- ggplot(
  treatment_summary_main %>%
    filter(!is.na(centroid_distance_to_convergence_region)),
  aes(x = treatment_mmol, y = centroid_distance_to_convergence_region)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  labs(
    x = "Nitrogen addition (mM)",
    y = "Centroid distance to convergence region (unitless)"
  ) +
  custom_theme +
  theme(aspect.ratio = 1)

convergence_plot <- ggarrange(
  plot_mean_dist,
  plot_centroid_dist,
  nrow = 1,
  labels = c("a", "b")
)

print(convergence_plot)

# Save directional convergence plot as PNG
ggsave(
  filename = "./figures/convergence_plot.png",
  plot = plot_centroid_dist,
  width = 5, height = 5,
  bg = "white"
)

# Save directional convergence plot as PNG
library(svglite)
ggsave(
  filename = "./figures/convergence_plot.svg",
  plot = plot_centroid_dist,
  width = 8, height = 6,
  bg = "white"
)

# Plot Figure 2 ----
# PCA on left, directional convergence on right

plot_centroid_dist_pretty <- ggplot(
  treatment_summary_main %>%
    filter(!is.na(centroid_distance_to_convergence_region)),
  aes(
    x = treatment_mmol,
    y = centroid_distance_to_convergence_region,
    color = factor(treatment_mmol),
    shape = factor(treatment_mmol)
  )
) +
  geom_point(size = 2.8) +
  geom_smooth(
    data = treatment_summary_main %>%
      filter(!is.na(centroid_distance_to_convergence_region)),
    aes(
      x = treatment_mmol,
      y = centroid_distance_to_convergence_region
    ),
    method = "lm",
    formula = y ~ x,
    color = "black",
    se = FALSE,
    inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  scale_color_manual(values = custom_colors1) +
  scale_shape_manual(values = custom_shapes) +
  labs(
    x = "Nitrogen addition (mM)",
    y = "Centroid distance to convergence region (unitless)"
  ) +
  custom_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )


figure2 <- ggarrange(
  pca_plot_treatment_ellipse,
  plot_centroid_dist_pretty,
  nrow = 1,
  labels = c("a", "b"),
  widths = c(1.2, 1),
  align = "hv"
)

print(figure2)

# Save combined figure as PNG
ggsave(
  filename = "./figures/Fig_2.png",
  plot = figure2,
  width = 10,
  height = 5,
  bg = "white"
)

# Save combined figure as SVG
library(svglite)
ggsave(
  filename = "./figures/Fig_2.svg",
  plot = figure2,
  width = 16,
  height = 6,
  bg = "white"
)

# Boxplots of aboveground biomass per species-treatment combination (Fig. S5)
# and median aboveground biomass-to-volume ratios (Table S2)
# for Malamud et al. Nitrogen trait convergence study
# Sean Michaletz (sean.michaletz@ubc.ca), 13 May 2026
# Revised 8 June 2026

# 0. Initialize ----

# Load data
df <- read.csv("./data/traits.csv",header=T)

# Libraries
library(ggplot2)

# Custom colors
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# Custom theme
custom_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 9, face = "italic"),
    strip.background = element_blank()
  )

# Treat treatment as a discrete factor for boxplot x-axis
df$treatment_mmol_f <- factor(df$treatment_mmol)

# 1. Build plot ----
p_S5 <- ggplot(df, aes(x = treatment_mmol_f, y = dry_whole_g,
                       color = species, fill = species)) +
  geom_boxplot(alpha = 0.25, linewidth = 0.6,
               outlier.shape = 16, outlier.size = 1.2) +
  facet_wrap(~ species, nrow = 1) +
  labs(
    x = "Nitrogen addition (mM)",
    y = "Aboveground biomass (g)",
    color = "Species", fill = "Species"
  ) +
  scale_color_manual(values = josef_colors) +
  scale_fill_manual(values = josef_colors) +
  custom_theme +
  theme(
    legend.position = "none",          # species identity carried by facet strips
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_S5)

# 2. Save ----
ggsave("Figure_S5.jpg", p_S5,
       width = 8.5, height = 3.75, dpi = 300, bg = "white")

library(svglite)
ggsave("Figure_S5.svg", p_S5,
       width = 8.5, height = 3.75, bg = "white")

# 3. Median biomass-to-volume ratios (Table S2) ----

pot_volume_L <- 0.7
df$bvr_g_per_L <- df$dry_whole_g / pot_volume_L

med_biomass <- aggregate(dry_whole_g ~ species + treatment_mmol, data = df,
                         FUN = function(x) round(median(x, na.rm = TRUE), 2))
med_bvr     <- aggregate(bvr_g_per_L ~ species + treatment_mmol, data = df,
                         FUN = function(x) round(median(x, na.rm = TRUE), 2))
n_tab       <- aggregate(dry_whole_g ~ species + treatment_mmol, data = df,
                         FUN = length)

bvr_summary <- merge(med_biomass, med_bvr, by = c("species", "treatment_mmol"))
bvr_summary <- merge(bvr_summary, n_tab,   by = c("species", "treatment_mmol"))
names(bvr_summary) <- c("species", "treatment_mmol",
                        "median_biomass_g", "median_bvr_g_per_L", "n")
bvr_summary <- bvr_summary[order(bvr_summary$species, bvr_summary$treatment_mmol), ]
print(bvr_summary)

# Wide matrix for the table; NA when a species-treatment combination is 
# absent (i.e. H. vulgare at 0 mM)
bvr_wide <- tapply(df$bvr_g_per_L,
                   list(treatment_mmol = df$treatment_mmol, species = df$species),
                   function(x) round(median(x, na.rm = TRUE), 2))
# reorder columns to manuscript species order
sp_order <- c("B. officinalis", "H. vulgare", "R. sativus")
bvr_wide <- bvr_wide[, sp_order[sp_order %in% colnames(bvr_wide)]]
print(bvr_wide)


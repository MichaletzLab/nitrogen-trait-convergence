# Boxplots of aboveground biomass per species-treatment combination (Fig. S5)
# for Malamud et al. Nitrogen trait convergence study
# Sean Michaletz (sean.michaletz@ubc.ca), 13 May 2026

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

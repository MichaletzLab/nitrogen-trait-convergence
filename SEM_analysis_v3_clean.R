# SEM analysis for Malamud et al. Nitrogen trait convergence study
# Sean Michaletz (sean.michaletz@ubc.ca), 2 Feb 2025
# Updated 20 Feb 2025, 18 March 2026

# 0. Initialize ----

# Load data
df <- read.csv("./data/traits.csv",header=T)

# Calculate aboveground growth rate
df$grt_g_d <- df$dry_whole_g / 6*7 #days in experiment


# 1. Summarize data ----
library(dplyr)

df_summary1 <- df %>%
  group_by(species, treatment_mmol) %>%
  summarise(count = n(), .groups = "drop")
print(df_summary1)

df_summary2 <- df %>%
  group_by(species, treatment_mmol) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(species) %>%
  summarise(avg_count = median(count), .groups = "drop")
print(df_summary2)


# 2. Trait-nitrogen relationships ----

### Model selection ----
# Compete models for log(trait) ~ treatment_mmol using:
# - Fixed effects only (lm)
# - Random intercept, random slope, and random slope+intercept by species (lmer)

# Load required packages
library(lmerTest)    # for lmer() with p-values
library(MuMIn)       # for r.squaredGLMM()
library(dplyr)

# Define the traits and predictor (only treatment_mmol)
traits <- c("LMA", "LDMC", "CHL")
predictors <- c("treatment_mmol")
predictor_labels <- c("treatment_mmol")

# Define random-effects types, including fixed-effects-only model
random_effects_types <- c("fixed", "intercept", "slope", "intercept+slope")

# Initialize a list to store results
results_list <- list()

# Loop over each trait, the predictor, and each random-effects structure
for (trait in traits) {
  for (i in seq_along(predictors)) {
    pred <- predictors[i]
    pred_label <- predictor_labels[i]
    
    for (re_type in random_effects_types) {
      # Define the random-effects part of the formula based on the type
      if (re_type == "fixed") {
        re_str <- ""  # No random effects
      } else if (re_type == "intercept") {
        re_str <- "(1 | species)"
      } else if (re_type == "slope") {
        re_str <- paste0("(0 + ", pred, " | species)")
      } else if (re_type == "intercept+slope") {
        re_str <- paste0("(", pred, " | species)")
      }
      
      # Construct formula differently for fixed-effects model
      if (re_type == "fixed") {
        formula_str <- paste0("log10(", trait, ") ~ ", pred)
      } else {
        formula_str <- paste0("log10(", trait, ") ~ ", pred, " + ", re_str)
      }
      
      formula_obj <- as.formula(formula_str)
      
      # Fit the model: use lm() for fixed-effects only, otherwise lmer()
      model_fit <- tryCatch(
        if (re_type == "fixed") {
          lm(formula_obj, data = df)
        } else {
          lmer(formula_obj, data = df)
        },
        error = function(e) NULL
      )
      
      # If the model failed to fit, record NAs and continue
      if (is.null(model_fit)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          trait = trait,
          predictor = pred_label,
          random_effects = re_type,
          singular = NA,
          AIC = NA,
          p_value = NA,
          r2_marginal = NA,
          r2_conditional = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # Check whether the fit was singular (only for mixed models)
      sing <- if (re_type == "fixed") FALSE else isSingular(model_fit, tol = 1e-4)
      
      # Get the AIC
      aic_val <- AIC(model_fit)
      
      # Extract the p-value for the fixed effect treatment_mmol
      summ <- summary(model_fit)
      p_val <- summ$coefficients[2, "Pr(>|t|)"]
      
      # Compute marginal and conditional R2
      if (re_type == "fixed") {
        r2_marg <- summary(model_fit)$r.squared
        r2_cond <- NA  # No conditional R² for fixed models
      } else {
        r2_vals <- r.squaredGLMM(model_fit)
        r2_marg <- r2_vals[1]
        r2_cond <- r2_vals[2]
      }
      
      # Store the results
      results_list[[length(results_list) + 1]] <- data.frame(
        trait = trait,
        predictor = pred_label,
        random_effects = re_type,
        singular = sing,
        AIC = aic_val,
        p_value = p_val,
        r2_marginal = r2_marg,
        r2_conditional = r2_cond,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine all results into one data frame
results_df <- do.call(rbind, results_list)

# Subset only non-singular fits (for mixed models)
results_df <- subset(results_df, is.na(singular) | singular == FALSE)

# Print results
print(results_df)
# AIC prefers only fixed effect for LMA and CHL, and random intercept for LDMC.


### Plots ----
# Make plots, fitting regressions fixed effects only for LMA and CHL, and
# with species as random intercept for LDMC

# Load required packages
library(ggplot2)
library(lme4)
library(dplyr)
library(rlang)

# Define a function to plot a given trait against treatment_mmol
plot_trait <- function(trait) {
  # Determine which model to use based on AIC results
  if (trait %in% c("LMA", "CHL")) {
    # Fixed-effects model (lm)
    model_formula <- as.formula(paste0("log10(", trait, ") ~ treatment_mmol"))
    model <- lm(model_formula, data = df)
    
    # Create new data for prediction
    pred_data <- data.frame(
      treatment_mmol = seq(min(df$treatment_mmol), max(df$treatment_mmol), length.out = 100)
    )
    
    # Get predictions
    pred_data$predicted_log_trait <- predict(model, newdata = pred_data)
    
    # Generate plot (black regression line for fixed effects)
    p <- ggplot(df, aes(x = treatment_mmol, y = !!sym(trait), color = species)) +
      geom_point() +
      geom_line(data = pred_data, aes(x = treatment_mmol, y = 10^(predicted_log_trait)), 
                color = "black", linewidth = 1) + 
      scale_y_log10() +
      labs(x = "Nitrogen treatment (mmol)",
           y = trait,
           title = paste(trait, "vs. Nitrogen Treatment"),
           color = "Species") +
      theme_minimal()
    
  } else if (trait == "LDMC") {
    # Mixed-effects model (lmer) with random intercept for species
    model_formula <- as.formula(paste0("log10(", trait, ") ~ treatment_mmol + (1 | species)"))
    model <- lmer(model_formula, data = df)
    
    # Create new data for prediction: a grid of treatment_mmol values for each species
    pred_data <- expand.grid(
      treatment_mmol = seq(min(df$treatment_mmol), max(df$treatment_mmol), length.out = 100),
      species = unique(df$species)
    )
    
    # Get predictions (random intercept model)
    pred_data$predicted_log_trait <- predict(model, newdata = pred_data, re.form = NULL)
    
    # Generate plot (colored regression lines per species)
    p <- ggplot(df, aes(x = treatment_mmol, y = !!sym(trait), color = species)) +
      geom_point() +
      geom_line(data = pred_data, aes(x = treatment_mmol, y = 10^(predicted_log_trait), color = species), 
                linewidth = 1) +  
      scale_y_log10() +
      labs(x = "Nitrogen treatment (mmol)",
           y = trait,
           title = paste(trait, "vs. Nitrogen Treatment"),
           color = "Species") +
      theme_minimal()
  } else {
    stop("Trait not recognized!")
  }
  
  return(p)
}

# Generate plots for each of the three traits.
p_LMA  <- plot_trait("LMA")
p_LDMC <- plot_trait("LDMC")
p_CHL  <- plot_trait("CHL")

# Display the plots
print(p_LMA)
print(p_LDMC)
print(p_CHL)


# 3. Growth-Nitrogen Relationships ----

### Model Selection ----
library(lme4)        # For lmer and singular fit checking
library(MuMIn)       # For r.squaredGLMM()
library(performance) # For alternative R^2 calculation

# Define candidate models: fixed effects, random intercept, random slope, and random intercept+slope
models <- list(
  FE = function(reml) {
    lm(log10(grt_g_d) ~ treatment_mmol, data = df)  # Fixed effects only
  },
  RI = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (1 | species),
         data = df, REML = reml)
  },
  RS = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (0 + treatment_mmol | species),
         data = df, REML = reml)
  },
  RIS = function(reml) {
    lmer(log10(grt_g_d) ~ treatment_mmol + (1 + treatment_mmol | species),
         data = df, REML = reml)
  }
)

# Prepare an empty results data frame
results_df <- data.frame(
  Model          = character(),
  AIC            = numeric(),
  p_value        = numeric(),
  R2_marginal    = numeric(),
  R2_conditional = numeric(),
  Singular_Fit   = logical(),
  stringsAsFactors = FALSE
)

# Loop over each candidate model
for(model_name in names(models)){
  
  # Track if any step emits a singular-fit warning (only for mixed models)
  singular_warning <- FALSE
  
  # Fit the model (handling fixed-effects model separately)
  if (model_name == "FE") {
    fit_ml <- models[[model_name]](reml = FALSE)
    
    # Extract AIC
    model_aic <- AIC(fit_ml)
    
    # Extract p-value for treatment_mmol
    model_summary <- summary(fit_ml)
    p_val <- model_summary$coefficients["treatment_mmol", "Pr(>|t|)"]
    
    # Compute R² for fixed-effects model using summary()
    r2_marg <- model_summary$r.squared
    r2_cond <- NA  # Conditional R² does not apply to fixed-effects models
    
  } else {
    # Mixed-effects models: Fit model with ML (REML = FALSE) for AIC and p-values
    fit_ml <- withCallingHandlers(
      models[[model_name]](reml = FALSE),
      warning = function(w) {
        if(grepl("singular", w$message, ignore.case = TRUE)) {
          singular_warning <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    
    # Extract AIC
    model_aic <- AIC(fit_ml)
    
    # Extract p-value for treatment_mmol
    model_summary <- summary(fit_ml)
    p_val <- NA
    if("treatment_mmol" %in% rownames(model_summary$coefficients)){
      p_val <- model_summary$coefficients["treatment_mmol", "Pr(>|t|)"]
    }
    
    # Refit with REML for R² calculation
    fit_reml <- withCallingHandlers(
      update(fit_ml, REML = TRUE),
      warning = function(w) {
        if(grepl("singular", w$message, ignore.case = TRUE)) {
          singular_warning <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    
    # Compute R² using MuMIn
    r2_vals <- withCallingHandlers(
      r.squaredGLMM(fit_reml),
      warning = function(w) {
        if(grepl("singular", w$message, ignore.case = TRUE)) {
          singular_warning <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    )
    
    r2_marg <- r2_vals["R2m"]
    r2_cond <- r2_vals["R2c"]
    
    # If R² is missing, try performance::r2()
    if(is.na(r2_marg) || is.na(r2_cond)){
      r2_perf <- withCallingHandlers(
        r2(fit_reml),
        warning = function(w) {
          if(grepl("singular", w$message, ignore.case = TRUE)) {
            singular_warning <<- TRUE
          }
          invokeRestart("muffleWarning")
        }
      )
      r2_marg <- r2_perf$R2_marginal
      r2_cond <- r2_perf$R2_conditional
    }
  }
  
  # Append results to the data frame
  results_df <- rbind(
    results_df,
    data.frame(
      Model            = model_name,
      AIC              = model_aic,
      p_value          = p_val,
      R2_marginal      = r2_marg,
      R2_conditional   = r2_cond,
      Singular_Fit     = ifelse(model_name == "FE", FALSE, singular_warning), # Only applies to mixed models
      stringsAsFactors = FALSE
    )
  )
}

# Subset only non-singular fits for mixed models
results_df <- subset(results_df, is.na(Singular_Fit) | Singular_Fit == FALSE)

# Print the summary data frame
print(results_df)
# AIC prefers random intercept by species


### Plot ----
# Use random intercept by species
library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)
library(purrr)

# Fit the mixed-effects model on the log-transformed response
model <- lmer(log10(grt_g_d) ~ treatment_mmol + (1 | species), data = df)

# Create new data for predictions for each species over the range of treatment_mmol
newdata <- df %>%
  group_by(species) %>%
  summarise(min_t = min(treatment_mmol), max_t = max(treatment_mmol)) %>%
  mutate(treatment_mmol = map2(min_t, max_t, ~ seq(.x, .y, length.out = 100))) %>%
  unnest(treatment_mmol) %>%
  select(species, treatment_mmol)

# Get predictions from the model (predictions are on the log10 scale)
newdata$pred <- predict(model, newdata = newdata)

# Add a new column to df for the log-transformed growth rate
df <- df %>% mutate(log_grt = log10(grt_g_d))

# Plot the data and the regression lines
p_grt <- ggplot(df, aes(x = treatment_mmol, y = log_grt, color = species)) +
  geom_point() +
  geom_line(data = newdata, aes(x = treatment_mmol, y = pred, color = species)) +
  labs(x = "Nitrogen treatment (mmol)",
       y = "Log10(Growth rate) (g/day)",
       title = "Growth rate vs. Nitrogen Treatment (Random Intercept Model)",
       color = "Species") +
  theme_bw()

print(p_grt)


# 4. Trait-trait relationships ----

### SMA regressions ----
library(smatr)
sma1 <- sma(LDMC ~ LMA * species, data = df, log = "xy")
summary(sma1) # all significant
sma2 <- sma(CHL ~ LMA * species, data = df, log = "xy")
summary(sma2) # B. officinalis and H. vulgare not significant
sma3 <- sma(CHL ~ LDMC * species, data = df, log = "xy")
summary(sma3) # B. officinalis and H. vulgare not significant (R. sativus is marginally significant)

### Plots ----
# Make plots, including regression lines only when significant
library(ggplot2)
library(ggpmisc) # For stat_ma_line()
library(dplyr)
library(gridExtra)

# A function that creates a pairwise plot with points and, conditionally, SMA regression lines by species.
plot_pair <- function(xvar, yvar) {
  p <- ggplot(df, aes_string(x = xvar, y = yvar, color = "species")) +
    geom_point() +
    # Plot raw data with log10-transformed axes.
    scale_x_log10() +
    scale_y_log10() +
    labs(x = xvar, y = yvar) +
    theme_bw()
  
  # Add regression lines conditionally based on xvar and yvar.
  if(xvar == "LMA" && yvar == "CHL") {
    # For CHL ~ LMA, add regression lines only for species other than B. officinalis and H. vulgare.
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = xvar, y = yvar, group = "species"),
      method = "SMA", se = FALSE, linewidth = 1
    )
  } else if(xvar == "LDMC" && yvar == "CHL") {
    # For CHL ~ LDMC, add regression lines only for species other than B. officinalis and H. vulgare.
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = xvar, y = yvar, group = "species"),
      method = "SMA", se = FALSE, linewidth = 1
    )
  } else {
    # For other pairs (e.g., LMA vs. LDMC), add regression lines for all species.
    p <- p + stat_ma_line(aes_string(group = "species"), method = "SMA", se = FALSE, linewidth = 1)
  }
  
  return(p)
}

# Create the three pairwise plots:
p1 <- plot_pair("LMA", "LDMC")   # Regression lines for all species.
p2 <- plot_pair("LMA", "CHL")    # Regression lines only for species other than B. officinalis and H. vulgare.
p3 <- plot_pair("LDMC", "CHL")   # Regression lines only for species other than B. officinalis and H. vulgare.

# Arrange plots
gridExtra::grid.arrange(p1, p2, p3, ncol = 2)


# 5. Growth-trait relationships ----

### SMA regressions ----
library(smatr)
sma4 <- sma(grt_g_d ~ LMA * species, data = df, log = "xy")
summary(sma4) # B. officinalis and H. vulgare not significant
sma5 <- sma(grt_g_d ~ LDMC * species, data = df, log = "xy")
summary(sma5) # B. officinalis and H. vulgare not significant
sma6 <- sma(grt_g_d ~ CHL * species, data = df, log = "xy")
summary(sma6) # All significant

### Plots ----
# Make plots, including regression lines only when significant
# Load required packages
library(ggplot2)
library(ggpmisc)   # Provides stat_ma_line() for Model II SMA lines
library(dplyr)
library(gridExtra)

# A function to create a plot for growth rate vs. a given trait.
# Both axes are log-scaled (without transforming the underlying data).
# A separate SMA regression line is fit for each species,
# but for grt_g_d ~ LMA, the line is omitted for B. officinalis only;
# for grt_g_d ~ LDMC, lines are omitted for both B. officinalis and H. vulgare.
plot_growth_vs_trait <- function(trait) {
  p <- ggplot(df, aes_string(x = trait, y = "grt_g_d", color = "species")) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = trait,
         y = "Growth rate (g/d)") +
    theme_bw()
  
  # For LMA, include a regression line for H. vulgare but not B. officinalis;
  # For LDMC, omit regression lines for both B. officinalis and H. vulgare;
  # For other traits (e.g., CHL), add regression lines for all species.
  if (trait == "LMA") {
    p <- p + stat_ma_line(
      data = subset(df, species != "B. officinalis"),  # H. vulgare significant, B. officinalis not
      mapping = aes_string(x = trait, y = "grt_g_d", group = "species"),
      method = "SMA",
      se = FALSE,
      linewidth = 1
    )
  } else if (trait == "LDMC") {
    p <- p + stat_ma_line(
      data = subset(df, !(species %in% c("B. officinalis", "H. vulgare"))),
      mapping = aes_string(x = trait, y = "grt_g_d", group = "species"),
      method = "SMA",
      se = FALSE,
      linewidth = 1
    )
  } else {
    p <- p + stat_ma_line(aes_string(group = "species"),
                          method = "SMA",
                          se = FALSE,
                          linewidth = 1)
  }
  
  return(p)
}

# Create the three plots:
p_LMA  <- plot_growth_vs_trait("LMA")   
p_LDMC <- plot_growth_vs_trait("LDMC")  # grt_g_d vs. LDMC: regression lines omitted for B. officinalis and H. vulgare
p_CHL  <- plot_growth_vs_trait("CHL")   # grt_g_d vs. CHL: regression lines for all species

# Arrange the plots in a matrix-style layout (e.g., 2 columns)
gridExtra::grid.arrange(p_LMA, p_LDMC, p_CHL, ncol = 2)


### Model selection ----

#### Simple Bivariate Models ----

# Load required packages
library(lmerTest)    # for lmer() with p-values
library(MuMIn)       # for r.squaredGLMM()
library(dplyr)

# The three traits to test.
traits <- c("LMA", "LDMC", "CHL")

# Four model types:
# 1) Fixed-effects only
# 2) Random intercept
# 3) Random slope
# 4) Random intercept + slope
model_types <- c("fixed", "intercept", "slope", "intercept+slope")

# Initialize a list to store model-fitting results.
results_list <- list()

# Loop over each trait and each model type.
for (trait in traits) {
  for (model_type in model_types) {
    
    # Define model formula based on model type
    if (model_type == "fixed") {
      # Fixed-effects only model
      formula_str <- paste0("log10(grt_g_d) ~ log10(", trait, ")")
    } else if (model_type == "intercept") {
      # Random intercept model
      re_str <- "(1 | species)"
      formula_str <- paste0("log10(grt_g_d) ~ log10(", trait, ") + ", re_str)
    } else if (model_type == "slope") {
      # Random slope only model (no random intercept)
      re_str <- paste0("(0 + log10(", trait, ") | species)")
      formula_str <- paste0("log10(grt_g_d) ~ log10(", trait, ") + ", re_str)
    } else if (model_type == "intercept+slope") {
      # Random intercept + random slope model
      re_str <- paste0("(log10(", trait, ") | species)")
      formula_str <- paste0("log10(grt_g_d) ~ log10(", trait, ") + ", re_str)
    }
    
    formula_obj <- as.formula(formula_str)
    
    # Fit the model (use `lm()` for fixed effects, `lmer()` for mixed effects)
    model_fit <- tryCatch(
      {
        if (model_type == "fixed") {
          lm(formula_obj, data = df)
        } else {
          lmer(formula_obj, data = df)
        }
      },
      error = function(e) {
        NULL  # If model fitting fails, return NULL
      }
    )
    
    # If the model failed, store NA values and continue.
    if (is.null(model_fit)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        trait           = trait,
        model_type      = model_type,
        singular        = NA,
        AIC             = NA,
        p_value         = NA,
        r2_marginal     = NA,
        r2_conditional  = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Check for singular fit (only applies to mixed models)
    sing <- if (model_type == "fixed") FALSE else isSingular(model_fit, tol = 1e-4)
    
    # Extract AIC
    aic_val <- AIC(model_fit)
    
    # Extract p-value for the fixed effect "log10(trait)"
    coefs <- summary(model_fit)$coefficients
    row_name <- paste0("log10(", trait, ")")
    p_val <- if (row_name %in% rownames(coefs)) coefs[row_name, "Pr(>|t|)"] else NA
    
    # Compute R²
    if (model_type == "fixed") {
      # For fixed-effects model (lm), use summary() for R²
      r2_marg <- summary(model_fit)$r.squared
      r2_cond <- NA  # No conditional R² for fixed-effects models
    } else {
      # For mixed-effects models (lmer), use MuMIn::r.squaredGLMM
      r2_vals <- r.squaredGLMM(model_fit)
      r2_marg <- r2_vals[1]
      r2_cond <- r2_vals[2]
    }
    
    # Store the results
    results_list[[length(results_list) + 1]] <- data.frame(
      trait           = trait,
      model_type      = model_type,
      singular        = sing,
      AIC             = aic_val,
      p_value         = p_val,
      r2_marginal     = r2_marg,
      r2_conditional  = r2_cond,
      stringsAsFactors = FALSE
    )
  }
}

# Combine all results into one final data frame.
results_df <- do.call(rbind, results_list)

# Subset only non-singular fits (for mixed models)
results_df <- subset(results_df, is.na(singular) | singular == FALSE)

# Print the results
print(results_df)

# AIC prefers random slope for LMA, random slope (or secondarily intercept) for 
# LDMC, and random intercept for CHL.


#### Multivariate Models ----

# Load required packages
library(lmerTest)   # For lmer() with p-values
library(MuMIn)      # For r.squaredGLMM()
library(dplyr)

# Define the four model options: fixed effects, random intercept, random slope, and random intercept+slope
re_options <- c("fixed", "intercept", "slope", "intercept+slope")

# Create a dataframe with all combinations for LMA, LDMC, CHL.
combo_df <- expand.grid(LMA_opt = re_options,
                        LDMC_opt = re_options,
                        CHL_opt = re_options,
                        stringsAsFactors = FALSE)

# Initialize a list to store results.
results_list <- list()

# Loop over each combination (64 total).
for(i in seq_len(nrow(combo_df))) {
  # Extract the options for each predictor
  opt_LMA   <- combo_df$LMA_opt[i]
  opt_LDMC  <- combo_df$LDMC_opt[i]
  opt_CHL   <- combo_df$CHL_opt[i]
  
  # Build the random effects terms.
  re_terms <- c()
  intercept_flag <- FALSE
  use_lm <- FALSE  # Flag for fixed-effects model
  
  # For LMA:
  if(opt_LMA == "fixed") {
    use_lm <- TRUE
  } else if(opt_LMA == "slope") {
    re_terms <- c(re_terms, "(0 + log10(LMA)|species)")
  } else if(opt_LMA == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(LMA)|species)")
  } else if(opt_LMA == "intercept") {
    intercept_flag <- TRUE
  }
  
  # For LDMC:
  if(opt_LDMC == "fixed") {
    use_lm <- TRUE
  } else if(opt_LDMC == "slope") {
    re_terms <- c(re_terms, "(0 + log10(LDMC)|species)")
  } else if(opt_LDMC == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(LDMC)|species)")
  } else if(opt_LDMC == "intercept") {
    intercept_flag <- TRUE
  }
  
  # For CHL:
  if(opt_CHL == "fixed") {
    use_lm <- TRUE
  } else if(opt_CHL == "slope") {
    re_terms <- c(re_terms, "(0 + log10(CHL)|species)")
  } else if(opt_CHL == "intercept+slope") {
    re_terms <- c(re_terms, "(log10(CHL)|species)")
  } else if(opt_CHL == "intercept") {
    intercept_flag <- TRUE
  }
  
  # If any predictor was set to "intercept", add a single random intercept.
  if(intercept_flag) {
    re_terms <- c(re_terms, "(1|species)")
  }
  
  # Remove duplicate intercept terms.
  re_terms <- unique(re_terms)
  
  # Build the full random effects string and the full model formula.
  if(length(re_terms) > 0 & !use_lm) {
    re_str <- paste(re_terms, collapse = " + ")
    full_formula_str <- paste0("log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL) + ", re_str)
  } else {
    full_formula_str <- "log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL)"
  }
  
  # Convert the string to a formula.
  model_formula <- as.formula(full_formula_str)
  
  # Fit the model
  fit <- tryCatch(
    {
      if (use_lm) {
        lm(model_formula, data = df)  # Fixed-effects model
      } else {
        suppressWarnings(lmer(model_formula, data = df))  # Mixed-effects model
      }
    },
    error = function(e) NULL  # If model fitting fails, return NULL
  )
  
  # If the model failed, record NA's.
  if(is.null(fit)) {
    results_list[[length(results_list) + 1]] <- data.frame(
      LMA_opt = opt_LMA,
      LDMC_opt = opt_LDMC,
      CHL_opt = opt_CHL,
      re_formula = full_formula_str,
      AIC = NA,
      singular = NA,
      p_logLMA = NA,
      p_logLDMC = NA,
      p_logCHL = NA,
      r2_marginal = NA,
      r2_conditional = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # Check for singular fit (only for mixed models)
  is_sing <- if (use_lm) FALSE else isSingular(fit, tol = 1e-4)
  
  # Extract AIC
  aic_val <- AIC(fit)
  
  # Get summary and extract fixed-effect coefficients.
  summ <- summary(fit)
  coefs <- summ$coefficients
  p_logLMA  <- if("log10(LMA)" %in% rownames(coefs))  coefs["log10(LMA)", "Pr(>|t|)"]  else NA
  p_logLDMC <- if("log10(LDMC)" %in% rownames(coefs)) coefs["log10(LDMC)", "Pr(>|t|)"] else NA
  p_logCHL  <- if("log10(CHL)" %in% rownames(coefs))  coefs["log10(CHL)", "Pr(>|t|)"]  else NA
  
  # Compute R²
  if (use_lm) {
    r2_marg <- summ$r.squared
    r2_cond <- NA  # No conditional R² for fixed-effects models
  } else {
    r2_vals <- try(r.squaredGLMM(fit), silent = TRUE)
    if(inherits(r2_vals, "try-error")) {
      r2_marg <- NA
      r2_cond <- NA
    } else {
      r2_marg <- r2_vals[1]
      r2_cond <- r2_vals[2]
    }
  }
  
  # Store results in a data frame row.
  results_list[[length(results_list) + 1]] <- data.frame(
    LMA_opt = opt_LMA,
    LDMC_opt = opt_LDMC,
    CHL_opt = opt_CHL,
    re_formula = full_formula_str,
    AIC = aic_val,
    singular = is_sing,
    p_logLMA = p_logLMA,
    p_logLDMC = p_logLDMC,
    p_logCHL = p_logCHL,
    r2_marginal = r2_marg,
    r2_conditional = r2_cond,
    stringsAsFactors = FALSE
  )
}

# Combine all results into one dataframe.
results_df <- do.call(rbind, results_list)

# Subset only non-singular fits (for mixed models)
results_df <- subset(results_df, is.na(singular) | singular == FALSE)

# Print the results
print(results_df)
# AIC equally prefers:
# 1) slope for LMA, intercept for LDMC and CHL
# 2) intercept for LMA, slope for LDMC and CHL
# 3) intercept+slope for LMA, intercept for LDMC and CHL
# To be most consistent with AIC results for simple bivariate, use slope for LMA
# and intercept for LDMC and CHL


# 6. SEM analysis ----

library(lme4)
library(piecewiseSEM)

### SEM v0 (starting model) ----
# All hypothesized paths and correlations from the MPF framework (Fig. 1b).
# For direct paths from nitrogen to traits, LMA and CHL use fixed effects only,
# while LDMC uses random intercepts by species
m0.1 <- lm(log10(LMA)  ~ treatment_mmol, data = df)
m0.2 <- lmer(log10(LDMC)  ~ treatment_mmol + (1 | species), data = df)
m0.3 <- lm(log10(CHL) ~ treatment_mmol, data = df)

# Growth rate model with a random intercept by species (the structure preferred
# by AIC for LDMC and CHL in bivariate selection) and a random slope for LMA
# by species (preferred by AIC for LMA in bivariate selection).
m0.4 <- lmer(
  log10(grt_g_d) ~ log10(LMA) + log10(LDMC) + log10(CHL) + treatment_mmol +
    (1 | species) +             # Random intercept for species
    (0 + log10(LMA) | species), # Random slope for LMA by species (no intercept)
  data = df
)


# Combine models into a piecewise SEM
sem0 <- psem(
  m0.1,
  m0.2,
  m0.3,
  m0.4,
  log10(LMA) %~~% log10(LDMC),
  log10(LMA) %~~% log10(CHL),
  log10(LDMC) %~~% log10(CHL)
)

# View summary
summary(sem0)
# chi-squared = 0, p = 1: model is saturated (no independence claims to test).
# AIC = -532.964. There are 3 non-significant paths.


### SEM v1 ----
# Remove most non-significant relationship: log10(LMA) ~~ log10(CHL)

# Combine models into a piecewise SEM
sem1 <- psem(
  m0.1,
  m0.2,
  m0.3,
  m0.4,
  log10(LMA) %~~% log10(LDMC),
  log10(LDMC) %~~% log10(CHL)
)

# View summary
summary(sem1)
# chi-squared = 0.004, p = 0.952. AIC = -532.964. Two non-significant paths remain.


### SEM v2 (final model) ----
# Remove most non-significant relationship: LDMC ~~ CHL.

# Combine models into a piecewise SEM
sem2 <- psem(
  m0.1,
  m0.2,
  m0.3,
  m0.4,
  log10(LMA) %~~% log10(LDMC)
)

# View summary
summary(sem2)
# chi-squared = 2.818, p = 0.244. AIC = -532.964.
# One non-significant path remains (LMA -> growth rate); evaluated below.


### SEM v3 (test of LMA -> growth rate path) ----
# Remove most non-significant relationship: log10(grt_g_d) ~ log10(LMA)
# This model tests whether that path should be retained via d-separation.

m3.4 <- lmer(
  log10(grt_g_d) ~ log10(LDMC) + log10(CHL) + treatment_mmol +
    (1 | species),             # Random intercept for species
  data = df
)

# Combine models into a piecewise SEM
sem3 <- psem(
  m0.1,
  m0.2,
  m0.3,
  m3.4,
  log10(LMA) %~~% log10(LDMC)
)

# View summary
summary(sem3)
# chi-squared = 9.95, p = 0.019. AIC = -517.393 (DAIC = +15.57 vs. sem2).
# D-separation test indicates missing path from LMA -> growth rate (p = 0.01).
# LMA -> growth rate is therefore retained; sem2 is the final model.


# 7. Multivariate normality test (Fig. S2) ----

# Chi-square Q-Q plot and Henze-Zirkler test on log-transformed SEM variables.
df2 <- df[, c("treatment_mmol", "dry_whole_g", "LMA", "LDMC", "CHL")]

# Q-Q plot (data log-transformed)
df3 <- df2
df3[, c("dry_whole_g", "LMA", "LDMC", "CHL")] <- log(df3[, c("dry_whole_g", "LMA", "LDMC", "CHL")])
library(MVN)
mvn(df3, mvnTest = "hz", multivariatePlot = "qq")



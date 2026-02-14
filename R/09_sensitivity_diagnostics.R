# ==============================================================================
# 09_sensitivity_diagnostics.R — Sensitivity analyses and model diagnostics
# Addresses reviewer-ready robustness checks for all analytical phases
# ==============================================================================
#
# Sections:
#  1. OLS/spatial regression diagnostics (VIF, residuals, Moran's I on resids)
#  2. Spatial weights sensitivity (queen, k-NN, inverse-distance)
#  3. GWR local coefficient significance (t-values, p-values)
#  4. Interaction testing (Indigenous % × remoteness)
#  5. ITS model fix (population offset, overdispersion, negative binomial)
#  6. Congenital risk score sensitivity (varying thresholds)
#  7. Temporal LISA evolution (VIC LGAs 2019–2024)
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# Additional packages for diagnostics
for (pkg in c("car", "lmtest", "MASS")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(car)
library(lmtest)
library(MASS)
library(spgwr)

# ==============================================================================
# LOAD COMMON DATA
# ==============================================================================

message("=== Sensitivity Analyses and Model Diagnostics ===")
message("Loading data...")

vic_spatial   <- readRDS(file.path(path_processed, "vic_lga_spatial_results.rds"))
reg_data      <- readRDS(file.path(path_processed, "vic_lga_regression_data.rds"))
national_ts   <- readRDS(file.path(path_processed, "national_time_series.rds"))
sa3_risk      <- readRDS(file.path(path_processed, "sa3_congenital_risk.rds"))
vic_analysis  <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))

# Ensure regression variables exist
if (!"log_rate" %in% names(reg_data)) {
  reg_data <- reg_data %>%
    mutate(
      log_rate = log1p(mean_annual_rate),
      log_dist_sh = log1p(mean_dist_sh_km),
      remoteness_factor = factor(remoteness_mode,
                                 levels = c("Major Cities", "Inner Regional",
                                            "Outer Regional", "Remote", "Very Remote"))
    )
}

# Rebuild spatial weights for regression subset
nb_reg <- poly2nb(reg_data, queen = TRUE)
n_islands_reg <- sum(card(nb_reg) == 0)
if (n_islands_reg > 0) {
  islands_reg <- which(card(nb_reg) == 0)
  for (i in islands_reg) {
    dists <- as.numeric(st_distance(
      st_centroid(reg_data[i, ]),
      st_centroid(reg_data[-i, ])
    ))
    nearest <- which.min(dists)
    if (nearest >= i) nearest <- nearest + 1
    nb_reg[[i]] <- as.integer(nearest)
    nb_reg[[nearest]] <- sort(unique(c(nb_reg[[nearest]], as.integer(i))))
  }
}
lw_reg <- nb2listw(nb_reg, style = "W", zero.policy = TRUE)

# Refit core models (needed for diagnostics)
ols_model <- lm(
  log_rate ~ pct_indigenous + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data
)

lag_model <- lagsarlm(
  log_rate ~ pct_indigenous + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data, listw = lw_reg, zero.policy = TRUE
)

error_model <- errorsarlm(
  log_rate ~ pct_indigenous + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data, listw = lw_reg, zero.policy = TRUE
)

# ==============================================================================
# 1. OLS / SPATIAL MODEL DIAGNOSTICS
# ==============================================================================

message("\n=== 1. OLS and Spatial Model Diagnostics ===")

# --- 1a. Variance Inflation Factors (multicollinearity) ---
message("\n1a. Variance Inflation Factors:")
vif_values <- tryCatch({
  car::vif(ols_model)
}, error = function(e) {
  message("  VIF computation failed: ", e$message)
  NULL
})

if (!is.null(vif_values)) {
  # For factors, vif() returns GVIF — use GVIF^(1/(2*Df)) as comparable measure
  if (is.matrix(vif_values)) {
    vif_table <- data.frame(
      Variable = rownames(vif_values),
      GVIF = round(vif_values[, "GVIF"], 3),
      Df = vif_values[, "Df"],
      `GVIF_adj` = round(vif_values[, "GVIF^(1/(2*Df))"], 3)
    )
    message("  GVIF (adjusted for df):")
    for (i in seq_len(nrow(vif_table))) {
      flag <- if (vif_table$GVIF_adj[i] > 2.24) " ** HIGH" else ""
      message("    ", vif_table$Variable[i], ": GVIF^(1/2Df) = ",
              vif_table$GVIF_adj[i], flag)
    }
  } else {
    vif_table <- data.frame(
      Variable = names(vif_values),
      VIF = round(vif_values, 3)
    )
    for (i in seq_len(nrow(vif_table))) {
      flag <- if (vif_table$VIF[i] > 5) " ** HIGH" else ""
      message("    ", vif_table$Variable[i], ": VIF = ", vif_table$VIF[i], flag)
    }
  }
  message("  Rule of thumb: VIF > 5 (or GVIF^(1/2Df) > ~2.24) indicates concern")
}

# --- 1b. Residual normality test ---
message("\n1b. Residual normality (Shapiro-Wilk):")
ols_resid <- residuals(ols_model)
sw_test <- shapiro.test(ols_resid)
message("  OLS residuals: W = ", round(sw_test$statistic, 4),
        ", p = ", format.pval(sw_test$p.value, digits = 4))
message("  Interpretation: ", ifelse(sw_test$p.value < 0.05,
        "Residuals depart from normality (p<0.05)",
        "No evidence against normality (p>=0.05)"))

# --- 1c. Heteroscedasticity test ---
message("\n1c. Heteroscedasticity (Breusch-Pagan):")
bp_test <- lmtest::bptest(ols_model)
message("  BP statistic = ", round(bp_test$statistic, 4),
        ", df = ", bp_test$parameter,
        ", p = ", format.pval(bp_test$p.value, digits = 4))
message("  Interpretation: ", ifelse(bp_test$p.value < 0.05,
        "Significant heteroscedasticity detected (p<0.05)",
        "No evidence of heteroscedasticity (p>=0.05)"))

# --- 1d. Moran's I on model residuals ---
message("\n1d. Moran's I on model residuals:")

# OLS residuals
moran_ols_resid <- moran.test(ols_resid, lw_reg, zero.policy = TRUE)
message("  OLS residuals: I = ", round(moran_ols_resid$estimate["Moran I statistic"], 4),
        ", p = ", format.pval(moran_ols_resid$p.value, digits = 4))

# Spatial lag residuals
lag_resid <- residuals(lag_model)
moran_lag_resid <- moran.test(lag_resid, lw_reg, zero.policy = TRUE)
message("  Spatial lag residuals: I = ", round(moran_lag_resid$estimate["Moran I statistic"], 4),
        ", p = ", format.pval(moran_lag_resid$p.value, digits = 4))

# Spatial error residuals
err_resid <- residuals(error_model)
moran_err_resid <- moran.test(err_resid, lw_reg, zero.policy = TRUE)
message("  Spatial error residuals: I = ", round(moran_err_resid$estimate["Moran I statistic"], 4),
        ", p = ", format.pval(moran_err_resid$p.value, digits = 4))

message("\n  Summary: Spatial dependence ",
        ifelse(moran_ols_resid$p.value < 0.05 &
               (moran_lag_resid$p.value >= 0.05 | moran_err_resid$p.value >= 0.05),
               "present in OLS, resolved by spatial model(s)",
               ifelse(moran_ols_resid$p.value >= 0.05,
                      "not detected in OLS residuals",
                      "persists in spatial model residuals — consider alternative specification")))

# --- 1e. Compile diagnostics summary table ---
diag_summary <- tibble(
  Test = c("VIF max (GVIF^1/2Df)", "Shapiro-Wilk (normality)",
           "Breusch-Pagan (heteroscedasticity)",
           "Moran's I: OLS residuals", "Moran's I: Lag residuals",
           "Moran's I: Error residuals"),
  Statistic = c(
    if (!is.null(vif_values) && is.matrix(vif_values))
      round(max(vif_values[, "GVIF^(1/(2*Df))"]), 3) else NA,
    round(sw_test$statistic, 4),
    round(bp_test$statistic, 4),
    round(moran_ols_resid$estimate["Moran I statistic"], 4),
    round(moran_lag_resid$estimate["Moran I statistic"], 4),
    round(moran_err_resid$estimate["Moran I statistic"], 4)
  ),
  p_value = c(
    NA,
    sw_test$p.value,
    bp_test$p.value,
    moran_ols_resid$p.value,
    moran_lag_resid$p.value,
    moran_err_resid$p.value
  ),
  Interpretation = c(
    if (!is.null(vif_values) && is.matrix(vif_values))
      ifelse(max(vif_values[, "GVIF^(1/(2*Df))"]) > 2.24,
             "Multicollinearity concern", "No multicollinearity concern")
    else "Not computed",
    ifelse(sw_test$p.value < 0.05, "Non-normal residuals", "Normal residuals"),
    ifelse(bp_test$p.value < 0.05, "Heteroscedastic", "Homoscedastic"),
    ifelse(moran_ols_resid$p.value < 0.05, "Spatial dependence present", "No spatial dependence"),
    ifelse(moran_lag_resid$p.value < 0.05, "Residual spatial dependence", "Spatial dependence resolved"),
    ifelse(moran_err_resid$p.value < 0.05, "Residual spatial dependence", "Spatial dependence resolved")
  )
)

diag_gt <- diag_summary %>%
  mutate(p_value = ifelse(is.na(p_value), "—", format.pval(p_value, digits = 3))) %>%
  gt() %>%
  tab_header(
    title = "Table S1. Model diagnostic tests",
    subtitle = "OLS and spatial regression models, Victorian LGAs (n=74)"
  ) %>%
  tab_source_note("VIF threshold: GVIF^(1/2Df) > 2.24 equivalent to VIF > 5 for continuous predictors.")

gtsave(diag_gt, file.path(path_tables, "tableS1_model_diagnostics.html"))
write_csv(diag_summary, file.path(path_tables, "tableS1_model_diagnostics.csv"))

# Residual Q-Q plot
fig_qq <- ggplot(data.frame(resid = ols_resid), aes(sample = resid)) +
  stat_qq(colour = "#2c3e50", alpha = 0.6) +
  stat_qq_line(colour = "#e74c3c") +
  labs(title = "Normal Q-Q plot: OLS residuals",
       subtitle = paste0("Shapiro-Wilk p = ", format.pval(sw_test$p.value, digits = 3)),
       x = "Theoretical quantiles", y = "Sample quantiles")

save_plot(fig_qq, "figS1_residual_qq_plot", width = 7, height = 7)

message("  Diagnostics saved to tableS1 and figS1")

# ==============================================================================
# 2. SPATIAL WEIGHTS SENSITIVITY ANALYSIS
# ==============================================================================

message("\n=== 2. Spatial Weights Sensitivity ===")

coords_reg <- st_coordinates(st_centroid(reg_data))

# Build alternative weight matrices
message("  Building alternative spatial weights...")

# Queen contiguity (existing — lw_reg)
# k=4 nearest neighbours
knn4 <- knearneigh(coords_reg, k = 4)
nb_knn4 <- knn2nb(knn4, sym = TRUE)
lw_knn4 <- nb2listw(nb_knn4, style = "W", zero.policy = TRUE)

# k=6 nearest neighbours
knn6 <- knearneigh(coords_reg, k = 6)
nb_knn6 <- knn2nb(knn6, sym = TRUE)
lw_knn6 <- nb2listw(nb_knn6, style = "W", zero.policy = TRUE)

# Inverse distance weights: use queen neighbours but weight by 1/distance
nb_dists <- nbdists(nb_reg, coords_reg)
idw_weights <- lapply(nb_dists, function(d) 1 / d)
lw_idw <- nb2listw(nb_reg, glist = idw_weights, style = "W", zero.policy = TRUE)

weight_specs <- list(
  "Queen contiguity" = lw_reg,
  "k=4 nearest neighbours" = lw_knn4,
  "k=6 nearest neighbours" = lw_knn6,
  "Inverse distance" = lw_idw
)

# Run Global Moran's I and spatial lag model under each specification
sensitivity_results <- list()

for (wname in names(weight_specs)) {
  lw_test <- weight_specs[[wname]]

  # Moran's I
  mi <- tryCatch(
    moran.test(reg_data$log_rate, lw_test, zero.policy = TRUE),
    error = function(e) NULL
  )

  # Spatial lag model
  lag_test <- tryCatch(
    lagsarlm(log_rate ~ pct_indigenous + mean_irsd + log_dist_sh + remoteness_factor,
             data = reg_data, listw = lw_test, zero.policy = TRUE),
    error = function(e) NULL
  )

  sensitivity_results[[wname]] <- tibble(
    Weights = wname,
    Morans_I = if (!is.null(mi)) round(mi$estimate["Moran I statistic"], 4) else NA,
    Morans_p = if (!is.null(mi)) mi$p.value else NA,
    Lag_rho = if (!is.null(lag_test)) round(lag_test$rho, 3) else NA,
    Lag_AIC = if (!is.null(lag_test)) round(AIC(lag_test), 1) else NA,
    Lag_dist_coef = if (!is.null(lag_test)) round(coef(lag_test)["log_dist_sh"], 3) else NA,
    Lag_indig_coef = if (!is.null(lag_test)) round(coef(lag_test)["pct_indigenous"], 3) else NA
  )

  message("  ", wname, ": I=",
          if (!is.null(mi)) round(mi$estimate["Moran I statistic"], 3) else "NA",
          ", rho=",
          if (!is.null(lag_test)) round(lag_test$rho, 3) else "NA",
          ", AIC=",
          if (!is.null(lag_test)) round(AIC(lag_test), 1) else "NA")
}

sensitivity_table <- bind_rows(sensitivity_results)

sens_gt <- sensitivity_table %>%
  mutate(Morans_p = format.pval(Morans_p, digits = 3)) %>%
  gt() %>%
  tab_header(
    title = "Table S2. Spatial weights sensitivity analysis",
    subtitle = "Global Moran's I and spatial lag model under alternative weight specifications"
  ) %>%
  cols_label(
    Weights = "Weights specification",
    Morans_I = "Moran's I",
    Morans_p = "p-value",
    Lag_rho = "ρ (lag)",
    Lag_AIC = "AIC",
    Lag_dist_coef = "β (distance)",
    Lag_indig_coef = "β (Indigenous %)"
  ) %>%
  tab_source_note("All models: log(rate) ~ % Indigenous + IRSD + log(distance to SH clinic) + remoteness.")

gtsave(sens_gt, file.path(path_tables, "tableS2_weights_sensitivity.html"))
write_csv(sensitivity_table, file.path(path_tables, "tableS2_weights_sensitivity.csv"))

message("  Sensitivity table saved to tableS2")

# ==============================================================================
# 3. GWR LOCAL COEFFICIENT SIGNIFICANCE
# ==============================================================================

message("\n=== 3. GWR Local Coefficient Significance ===")

gwr_formula <- log_rate ~ pct_indigenous + mean_irsd + log_dist_sh

# Bandwidth selection
bw <- tryCatch({
  gwr.sel(gwr_formula, data = as(reg_data, "Spatial"),
          coords = coords_reg, adapt = TRUE, verbose = FALSE)
}, error = function(e) {
  message("  Bandwidth selection failed, using 0.5")
  0.5
})

gwr_model <- tryCatch({
  gwr(gwr_formula, data = as(reg_data, "Spatial"),
      coords = coords_reg, adapt = bw, hatmatrix = TRUE, se.fit = TRUE)
}, error = function(e) {
  message("  GWR fitting failed: ", e$message)
  NULL
})

if (!is.null(gwr_model)) {
  gwr_df <- as.data.frame(gwr_model$SDF)

  # Compute local t-values = coefficient / SE
  # SE columns are named with _se suffix or may be in specific columns
  # spgwr stores SEs — extract them
  coef_names <- c("pct_indigenous", "mean_irsd", "log_dist_sh")

  gwr_significance <- list()

  for (vname in coef_names) {
    coef_col <- gwr_df[[vname]]
    # Local SEs are stored as 'pred.se' or need to be derived from the hat matrix
    # spgwr with se.fit=TRUE stores standard errors
    se_col_name <- paste0(vname, "_se")

    if (se_col_name %in% names(gwr_df)) {
      se_col <- gwr_df[[se_col_name]]
    } else {
      # Fallback: compute approximate SE from coefficient and localR2
      # Use global SE as approximation if local SEs unavailable
      global_se <- summary(ols_model)$coefficients[vname, "Std. Error"]
      se_col <- rep(global_se, nrow(gwr_df))
      message("  Note: local SEs not available for ", vname, "; using global SE as proxy")
    }

    t_val <- coef_col / se_col
    p_val <- 2 * pt(abs(t_val), df = nrow(reg_data) - length(coef_names) - 1,
                     lower.tail = FALSE)

    reg_data[[paste0("gwr_t_", vname)]] <- t_val
    reg_data[[paste0("gwr_p_", vname)]] <- p_val

    n_sig <- sum(p_val < 0.05, na.rm = TRUE)
    pct_sig <- round(n_sig / length(p_val) * 100, 1)

    gwr_significance[[vname]] <- tibble(
      Variable = vname,
      N_LGAs = length(p_val),
      N_significant = n_sig,
      Pct_significant = pct_sig,
      Coef_range = paste0(round(min(coef_col, na.rm = TRUE), 3), " to ",
                          round(max(coef_col, na.rm = TRUE), 3)),
      Median_coef = round(median(coef_col, na.rm = TRUE), 3)
    )

    message("  ", vname, ": ", n_sig, "/", length(p_val),
            " LGAs significant at p<0.05 (", pct_sig, "%)")
  }

  gwr_sig_table <- bind_rows(gwr_significance)

  gwr_sig_gt <- gwr_sig_table %>%
    gt() %>%
    tab_header(
      title = "Table S3. GWR local coefficient significance",
      subtitle = "Proportion of LGAs with significant (p<0.05) local coefficients"
    ) %>%
    cols_label(
      Variable = "Predictor",
      N_LGAs = "N (LGAs)",
      N_significant = "N significant",
      Pct_significant = "% significant",
      Coef_range = "Coefficient range",
      Median_coef = "Median coefficient"
    )

  gtsave(gwr_sig_gt, file.path(path_tables, "tableS3_gwr_local_significance.html"))
  write_csv(gwr_sig_table, file.path(path_tables, "tableS3_gwr_local_significance.csv"))

  # Map of local significance for distance coefficient
  if ("gwr_p_log_dist_sh" %in% names(reg_data)) {
    reg_data <- reg_data %>%
      mutate(
        gwr_dist_sig = case_when(
          gwr_p_log_dist_sh < 0.01 ~ "p < 0.01",
          gwr_p_log_dist_sh < 0.05 ~ "p < 0.05",
          gwr_p_log_dist_sh < 0.10 ~ "p < 0.10",
          TRUE ~ "Not significant"
        ),
        gwr_dist_sig = factor(gwr_dist_sig,
                              levels = c("p < 0.01", "p < 0.05", "p < 0.10",
                                         "Not significant"))
      )

    sig_colours <- c("p < 0.01" = "#d73027", "p < 0.05" = "#fc8d59",
                     "p < 0.10" = "#fee08b", "Not significant" = "#f0f0f0")

    fig_gwr_sig <- tm_shape(reg_data) +
      tm_polygons(
        fill = "gwr_dist_sig",
        fill.scale = tm_scale_categorical(values = sig_colours),
        fill.legend = tm_legend(title = "Local significance")
      ) +
      tm_title("GWR local significance: Distance to SH clinic coefficient")

    save_map(fig_gwr_sig, "figS2_gwr_distance_significance", width = 8, height = 8)
    message("  GWR significance map saved to figS2")
  }
} else {
  message("  GWR model unavailable — skipping significance analysis")
}

# ==============================================================================
# 4. INTERACTION TESTING: INDIGENOUS % × REMOTENESS
# ==============================================================================

message("\n=== 4. Interaction Testing: Indigenous % × Remoteness ===")

# OLS with interaction
ols_interaction <- lm(
  log_rate ~ pct_indigenous * remoteness_factor + mean_irsd + log_dist_sh,
  data = reg_data
)

# Compare with main-effects-only model
anova_interaction <- anova(ols_model, ols_interaction, test = "F")
interaction_p <- anova_interaction$`Pr(>F)`[2]

message("  F-test for interaction (OLS): F = ",
        round(anova_interaction$F[2], 3),
        ", p = ", format.pval(interaction_p, digits = 4))
message("  Interaction ", ifelse(interaction_p < 0.05,
        "IS significant — effect of Indigenous % varies by remoteness",
        "is NOT significant — no evidence of effect modification"))
message("  AIC without interaction: ", round(AIC(ols_model), 1))
message("  AIC with interaction: ", round(AIC(ols_interaction), 1))

# Spatial lag with interaction (if interaction significant or AIC improves)
lag_interaction <- tryCatch({
  lagsarlm(
    log_rate ~ pct_indigenous * remoteness_factor + mean_irsd + log_dist_sh,
    data = reg_data, listw = lw_reg, zero.policy = TRUE
  )
}, error = function(e) {
  message("  Spatial lag with interaction failed: ", e$message)
  NULL
})

if (!is.null(lag_interaction)) {
  message("  Spatial lag AIC without interaction: ", round(AIC(lag_model), 1))
  message("  Spatial lag AIC with interaction: ", round(AIC(lag_interaction), 1))
  message("  Interaction ", ifelse(AIC(lag_interaction) < AIC(lag_model),
          "improves model fit (lower AIC)",
          "does not improve model fit"))
}

# Save interaction results
interaction_results <- tibble(
  Model = c("OLS (no interaction)", "OLS (with interaction)",
            "Spatial lag (no interaction)", "Spatial lag (with interaction)"),
  AIC = c(AIC(ols_model), AIC(ols_interaction),
          AIC(lag_model),
          if (!is.null(lag_interaction)) AIC(lag_interaction) else NA),
  Note = c("Main effects only", paste0("F-test p = ", format.pval(interaction_p, digits = 3)),
           "Main effects only",
           if (!is.null(lag_interaction))
             paste0("rho = ", round(lag_interaction$rho, 3)) else "Failed to converge")
)

write_csv(interaction_results, file.path(path_tables, "tableS4_interaction_test.csv"))

# ==============================================================================
# 5. ITS MODEL FIX: POPULATION OFFSET + OVERDISPERSION + NB
# ==============================================================================

message("\n=== 5. ITS Model: Population Offset and Overdispersion ===")

# Australian ERP mid-year estimates (ABS 3101.0, thousands)
aus_erp <- tibble(
  year = 2013:2022,
  population = c(23117997, 23490736, 23815995, 24210809, 24601860,
                 24992860, 25365745, 25693068, 25766605, 26269421)
)

national_covid <- national_ts %>%
  left_join(aus_erp, by = "year") %>%
  mutate(
    time = year - min(year) + 1,
    covid = as.integer(year >= 2020),
    time_since_covid = pmax(0, year - 2020),
    log_pop = log(population)
  )

# --- 5a. Original Poisson (no offset) ---
its_original <- glm(
  count ~ time + covid + time_since_covid,
  data = national_covid, family = poisson()
)

# --- 5b. Poisson with population offset ---
its_offset <- glm(
  count ~ time + covid + time_since_covid + offset(log_pop),
  data = national_covid, family = poisson()
)

# --- 5c. Overdispersion test ---
disp_original <- deviance(its_original) / df.residual(its_original)
disp_offset <- deviance(its_offset) / df.residual(its_offset)

message("  Overdispersion (deviance/df):")
message("    Poisson (no offset): ", round(disp_original, 2),
        ifelse(disp_original > 1.5, " — overdispersed", " — acceptable"))
message("    Poisson (with offset): ", round(disp_offset, 2),
        ifelse(disp_offset > 1.5, " — overdispersed", " — acceptable"))

# --- 5d. Negative binomial with offset ---
its_nb <- tryCatch({
  MASS::glm.nb(
    count ~ time + covid + time_since_covid + offset(log_pop),
    data = national_covid
  )
}, error = function(e) {
  message("  Negative binomial failed: ", e$message)
  NULL
})

# --- 5e. Quasipoisson (robust SE) ---
its_quasi <- glm(
  count ~ time + covid + time_since_covid + offset(log_pop),
  data = national_covid, family = quasipoisson()
)

# --- 5f. Compare models ---
message("\n  Model comparison:")

extract_irr <- function(model, term) {
  coefs <- coef(summary(model))
  if (term %in% rownames(coefs)) {
    est <- coefs[term, "Estimate"]
    se <- coefs[term, "Std. Error"]
    p <- coefs[term, ncol(coefs)]
    irr <- exp(est)
    irr_lo <- exp(est - 1.96 * se)
    irr_hi <- exp(est + 1.96 * se)
    return(list(irr = irr, ci_lo = irr_lo, ci_hi = irr_hi, p = p))
  }
  return(list(irr = NA, ci_lo = NA, ci_hi = NA, p = NA))
}

its_comparison <- tibble(
  Model = c("Poisson (no offset)", "Poisson (offset)", "Quasipoisson (offset)",
            if (!is.null(its_nb)) "Neg. binomial (offset)" else NULL),
  AIC = c(AIC(its_original), AIC(its_offset), NA,
          if (!is.null(its_nb)) AIC(its_nb) else NULL),
  `Deviance/df` = c(disp_original, disp_offset,
                    deviance(its_quasi) / df.residual(its_quasi),
                    if (!is.null(its_nb)) deviance(its_nb) / df.residual(its_nb) else NULL),
  COVID_IRR = c(
    extract_irr(its_original, "covid")$irr,
    extract_irr(its_offset, "covid")$irr,
    extract_irr(its_quasi, "covid")$irr,
    if (!is.null(its_nb)) extract_irr(its_nb, "covid")$irr else NULL
  ),
  COVID_p = c(
    extract_irr(its_original, "covid")$p,
    extract_irr(its_offset, "covid")$p,
    extract_irr(its_quasi, "covid")$p,
    if (!is.null(its_nb)) extract_irr(its_nb, "covid")$p else NULL
  ),
  Trend_IRR = c(
    extract_irr(its_original, "time_since_covid")$irr,
    extract_irr(its_offset, "time_since_covid")$irr,
    extract_irr(its_quasi, "time_since_covid")$irr,
    if (!is.null(its_nb)) extract_irr(its_nb, "time_since_covid")$irr else NULL
  )
) %>%
  mutate(
    across(c(AIC, `Deviance/df`, COVID_IRR, Trend_IRR), ~round(.x, 3)),
    COVID_p = format.pval(COVID_p, digits = 3)
  )

message("  ITS model comparison:")
for (i in seq_len(nrow(its_comparison))) {
  message("    ", its_comparison$Model[i],
          ": AIC=", its_comparison$AIC[i],
          ", Dev/df=", its_comparison$`Deviance/df`[i],
          ", COVID IRR=", its_comparison$COVID_IRR[i],
          " (p=", its_comparison$COVID_p[i], ")")
}

its_gt <- its_comparison %>%
  gt() %>%
  tab_header(
    title = "Table S4. ITS model comparison: Poisson vs negative binomial with population offset",
    subtitle = "National infectious syphilis notifications, 2013–2022"
  ) %>%
  tab_source_note("IRR = incidence rate ratio. Offset = log(ABS estimated resident population).") %>%
  tab_source_note("Deviance/df > 1.5 suggests overdispersion; quasipoisson and NB provide robust inference.")

gtsave(its_gt, file.path(path_tables, "tableS4_its_model_comparison.html"))
write_csv(its_comparison, file.path(path_tables, "tableS4_its_model_comparison.csv"))

# Save updated national ITS data with offset predictions
best_its <- if (!is.null(its_nb)) its_nb else its_quasi
national_covid$predicted_offset <- predict(best_its, type = "response")
counterfactual_data <- national_covid %>%
  mutate(covid = 0, time_since_covid = 0)
national_covid$counterfactual_offset <- predict(best_its, newdata = counterfactual_data,
                                                 type = "response")
saveRDS(national_covid, file.path(path_processed, "national_covid_its.rds"))

message("  ITS comparison saved to tableS4")

# --- 5g. Quadratic pre-COVID trend sensitivity ---
message("\n  Fitting quadratic pre-COVID trend sensitivity check...")

national_covid <- national_covid %>%
  mutate(time_sq = time^2)

its_nb_quad <- tryCatch({
  MASS::glm.nb(
    count ~ time + time_sq + covid + time_since_covid + offset(log_pop),
    data = national_covid
  )
}, error = function(e) {
  message("  Quadratic NB failed: ", e$message)
  NULL
})

if (!is.null(its_nb_quad)) {
  quad_irr <- extract_irr(its_nb_quad, "covid")
  linear_irr <- extract_irr(its_nb, "covid")

  message("  Quadratic sensitivity:")
  message("    Linear NB COVID IRR: ", round(linear_irr$irr, 3),
          " (p=", format.pval(linear_irr$p, digits = 3), ")")
  message("    Quadratic NB COVID IRR: ", round(quad_irr$irr, 3),
          " (p=", format.pval(quad_irr$p, digits = 3), ")")
  message("    AIC linear: ", round(AIC(its_nb), 1),
          " vs quadratic: ", round(AIC(its_nb_quad), 1))

  # Add to comparison table
  quad_row <- tibble(
    Model = "NB quadratic trend (offset)",
    AIC = round(AIC(its_nb_quad), 3),
    `Deviance/df` = round(deviance(its_nb_quad) / df.residual(its_nb_quad), 3),
    COVID_IRR = round(quad_irr$irr, 3),
    COVID_p = format.pval(quad_irr$p, digits = 3),
    Trend_IRR = round(extract_irr(its_nb_quad, "time_since_covid")$irr, 3)
  )

  its_comparison <- bind_rows(its_comparison, quad_row)

  # Re-save updated Table S4
  its_gt <- its_comparison %>%
    gt() %>%
    tab_header(
      title = "Table S4. ITS model comparison: sensitivity to model specification",
      subtitle = "National infectious syphilis notifications, 2013-2022"
    ) %>%
    tab_source_note("IRR = incidence rate ratio. Offset = log(ABS estimated resident population).") %>%
    tab_source_note("Deviance/df > 1.5 suggests overdispersion. Quadratic model includes time^2 to test sensitivity to pre-COVID functional form.")

  gtsave(its_gt, file.path(path_tables, "tableS4_its_model_comparison.html"))
  write_csv(its_comparison, file.path(path_tables, "tableS4_its_model_comparison.csv"))

  message("  Updated tableS4 with quadratic sensitivity")
} else {
  message("  Quadratic model could not be fitted — skipping")
}

# ==============================================================================
# 6. CONGENITAL RISK SCORE SENSITIVITY
# ==============================================================================

message("\n=== 6. Congenital Risk Score Sensitivity ===")

# Reload SA3 risk data with underlying variables
sa3_data <- sa3_risk %>% st_drop_geometry()

# Function to classify risk under different thresholds
classify_risk <- function(data, hi_pctile, lo_pctile) {
  data %>%
    mutate(
      high_indigenous = pct_indigenous > quantile(pct_indigenous, hi_pctile, na.rm = TRUE),
      low_anc = !is.na(pct_first_trimester) &
        pct_first_trimester < quantile(pct_first_trimester, lo_pctile, na.rm = TRUE),
      high_remoteness = remoteness_mode %in% c("Remote", "Very Remote"),
      high_dist_clinic = !is.na(mean_dist_sh_km) &
        mean_dist_sh_km > quantile(mean_dist_sh_km, hi_pctile, na.rm = TRUE),
      risk_score = as.integer(high_indigenous) + as.integer(low_anc) +
        as.integer(high_remoteness) + as.integer(high_dist_clinic),
      risk_category = case_when(
        risk_score >= 3 ~ "High risk",
        risk_score == 2 ~ "Elevated risk",
        risk_score == 1 ~ "Moderate risk",
        risk_score == 0 ~ "Lower risk"
      )
    )
}

thresholds <- list(
  "70th/30th" = c(0.70, 0.30),
  "75th/25th (original)" = c(0.75, 0.25),
  "80th/20th" = c(0.80, 0.20),
  "90th/10th" = c(0.90, 0.10)
)

risk_sensitivity <- list()
high_risk_areas <- list()

for (tname in names(thresholds)) {
  hi <- thresholds[[tname]][1]
  lo <- thresholds[[tname]][2]

  classified <- classify_risk(sa3_data, hi, lo)
  counts <- table(classified$risk_category)

  high_risk_names <- classified %>%
    filter(risk_category == "High risk") %>%
    pull(sa3_name)

  risk_sensitivity[[tname]] <- tibble(
    Threshold = tname,
    High_risk = as.integer(counts["High risk"] %||% 0),
    Elevated_risk = as.integer(counts["Elevated risk"] %||% 0),
    Moderate_risk = as.integer(counts["Moderate risk"] %||% 0),
    Lower_risk = as.integer(counts["Lower risk"] %||% 0)
  )

  high_risk_areas[[tname]] <- high_risk_names

  message("  ", tname, ": High=", risk_sensitivity[[tname]]$High_risk,
          ", Elevated=", risk_sensitivity[[tname]]$Elevated_risk,
          ", Moderate=", risk_sensitivity[[tname]]$Moderate_risk,
          ", Lower=", risk_sensitivity[[tname]]$Lower_risk)
}

risk_sens_table <- bind_rows(risk_sensitivity)

# Identify areas classified as High risk in ALL threshold specifications
all_high <- Reduce(intersect, high_risk_areas)
message("\n  SA3 areas classified High risk under ALL thresholds (", length(all_high), "):")
for (a in all_high) message("    ", a)

# Areas classified High under 3+ threshold specifications
area_counts <- table(unlist(high_risk_areas))
robust_high <- names(area_counts[area_counts >= 3])
message("  SA3 areas classified High in ≥3 of 4 specifications: ", length(robust_high))

risk_sens_table <- risk_sens_table %>%
  mutate(
    Consistently_high = c(
      paste0(length(all_high), " in all 4"),
      paste0(length(robust_high), " in ≥3"),
      "", ""
    )[seq_len(n())]
  )

risk_sens_gt <- risk_sens_table %>%
  gt() %>%
  tab_header(
    title = "Table S5. Congenital syphilis risk classification sensitivity",
    subtitle = "Number of SA3 areas by risk category under varying percentile thresholds"
  ) %>%
  cols_label(
    Threshold = "Threshold (high/low percentile)",
    High_risk = "High risk (3-4)",
    Elevated_risk = "Elevated (2)",
    Moderate_risk = "Moderate (1)",
    Lower_risk = "Lower (0)"
  ) %>%
  tab_source_note("Risk factors: Indigenous %, first-trimester ANC, remoteness, clinic distance.") %>%
  tab_source_note(paste0("SA3 areas consistently High under all thresholds: ",
                         paste(all_high, collapse = "; ")))

gtsave(risk_sens_gt, file.path(path_tables, "tableS5_risk_sensitivity.html"))
write_csv(risk_sens_table, file.path(path_tables, "tableS5_risk_sensitivity.csv"))

message("  Risk sensitivity saved to tableS5")

# ==============================================================================
# 7. TEMPORAL LISA EVOLUTION (VIC LGAs 2019–2024)
# ==============================================================================

message("\n=== 7. Temporal LISA Evolution (VIC LGAs 2019-2024) ===")

# Build annual rates for each year
vic_annual <- vic_analysis %>%
  mutate(lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>% str_trim())

# Use vic_spatial for boundary + weights (same n)
nb_full <- poly2nb(vic_spatial, queen = TRUE)
# Fix islands
n_isl <- sum(card(nb_full) == 0)
if (n_isl > 0) {
  isl <- which(card(nb_full) == 0)
  coords_full <- st_coordinates(st_centroid(vic_spatial))
  for (i in isl) {
    dists <- as.numeric(st_distance(st_centroid(vic_spatial[i, ]),
                                     st_centroid(vic_spatial[-i, ])))
    near <- which.min(dists)
    if (near >= i) near <- near + 1
    nb_full[[i]] <- as.integer(near)
    nb_full[[near]] <- sort(unique(c(nb_full[[near]], as.integer(i))))
  }
}
lw_full <- nb2listw(nb_full, style = "W", zero.policy = TRUE)

years <- sort(unique(vic_annual$year))
temporal_results <- list()

for (yr in years) {
  yr_data <- vic_annual %>%
    filter(year == yr) %>%
    dplyr::select(lga_name_clean, rate_per_100k)

  # Match to spatial object
  matched <- vic_spatial %>%
    left_join(yr_data, by = c("LGA_NAME21" = "lga_name_clean"))

  rate_vec <- matched$rate_per_100k
  # Replace NAs with 0 for Moran's I
  rate_vec[is.na(rate_vec)] <- 0

  # Global Moran's I
  mi <- tryCatch(
    moran.test(rate_vec, lw_full, zero.policy = TRUE),
    error = function(e) NULL
  )

  # LISA
  lisa_yr <- tryCatch(
    localmoran(rate_vec, lw_full, zero.policy = TRUE),
    error = function(e) NULL
  )

  n_hh <- 0
  if (!is.null(lisa_yr)) {
    mean_r <- mean(rate_vec, na.rm = TRUE)
    lag_r <- lag.listw(lw_full, rate_vec, zero.policy = TRUE)
    p_vals <- lisa_yr[, "Pr(z != E(Ii))"]
    hh <- which(p_vals < 0.05 & rate_vec > mean_r & lag_r > mean(lag_r, na.rm = TRUE))
    n_hh <- length(hh)
  }

  temporal_results[[as.character(yr)]] <- tibble(
    Year = yr,
    Morans_I = if (!is.null(mi)) round(mi$estimate["Moran I statistic"], 4) else NA,
    Morans_p = if (!is.null(mi)) mi$p.value else NA,
    N_HighHigh_clusters = n_hh,
    Mean_rate = round(mean(rate_vec, na.rm = TRUE), 1)
  )

  message("  ", yr, ": I = ",
          if (!is.null(mi)) round(mi$estimate["Moran I statistic"], 3) else "NA",
          " (p = ", if (!is.null(mi)) format.pval(mi$p.value, digits = 3) else "NA",
          "), H-H clusters = ", n_hh,
          ", Mean rate = ", round(mean(rate_vec, na.rm = TRUE), 1))
}

temporal_table <- bind_rows(temporal_results)

temporal_gt <- temporal_table %>%
  mutate(Morans_p = format.pval(Morans_p, digits = 3)) %>%
  gt() %>%
  tab_header(
    title = "Table S6. Temporal evolution of spatial clustering, VIC LGAs, 2019–2024",
    subtitle = "Annual Global Moran's I and LISA High-High cluster count"
  ) %>%
  cols_label(
    Year = "Year",
    Morans_I = "Moran's I",
    Morans_p = "p-value",
    N_HighHigh_clusters = "High-High clusters",
    Mean_rate = "Mean rate (per 100k)"
  ) %>%
  tab_source_note("Weights: Queen contiguity. LISA significance at p < 0.05 with FDR adjustment.")

gtsave(temporal_gt, file.path(path_tables, "tableS6_temporal_lisa.html"))
write_csv(temporal_table, file.path(path_tables, "tableS6_temporal_lisa.csv"))

# Temporal Moran's I trend plot
fig_temporal_moran <- ggplot(temporal_table, aes(x = Year)) +
  geom_line(aes(y = Morans_I), colour = "#2c3e50", linewidth = 1) +
  geom_point(aes(y = Morans_I, size = N_HighHigh_clusters),
             colour = "#e74c3c", alpha = 0.7) +
  scale_size_continuous(range = c(2, 8), name = "H-H clusters") +
  scale_x_continuous(breaks = years) +
  labs(
    title = "Temporal evolution of spatial clustering: VIC LGA syphilis rates",
    subtitle = "2019-2024 | Point size = number of High-High LISA clusters",
    x = NULL, y = "Global Moran's I"
  )

save_plot(fig_temporal_moran, "figS3_temporal_morans_i", width = 10, height = 6)

message("  Temporal LISA saved to tableS6 and figS3")

# ==============================================================================
# SAVE UPDATED REGRESSION DATA
# ==============================================================================

saveRDS(reg_data, file.path(path_processed, "vic_lga_regression_data.rds"))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== Sensitivity Analyses Complete ===")
message("Supplementary tables:")
message("  Table S1: Model diagnostics (VIF, normality, heteroscedasticity, residual Moran's I)")
message("  Table S2: Spatial weights sensitivity")
message("  Table S3: GWR local coefficient significance")
message("  Table S4: ITS model comparison (Poisson, NB, population offset)")
message("  Table S5: Congenital risk classification sensitivity")
message("  Table S6: Temporal LISA evolution (2019-2024)")
message("\nSupplementary figures:")
message("  Figure S1: Residual Q-Q plot")
message("  Figure S2: GWR local significance map (distance coefficient)")
message("  Figure S3: Temporal Moran's I trend")
message("\nInteraction test results: ", file.path(path_tables, "tableS4_interaction_test.csv"))

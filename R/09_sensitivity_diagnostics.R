# ==============================================================================
# 09_sensitivity_diagnostics.R — Sensitivity analyses and model diagnostics
# MSM-focused rewrite: updated covariates (pct_male_20_44, log_pop_density)
# ==============================================================================
#
# Sections:
#  1. OLS/spatial regression diagnostics (VIF, residuals, Moran's I on resids)
#  2. Spatial weights sensitivity (queen, k-NN, inverse-distance)
#  3. GWR local coefficient significance (t-values, p-values)
#  4. ITS model sensitivity (population offset, overdispersion, NB)
#  5. Moran's I by year (2019-2024) — spatial autocorrelation persistence
#  6. Temporal LISA evolution + persistence metric
#  7. Victorian temporal trends by cluster status (hot spot vs other)
#  8. Leave-one-out influence analysis (inner Melbourne LGAs)
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

message("=== Sensitivity Analyses and Model Diagnostics (MSM rewrite) ===")
message("Loading data...")

vic_spatial   <- readRDS(file.path(path_processed, "vic_lga_spatial_results.rds"))
reg_data      <- readRDS(file.path(path_processed, "vic_lga_regression_data.rds"))
national_ts   <- readRDS(file.path(path_processed, "national_time_series.rds"))
vic_analysis  <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))

# Ensure regression variables exist
if (!"log_rate" %in% names(reg_data)) {
  reg_data <- reg_data %>%
    mutate(
      log_rate = log1p(mean_annual_rate),
      log_dist_sh = log1p(mean_dist_sh_km),
      log_pop_density = log1p(pop_density),
      remoteness_factor = factor(remoteness_mode,
                                 levels = c("Major Cities", "Inner Regional",
                                            "Outer Regional", "Remote", "Very Remote"))
    )
}

# New regression formula (MSM rewrite)
reg_formula <- log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor

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

# Refit core models with new covariates
ols_model <- lm(reg_formula, data = reg_data)
lag_model <- lagsarlm(reg_formula, data = reg_data, listw = lw_reg, zero.policy = TRUE)
error_model <- errorsarlm(reg_formula, data = reg_data, listw = lw_reg, zero.policy = TRUE)

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

moran_ols_resid <- moran.test(ols_resid, lw_reg, zero.policy = TRUE)
message("  OLS residuals: I = ", round(moran_ols_resid$estimate["Moran I statistic"], 4),
        ", p = ", format.pval(moran_ols_resid$p.value, digits = 4))

lag_resid <- residuals(lag_model)
moran_lag_resid <- moran.test(lag_resid, lw_reg, zero.policy = TRUE)
message("  Spatial lag residuals: I = ", round(moran_lag_resid$estimate["Moran I statistic"], 4),
        ", p = ", format.pval(moran_lag_resid$p.value, digits = 4))

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
  mutate(p_value = ifelse(is.na(p_value), "\u2014", format.pval(p_value, digits = 3))) %>%
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

save_plot(fig_qq, "figS_residual_qq_plot", width = 7, height = 7)

message("  Diagnostics saved to tableS1 and figS_residual_qq_plot")

# ==============================================================================
# 2. SPATIAL WEIGHTS SENSITIVITY ANALYSIS
# ==============================================================================

message("\n=== 2. Spatial Weights Sensitivity ===")

coords_reg <- st_coordinates(st_centroid(reg_data))

message("  Building alternative spatial weights...")

# k=4 nearest neighbours
knn4 <- knearneigh(coords_reg, k = 4)
nb_knn4 <- knn2nb(knn4, sym = TRUE)
lw_knn4 <- nb2listw(nb_knn4, style = "W", zero.policy = TRUE)

# k=6 nearest neighbours
knn6 <- knearneigh(coords_reg, k = 6)
nb_knn6 <- knn2nb(knn6, sym = TRUE)
lw_knn6 <- nb2listw(nb_knn6, style = "W", zero.policy = TRUE)

# Rook contiguity
nb_rook <- poly2nb(reg_data, queen = FALSE)
n_rook_islands <- sum(card(nb_rook) == 0)
if (n_rook_islands > 0) {
  rook_islands <- which(card(nb_rook) == 0)
  for (i in rook_islands) {
    dists <- as.numeric(st_distance(
      st_centroid(reg_data[i, ]),
      st_centroid(reg_data[-i, ])
    ))
    nearest <- which.min(dists)
    if (nearest >= i) nearest <- nearest + 1
    nb_rook[[i]] <- as.integer(nearest)
    nb_rook[[nearest]] <- sort(unique(c(nb_rook[[nearest]], as.integer(i))))
  }
}
lw_rook <- nb2listw(nb_rook, style = "W", zero.policy = TRUE)

# Inverse distance weights
nb_dists <- nbdists(nb_reg, coords_reg)
idw_weights <- lapply(nb_dists, function(d) 1 / d)
lw_idw <- nb2listw(nb_reg, glist = idw_weights, style = "W", zero.policy = TRUE)

weight_specs <- list(
  "Queen contiguity" = lw_reg,
  "Rook contiguity" = lw_rook,
  "k=4 nearest neighbours" = lw_knn4,
  "k=6 nearest neighbours" = lw_knn6,
  "Inverse distance" = lw_idw
)

sensitivity_results <- list()

for (wname in names(weight_specs)) {
  lw_test <- weight_specs[[wname]]

  mi <- tryCatch(
    moran.test(reg_data$log_rate, lw_test, zero.policy = TRUE),
    error = function(e) NULL
  )

  lag_test <- tryCatch(
    lagsarlm(reg_formula, data = reg_data, listw = lw_test, zero.policy = TRUE),
    error = function(e) NULL
  )

  sensitivity_results[[wname]] <- tibble(
    Weights = wname,
    Morans_I = if (!is.null(mi)) round(mi$estimate["Moran I statistic"], 4) else NA,
    Morans_p = if (!is.null(mi)) mi$p.value else NA,
    Lag_rho = if (!is.null(lag_test)) round(lag_test$rho, 3) else NA,
    Lag_AIC = if (!is.null(lag_test)) round(AIC(lag_test), 1) else NA,
    Lag_male2044_coef = if (!is.null(lag_test)) round(coef(lag_test)["pct_male_20_44"], 3) else NA,
    Lag_dist_coef = if (!is.null(lag_test)) round(coef(lag_test)["log_dist_sh"], 3) else NA
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
    Lag_rho = "\u03c1 (lag)",
    Lag_AIC = "AIC",
    Lag_male2044_coef = "\u03b2 (% Males 20-44)",
    Lag_dist_coef = "\u03b2 (distance)"
  ) %>%
  tab_source_note("All models: log(rate) ~ % Males 20-44 + log(pop density) + IRSD + log(distance to SH clinic) + remoteness.")

gtsave(sens_gt, file.path(path_tables, "tableS2_weights_sensitivity.html"))
write_csv(sensitivity_table, file.path(path_tables, "tableS2_weights_sensitivity.csv"))

message("  Sensitivity table saved to tableS2")

# ==============================================================================
# 3. GWR LOCAL COEFFICIENT SIGNIFICANCE
# ==============================================================================

message("\n=== 3. GWR Local Coefficient Significance ===")

gwr_formula <- log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh

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

  coef_names <- c("pct_male_20_44", "log_pop_density", "mean_irsd", "log_dist_sh")

  gwr_significance <- list()

  for (vname in coef_names) {
    coef_col <- gwr_df[[vname]]
    se_col_name <- paste0(vname, "_se")

    if (se_col_name %in% names(gwr_df)) {
      se_col <- gwr_df[[se_col_name]]
    } else {
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

    save_map(fig_gwr_sig, "figS_gwr_distance_significance", width = 8, height = 8)
    message("  GWR significance map saved")
  }
} else {
  message("  GWR model unavailable — skipping significance analysis")
}

# ==============================================================================
# 4. ITS MODEL SENSITIVITY (POPULATION OFFSET + OVERDISPERSION + NB)
# ==============================================================================

message("\n=== 4. ITS Model: Population Offset and Overdispersion ===")

# Australian ERP mid-year estimates (ABS 3101.0)
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

its_original <- glm(count ~ time + covid + time_since_covid,
                     data = national_covid, family = poisson())
its_offset <- glm(count ~ time + covid + time_since_covid + offset(log_pop),
                   data = national_covid, family = poisson())

disp_original <- deviance(its_original) / df.residual(its_original)
disp_offset <- deviance(its_offset) / df.residual(its_offset)

message("  Overdispersion (deviance/df):")
message("    Poisson (no offset): ", round(disp_original, 2))
message("    Poisson (with offset): ", round(disp_offset, 2))

its_nb <- tryCatch({
  MASS::glm.nb(count ~ time + covid + time_since_covid + offset(log_pop),
               data = national_covid)
}, error = function(e) { message("  NB failed: ", e$message); NULL })

its_quasi <- glm(count ~ time + covid + time_since_covid + offset(log_pop),
                  data = national_covid, family = quasipoisson())

extract_irr <- function(model, term) {
  coefs <- coef(summary(model))
  if (term %in% rownames(coefs)) {
    est <- coefs[term, "Estimate"]
    se <- coefs[term, "Std. Error"]
    p <- coefs[term, ncol(coefs)]
    return(list(irr = exp(est), ci_lo = exp(est - 1.96 * se),
                ci_hi = exp(est + 1.96 * se), p = p))
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
  COVID_IRR = c(extract_irr(its_original, "covid")$irr,
                extract_irr(its_offset, "covid")$irr,
                extract_irr(its_quasi, "covid")$irr,
                if (!is.null(its_nb)) extract_irr(its_nb, "covid")$irr else NULL),
  COVID_p = c(extract_irr(its_original, "covid")$p,
              extract_irr(its_offset, "covid")$p,
              extract_irr(its_quasi, "covid")$p,
              if (!is.null(its_nb)) extract_irr(its_nb, "covid")$p else NULL),
  Trend_IRR = c(extract_irr(its_original, "time_since_covid")$irr,
                extract_irr(its_offset, "time_since_covid")$irr,
                extract_irr(its_quasi, "time_since_covid")$irr,
                if (!is.null(its_nb)) extract_irr(its_nb, "time_since_covid")$irr else NULL)
) %>%
  mutate(across(c(AIC, `Deviance/df`, COVID_IRR, Trend_IRR), ~round(.x, 3)),
         COVID_p = format.pval(COVID_p, digits = 3))

its_gt <- its_comparison %>%
  gt() %>%
  tab_header(
    title = "Table S4. ITS model comparison: sensitivity to model specification",
    subtitle = "National infectious syphilis notifications, 2013-2022"
  ) %>%
  tab_source_note("IRR = incidence rate ratio. Offset = log(ABS estimated resident population).") %>%
  tab_source_note("Deviance/df > 1.5 suggests overdispersion.")

gtsave(its_gt, file.path(path_tables, "tableS4_its_model_comparison.html"))
write_csv(its_comparison, file.path(path_tables, "tableS4_its_model_comparison.csv"))

# Save ITS data for manuscript figures
best_its <- if (!is.null(its_nb)) its_nb else its_quasi
national_covid$predicted_offset <- predict(best_its, type = "response")
cf_data <- national_covid %>% mutate(covid = 0, time_since_covid = 0)
national_covid$counterfactual_offset <- predict(best_its, newdata = cf_data, type = "response")
saveRDS(national_covid, file.path(path_processed, "national_covid_its.rds"))

message("  ITS comparison saved to tableS4")

# ==============================================================================
# 5. MORAN'S I BY YEAR (2019-2024) — Spatial autocorrelation persistence
# ==============================================================================

message("\n=== 5. Moran's I by Year (2019-2024) ===")

# Build spatial weights on full VIC spatial object
nb_full <- poly2nb(vic_spatial, queen = TRUE)
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

vic_annual <- vic_analysis %>%
  mutate(lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>% str_trim())

years <- sort(unique(vic_annual$year))
temporal_results <- list()
annual_lisa_lists <- list()  # Store H-H LGA names for persistence analysis

for (yr in years) {
  yr_data <- vic_annual %>%
    filter(year == yr) %>%
    dplyr::select(lga_name_clean, rate_per_100k)

  matched <- vic_spatial %>%
    left_join(yr_data, by = c("LGA_NAME21" = "lga_name_clean"))

  rate_vec <- matched$rate_per_100k
  rate_vec[is.na(rate_vec)] <- 0

  mi <- tryCatch(
    moran.test(rate_vec, lw_full, zero.policy = TRUE),
    error = function(e) NULL
  )

  lisa_yr <- tryCatch(
    localmoran(rate_vec, lw_full, zero.policy = TRUE),
    error = function(e) NULL
  )

  hh_names <- character(0)
  n_hh <- 0
  if (!is.null(lisa_yr)) {
    mean_r <- mean(rate_vec, na.rm = TRUE)
    lag_r <- lag.listw(lw_full, rate_vec, zero.policy = TRUE)
    p_vals <- lisa_yr[, "Pr(z != E(Ii))"]
    hh_idx <- which(p_vals < 0.05 & rate_vec > mean_r & lag_r > mean(lag_r, na.rm = TRUE))
    n_hh <- length(hh_idx)
    hh_names <- matched$LGA_NAME21[hh_idx]
  }

  annual_lisa_lists[[as.character(yr)]] <- hh_names

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
    title = "Table S5. Temporal evolution of spatial clustering, VIC LGAs, 2019-2024",
    subtitle = "Annual Global Moran's I and LISA High-High cluster count"
  ) %>%
  cols_label(
    Year = "Year",
    Morans_I = "Moran's I",
    Morans_p = "p-value",
    N_HighHigh_clusters = "High-High clusters",
    Mean_rate = "Mean rate (per 100k)"
  ) %>%
  tab_source_note("Weights: Queen contiguity. LISA significance at p < 0.05.")

gtsave(temporal_gt, file.path(path_tables, "tableS5_temporal_lisa.html"))
write_csv(temporal_table, file.path(path_tables, "tableS5_temporal_lisa.csv"))

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

save_plot(fig_temporal_moran, "figS_temporal_morans_i", width = 10, height = 6)

message("  Temporal Moran's I saved to tableS5 and figS_temporal_morans_i")

# ==============================================================================
# 6. TEMPORAL LISA PERSISTENCE METRIC
# ==============================================================================

message("\n=== 6. LISA Cluster Persistence (2019-2024) ===")

# Which LGAs appeared in H-H clusters in how many years?
all_hh_lgas <- unique(unlist(annual_lisa_lists))
persistence <- tibble(LGA = all_hh_lgas) %>%
  mutate(
    years_in_cluster = sapply(LGA, function(lga) {
      sum(sapply(annual_lisa_lists, function(yr_list) lga %in% yr_list))
    }),
    pct_years = round(years_in_cluster / length(years) * 100, 1)
  ) %>%
  arrange(desc(years_in_cluster), LGA)

message("  LGAs ever in H-H cluster: ", nrow(persistence))
message("  LGAs in H-H cluster ALL ", length(years), " years: ",
        sum(persistence$years_in_cluster == length(years)))
message("  LGAs in H-H cluster >=4 of ", length(years), " years: ",
        sum(persistence$years_in_cluster >= 4))

message("\n  Persistence detail:")
for (i in seq_len(nrow(persistence))) {
  message("    ", persistence$LGA[i], ": ",
          persistence$years_in_cluster[i], "/", length(years),
          " years (", persistence$pct_years[i], "%)")
}

persistence_gt <- persistence %>%
  gt() %>%
  tab_header(
    title = "Table S6. LISA cluster persistence: LGAs in High-High clusters, 2019-2024",
    subtitle = paste0("Number of years (out of ", length(years),
                      ") each LGA appeared in a significant High-High LISA cluster")
  ) %>%
  cols_label(
    LGA = "LGA name",
    years_in_cluster = "Years in H-H cluster",
    pct_years = "% of study period"
  ) %>%
  tab_source_note("Queen contiguity weights, LISA significance at p < 0.05.")

gtsave(persistence_gt, file.path(path_tables, "tableS6_lisa_persistence.html"))
write_csv(persistence, file.path(path_tables, "tableS6_lisa_persistence.csv"))

message("  Persistence table saved to tableS6")

# ==============================================================================
# 7. VICTORIAN TEMPORAL TRENDS BY CLUSTER STATUS
# ==============================================================================

message("\n=== 7. Victorian Temporal Trends by Cluster Status ===")

# Classify LGAs as persistent hot spot (>=4 of 6 years) vs other
persistent_lgas <- persistence$LGA[persistence$years_in_cluster >= 4]

vic_temporal_cluster <- vic_annual %>%
  mutate(
    cluster_status = ifelse(lga_name_clean %in% persistent_lgas,
                            "Persistent hot spot", "Other LGAs")
  ) %>%
  group_by(year, cluster_status) %>%
  summarise(
    n_lgas = n_distinct(lga_name_clean),
    total_count = sum(count, na.rm = TRUE),
    mean_pop = mean(population, na.rm = TRUE),
    mean_rate = mean(rate_per_100k, na.rm = TRUE),
    median_rate = median(rate_per_100k, na.rm = TRUE),
    .groups = "drop"
  )

fig_cluster_trends <- ggplot(vic_temporal_cluster,
                              aes(x = year, y = mean_rate, colour = cluster_status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Persistent hot spot" = "#d73027",
                                  "Other LGAs" = "#4575b4")) +
  scale_x_continuous(breaks = years) +
  labs(
    title = "Mean syphilis notification rate by spatial cluster status",
    subtitle = paste0("Victoria, 2019-2024 | Persistent hot spots: ",
                      length(persistent_lgas), " LGAs in H-H cluster \u22654 of 6 years"),
    x = NULL, y = "Mean rate per 100,000",
    colour = NULL
  ) +
  theme(legend.position = "bottom")

save_plot(fig_cluster_trends, "figS_cluster_temporal_trends", width = 8, height = 6)

# Summary table
cluster_trend_gt <- vic_temporal_cluster %>%
  mutate(across(c(mean_rate, median_rate), ~round(.x, 1))) %>%
  gt(groupname_col = "cluster_status") %>%
  tab_header(
    title = "Table S7. Temporal trends by spatial cluster status, VIC LGAs, 2019-2024",
    subtitle = paste0("Persistent hot spots: LGAs in LISA H-H cluster \u22654 of 6 years (n=",
                      length(persistent_lgas), ")")
  ) %>%
  cols_label(
    year = "Year",
    n_lgas = "N LGAs",
    total_count = "Total notifications",
    mean_rate = "Mean rate",
    median_rate = "Median rate"
  ) %>%
  cols_hide(mean_pop)

gtsave(cluster_trend_gt, file.path(path_tables, "tableS7_cluster_temporal_trends.html"))
write_csv(vic_temporal_cluster, file.path(path_tables, "tableS7_cluster_temporal_trends.csv"))

message("  Cluster temporal trends saved to tableS7 and figS_cluster_temporal_trends")

# ==============================================================================
# 8. LEAVE-ONE-OUT INFLUENCE ANALYSIS (Inner Melbourne LGAs)
# ==============================================================================

message("\n=== 8. Leave-One-Out Influence Analysis ===")

# Identify inner Melbourne LGAs (those ever in H-H cluster)
inner_melb_lgas <- unique(unlist(annual_lisa_lists))
message("  Testing influence of ", length(inner_melb_lgas), " inner Melbourne LGAs")
message("  LGAs: ", paste(inner_melb_lgas, collapse = ", "))

loo_results <- list()

# Full model baseline
full_rho <- lag_model$rho
full_aic <- AIC(lag_model)
full_male2044 <- coef(lag_model)["pct_male_20_44"]
full_dist <- coef(lag_model)["log_dist_sh"]

for (lga_name in inner_melb_lgas) {
  # Identify row index
  idx <- which(reg_data$LGA_NAME21 == lga_name)
  if (length(idx) == 0) next

  loo_data <- reg_data[-idx, ]

  # Rebuild spatial weights without this LGA
  nb_loo <- tryCatch({
    nb <- poly2nb(loo_data, queen = TRUE)
    # Fix islands
    n_isl_loo <- sum(card(nb) == 0)
    if (n_isl_loo > 0) {
      for (j in which(card(nb) == 0)) {
        d <- as.numeric(st_distance(st_centroid(loo_data[j, ]),
                                     st_centroid(loo_data[-j, ])))
        nr <- which.min(d)
        if (nr >= j) nr <- nr + 1
        nb[[j]] <- as.integer(nr)
        nb[[nr]] <- sort(unique(c(nb[[nr]], as.integer(j))))
      }
    }
    nb
  }, error = function(e) NULL)

  if (is.null(nb_loo)) {
    message("  Skipping ", lga_name, " — weights construction failed")
    next
  }

  lw_loo <- nb2listw(nb_loo, style = "W", zero.policy = TRUE)

  loo_lag <- tryCatch(
    lagsarlm(reg_formula, data = loo_data, listw = lw_loo, zero.policy = TRUE),
    error = function(e) NULL
  )

  if (!is.null(loo_lag)) {
    loo_results[[lga_name]] <- tibble(
      LGA_removed = lga_name,
      rho = round(loo_lag$rho, 4),
      rho_change = round(loo_lag$rho - full_rho, 4),
      AIC = round(AIC(loo_lag), 1),
      AIC_change = round(AIC(loo_lag) - full_aic, 1),
      beta_male2044 = round(coef(loo_lag)["pct_male_20_44"], 4),
      beta_male2044_change = round(coef(loo_lag)["pct_male_20_44"] - full_male2044, 4),
      beta_dist = round(coef(loo_lag)["log_dist_sh"], 4),
      beta_dist_change = round(coef(loo_lag)["log_dist_sh"] - full_dist, 4)
    )
    message("  Dropped ", lga_name,
            ": rho=", round(loo_lag$rho, 3),
            " (", sprintf("%+.3f", loo_lag$rho - full_rho), ")",
            ", AIC=", round(AIC(loo_lag), 1),
            " (", sprintf("%+.1f", AIC(loo_lag) - full_aic), ")")
  } else {
    message("  Skipping ", lga_name, " — lag model failed to converge")
  }
}

loo_table <- bind_rows(loo_results)

# Add full model row for reference
full_row <- tibble(
  LGA_removed = "None (full model)",
  rho = round(full_rho, 4), rho_change = 0,
  AIC = round(full_aic, 1), AIC_change = 0,
  beta_male2044 = round(full_male2044, 4), beta_male2044_change = 0,
  beta_dist = round(full_dist, 4), beta_dist_change = 0
)

loo_table <- bind_rows(full_row, loo_table)

loo_gt <- loo_table %>%
  gt() %>%
  tab_header(
    title = "Table S8. Leave-one-out influence analysis: inner Melbourne LGAs",
    subtitle = paste0("Spatial lag model re-estimated dropping each H-H cluster LGA (n=",
                      length(inner_melb_lgas), ")")
  ) %>%
  cols_label(
    LGA_removed = "LGA removed",
    rho = "\u03c1",
    rho_change = "\u0394\u03c1",
    AIC = "AIC",
    AIC_change = "\u0394AIC",
    beta_male2044 = "\u03b2 (% Males 20-44)",
    beta_male2044_change = "\u0394\u03b2",
    beta_dist = "\u03b2 (distance)",
    beta_dist_change = "\u0394\u03b2"
  ) %>%
  tab_source_note("Change (\u0394) relative to full model. Large changes indicate influential observations.")

gtsave(loo_gt, file.path(path_tables, "tableS8_leave_one_out.html"))
write_csv(loo_table, file.path(path_tables, "tableS8_leave_one_out.csv"))

# Influence plot
if (nrow(loo_table) > 1) {
  loo_plot_data <- loo_table %>% filter(LGA_removed != "None (full model)")

  fig_influence <- ggplot(loo_plot_data, aes(x = reorder(LGA_removed, abs(rho_change)),
                                              y = rho_change)) +
    geom_col(fill = ifelse(loo_plot_data$rho_change > 0, "#d73027", "#4575b4"),
             alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(
      title = "Leave-one-out influence on spatial lag parameter (\u03c1)",
      subtitle = paste0("Full model \u03c1 = ", round(full_rho, 3),
                        " | Bars show change when each LGA is removed"),
      x = NULL, y = paste0("Change in \u03c1 (relative to full model)")
    )

  save_plot(fig_influence, "figS_leave_one_out_influence", width = 8, height = 6)
  message("  Influence plot saved")
}

message("  Leave-one-out analysis saved to tableS8")

# ==============================================================================
# 9. MSM PRECINCT VALIDATION (qualitative)
# ==============================================================================

message("\n9. MSM precinct validation...")

# Known MSM community precincts in inner Melbourne:
# - Commercial Road / Prahran → Stonnington
# - Fitzroy / Smith Street → Yarra
# - St Kilda → Port Phillip
# - Collingwood / Abbotsford → Yarra
# - Northcote / Brunswick → Darebin / Merri-bek (formerly Moreland)
# - South Melbourne / Albert Park → Port Phillip

msm_precincts <- tibble(
  precinct = c("Commercial Rd / Prahran", "Fitzroy / Smith St",
               "St Kilda", "Collingwood / Abbotsford",
               "Northcote / Brunswick", "CBD"),
  expected_lga = c("Stonnington", "Yarra", "Port Phillip",
                   "Yarra", "Moreland", "Melbourne")
)

# Check which expected LGAs are in persistent H-H clusters
persistent_hh <- vic_spatial %>%
  st_drop_geometry() %>%
  filter(lisa_cluster == "High-High") %>%
  pull(LGA_NAME21)

msm_validation <- msm_precincts %>%
  mutate(
    in_lisa_hh = expected_lga %in% persistent_hh,
    validation = ifelse(in_lisa_hh, "Confirmed", "Not in H-H cluster")
  )

message("  MSM precinct → LISA cluster concordance:")
for (i in seq_len(nrow(msm_validation))) {
  message("    ", msm_validation$precinct[i], " (", msm_validation$expected_lga[i],
          "): ", msm_validation$validation[i])
}
concordance_pct <- round(mean(msm_validation$in_lisa_hh) * 100, 0)
message("  Overall concordance: ", concordance_pct, "%")

write_csv(msm_validation, file.path(path_tables, "tableS9_msm_precinct_validation.csv"))

# ==============================================================================
# 10. DUAL-PROXY VALIDATION (if SSCF data available)
# ==============================================================================

message("\n10. Dual-proxy validation...")

dual_path <- file.path(path_processed, "dual_proxy_results.rds")
if (file.exists(dual_path)) {
  dual <- readRDS(dual_path)

  message("  Proxy correlation:")
  message("    Pearson r = ", round(dual$proxy_cor$estimate, 3),
          " (p = ", format.pval(dual$proxy_cor$p.value, digits = 3), ")")
  message("    Spearman rho = ", round(dual$proxy_cor_spearman$estimate, 3),
          " (p = ", format.pval(dual$proxy_cor_spearman$p.value, digits = 3), ")")

  # Create comprehensive proxy comparison table
  m1_ols <- summary(dual$ols_m1)
  m2_ols <- summary(dual$ols_m2)
  m3_ols <- summary(dual$ols_m3)

  proxy_comparison <- tibble(
    Metric = c("OLS R-squared", "OLS Adj R-squared", "OLS AIC",
               "Lag AIC", "Lag rho", "Lag rho p-value",
               "MSM proxy coefficient (OLS)", "MSM proxy p-value (OLS)"),
    `Model 1: % Male 20-44` = c(
      round(m1_ols$r.squared, 4), round(m1_ols$adj.r.squared, 4),
      round(AIC(dual$ols_m1), 1), round(AIC(dual$lag_m1), 1),
      round(dual$lag_m1$rho, 4),
      format.pval(summary(dual$lag_m1)$Wald1$p.value, digits = 3),
      round(coef(m1_ols)["pct_male_20_44", 1], 4),
      format.pval(coef(m1_ols)["pct_male_20_44", 4], digits = 3)
    ),
    `Model 2: SS couple density` = c(
      round(m2_ols$r.squared, 4), round(m2_ols$adj.r.squared, 4),
      round(AIC(dual$ols_m2), 1), round(AIC(dual$lag_m2), 1),
      round(dual$lag_m2$rho, 4),
      format.pval(summary(dual$lag_m2)$Wald1$p.value, digits = 3),
      round(coef(m2_ols)["male_ss_per_1000", 1], 4),
      format.pval(coef(m2_ols)["male_ss_per_1000", 4], digits = 3)
    ),
    `Model 3: Both proxies` = c(
      round(m3_ols$r.squared, 4), round(m3_ols$adj.r.squared, 4),
      round(AIC(dual$ols_m3), 1), round(AIC(dual$lag_m3), 1),
      round(dual$lag_m3$rho, 4),
      format.pval(summary(dual$lag_m3)$Wald1$p.value, digits = 3),
      "see individual coefficients", "—"
    )
  )

  message("\n  Proxy comparison:")
  print(as.data.frame(proxy_comparison))

  write_csv(proxy_comparison,
            file.path(path_tables, "tableS10_proxy_comparison.csv"))

  # Save proxy comparison as gt table
  proxy_gt <- proxy_comparison %>%
    gt() %>%
    tab_header(
      title = "Table S10. MSM proxy comparison: demographic vs Census same-sex couple density",
      subtitle = paste0("Victorian LGAs (n=", nrow(dual$reg_data_ss), ")")
    ) %>%
    tab_source_note("Model 1 uses % males aged 20-44 as demographic MSM proxy. Model 2 uses male-male couples per 1,000 couple families from Census 2021 SSCF. Model 3 includes both proxies simultaneously.")

  gtsave(proxy_gt, file.path(path_tables, "tableS10_proxy_comparison.html"))

  message("  Proxy comparison table saved (S10)")
} else {
  message("  Dual-proxy results not available — run 04_spatial_autocorrelation.R with SSCF data")
}

# ==============================================================================
# SAVE UPDATED REGRESSION DATA
# ==============================================================================

saveRDS(reg_data, file.path(path_processed, "vic_lga_regression_data.rds"))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== Sensitivity Analyses Complete (MSM rewrite) ===")
message("Supplementary tables:")
message("  Table S1: Model diagnostics (VIF, normality, heteroscedasticity, residual Moran's I)")
message("  Table S2: Spatial weights sensitivity (queen, rook, k-NN, IDW)")
message("  Table S3: GWR local coefficient significance")
message("  Table S4: ITS model comparison (Poisson, NB, population offset)")
message("  Table S5: Temporal Moran's I evolution (2019-2024)")
message("  Table S6: LISA cluster persistence by LGA")
message("  Table S7: Temporal trends by cluster status")
message("  Table S8: Leave-one-out influence (inner Melbourne LGAs)")
message("  Table S9: MSM precinct validation")
if (file.exists(dual_path)) {
  message("  Table S10: Dual-proxy model comparison")
}
message("\nSupplementary figures:")
message("  figS_residual_qq_plot: OLS residual Q-Q plot")
message("  figS_gwr_distance_significance: GWR local significance map")
message("  figS_temporal_morans_i: Temporal Moran's I trend")
message("  figS_cluster_temporal_trends: Rate trends by cluster status")
message("  figS_leave_one_out_influence: LOO influence on rho")
if (file.exists(dual_path)) {
  message("  figS_proxy_validation: Scatter plot of proxy correlation")
}

# ==============================================================================
# 11. TEMPORAL STABILITY: 2016 vs 2021 Census same-sex couple distribution
# ==============================================================================
# Requires: data/raw/demographics/ABS_SSCF_VIC_LGA_2016.csv
# Download from ABS TableBuilder: 2016 Census > Counting Families > SSCF x LGA (VIC)
# Same process as 2021 data download

message("\n11. Temporal stability of MSM proxy (2016 vs 2021)...")

sscf_2016_path <- file.path("data/raw/demographics", "ABS_SSCF_VIC_LGA_2016.csv")
if (file.exists(sscf_2016_path)) {
  sscf_2016 <- read_csv(sscf_2016_path, show_col_types = FALSE)
  sscf_2021 <- readRDS(file.path(path_processed, "sscf_vic_lga.rds"))
  
  # Match LGA names between years
  merged <- sscf_2021 %>%
    select(lga_name_clean, couples_2021 = male_ss_couples) %>%
    inner_join(
      sscf_2016 %>%
        transmute(lga_name_clean = lga, couples_2016 = male_ss),
      by = "lga_name_clean"
    ) %>%
    filter(!is.na(couples_2016), !is.na(couples_2021))
  
  cor_pearson <- cor.test(merged$couples_2016, merged$couples_2021, method = "pearson")
  cor_spearman <- cor.test(merged$couples_2016, merged$couples_2021, method = "spearman", exact = FALSE)
  
  message("  Matched LGAs: ", nrow(merged))
  message("  Pearson r = ", round(cor_pearson$estimate, 4), " (p = ", format.pval(cor_pearson$p.value, 3), ")")
  message("  Spearman rho = ", round(cor_spearman$estimate, 4), " (p = ", format.pval(cor_spearman$p.value, 3), ")")
  
  fig_temporal <- ggplot(merged, aes(x = couples_2016, y = couples_2021)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", se = TRUE, colour = "#e74c3c") +
    labs(title = "Temporal stability: 2016 vs 2021 Census male-male couples by LGA",
         subtitle = paste0("r = ", round(cor_pearson$estimate, 3)),
         x = "Male-male couples (2016 Census)", y = "Male-male couples (2021 Census)")
  
  ggsave(file.path(path_figures, "figS_temporal_stability_sscf.png"), fig_temporal, width = 8, height = 7, dpi = 300)
  message("  Temporal stability figure saved")
} else {
  message("  2016 Census SSCF data not available")
  message("  Download from TableBuilder: 2016 Census > Counting Families > SSCF x LGA (VIC)")
  message("  Save to: ", sscf_2016_path)
}

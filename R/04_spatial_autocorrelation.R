# ==============================================================================
# 04_spatial_autocorrelation.R — Phase 3: Spatial Statistical Analysis
# Spatial autocorrelation, hot spot analysis, spatial regression, GWR
# ==============================================================================
#
# NOTE: SA2-level notification rates are not available — notifications are
# reported at sub-state level (NSW=LHD, VIC=LGA, QLD=HHS). This script uses
# VIC LGA data (n≈80) as the primary spatial analysis unit because it has the
# finest sub-state geography with both counts and population denominators.
# SA2 covariates are aggregated to LGA level for regression modelling.
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. PREPARE SPATIAL ANALYSIS DATASET
# ==============================================================================

message("=== Phase 3: Spatial Statistical Analysis ===")
message("1. Preparing VIC LGA spatial analysis dataset...")

# Load data
vic_lga_boundaries <- readRDS(file.path(path_processed, "boundaries_vic_lga.rds"))
vic_analysis <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))
sa2_sf <- readRDS(file.path(path_processed, "sa2_analysis_sf.rds"))

# Aggregate VIC notifications to get average annual rate (2019-2024)
vic_rates <- vic_analysis %>%
  mutate(
    lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>% str_trim()
  ) %>%
  group_by(lga_name_clean) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    mean_population = mean(population, na.rm = TRUE),
    n_years = n_distinct(year),
    mean_annual_rate = (total_count / n_years) / mean_population * 100000,
    rate_2024 = {
      d <- pick(everything())
      r <- d$rate_per_100k[d$year == 2024]
      if (length(r) == 0) NA_real_ else r[1]
    },
    rate_2019 = {
      d <- pick(everything())
      r <- d$rate_per_100k[d$year == 2019]
      if (length(r) == 0) NA_real_ else r[1]
    },
    .groups = "drop"
  ) %>%
  mutate(
    rate_change = rate_2024 - rate_2019
  )

# Aggregate SA2 covariates to VIC LGA via spatial join
sa2_vic <- sa2_sf %>%
  filter(STE_CODE21 == "2", !is.na(pop_total), pop_total > 0)

centroids <- st_centroid(sa2_vic)

joined <- st_join(
  centroids,
  vic_lga_boundaries[, c("LGA_CODE21", "LGA_NAME21")],
  left = TRUE
) %>%
  st_drop_geometry() %>%
  filter(!is.na(LGA_NAME21))

# Weighted aggregation — handle NA weights properly
lga_covariates <- joined %>%
  group_by(LGA_CODE21, LGA_NAME21) %>%
  summarise(
    n_sa2 = n(),
    lga_pop = sum(pop_total, na.rm = TRUE),
    lga_pop_indigenous = sum(pop_indigenous, na.rm = TRUE),
    pct_indigenous = lga_pop_indigenous / lga_pop * 100,
    lga_pop_male_20_44 = sum(pop_male_20_44, na.rm = TRUE),
    pct_male_20_44 = lga_pop_male_20_44 / lga_pop * 100,
    lga_area_sqkm = sum(area_sqkm, na.rm = TRUE),
    pop_density = ifelse(lga_area_sqkm > 0, lga_pop / lga_area_sqkm, NA_real_),
    mean_irsd = {
      valid <- !is.na(irsd_score) & !is.na(pop_total)
      if (sum(valid) > 0) weighted.mean(irsd_score[valid], pop_total[valid])
      else NA_real_
    },
    mean_dist_sh_km = {
      valid <- !is.na(dist_sh_clinic_km) & !is.na(pop_total)
      if (sum(valid) > 0) weighted.mean(dist_sh_clinic_km[valid], pop_total[valid])
      else NA_real_
    },
    mean_dist_ah_km = {
      valid <- !is.na(dist_ah_service_km) & !is.na(pop_total)
      if (sum(valid) > 0) weighted.mean(dist_ah_service_km[valid], pop_total[valid])
      else NA_real_
    },
    remoteness_mode = names(which.max(table(remoteness[!is.na(remoteness)]))),
    .groups = "drop"
  )

# Join rates + covariates to boundaries
vic_spatial <- vic_lga_boundaries %>%
  left_join(vic_rates, by = c("LGA_NAME21" = "lga_name_clean")) %>%
  left_join(lga_covariates %>% select(-LGA_NAME21), by = "LGA_CODE21") %>%
  filter(!is.na(mean_annual_rate), mean_annual_rate > 0 | total_count >= 0)

# Remove LGAs with no notification match (e.g., Unincorporated VIC)
vic_spatial <- vic_spatial %>% filter(!is.na(total_count))

message("  VIC LGA spatial dataset: ", nrow(vic_spatial), " LGAs with rates + covariates")
message("  Mean annual rate range: ",
        round(min(vic_spatial$mean_annual_rate, na.rm = TRUE), 1), " - ",
        round(max(vic_spatial$mean_annual_rate, na.rm = TRUE), 1), " per 100k")

# ==============================================================================
# 2. SPATIAL WEIGHTS MATRIX
# ==============================================================================

message("\n2. Building spatial weights matrix...")

# Queen contiguity neighbours
nb_queen <- poly2nb(vic_spatial, queen = TRUE)

# Check for islands (LGAs with no neighbours)
n_islands <- sum(card(nb_queen) == 0)
if (n_islands > 0) {
  message("  Warning: ", n_islands, " island(s) with no neighbours detected")
  # Add nearest neighbour for islands to prevent errors
  islands <- which(card(nb_queen) == 0)
  coords <- st_coordinates(st_centroid(vic_spatial))
  for (i in islands) {
    dists <- as.numeric(st_distance(
      st_centroid(vic_spatial[i, ]),
      st_centroid(vic_spatial[-i, ])
    ))
    nearest <- which.min(dists)
    # Adjust index (since we excluded i)
    if (nearest >= i) nearest <- nearest + 1
    nb_queen[[i]] <- as.integer(nearest)
    # Make symmetric
    nb_queen[[nearest]] <- sort(unique(c(nb_queen[[nearest]], as.integer(i))))
  }
  message("  Fixed by adding nearest neighbour connections")
}

# Row-standardised weights
lw <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)

message("  Neighbours: mean = ", round(mean(card(nb_queen)), 1),
        ", min = ", min(card(nb_queen)),
        ", max = ", max(card(nb_queen)))

# ==============================================================================
# 3. GLOBAL MORAN'S I
# ==============================================================================

message("\n3. Global Moran's I test...")

# Test on mean annual rate
moran_result <- moran.test(vic_spatial$mean_annual_rate, lw, zero.policy = TRUE)

message("  Moran's I = ", round(moran_result$estimate["Moran I statistic"], 4))
message("  Expected I = ", round(moran_result$estimate["Expectation"], 4))
message("  Variance = ", round(moran_result$estimate["Variance"], 6))
message("  p-value = ", format.pval(moran_result$p.value, digits = 4))
message("  Interpretation: ",
        ifelse(moran_result$p.value < 0.05,
               "Significant spatial autocorrelation detected",
               "No significant spatial autocorrelation"))

# Moran's I on 2024 rate
moran_2024 <- moran.test(
  vic_spatial$rate_2024,
  lw, zero.policy = TRUE, na.action = na.exclude
)

message("\n  2024 rate Moran's I = ", round(moran_2024$estimate["Moran I statistic"], 4),
        " (p = ", format.pval(moran_2024$p.value, digits = 4), ")")

# Moran scatterplot
moran_mc <- moran.mc(vic_spatial$mean_annual_rate, lw,
                     nsim = 999, zero.policy = TRUE)

fig_moran <- ggplot(data.frame(
  rate = vic_spatial$mean_annual_rate,
  lag_rate = lag.listw(lw, vic_spatial$mean_annual_rate, zero.policy = TRUE)
), aes(x = rate, y = lag_rate)) +
  geom_point(alpha = 0.6, colour = "#2c3e50") +
  geom_smooth(method = "lm", se = FALSE, colour = "#e74c3c", linewidth = 0.8) +
  geom_hline(yintercept = mean(lag.listw(lw, vic_spatial$mean_annual_rate,
                                          zero.policy = TRUE)),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = mean(vic_spatial$mean_annual_rate),
             linetype = "dashed", colour = "grey50") +
  labs(
    title = "Moran's I scatterplot: Infectious syphilis rates",
    subtitle = paste0("VIC LGAs, mean annual rate 2019-2024 | I = ",
                      round(moran_result$estimate["Moran I statistic"], 3),
                      ", p = ", format.pval(moran_result$p.value, digits = 3)),
    x = "Mean annual rate per 100,000",
    y = "Spatially lagged rate"
  )

save_plot(fig_moran, "fig10_moran_scatterplot", width = 8, height = 7)

# ==============================================================================
# 4. LOCAL MORAN'S I (LISA) CLUSTER ANALYSIS
# ==============================================================================

message("\n4. Local Moran's I (LISA) cluster analysis...")

# Compute local Moran's I
lisa <- localmoran(vic_spatial$mean_annual_rate, lw, zero.policy = TRUE)

# Classify LISA clusters
mean_rate <- mean(vic_spatial$mean_annual_rate)
lag_rate <- lag.listw(lw, vic_spatial$mean_annual_rate, zero.policy = TRUE)

vic_spatial <- vic_spatial %>%
  mutate(
    lisa_ii = lisa[, "Ii"],
    lisa_p = lisa[, "Pr(z != E(Ii))"],
    lisa_cluster = case_when(
      lisa_p >= 0.05 ~ "Not significant",
      mean_annual_rate > mean_rate & lag_rate > mean(lag_rate) ~ "High-High",
      mean_annual_rate < mean_rate & lag_rate < mean(lag_rate) ~ "Low-Low",
      mean_annual_rate > mean_rate & lag_rate < mean(lag_rate) ~ "High-Low",
      mean_annual_rate < mean_rate & lag_rate > mean(lag_rate) ~ "Low-High",
      TRUE ~ "Not significant"
    ),
    lisa_cluster = factor(lisa_cluster,
                          levels = c("High-High", "High-Low", "Low-High",
                                     "Low-Low", "Not significant"))
  )

# LISA cluster counts
lisa_counts <- table(vic_spatial$lisa_cluster)
message("  LISA clusters:")
for (cl in names(lisa_counts)) {
  message("    ", cl, ": ", lisa_counts[cl])
}

# LISA cluster map
lisa_colours <- c(
  "High-High" = "#d73027",
  "High-Low" = "#fdae61",
  "Low-High" = "#abd9e9",
  "Low-Low" = "#4575b4",
  "Not significant" = "#f0f0f0"
)

fig_lisa <- tm_shape(vic_spatial) +
  tm_polygons(
    fill = "lisa_cluster",
    fill.scale = tm_scale_categorical(
      values = lisa_colours
    ),
    fill.legend = tm_legend(title = "LISA cluster")
  ) +
  tm_title("Local Moran's I clusters: Infectious syphilis rates, VIC LGAs (2019-2024)")

save_map(fig_lisa, "fig11_lisa_clusters_vic", width = 8, height = 8)

# ==============================================================================
# 5. GETIS-ORD Gi* HOT SPOT ANALYSIS
# ==============================================================================

message("\n5. Getis-Ord Gi* hot spot analysis...")

# Gi* requires binary weights (not row-standardised)
lw_binary <- nb2listw(nb_queen, style = "B", zero.policy = TRUE)

gi_star <- localG(vic_spatial$mean_annual_rate, lw_binary, zero.policy = TRUE)

vic_spatial <- vic_spatial %>%
  mutate(
    gi_z = as.numeric(gi_star),
    hotspot_class = case_when(
      gi_z > 2.58  ~ "Hot spot (99% CI)",
      gi_z > 1.96  ~ "Hot spot (95% CI)",
      gi_z > 1.65  ~ "Hot spot (90% CI)",
      gi_z < -2.58 ~ "Cold spot (99% CI)",
      gi_z < -1.96 ~ "Cold spot (95% CI)",
      gi_z < -1.65 ~ "Cold spot (90% CI)",
      TRUE ~ "Not significant"
    ),
    hotspot_class = factor(hotspot_class,
                           levels = c("Hot spot (99% CI)", "Hot spot (95% CI)",
                                      "Hot spot (90% CI)", "Not significant",
                                      "Cold spot (90% CI)", "Cold spot (95% CI)",
                                      "Cold spot (99% CI)"))
  )

hotspot_counts <- table(vic_spatial$hotspot_class)
message("  Gi* hot/cold spots:")
for (cl in names(hotspot_counts)) {
  if (hotspot_counts[cl] > 0) message("    ", cl, ": ", hotspot_counts[cl])
}

# Hot spot map
hotspot_colours <- c(
  "Hot spot (99% CI)" = "#d73027",
  "Hot spot (95% CI)" = "#fc8d59",
  "Hot spot (90% CI)" = "#fee08b",
  "Not significant"   = "#f0f0f0",
  "Cold spot (90% CI)" = "#d9ef8b",
  "Cold spot (95% CI)" = "#91bfdb",
  "Cold spot (99% CI)" = "#4575b4"
)

fig_hotspot <- tm_shape(vic_spatial) +
  tm_polygons(
    fill = "hotspot_class",
    fill.scale = tm_scale_categorical(
      values = hotspot_colours
    ),
    fill.legend = tm_legend(title = "Gi* classification")
  ) +
  tm_title("Getis-Ord Gi* hot spots: Infectious syphilis rates, VIC LGAs (2019-2024)")

save_map(fig_hotspot, "fig12_hotspots_gi_star_vic", width = 8, height = 8)

# ==============================================================================
# 6. OLS REGRESSION (BASELINE)
# ==============================================================================

message("\n6. OLS regression (baseline model)...")

# Prepare regression data — exclude LGAs with missing covariates
reg_data <- vic_spatial %>%
  filter(
    !is.na(mean_annual_rate),
    !is.na(pct_male_20_44),
    !is.na(pop_density),
    !is.na(mean_irsd),
    !is.na(mean_dist_sh_km),
    !is.na(remoteness_mode)
  )

message("  Regression dataset: ", nrow(reg_data), " LGAs (of ", nrow(vic_spatial), ")")

# Log-transform distance, density, and rate for better model fit
reg_data <- reg_data %>%
  mutate(
    log_rate = log1p(mean_annual_rate),
    log_dist_sh = log1p(mean_dist_sh_km),
    log_pop_density = log1p(pop_density),
    remoteness_factor = factor(remoteness_mode,
                               levels = c("Major Cities", "Inner Regional",
                                          "Outer Regional", "Remote", "Very Remote"))
  )

# OLS model
ols_model <- lm(
  log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data
)

message("\n  OLS Model Summary:")
ols_summary <- summary(ols_model)
message("  R-squared: ", round(ols_summary$r.squared, 4))
message("  Adj R-squared: ", round(ols_summary$adj.r.squared, 4))
message("  F-statistic: ", round(ols_summary$fstatistic[1], 2),
        " (p = ", format.pval(pf(ols_summary$fstatistic[1],
                                  ols_summary$fstatistic[2],
                                  ols_summary$fstatistic[3],
                                  lower.tail = FALSE), digits = 4), ")")

# Print coefficients
cat("\n  Coefficients:\n")
coefs <- coef(ols_summary)
for (i in seq_len(nrow(coefs))) {
  sig <- ""
  if (coefs[i, 4] < 0.001) sig <- " ***"
  else if (coefs[i, 4] < 0.01) sig <- " **"
  else if (coefs[i, 4] < 0.05) sig <- " *"
  else if (coefs[i, 4] < 0.1) sig <- " ."
  message("    ", rownames(coefs)[i], ": ",
          round(coefs[i, 1], 4), " (SE=", round(coefs[i, 2], 4),
          ", p=", format.pval(coefs[i, 4], digits = 3), ")", sig)
}

# ==============================================================================
# 7. LAGRANGE MULTIPLIER TESTS FOR SPATIAL DEPENDENCE
# ==============================================================================

message("\n7. Lagrange Multiplier tests for spatial dependence...")

# Rebuild weights for regression subset
nb_reg <- poly2nb(reg_data, queen = TRUE)

# Fix islands in regression subset
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

rs_tests <- lm.RStests(ols_model, lw_reg, test = "all", zero.policy = TRUE)

message("  RS error:      statistic = ",
        round(rs_tests$RSerr$statistic, 4),
        ", p = ", format.pval(rs_tests$RSerr$p.value, digits = 4))
message("  RS lag:        statistic = ",
        round(rs_tests$RSlag$statistic, 4),
        ", p = ", format.pval(rs_tests$RSlag$p.value, digits = 4))
message("  Adj RS error:  statistic = ",
        round(rs_tests$adjRSerr$statistic, 4),
        ", p = ", format.pval(rs_tests$adjRSerr$p.value, digits = 4))
message("  Adj RS lag:    statistic = ",
        round(rs_tests$adjRSlag$statistic, 4),
        ", p = ", format.pval(rs_tests$adjRSlag$p.value, digits = 4))
message("  SARMA:         statistic = ",
        round(rs_tests$SARMA$statistic, 4),
        ", p = ", format.pval(rs_tests$SARMA$p.value, digits = 4))

# Determine which spatial model to use
lm_lag_sig <- rs_tests$RSlag$p.value < 0.05
lm_err_sig <- rs_tests$RSerr$p.value < 0.05
rlm_lag_sig <- rs_tests$adjRSlag$p.value < 0.05
rlm_err_sig <- rs_tests$adjRSerr$p.value < 0.05

if (lm_lag_sig & lm_err_sig) {
  message("\n  Both RS tests significant - compare adjusted versions:")
  if (rlm_lag_sig & !rlm_err_sig) {
    message("  -> Adj RS lag significant: spatial LAG model preferred")
  } else if (!rlm_lag_sig & rlm_err_sig) {
    message("  -> Adj RS error significant: spatial ERROR model preferred")
  } else {
    message("  -> Both/neither adjusted tests significant: fit both models, compare AIC")
  }
} else if (lm_lag_sig) {
  message("\n  -> RS lag significant: spatial LAG model preferred")
} else if (lm_err_sig) {
  message("\n  -> RS error significant: spatial ERROR model preferred")
} else {
  message("\n  -> Neither RS test significant: OLS may be adequate (no spatial dependence)")
}

# ==============================================================================
# 8. SPATIAL LAG MODEL
# ==============================================================================

message("\n8. Spatial lag model...")

lag_model <- lagsarlm(
  log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data,
  listw = lw_reg,
  zero.policy = TRUE
)

lag_summary <- summary(lag_model)
message("  Rho (spatial lag parameter): ", round(lag_model$rho, 4),
        " (p = ", format.pval(lag_summary$Wald1$p.value, digits = 4), ")")
message("  AIC: ", round(AIC(lag_model), 2))
message("  Log-likelihood: ", round(lag_model$LL, 4))

cat("\n  Lag model coefficients:\n")
lag_coefs <- lag_summary$Coef
for (i in seq_len(nrow(lag_coefs))) {
  sig <- ""
  if (lag_coefs[i, 4] < 0.001) sig <- " ***"
  else if (lag_coefs[i, 4] < 0.01) sig <- " **"
  else if (lag_coefs[i, 4] < 0.05) sig <- " *"
  else if (lag_coefs[i, 4] < 0.1) sig <- " ."
  message("    ", rownames(lag_coefs)[i], ": ",
          round(lag_coefs[i, 1], 4), " (SE=", round(lag_coefs[i, 2], 4),
          ", p=", format.pval(lag_coefs[i, 4], digits = 3), ")", sig)
}

# ==============================================================================
# 9. SPATIAL ERROR MODEL
# ==============================================================================

message("\n9. Spatial error model...")

error_model <- errorsarlm(
  log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
  data = reg_data,
  listw = lw_reg,
  zero.policy = TRUE
)

error_summary <- summary(error_model)
message("  Lambda (spatial error parameter): ", round(error_model$lambda, 4),
        " (p = ", format.pval(error_summary$Wald1$p.value, digits = 4), ")")
message("  AIC: ", round(AIC(error_model), 2))
message("  Log-likelihood: ", round(error_model$LL, 4))

cat("\n  Error model coefficients:\n")
err_coefs <- error_summary$Coef
for (i in seq_len(nrow(err_coefs))) {
  sig <- ""
  if (err_coefs[i, 4] < 0.001) sig <- " ***"
  else if (err_coefs[i, 4] < 0.01) sig <- " **"
  else if (err_coefs[i, 4] < 0.05) sig <- " *"
  else if (err_coefs[i, 4] < 0.1) sig <- " ."
  message("    ", rownames(err_coefs)[i], ": ",
          round(err_coefs[i, 1], 4), " (SE=", round(err_coefs[i, 2], 4),
          ", p=", format.pval(err_coefs[i, 4], digits = 3), ")", sig)
}

# ==============================================================================
# 10. MODEL COMPARISON
# ==============================================================================

message("\n10. Model comparison...")

model_comparison <- tibble(
  model = c("OLS", "Spatial Lag", "Spatial Error"),
  aic = c(AIC(ols_model), AIC(lag_model), AIC(error_model)),
  log_lik = c(logLik(ols_model), lag_model$LL, error_model$LL),
  spatial_param = c(NA, lag_model$rho, error_model$lambda),
  spatial_p = c(NA, lag_summary$Wald1$p.value, error_summary$Wald1$p.value)
)

message("  Model comparison:")
message("  ", sprintf("%-15s AIC=%-8.2f LL=%-8.2f Spatial=%-6s p=%s",
                       model_comparison$model,
                       model_comparison$aic,
                       model_comparison$log_lik,
                       ifelse(is.na(model_comparison$spatial_param), "-",
                              round(model_comparison$spatial_param, 3)),
                       ifelse(is.na(model_comparison$spatial_p), "-",
                              format.pval(model_comparison$spatial_p, digits = 3))))

best_model <- model_comparison$model[which.min(model_comparison$aic)]
message("\n  Best model (lowest AIC): ", best_model)

# ==============================================================================
# 11. REGRESSION RESULTS TABLE
# ==============================================================================

message("\n11. Creating regression results table...")

# Build comparison table
build_coef_table <- function(ols_mod, lag_mod, err_mod) {
  ols_c <- coef(summary(ols_mod))
  lag_c <- summary(lag_mod)$Coef
  err_c <- summary(err_mod)$Coef

  vars <- rownames(ols_c)

  tibble(
    variable = vars,
    ols_estimate = ols_c[, 1],
    ols_se = ols_c[, 2],
    ols_p = ols_c[, 4],
    lag_estimate = lag_c[vars, 1],
    lag_se = lag_c[vars, 2],
    lag_p = lag_c[vars, 4],
    error_estimate = err_c[vars, 1],
    error_se = err_c[vars, 2],
    error_p = err_c[vars, 4]
  )
}

coef_table <- build_coef_table(ols_model, lag_model, error_model)

# Format for display
format_coef <- function(est, se, p) {
  sig <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                   p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
  paste0(round(est, 3), " (", round(se, 3), ")", sig)
}

table2_display <- coef_table %>%
  transmute(
    Variable = case_when(
      variable == "(Intercept)" ~ "Intercept",
      variable == "pct_male_20_44" ~ "% Males aged 20-44",
      variable == "log_pop_density" ~ "log(Population density, per km2)",
      variable == "mean_irsd" ~ "IRSD score (socioeconomic)",
      variable == "log_dist_sh" ~ "log(Distance to SH clinic, km)",
      variable == "remoteness_factorInner Regional" ~ "Inner Regional (ref: Major Cities)",
      variable == "remoteness_factorOuter Regional" ~ "Outer Regional",
      variable == "remoteness_factorRemote" ~ "Remote",
      variable == "remoteness_factorVery Remote" ~ "Very Remote",
      TRUE ~ variable
    ),
    OLS = format_coef(ols_estimate, ols_se, ols_p),
    `Spatial Lag` = format_coef(lag_estimate, lag_se, lag_p),
    `Spatial Error` = format_coef(error_estimate, error_se, error_p)
  )

# Add model fit stats
fit_row <- tibble(
  Variable = c("AIC", "Log-likelihood", "Spatial parameter", "Spatial p-value",
               "N (LGAs)"),
  OLS = c(round(AIC(ols_model), 1),
          round(as.numeric(logLik(ols_model)), 1),
          "-", "-", nrow(reg_data)),
  `Spatial Lag` = c(round(AIC(lag_model), 1), round(lag_model$LL, 1),
                    paste0("rho=", round(lag_model$rho, 3)),
                    format.pval(lag_summary$Wald1$p.value, digits = 3),
                    nrow(reg_data)),
  `Spatial Error` = c(round(AIC(error_model), 1), round(error_model$LL, 1),
                      paste0("lambda=", round(error_model$lambda, 3)),
                      format.pval(error_summary$Wald1$p.value, digits = 3),
                      nrow(reg_data))
)

table2_full <- bind_rows(table2_display, fit_row)

table2_gt <- table2_full %>%
  gt() %>%
  tab_header(
    title = "Table 2. Spatial regression models: Predictors of infectious syphilis rates",
    subtitle = "Victorian LGAs, mean annual rate 2019-2024 (log-transformed)"
  ) %>%
  tab_source_note("Coefficients shown as estimate (SE). Significance: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1") %>%
  tab_source_note("Source: VIC Health notifications, ABS Census 2021, NHSD 2025")

gtsave(table2_gt, file.path(path_tables, "table2_spatial_regression.html"))
write_csv(table2_full, file.path(path_tables, "table2_spatial_regression.csv"))

# ==============================================================================
# 11b. SAME-SEX COUPLE PROXY REGRESSION (if data available)
# ==============================================================================

message("\n11b. Same-sex couple proxy regression models...")

sscf_density_path <- file.path(path_processed, "sscf_vic_lga_density.rds")

if (file.exists(sscf_density_path)) {
  sscf_density <- readRDS(sscf_density_path)

  # Join SSCF density to regression data
  reg_data_ss <- reg_data %>%
    left_join(sscf_density, by = c("LGA_CODE21")) %>%
    filter(!is.na(male_ss_per_1000))

  message("  LGAs with SSCF data: ", nrow(reg_data_ss), " of ", nrow(reg_data))

  if (nrow(reg_data_ss) >= 30) {
    # ---- Proxy validation: correlate pct_male_20_44 with male_ss_per_1000 ----
    proxy_cor <- cor.test(reg_data_ss$pct_male_20_44, reg_data_ss$male_ss_per_1000,
                          method = "pearson")
    proxy_cor_spearman <- cor.test(reg_data_ss$pct_male_20_44, reg_data_ss$male_ss_per_1000,
                                   method = "spearman")

    message("  Proxy correlation (pct_male_20_44 vs male_ss_per_1000):")
    message("    Pearson r = ", round(proxy_cor$estimate, 3),
            " (p = ", format.pval(proxy_cor$p.value, digits = 3), ")")
    message("    Spearman rho = ", round(proxy_cor_spearman$estimate, 3),
            " (p = ", format.pval(proxy_cor_spearman$p.value, digits = 3), ")")

    # Scatter plot: proxy validation
    fig_proxy <- ggplot(reg_data_ss, aes(x = pct_male_20_44, y = male_ss_per_1000)) +
      geom_point(alpha = 0.6, colour = "#2c3e50") +
      geom_smooth(method = "lm", se = TRUE, colour = "#e74c3c", linewidth = 0.8) +
      labs(
        title = "Proxy validation: % Males 20-44 vs Same-sex couple density",
        subtitle = paste0("Victorian LGAs (n=", nrow(reg_data_ss),
                          ") | r = ", round(proxy_cor$estimate, 3),
                          ", p = ", format.pval(proxy_cor$p.value, digits = 3)),
        x = "% Males aged 20-44 (Census 2021)",
        y = "Male-male couples per 1,000 couple families"
      )
    save_plot(fig_proxy, "figS_proxy_validation", width = 8, height = 7)

    # ---- Model 2: Same-sex couple density as MSM proxy ----
    # Rebuild spatial weights for SSCF subset
    nb_ss <- poly2nb(reg_data_ss, queen = TRUE)
    n_islands_ss <- sum(card(nb_ss) == 0)
    if (n_islands_ss > 0) {
      islands_ss <- which(card(nb_ss) == 0)
      for (i in islands_ss) {
        dists <- as.numeric(st_distance(st_centroid(reg_data_ss[i, ]),
                                         st_centroid(reg_data_ss[-i, ])))
        nearest <- which.min(dists)
        if (nearest >= i) nearest <- nearest + 1
        nb_ss[[i]] <- as.integer(nearest)
        nb_ss[[nearest]] <- sort(unique(c(nb_ss[[nearest]], as.integer(i))))
      }
    }
    lw_ss <- nb2listw(nb_ss, style = "W", zero.policy = TRUE)

    # OLS Model 2: same-sex proxy
    ols_ss <- lm(
      log_rate ~ male_ss_per_1000 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
      data = reg_data_ss
    )
    ols_ss_summ <- summary(ols_ss)

    message("\n  Model 2 (same-sex proxy) OLS:")
    message("    R-squared: ", round(ols_ss_summ$r.squared, 4))
    message("    Adj R-squared: ", round(ols_ss_summ$adj.r.squared, 4))
    cat("    Coefficients:\n")
    coefs_ss <- coef(ols_ss_summ)
    for (i in seq_len(nrow(coefs_ss))) {
      sig <- ifelse(coefs_ss[i, 4] < 0.001, "***",
             ifelse(coefs_ss[i, 4] < 0.01, "**",
             ifelse(coefs_ss[i, 4] < 0.05, "*", "")))
      message("      ", rownames(coefs_ss)[i], ": ",
              round(coefs_ss[i, 1], 4), " (p=", format.pval(coefs_ss[i, 4], digits = 3), ")", sig)
    }

    # Spatial lag Model 2
    lag_ss <- lagsarlm(
      log_rate ~ male_ss_per_1000 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
      data = reg_data_ss, listw = lw_ss, zero.policy = TRUE
    )
    lag_ss_summ <- summary(lag_ss)

    message("\n  Model 2 spatial lag:")
    message("    Rho: ", round(lag_ss$rho, 4),
            " (p = ", format.pval(lag_ss_summ$Wald1$p.value, digits = 4), ")")
    message("    AIC: ", round(AIC(lag_ss), 2))

    # Spatial error Model 2
    err_ss <- errorsarlm(
      log_rate ~ male_ss_per_1000 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
      data = reg_data_ss, listw = lw_ss, zero.policy = TRUE
    )

    # ---- Model 3: Both proxies combined ----
    ols_both <- lm(
      log_rate ~ pct_male_20_44 + male_ss_per_1000 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
      data = reg_data_ss
    )
    ols_both_summ <- summary(ols_both)

    message("\n  Model 3 (combined proxies) OLS:")
    message("    R-squared: ", round(ols_both_summ$r.squared, 4))
    coefs_both <- coef(ols_both_summ)
    for (i in seq_len(nrow(coefs_both))) {
      sig <- ifelse(coefs_both[i, 4] < 0.001, "***",
             ifelse(coefs_both[i, 4] < 0.01, "**",
             ifelse(coefs_both[i, 4] < 0.05, "*", "")))
      message("      ", rownames(coefs_both)[i], ": ",
              round(coefs_both[i, 1], 4), " (p=", format.pval(coefs_both[i, 4], digits = 3), ")", sig)
    }

    lag_both <- lagsarlm(
      log_rate ~ pct_male_20_44 + male_ss_per_1000 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
      data = reg_data_ss, listw = lw_ss, zero.policy = TRUE
    )

    # ---- Model comparison table (all specifications) ----
    model_comparison_dual <- tibble(
      Model = c("Model 1: % Male 20-44 (OLS)", "Model 1: % Male 20-44 (Lag)",
                "Model 2: SS couple density (OLS)", "Model 2: SS couple density (Lag)",
                "Model 3: Both proxies (OLS)", "Model 3: Both proxies (Lag)"),
      N = nrow(reg_data_ss),
      AIC = c(AIC(lm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor, data = reg_data_ss)),
              AIC(lagsarlm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
                           data = reg_data_ss, listw = lw_ss, zero.policy = TRUE)),
              AIC(ols_ss), AIC(lag_ss),
              AIC(ols_both), AIC(lag_both)),
      R2_or_rho = c(
        round(summary(lm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor, data = reg_data_ss))$r.squared, 3),
        round(lagsarlm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
                        data = reg_data_ss, listw = lw_ss, zero.policy = TRUE)$rho, 3),
        round(ols_ss_summ$r.squared, 3), round(lag_ss$rho, 3),
        round(ols_both_summ$r.squared, 3), round(lag_both$rho, 3)
      )
    )

    message("\n  Dual-proxy model comparison:")
    print(as.data.frame(model_comparison_dual))

    # Save comparison table
    write_csv(model_comparison_dual,
              file.path(path_tables, "tableS_dual_proxy_comparison.csv"))

    # Save results for manuscript figures
    saveRDS(list(
      proxy_cor = proxy_cor,
      proxy_cor_spearman = proxy_cor_spearman,
      ols_m1 = lm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor, data = reg_data_ss),
      ols_m2 = ols_ss,
      ols_m3 = ols_both,
      lag_m1 = lagsarlm(log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh + remoteness_factor,
                         data = reg_data_ss, listw = lw_ss, zero.policy = TRUE),
      lag_m2 = lag_ss,
      lag_m3 = lag_both,
      err_m2 = err_ss,
      comparison = model_comparison_dual,
      reg_data_ss = reg_data_ss
    ), file.path(path_processed, "dual_proxy_results.rds"))

    message("  Dual-proxy results saved")
  } else {
    message("  Insufficient LGAs with SSCF data for regression (need >= 30, have ", nrow(reg_data_ss), ")")
  }
} else {
  message("  SSCF density data not available — skipping dual-proxy analysis")
  message("  Run 01_data_acquisition.R after placing ABS_SSCF_VIC_LGA_2021.csv in data/raw/demographics/")
}

# ==============================================================================
# 12. GEOGRAPHICALLY WEIGHTED REGRESSION (GWR)
# ==============================================================================

message("\n12. Geographically Weighted Regression...")

library(spgwr)

# Convert to Spatial for spgwr compatibility
coords <- st_coordinates(st_centroid(reg_data))

# Simpler model for GWR (fewer parameters given n~80)
# GWR needs n >> parameters at each local fit
# Uses pct_male_20_44 and log_pop_density instead of pct_indigenous
gwr_formula <- log_rate ~ pct_male_20_44 + log_pop_density + mean_irsd + log_dist_sh

# Bandwidth selection via cross-validation
message("  Selecting optimal bandwidth (cross-validation)...")
bw <- tryCatch({
  gwr.sel(
    gwr_formula,
    data = as(reg_data, "Spatial"),
    coords = coords,
    adapt = TRUE,  # Adaptive bandwidth (proportion of data)
    verbose = FALSE
  )
}, error = function(e) {
  message("  GWR bandwidth selection failed: ", e$message)
  message("  Using fixed adaptive bandwidth of 0.5")
  0.5
})

message("  Optimal adaptive bandwidth: ", round(bw, 4))

# Fit GWR model
gwr_model <- tryCatch({
  gwr(
    gwr_formula,
    data = as(reg_data, "Spatial"),
    coords = coords,
    adapt = bw,
    hatmatrix = TRUE
  )
}, error = function(e) {
  message("  GWR model fitting failed: ", e$message)
  NULL
})

if (!is.null(gwr_model)) {
  message("  GWR fitted successfully")
  message("  Global R-squared: ", round(ols_summary$r.squared, 4))

  # Extract local results
  gwr_results <- as.data.frame(gwr_model$SDF)

  reg_data <- reg_data %>%
    mutate(
      gwr_localR2 = gwr_results$localR2,
      gwr_coef_male2044 = gwr_results$pct_male_20_44,
      gwr_coef_density = gwr_results$log_pop_density,
      gwr_coef_irsd = gwr_results$mean_irsd,
      gwr_coef_dist = gwr_results$log_dist_sh
    )

  message("  Local R-squared range: ",
          round(min(reg_data$gwr_localR2, na.rm = TRUE), 3), " - ",
          round(max(reg_data$gwr_localR2, na.rm = TRUE), 3))

  # Map local R-squared
  fig_gwr_r2 <- tm_shape(reg_data) +
    tm_polygons(
      fill = "gwr_localR2",
      fill.scale = tm_scale_continuous(
        values = "brewer.yl_gn_bu",
        midpoint = NA
      ),
      fill.legend = tm_legend(title = "Local R-squared")
    ) +
    tm_title("GWR local R-squared: Syphilis rate model, VIC LGAs")

  save_map(fig_gwr_r2, "fig13_gwr_local_r2", width = 8, height = 8)

  # Map local coefficient for distance to SH clinic
  fig_gwr_dist <- tm_shape(reg_data) +
    tm_polygons(
      fill = "gwr_coef_dist",
      fill.scale = tm_scale_continuous(
        values = "brewer.rd_yl_bu",
        midpoint = 0
      ),
      fill.legend = tm_legend(title = "Coefficient:\nlog(Distance to\nSH clinic)")
    ) +
    tm_title("GWR local coefficient: Distance to sexual health clinic, VIC LGAs")

  save_map(fig_gwr_dist, "fig14_gwr_coef_distance", width = 8, height = 8)

  # Map local coefficient for % males 20-44
  fig_gwr_male2044 <- tm_shape(reg_data) +
    tm_polygons(
      fill = "gwr_coef_male2044",
      fill.scale = tm_scale_continuous(
        values = "brewer.rd_yl_bu",
        midpoint = 0
      ),
      fill.legend = tm_legend(title = "Coefficient:\n% Males 20-44")
    ) +
    tm_title("GWR local coefficient: % Males aged 20-44, VIC LGAs")

  save_map(fig_gwr_male2044, "fig15_gwr_coef_male2044", width = 8, height = 8)

  # Map local coefficient for population density
  fig_gwr_density <- tm_shape(reg_data) +
    tm_polygons(
      fill = "gwr_coef_density",
      fill.scale = tm_scale_continuous(
        values = "brewer.rd_yl_bu",
        midpoint = 0
      ),
      fill.legend = tm_legend(title = "Coefficient:\nlog(Pop density)")
    ) +
    tm_title("GWR local coefficient: Population density, VIC LGAs")

  save_map(fig_gwr_density, "fig16_gwr_coef_density", width = 8, height = 8)
} else {
  message("  GWR skipped due to fitting error")
}

# ==============================================================================
# 13. QLD HHS SPATIAL AUTOCORRELATION (SUPPLEMENTARY) — REMOVED in MSM rewrite
# ==============================================================================
# QLD HHS analysis removed: paper now focuses on Victoria only.
# Original code retained as comments for reference.

if (FALSE) {
message("\n13. Supplementary: QLD HHS spatial autocorrelation...")

qld_hhs_boundaries <- readRDS(file.path(path_processed, "boundaries_qld_hhs.rds"))
qld_hhs_notif <- readRDS(file.path(path_processed, "qld_hhs_notifications.rds"))

# Average rate 2020-2024
qld_rates <- qld_hhs_notif %>%
  filter(year >= 2020, area_name != "Queensland") %>%
  group_by(area_name) %>%
  summarise(mean_rate = mean(rate, na.rm = TRUE), .groups = "drop")

qld_spatial <- qld_hhs_boundaries %>%
  left_join(qld_rates, by = c("hhs" = "area_name")) %>%
  filter(!is.na(mean_rate))

if (nrow(qld_spatial) >= 5) {
  nb_qld <- poly2nb(qld_spatial, queen = TRUE)

  # Fix islands
  n_islands_q <- sum(card(nb_qld) == 0)
  if (n_islands_q > 0) {
    islands_q <- which(card(nb_qld) == 0)
    for (i in islands_q) {
      dists <- as.numeric(st_distance(
        st_centroid(qld_spatial[i, ]),
        st_centroid(qld_spatial[-i, ])
      ))
      nearest <- which.min(dists)
      if (nearest >= i) nearest <- nearest + 1
      nb_qld[[i]] <- as.integer(nearest)
      nb_qld[[nearest]] <- sort(unique(c(nb_qld[[nearest]], as.integer(i))))
    }
  }

  lw_qld <- nb2listw(nb_qld, style = "W", zero.policy = TRUE)

  moran_qld <- moran.test(qld_spatial$mean_rate, lw_qld, zero.policy = TRUE)
  message("  QLD HHS Moran's I = ", round(moran_qld$estimate["Moran I statistic"], 4),
          " (p = ", format.pval(moran_qld$p.value, digits = 4), ")")
  message("  N = ", nrow(qld_spatial), " HHS areas (limited for spatial statistics)")

  # Gi* for QLD
  lw_qld_b <- nb2listw(nb_qld, style = "B", zero.policy = TRUE)
  gi_qld <- localG(qld_spatial$mean_rate, lw_qld_b, zero.policy = TRUE)

  qld_spatial <- qld_spatial %>%
    mutate(
      gi_z = as.numeric(gi_qld),
      hotspot_class = case_when(
        gi_z > 1.96  ~ "Hot spot",
        gi_z < -1.96 ~ "Cold spot",
        TRUE ~ "Not significant"
      )
    )

  hotspot_colours_qld <- c(
    "Hot spot" = "#d73027",
    "Not significant" = "#f0f0f0",
    "Cold spot" = "#4575b4"
  )

  fig_qld_hotspot <- tm_shape(qld_spatial) +
    tm_polygons(
      fill = "hotspot_class",
      fill.scale = tm_scale_categorical(values = hotspot_colours_qld),
      fill.legend = tm_legend(title = "Gi* classification")
    ) +
    tm_title("Getis-Ord Gi* hot spots: Infectious syphilis rates, QLD HHS (2020-2024)")

  save_map(fig_qld_hotspot, "fig16_hotspots_qld_hhs", width = 8, height = 10)
} else {
  message("  Too few QLD HHS areas for spatial analysis")
}
} # end if(FALSE) — QLD HHS section removed

# ==============================================================================
# 14. SAVE SPATIAL ANALYSIS DATASETS
# ==============================================================================

message("\n14. Saving spatial analysis outputs...")

saveRDS(vic_spatial, file.path(path_processed, "vic_lga_spatial_results.rds"))
if (exists("reg_data")) {
  saveRDS(reg_data, file.path(path_processed, "vic_lga_regression_data.rds"))
}
saveRDS(model_comparison, file.path(path_processed, "model_comparison.rds"))

# ==============================================================================
# 15. SUMMARY
# ==============================================================================

message("\n=== Phase 3 Output Summary ===")
message("Figures saved to: ", path_figures)
message("Maps saved to: ", path_maps)
message("Tables saved to: ", path_tables)

fig_files <- list.files(path_figures, pattern = "fig1[0-6]")
map_files <- list.files(path_maps, pattern = "fig1[1-6]")
table_files <- list.files(path_tables, pattern = "table2")

message("  New figures: ", length(fig_files))
for (f in fig_files) message("    ", f)
message("  New maps: ", length(map_files))
for (f in map_files) message("    ", f)
message("  New tables: ", length(table_files))
for (f in table_files) message("    ", f)

message("\n  Key results:")
message("  - Global Moran's I (VIC mean rate): ",
        round(moran_result$estimate["Moran I statistic"], 4),
        " (p=", format.pval(moran_result$p.value, digits = 3), ")")
message("  - Best regression model: ", best_model,
        " (AIC=", round(min(model_comparison$aic), 1), ")")
# QLD HHS analysis removed in MSM-focused rewrite

message("\nPhase 3 complete.")

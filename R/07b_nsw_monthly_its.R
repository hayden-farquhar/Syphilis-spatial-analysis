# ==============================================================================
# 07b_nsw_monthly_its.R — Sensitivity Analysis: NSW Monthly ITS
# Monthly interrupted time series for NSW infectious syphilis (2009-2025)
# ==============================================================================
#
# Reviewer requested a sub-analysis with finer temporal resolution to
# strengthen the national annual ITS (which has only 10 data points).
# NSW provides monthly notification data (204 observations), dramatically
# increasing degrees of freedom and enabling seasonal adjustment.
#
# Data source: NSW Health Notifiable Conditions Information Management System
# (NCIMS), extracted from NSW_syphilis_by_month_2009-2026.xls
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

library(MASS)  # for glm.nb

# ==============================================================================
# 1. LOAD NSW MONTHLY DATA
# ==============================================================================

message("=== NSW Monthly ITS Sensitivity Analysis ===")
message("1. Loading NSW monthly infectious syphilis data...")

nsw_monthly <- read_csv(
  file.path(path_notifications, "NSW_infectious_syphilis_monthly_2009-2025.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    date = as.Date(paste(year, month_num, 1, sep = "-")),
    time_index = row_number(),
    # COVID periods for NSW
    # First case: 25 Jan 2020; first lockdown: 23 March 2020
    # Major lockdowns: Mar-Jun 2020, Jun-Oct 2021
    # Border restrictions eased: Dec 2021
    covid = as.integer(date >= as.Date("2020-03-01") & date <= as.Date("2021-12-01")),
    post_covid = as.integer(date > as.Date("2021-12-01")),
    time_since_covid = pmax(0, time_index - min(time_index[covid == 1], na.rm = TRUE)),
    time_since_post = pmax(0, time_index - min(time_index[post_covid == 1], na.rm = TRUE)),
    # Seasonal terms (harmonic)
    sin12 = sin(2 * pi * month_num / 12),
    cos12 = cos(2 * pi * month_num / 12),
    sin6 = sin(4 * pi * month_num / 12),
    cos6 = cos(4 * pi * month_num / 12)
  )

message("  Loaded: ", nrow(nsw_monthly), " monthly observations (",
        min(nsw_monthly$year), "-", max(nsw_monthly$year), ")")
message("  Pre-COVID months: ", sum(nsw_monthly$covid == 0 & nsw_monthly$post_covid == 0))
message("  COVID restriction months: ", sum(nsw_monthly$covid == 1))
message("  Post-COVID months: ", sum(nsw_monthly$post_covid == 1))
message("  Total infectious syphilis notifications: ", format(sum(nsw_monthly$count), big.mark = ","))

# Save processed monthly data
saveRDS(nsw_monthly, file.path(path_processed, "nsw_monthly_infectious_syphilis.rds"))

# ==============================================================================
# 2. DESCRIPTIVE OVERVIEW
# ==============================================================================

message("\n2. Descriptive statistics...")

annual_nsw <- nsw_monthly %>%
  group_by(year) %>%
  summarise(
    total = sum(count),
    mean_monthly = mean(count),
    .groups = "drop"
  )

message("  Annual totals:")
for (i in seq_len(nrow(annual_nsw))) {
  message("    ", annual_nsw$year[i], ": ", annual_nsw$total[i],
          " (mean ", round(annual_nsw$mean_monthly[i], 1), "/month)")
}

# ==============================================================================
# 3. ITS MODEL: NEGATIVE BINOMIAL WITH SEASONAL TERMS
# ==============================================================================

message("\n3. Fitting ITS models...")

# Model 1: Poisson ITS with seasonal adjustment
its_pois <- glm(
  count ~ time_index + covid + time_since_covid +
    post_covid + time_since_post +
    sin12 + cos12 + sin6 + cos6,
  data = nsw_monthly,
  family = poisson()
)

# Check overdispersion
disp_ratio <- sum(residuals(its_pois, type = "pearson")^2) / its_pois$df.residual
message("  Poisson dispersion ratio: ", round(disp_ratio, 2),
        ifelse(disp_ratio > 1.5, " -> overdispersed, using negative binomial", " -> acceptable"))

# Model 2: Negative binomial ITS with seasonal adjustment
its_nb <- glm.nb(
  count ~ time_index + covid + time_since_covid +
    post_covid + time_since_post +
    sin12 + cos12 + sin6 + cos6,
  data = nsw_monthly
)

its_nb_summary <- summary(its_nb)

message("\n  Negative binomial ITS model (NSW monthly):")
nb_coefs <- coef(its_nb_summary)
key_vars <- c("time_index", "covid", "time_since_covid", "post_covid", "time_since_post")
for (v in key_vars) {
  if (v %in% rownames(nb_coefs)) {
    irr <- exp(nb_coefs[v, 1])
    ci_lo <- exp(nb_coefs[v, 1] - 1.96 * nb_coefs[v, 2])
    ci_hi <- exp(nb_coefs[v, 1] + 1.96 * nb_coefs[v, 2])
    sig <- ""
    if (nb_coefs[v, 4] < 0.001) sig <- " ***"
    else if (nb_coefs[v, 4] < 0.01) sig <- " **"
    else if (nb_coefs[v, 4] < 0.05) sig <- " *"
    message("    ", v, ": IRR=", round(irr, 3),
            " (95% CI: ", round(ci_lo, 3), "-", round(ci_hi, 3),
            ", p=", format.pval(nb_coefs[v, 4], digits = 3), ")", sig)
  }
}

message("  AIC: ", round(AIC(its_nb), 1))
message("  Theta (NB dispersion): ", round(its_nb$theta, 2))

# ==============================================================================
# 4. SIMPLIFIED MODEL (MATCHING NATIONAL SPECIFICATION)
# ==============================================================================

message("\n4. Simplified ITS model (matching national specification)...")

# Simpler model: just level change + trend change at COVID onset
its_nb_simple <- glm.nb(
  count ~ time_index + covid + time_since_covid +
    sin12 + cos12,
  data = nsw_monthly
)

its_simple_summary <- summary(its_nb_simple)

message("  Simplified NB ITS model:")
simple_coefs <- coef(its_simple_summary)
for (v in c("time_index", "covid", "time_since_covid")) {
  if (v %in% rownames(simple_coefs)) {
    irr <- exp(simple_coefs[v, 1])
    ci_lo <- exp(simple_coefs[v, 1] - 1.96 * simple_coefs[v, 2])
    ci_hi <- exp(simple_coefs[v, 1] + 1.96 * simple_coefs[v, 2])
    message("    ", v, ": IRR=", round(irr, 3),
            " (95% CI: ", round(ci_lo, 3), "-", round(ci_hi, 3),
            ", p=", format.pval(simple_coefs[v, 4], digits = 3), ")")
  }
}

# Extract key COVID IRR for manuscript
covid_irr <- exp(simple_coefs["covid", 1])
covid_ci_lo <- exp(simple_coefs["covid", 1] - 1.96 * simple_coefs["covid", 2])
covid_ci_hi <- exp(simple_coefs["covid", 1] + 1.96 * simple_coefs["covid", 2])
covid_p <- simple_coefs["covid", 4]

message("\n  Key result for manuscript:")
message("  COVID level change IRR = ", round(covid_irr, 2),
        " (95% CI: ", round(covid_ci_lo, 2), "-", round(covid_ci_hi, 2),
        ", p = ", format.pval(covid_p, digits = 3), ")")
pct_change <- round((covid_irr - 1) * 100, 0)
message("  Interpretation: ~", abs(pct_change), "% ",
        ifelse(pct_change < 0, "reduction", "increase"),
        " in monthly notification rate during COVID restrictions")

# ==============================================================================
# 5. PREDICTED VALUES AND COUNTERFACTUAL
# ==============================================================================

message("\n5. Generating predicted values and counterfactual...")

nsw_monthly$predicted <- predict(its_nb_simple, type = "response")

# Counterfactual: what would have happened without COVID
counterfactual_data <- nsw_monthly %>%
  mutate(covid = 0, time_since_covid = 0)
nsw_monthly$counterfactual <- predict(its_nb_simple,
                                       newdata = counterfactual_data,
                                       type = "response")

# ==============================================================================
# 6. ITS FIGURE
# ==============================================================================

message("\n6. Creating ITS figure...")

fig_nsw_its <- ggplot(nsw_monthly, aes(x = date)) +
  # COVID shading
  annotate("rect",
           xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.1) +
  annotate("text",
           x = as.Date("2020-12-01"),
           y = max(nsw_monthly$count) * 0.95,
           label = "COVID\nrestrictions", size = 3, colour = "#e74c3c") +
  # Observed counts
  geom_col(aes(y = count), fill = "#2c3e50", alpha = 0.4, width = 25) +
  # Counterfactual (no COVID)
  geom_line(aes(y = counterfactual), colour = "#e74c3c",
            linetype = "dashed", linewidth = 0.7) +
  # ITS predicted
  geom_line(aes(y = predicted), colour = "#3498db", linewidth = 0.8) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "NSW monthly infectious syphilis: Interrupted time series",
    subtitle = paste0("2009-2025 (n=204 months) | NB regression with seasonal adjustment\n",
                      "COVID IRR = ", round(covid_irr, 2),
                      " (95% CI: ", round(covid_ci_lo, 2), "-", round(covid_ci_hi, 2),
                      ", p ", ifelse(covid_p < 0.001, "< 0.001",
                                     paste0("= ", round(covid_p, 3))), ")"),
    x = NULL,
    y = "Monthly notifications (n)",
    caption = "Source: NSW Health NCIMS | Blue line = ITS model, red dashed = counterfactual (no COVID)"
  )

save_plot(fig_nsw_its, "fig33_nsw_monthly_its", width = 14, height = 7)

# ==============================================================================
# 7. SEASONAL PATTERN FIGURE
# ==============================================================================

message("\n7. Creating seasonal pattern figure...")

seasonal_summary <- nsw_monthly %>%
  mutate(
    period = case_when(
      year < 2020 ~ "Pre-COVID\n(2009-2019)",
      year <= 2021 ~ "COVID\n(2020-2021)",
      TRUE ~ "Post-COVID\n(2022-2025)"
    ),
    period = factor(period, levels = c("Pre-COVID\n(2009-2019)",
                                        "COVID\n(2020-2021)",
                                        "Post-COVID\n(2022-2025)"))
  ) %>%
  group_by(period, month_num) %>%
  summarise(
    mean_count = mean(count),
    se_count = sd(count) / sqrt(n()),
    .groups = "drop"
  )

fig_seasonal <- ggplot(seasonal_summary,
                        aes(x = month_num, y = mean_count, colour = period)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_count - se_count,
                  ymax = mean_count + se_count,
                  fill = period),
              alpha = 0.1, colour = NA) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb) +
  scale_colour_manual(values = c("#3498db", "#e74c3c", "#2ecc71")) +
  scale_fill_manual(values = c("#3498db", "#e74c3c", "#2ecc71")) +
  labs(
    title = "Seasonal patterns in NSW infectious syphilis notifications",
    subtitle = "Mean monthly notifications by period (with standard error)",
    x = "Month", y = "Mean monthly notifications",
    colour = "Period", fill = "Period",
    caption = "Source: NSW Health NCIMS"
  )

save_plot(fig_seasonal, "fig34_nsw_seasonal_pattern", width = 10, height = 6)

# ==============================================================================
# 8. MODEL COMPARISON TABLE
# ==============================================================================

message("\n8. Creating model comparison table...")

# Compare Poisson vs NB, simple vs full
model_names <- c("Poisson (full)", "NB (full)", "NB (simple)")
model_list <- list(its_pois, its_nb, its_nb_simple)

its_comparison <- tibble(
  model = model_names,
  aic = sapply(model_list, AIC),
  covid_irr = sapply(model_list, function(m) {
    round(exp(coef(m)["covid"]), 3)
  }),
  covid_ci = sapply(model_list, function(m) {
    se <- summary(m)$coefficients["covid", 2]
    lo <- exp(coef(m)["covid"] - 1.96 * se)
    hi <- exp(coef(m)["covid"] + 1.96 * se)
    paste0(round(lo, 3), "-", round(hi, 3))
  }),
  covid_p = sapply(model_list, function(m) {
    format.pval(summary(m)$coefficients["covid", 4], digits = 3)
  }),
  n_params = sapply(model_list, function(m) length(coef(m))),
  n_obs = nrow(nsw_monthly)
)

message("  Model comparison:")
for (i in seq_len(nrow(its_comparison))) {
  message("    ", its_comparison$model[i],
          ": AIC=", round(its_comparison$aic[i], 1),
          ", COVID IRR=", its_comparison$covid_irr[i],
          " (", its_comparison$covid_ci[i], ")",
          ", p=", its_comparison$covid_p[i])
}

# Save as supplementary table
its_comp_gt <- its_comparison %>%
  gt() %>%
  tab_header(
    title = "Supplementary Table S7. NSW monthly ITS model comparison",
    subtitle = "Sensitivity analysis: Model specification and COVID-19 impact estimates"
  ) %>%
  cols_label(
    model = "Model",
    aic = "AIC",
    covid_irr = "COVID IRR",
    covid_ci = "95% CI",
    covid_p = "p-value",
    n_params = "Parameters",
    n_obs = "Observations"
  ) %>%
  fmt_number(columns = aic, decimals = 1) %>%
  tab_source_note("Source: NSW Health NCIMS (2009-2025). NB = negative binomial.") %>%
  tab_source_note("COVID period: March 2020 - December 2021.")

gtsave(its_comp_gt, file.path(path_tables, "tableS7_nsw_monthly_its_comparison.html"))
write_csv(its_comparison, file.path(path_tables, "tableS7_nsw_monthly_its_comparison.csv"))

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================

message("\n9. Saving results...")

saveRDS(nsw_monthly, file.path(path_processed, "nsw_monthly_its_results.rds"))
saveRDS(its_comparison, file.path(path_processed, "nsw_its_model_comparison.rds"))

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

message("\n=== NSW Monthly ITS Summary ===")
message("  Data: ", nrow(nsw_monthly), " monthly observations (2009-2025)")
message("  Pre-COVID: ", sum(nsw_monthly$covid == 0 & nsw_monthly$post_covid == 0), " months")
message("  COVID: ", sum(nsw_monthly$covid == 1), " months")
message("  Post-COVID: ", sum(nsw_monthly$post_covid == 1), " months")
message("  Overdispersion ratio: ", round(disp_ratio, 2))
message("  COVID IRR (NB simple): ", round(covid_irr, 2),
        " (95% CI: ", round(covid_ci_lo, 2), "-", round(covid_ci_hi, 2), ")")
message("  Interpretation: ", abs(pct_change), "% ",
        ifelse(pct_change < 0, "reduction", "increase"),
        " in monthly rate during COVID restrictions")
message("\nFigures: fig33_nsw_monthly_its, fig34_nsw_seasonal_pattern")
message("Tables: tableS7_nsw_monthly_its_comparison")
message("\nNSW Monthly ITS analysis complete.")

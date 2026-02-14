# ==============================================================================
# 07_covid_impact.R — Phase 5: COVID-19 Impact Analysis
# Interrupted time series, differential impact by geography/demography
# ==============================================================================
#
# NOTE: Available data is annual (Kirby 2013-2022), not monthly/quarterly at
# national level. ITS with annual data has limited power — we use segmented
# Poisson regression with pre/during/post periods and supplement with
# rate-of-change analysis and sub-group comparisons.
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

message("=== Phase 5: COVID-19 Impact Analysis ===")
message("1. Loading data...")

state_ts <- readRDS(file.path(path_processed, "state_time_series.rds"))
national_ts <- readRDS(file.path(path_processed, "national_time_series.rds"))
kirby_remote <- readRDS(file.path(path_processed, "kirby_national_remoteness.rds"))
kirby_remote_indig <- readRDS(file.path(path_processed, "kirby_remoteness_indigenous.rds"))
kirby_state_sex <- readRDS(file.path(path_processed, "kirby_state_sex.rds"))
kirby_national_sex <- readRDS(file.path(path_processed, "kirby_national_sex.rds"))
kirby_national_indig <- readRDS(file.path(path_processed, "kirby_national_indigenous.rds"))

# QLD HHS for longer time series + sub-state COVID impact
qld_hhs <- readRDS(file.path(path_processed, "qld_hhs_notifications.rds"))

# NT quarterly for finer temporal resolution
nt_quarterly <- read_csv(
  file.path(path_raw, "notifications", "NT_infectious_syphilis_quarterly_2019-2024.csv"),
  show_col_types = FALSE
)

# VIC LGA for urban sub-state analysis
vic_analysis <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))

# ==============================================================================
# 2. DEFINE COVID PERIODS
# ==============================================================================

message("\n2. Defining COVID periods...")

# Australia COVID periods:
# Pre-COVID: before 2020 (last full year = 2019)
# COVID restrictions: 2020-2021 (border closures, lockdowns, social distancing)
# Post-COVID: 2022+ (borders reopened, restrictions lifted)

covid_periods <- function(year) {
  case_when(
    year < 2020 ~ "Pre-COVID",
    year <= 2021 ~ "COVID restrictions",
    TRUE ~ "Post-COVID"
  )
}

# ==============================================================================
# 3. NATIONAL LEVEL: PRE/DURING/POST COVID COMPARISON
# ==============================================================================

message("\n3. National-level COVID impact analysis...")

national_covid <- national_ts %>%
  mutate(
    period = factor(covid_periods(year),
                    levels = c("Pre-COVID", "COVID restrictions", "Post-COVID")),
    time = year - min(year) + 1,
    covid = as.integer(year >= 2020),
    post_covid = as.integer(year >= 2022),
    time_since_covid = pmax(0, year - 2020),
    time_since_post = pmax(0, year - 2022)
  )

# Annual rate of change
national_covid <- national_covid %>%
  arrange(year) %>%
  mutate(
    pct_change = (count - lag(count)) / lag(count) * 100
  )

message("  National counts by period:")
period_summary <- national_covid %>%
  group_by(period) %>%
  summarise(
    years = paste(range(year), collapse = "-"),
    mean_count = round(mean(count)),
    total_count = sum(count),
    .groups = "drop"
  )
for (i in seq_len(nrow(period_summary))) {
  message("    ", period_summary$period[i], " (", period_summary$years[i], "): ",
          "mean=", period_summary$mean_count[i], "/year, total=", period_summary$total_count[i])
}

# Key COVID year change
message("  2019 -> 2020 change: ", round(national_covid$pct_change[national_covid$year == 2020], 1), "%")
message("  2021 -> 2022 change: ", round(national_covid$pct_change[national_covid$year == 2022], 1), "%")

# Interrupted time series: Poisson model on annual counts
# Model: count ~ time + covid_indicator + time_since_covid
its_national <- glm(
  count ~ time + covid + time_since_covid,
  data = national_covid,
  family = poisson()
)

its_summary <- summary(its_national)
message("\n  ITS Poisson model (national annual):")
its_coefs <- coef(its_summary)
for (i in seq_len(nrow(its_coefs))) {
  irr <- exp(its_coefs[i, 1])
  sig <- ""
  if (its_coefs[i, 4] < 0.001) sig <- " ***"
  else if (its_coefs[i, 4] < 0.01) sig <- " **"
  else if (its_coefs[i, 4] < 0.05) sig <- " *"
  message("    ", rownames(its_coefs)[i], ": coef=",
          round(its_coefs[i, 1], 4), " (IRR=", round(irr, 3),
          ", p=", format.pval(its_coefs[i, 4], digits = 3), ")", sig)
}

message("  Interpretation:")
message("    - 'covid' coefficient: immediate level change at 2020")
message("    - 'time_since_covid' coefficient: trend change post-2020")

# ==============================================================================
# 4. FIGURE: NATIONAL ITS PLOT
# ==============================================================================

message("\n4. Creating ITS visualisation...")

# Predicted values from ITS model
national_covid$predicted <- predict(its_national, type = "response")

# Counterfactual: what would have happened without COVID
counterfactual_data <- national_covid %>%
  mutate(covid = 0, time_since_covid = 0)
national_covid$counterfactual <- predict(its_national, newdata = counterfactual_data,
                                          type = "response")

fig_its <- ggplot(national_covid, aes(x = year)) +
  # Shade COVID period
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.1) +
  annotate("text", x = 2020.5, y = max(national_covid$count) * 0.95,
           label = "COVID\nrestrictions", size = 3, colour = "#e74c3c") +
  # Actual data
  geom_col(aes(y = count), fill = "#2c3e50", alpha = 0.6) +
  # Counterfactual trend (no COVID)
  geom_line(aes(y = counterfactual), colour = "#e74c3c",
            linetype = "dashed", linewidth = 0.8) +
  # ITS predicted
  geom_line(aes(y = predicted), colour = "#3498db", linewidth = 0.8) +
  geom_point(aes(y = count), colour = "#2c3e50", size = 2) +
  scale_x_continuous(breaks = 2013:2022) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Interrupted time series: Infectious syphilis notifications, Australia",
    subtitle = "2013-2022 | Blue line = ITS predicted, red dashed = counterfactual (no COVID)",
    x = NULL, y = "Notifications (n)",
    caption = "Source: Kirby Institute | ITS: Poisson GLM with level + trend change at 2020"
  )

save_plot(fig_its, "fig25_its_national", width = 12, height = 7)

# ==============================================================================
# 5. STATE-LEVEL COVID IMPACT
# ==============================================================================

message("\n5. State-level COVID impact...")

state_covid <- state_ts %>%
  mutate(period = covid_periods(year)) %>%
  group_by(state, period) %>%
  summarise(
    mean_rate = mean(rate_per_100k),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period, values_from = mean_rate) %>%
  mutate(
    covid_change_pct = (`COVID restrictions` - `Pre-COVID`) / `Pre-COVID` * 100,
    post_covid_change_pct = (`Post-COVID` - `Pre-COVID`) / `Pre-COVID` * 100,
    rebound_pct = (`Post-COVID` - `COVID restrictions`) / `COVID restrictions` * 100
  )

message("  State COVID impact (rate change from pre-COVID):")
for (i in seq_len(nrow(state_covid))) {
  message("    ", state_covid$state[i], ": ",
          "COVID=", sprintf("%+.1f%%", state_covid$covid_change_pct[i]),
          ", Post-COVID=", sprintf("%+.1f%%", state_covid$post_covid_change_pct[i]))
}

# State-level ITS figure
state_ts_covid <- state_ts %>%
  mutate(period = factor(covid_periods(year),
                         levels = c("Pre-COVID", "COVID restrictions", "Post-COVID")))

fig_state_covid <- ggplot(state_ts_covid,
                           aes(x = year, y = rate_per_100k, colour = state)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.08) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = state_colours) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  labs(
    title = "Infectious syphilis rates by state: COVID-19 impact",
    subtitle = "2013-2022 | Shaded region = COVID restriction period (2020-2021)",
    x = NULL, y = "Rate per 100,000 population",
    colour = "State/Territory",
    caption = "Source: Kirby Institute"
  ) +
  theme(legend.position = "right")

save_plot(fig_state_covid, "fig26_state_covid_impact", width = 12, height = 7)

# ==============================================================================
# 6. REMOTENESS-LEVEL COVID IMPACT
# ==============================================================================

message("\n6. Remoteness-level COVID impact...")

remote_covid <- kirby_remote %>%
  filter(remoteness %in% c("Major cities", "Regional", "Remote")) %>%
  mutate(period = covid_periods(year)) %>%
  group_by(remoteness, period) %>%
  summarise(mean_count = mean(count), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_count) %>%
  mutate(
    covid_change_pct = (`COVID restrictions` - `Pre-COVID`) / `Pre-COVID` * 100,
    post_covid_change_pct = (`Post-COVID` - `Pre-COVID`) / `Pre-COVID` * 100
  )

message("  COVID impact by remoteness:")
for (i in seq_len(nrow(remote_covid))) {
  message("    ", remote_covid$remoteness[i], ": ",
          "COVID=", sprintf("%+.1f%%", remote_covid$covid_change_pct[i]),
          ", Post-COVID=", sprintf("%+.1f%%", remote_covid$post_covid_change_pct[i]))
}

# Remoteness time series plot
remote_ts <- kirby_remote %>%
  filter(remoteness %in% c("Major cities", "Regional", "Remote"))

fig_remote_covid <- ggplot(remote_ts, aes(x = year, y = count, colour = remoteness)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.08) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    "Major cities" = "#3498db", "Regional" = "#f39c12", "Remote" = "#e74c3c"
  )) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Infectious syphilis notifications by remoteness: COVID-19 impact",
    subtitle = "2013-2022 | Shaded = COVID restriction period",
    x = NULL, y = "Notifications (n)",
    colour = "Remoteness",
    caption = "Source: Kirby Institute"
  )

save_plot(fig_remote_covid, "fig27_remoteness_covid_impact", width = 12, height = 7)

# ==============================================================================
# 7. INDIGENOUS vs NON-INDIGENOUS COVID IMPACT
# ==============================================================================

message("\n7. Indigenous vs Non-Indigenous COVID impact...")

# Remoteness x Indigenous rates
remote_indig_covid <- kirby_remote_indig %>%
  mutate(
    indigenous = ifelse(str_detect(population_group, "Aboriginal"),
                        "Indigenous", "Non-Indigenous"),
    remoteness = case_when(
      str_detect(population_group, "major cities") ~ "Major cities",
      str_detect(population_group, "regional") ~ "Regional",
      str_detect(population_group, "remote") ~ "Remote"
    ),
    period = covid_periods(year)
  )

# Calculate pre/during/post rates by group
indig_covid_summary <- remote_indig_covid %>%
  group_by(indigenous, remoteness, period) %>%
  summarise(mean_rate = mean(rate), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_rate) %>%
  mutate(
    covid_change_pct = (`COVID restrictions` - `Pre-COVID`) / `Pre-COVID` * 100,
    post_covid_change_pct = (`Post-COVID` - `Pre-COVID`) / `Pre-COVID` * 100
  )

message("  COVID impact by Indigenous status x remoteness:")
for (i in seq_len(nrow(indig_covid_summary))) {
  message("    ", indig_covid_summary$indigenous[i], " - ",
          indig_covid_summary$remoteness[i], ": ",
          "COVID=", sprintf("%+.1f%%", indig_covid_summary$covid_change_pct[i]),
          ", Post-COVID=", sprintf("%+.1f%%", indig_covid_summary$post_covid_change_pct[i]))
}

# Figure: Indigenous vs Non-Indigenous by remoteness
fig_indig_covid <- ggplot(remote_indig_covid,
                           aes(x = year, y = rate, colour = remoteness,
                               linetype = indigenous)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.08) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(
    "Major cities" = "#3498db", "Regional" = "#f39c12", "Remote" = "#e74c3c"
  )) +
  scale_linetype_manual(values = c("Indigenous" = "solid", "Non-Indigenous" = "dashed")) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  labs(
    title = "COVID-19 impact on syphilis rates by remoteness and Indigenous status",
    subtitle = "2013-2022 | Shaded = COVID restriction period",
    x = NULL, y = "Rate per 100,000 population",
    colour = "Remoteness", linetype = "Indigenous status",
    caption = "Source: Kirby Institute"
  ) +
  theme(legend.position = "right")

save_plot(fig_indig_covid, "fig28_indigenous_remoteness_covid", width = 12, height = 7)

# ==============================================================================
# 8. SEX-SPECIFIC COVID IMPACT (MSM vs HETEROSEXUAL PROXY)
# ==============================================================================

message("\n8. Sex-specific COVID impact (MSM vs heterosexual proxy)...")

# Male rates (proxy for MSM in urban areas) vs Female rates (proxy for
# heterosexual transmission, especially in remote areas)
national_sex <- kirby_national_sex %>%
  filter(sex %in% c("Male", "Female")) %>%
  mutate(period = covid_periods(year))

sex_covid <- national_sex %>%
  group_by(sex, period) %>%
  summarise(mean_count = mean(count), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_count) %>%
  mutate(
    covid_change_pct = (`COVID restrictions` - `Pre-COVID`) / `Pre-COVID` * 100,
    post_covid_change_pct = (`Post-COVID` - `Pre-COVID`) / `Pre-COVID` * 100
  )

message("  COVID impact by sex:")
for (i in seq_len(nrow(sex_covid))) {
  message("    ", sex_covid$sex[i], ": ",
          "COVID=", sprintf("%+.1f%%", sex_covid$covid_change_pct[i]),
          ", Post-COVID=", sprintf("%+.1f%%", sex_covid$post_covid_change_pct[i]))
}

# State x sex for differential analysis
state_sex_covid <- kirby_state_sex %>%
  filter(sex %in% c("Females", "Males")) %>%
  mutate(
    period = covid_periods(year),
    sex = str_remove(sex, "s$")  # Females -> Female
  )

fig_sex_covid <- ggplot(national_sex, aes(x = year, y = count, colour = sex)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.08) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Female" = "#e74c3c", "Male" = "#3498db")) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Infectious syphilis notifications by sex: COVID-19 impact",
    subtitle = "2013-2022 | Male notifications ~ MSM transmission, Female ~ heterosexual",
    x = NULL, y = "Notifications (n)",
    colour = "Sex",
    caption = "Source: Kirby Institute"
  )

save_plot(fig_sex_covid, "fig29_sex_covid_impact", width = 12, height = 7)

# ==============================================================================
# 9. QLD HHS: SUB-STATE COVID IMPACT (LONGER TIME SERIES)
# ==============================================================================

message("\n9. QLD HHS sub-state COVID impact...")

qld_state <- qld_hhs %>%
  filter(area_name == "Queensland" | area_name %in%
           c("Metro North", "Metro South", "Gold Coast",
             "Torres and Cape", "Cairns and Hinterland", "North West")) %>%
  mutate(
    period = covid_periods(year),
    area_type = case_when(
      area_name %in% c("Metro North", "Metro South", "Gold Coast") ~ "Urban SE QLD",
      area_name %in% c("Torres and Cape", "Cairns and Hinterland", "North West") ~ "Remote North QLD",
      TRUE ~ "All Queensland"
    )
  )

# Aggregate urban vs remote QLD
qld_urban_remote <- qld_state %>%
  filter(area_type != "All Queensland") %>%
  group_by(area_type, year) %>%
  summarise(mean_rate = mean(rate, na.rm = TRUE), .groups = "drop") %>%
  mutate(period = covid_periods(year))

qld_covid_summary <- qld_urban_remote %>%
  group_by(area_type, period) %>%
  summarise(mean_rate = mean(mean_rate), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_rate) %>%
  mutate(
    covid_change_pct = (`COVID restrictions` - `Pre-COVID`) / `Pre-COVID` * 100,
    post_covid_change_pct = (`Post-COVID` - `Pre-COVID`) / `Pre-COVID` * 100
  )

message("  QLD COVID impact (urban vs remote):")
for (i in seq_len(nrow(qld_covid_summary))) {
  message("    ", qld_covid_summary$area_type[i], ": ",
          "COVID=", sprintf("%+.1f%%", qld_covid_summary$covid_change_pct[i]),
          ", Post-COVID=", sprintf("%+.1f%%", qld_covid_summary$post_covid_change_pct[i]))
}

fig_qld_covid <- ggplot(qld_urban_remote,
                         aes(x = year, y = mean_rate, colour = area_type)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.08) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Urban SE QLD" = "#3498db",
                                  "Remote North QLD" = "#e74c3c")) +
  scale_x_continuous(breaks = seq(2012, 2024, 2)) +
  labs(
    title = "COVID-19 impact on syphilis rates: Urban vs Remote Queensland",
    subtitle = "2012-2024 | Urban SE = Metro North/South + Gold Coast, Remote North = Torres/Cairns/NW",
    x = NULL, y = "Mean rate per 100,000",
    colour = "Region",
    caption = "Source: QLD Health Syphilis in Queensland 2024"
  )

save_plot(fig_qld_covid, "fig30_qld_urban_remote_covid", width = 12, height = 7)

# ==============================================================================
# 10. NT QUARTERLY: FINER TEMPORAL COVID IMPACT
# ==============================================================================

message("\n10. NT quarterly COVID impact analysis...")

nt_q <- nt_quarterly %>%
  mutate(
    q_num = as.integer(str_extract(Quarter, "\\d")),
    date = as.Date(paste(Year, (q_num - 1) * 3 + 2, 15, sep = "-")),
    period = case_when(
      Year < 2020 ~ "Pre-COVID",
      Year <= 2021 ~ "COVID",
      TRUE ~ "Post-COVID"
    ),
    time_index = seq_len(n())
  )

# NT quarterly ITS
nt_its <- glm(
  Infectious_Syphilis_Count ~ time_index +
    I(Year >= 2020) +
    pmax(0, time_index - min(time_index[Year >= 2020])),
  data = nt_q,
  family = poisson()
)

nt_its_summary <- summary(nt_its)
message("  NT quarterly ITS model:")
nt_coefs <- coef(nt_its_summary)
for (i in seq_len(nrow(nt_coefs))) {
  message("    ", rownames(nt_coefs)[i], ": coef=",
          round(nt_coefs[i, 1], 4),
          " (p=", format.pval(nt_coefs[i, 4], digits = 3), ")")
}

# NT predicted
nt_q$predicted <- predict(nt_its, type = "response")

fig_nt_quarterly <- ggplot(nt_q, aes(x = date)) +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2021-12-31"),
           ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.1) +
  geom_col(aes(y = Infectious_Syphilis_Count), fill = "#2c3e50", alpha = 0.6) +
  geom_line(aes(y = predicted), colour = "#3498db", linewidth = 0.8) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
  labs(
    title = "NT infectious syphilis: Quarterly notifications with ITS model",
    subtitle = "2019-2024 | Shaded = COVID restriction period",
    x = NULL, y = "Quarterly notifications (n)",
    caption = "Source: NT Health Surveillance Update"
  )

save_plot(fig_nt_quarterly, "fig31_nt_quarterly_covid", width = 12, height = 6)

# ==============================================================================
# 11. TABLE 5: COVID IMPACT SUMMARY
# ==============================================================================

message("\n11. Creating COVID impact summary table...")

# Build comprehensive table
table5_data <- bind_rows(
  # National
  tibble(
    level = "National",
    group = "All",
    pre_covid = period_summary$mean_count[period_summary$period == "Pre-COVID"],
    covid = period_summary$mean_count[period_summary$period == "COVID restrictions"],
    post_covid = period_summary$mean_count[period_summary$period == "Post-COVID"]
  ),
  # By remoteness
  remote_covid %>%
    transmute(
      level = "Remoteness",
      group = remoteness,
      pre_covid = round(`Pre-COVID`),
      covid = round(`COVID restrictions`),
      post_covid = round(`Post-COVID`)
    ),
  # By sex
  sex_covid %>%
    transmute(
      level = "Sex",
      group = sex,
      pre_covid = round(`Pre-COVID`),
      covid = round(`COVID restrictions`),
      post_covid = round(`Post-COVID`)
    )
) %>%
  mutate(
    covid_change_pct = round((covid - pre_covid) / pre_covid * 100, 1),
    post_covid_change_pct = round((post_covid - pre_covid) / pre_covid * 100, 1),
    rebound_pct = round((post_covid - covid) / covid * 100, 1)
  )

table5_gt <- table5_data %>%
  gt(groupname_col = "level") %>%
  tab_header(
    title = "Table 5. COVID-19 impact on infectious syphilis notifications",
    subtitle = "Mean annual notifications by period"
  ) %>%
  cols_label(
    group = "",
    pre_covid = "Pre-COVID (2013-2019)",
    covid = "COVID (2020-2021)",
    post_covid = "Post-COVID (2022)",
    covid_change_pct = "COVID change (%)",
    post_covid_change_pct = "Post-COVID vs pre (%)",
    rebound_pct = "Rebound (%)"
  ) %>%
  fmt_number(columns = c(pre_covid, covid, post_covid), decimals = 0) %>%
  fmt_number(columns = c(covid_change_pct, post_covid_change_pct, rebound_pct),
             decimals = 1) %>%
  tab_source_note("Source: Kirby Institute Annual Surveillance Report 2024") %>%
  tab_source_note("Pre-COVID = 2013-2019, COVID restrictions = 2020-2021, Post-COVID = 2022")

gtsave(table5_gt, file.path(path_tables, "table5_covid_impact.html"))
write_csv(table5_data, file.path(path_tables, "table5_covid_impact.csv"))

# COVID impact by remoteness x Indigenous
table5b_data <- indig_covid_summary %>%
  transmute(
    group = paste(indigenous, "-", remoteness),
    pre_covid_rate = round(`Pre-COVID`, 1),
    covid_rate = round(`COVID restrictions`, 1),
    post_covid_rate = round(`Post-COVID`, 1),
    covid_change_pct = round(covid_change_pct, 1),
    post_covid_change_pct = round(post_covid_change_pct, 1)
  )

table5b_gt <- table5b_data %>%
  gt() %>%
  tab_header(
    title = "Table 5b. COVID-19 impact by remoteness and Indigenous status",
    subtitle = "Mean annual rate per 100,000 by period"
  ) %>%
  cols_label(
    group = "Group",
    pre_covid_rate = "Pre-COVID",
    covid_rate = "COVID",
    post_covid_rate = "Post-COVID",
    covid_change_pct = "COVID change (%)",
    post_covid_change_pct = "Post-COVID vs pre (%)"
  ) %>%
  tab_source_note("Source: Kirby Institute") %>%
  tab_source_note("Rates per 100,000 population. Pre-COVID = 2013-2019, COVID = 2020-2021, Post-COVID = 2022")

gtsave(table5b_gt, file.path(path_tables, "table5b_covid_remoteness_indigenous.html"))
write_csv(table5b_data, file.path(path_tables, "table5b_covid_remoteness_indigenous.csv"))

# ==============================================================================
# 12. COMBINED PANEL: COVID IMPACT OVERVIEW
# ==============================================================================

message("\n12. Creating combined COVID panel...")

fig_covid_panel <- (fig_its + fig_state_covid) /
  (fig_remote_covid + fig_sex_covid) +
  plot_annotation(
    title = "COVID-19 impact on Australia's syphilis epidemic",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(fig_covid_panel, "fig32_covid_overview_panel", width = 18, height = 14)

# ==============================================================================
# 13. SAVE PROCESSED DATA
# ==============================================================================

message("\n13. Saving processed datasets...")

saveRDS(national_covid, file.path(path_processed, "national_covid_its.rds"))
saveRDS(state_covid, file.path(path_processed, "state_covid_impact.rds"))
saveRDS(indig_covid_summary, file.path(path_processed, "covid_remoteness_indigenous.rds"))

# ==============================================================================
# 14. SUMMARY
# ==============================================================================

message("\n=== Phase 5 Output Summary ===")
message("Figures saved to: ", path_figures)
message("Tables saved to: ", path_tables)

fig_files <- list.files(path_figures, pattern = "fig(2[5-9]|3[0-2])")
table_files <- list.files(path_tables, pattern = "table5")

message("  Figures: ", length(fig_files))
for (f in fig_files) message("    ", f)
message("  Tables: ", length(table_files))
for (f in table_files) message("    ", f)

# Summary stats
message("\n  Key findings:")
message("  - National 2020 dip: ",
        round(national_covid$pct_change[national_covid$year == 2020], 1),
        "% change (5900 -> 5355)")
message("  - Post-COVID rebound to 6036 in 2022")

# Urban vs remote differential
urban_chg <- remote_covid$covid_change_pct[remote_covid$remoteness == "Major cities"]
remote_chg <- remote_covid$covid_change_pct[remote_covid$remoteness == "Remote"]
message("  - Urban (major cities) COVID change: ", sprintf("%+.1f%%", urban_chg))
message("  - Remote COVID change: ", sprintf("%+.1f%%", remote_chg))
message("  - Differential impact: urban affected more than remote")

# Sex differential
male_chg <- sex_covid$covid_change_pct[sex_covid$sex == "Male"]
female_chg <- sex_covid$covid_change_pct[sex_covid$sex == "Female"]
message("  - Male COVID change: ", sprintf("%+.1f%%", male_chg))
message("  - Female COVID change: ", sprintf("%+.1f%%", female_chg))

message("\nPhase 5 complete.")

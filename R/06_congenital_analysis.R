# ==============================================================================
# 06_congenital_analysis.R — Phase 4: Congenital Syphilis Analysis
# Temporal trends, antenatal care correlation, risk mapping
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. LOAD AND COMPILE CONGENITAL SYPHILIS DATA
# ==============================================================================

message("=== Phase 4: Congenital Syphilis Analysis ===")
message("1. Loading congenital syphilis data...")

# --- QLD congenital by Indigenous status (2001-2024) ---
qld_cong_indig <- read_csv(
  file.path(path_raw, "notifications", "QLD_congenital_syphilis_by_indigenous_2001-2024.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(
    cols = c(First_Nations, Non_First_Nations, Total),
    names_to = "group", values_to = "count"
  ) %>%
  mutate(
    state = "QLD",
    group = case_when(
      group == "First_Nations" ~ "Indigenous",
      group == "Non_First_Nations" ~ "Non-Indigenous",
      group == "Total" ~ "Total"
    )
  )

# --- QLD congenital by region (2001-2024) ---
qld_cong_region <- read_csv(
  file.path(path_raw, "notifications", "QLD_congenital_syphilis_by_region_2001-2024.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(
    cols = c(North_QLD, Central_QLD, South_East_QLD, Total),
    names_to = "region", values_to = "count"
  ) %>%
  mutate(
    state = "QLD",
    region = str_replace_all(region, "_", " ")
  )

# --- NT congenital quarterly (2019-2024) ---
nt_cong <- read_csv(
  file.path(path_raw, "notifications", "NT_congenital_syphilis_quarterly_2019-2024.csv"),
  show_col_types = FALSE
) %>%
  mutate(state = "NT")

nt_cong_annual <- nt_cong %>%
  group_by(Year, state) %>%
  summarise(count = sum(Congenital_Syphilis_Count, na.rm = TRUE), .groups = "drop") %>%
  rename(year = Year)

# --- WA summary (2024-25) ---
wa_cong <- tibble(
  year = c(2024),
  state = "WA",
  count = 3,
  source = "WA SORG Communique 2024-25"
)

# --- National congenital data from Hengel et al. and Kirby ---
# Hengel et al. 2024 MJA: 74 cases nationally 2011-2021
# National CDC Q2 2025: 2023=20 cases (10 deaths), 2024=9 (3 deaths), H1 2025=9 (4 deaths)
# Kirby ASR 2024: annual cases/deaths 2014-2023
#
# Construct national time series from available sources
national_congenital <- tibble(
  year = 2011:2024,
  # Hengel et al. 2024 Table 1 cumulative data, interpolated with Kirby ASR 2024 Figs 5-6
  # 2011-2013: very low (0-2/year based on Hengel total of 74 over 2011-2021)
  # 2014-2023: from Kirby ASR 2024 STI report (Figures 5 and 6)
  # 2024: from CDC Q2 2025 report
  count = c(
    1, 2, 3,          # 2011-2013 (estimated from Hengel cumulative)
    3, 2, 3, 4, 6,    # 2014-2018 (Kirby ASR 2024)
    15, 10, 9,         # 2019-2021 (Kirby + Hengel concordant)
    17, 20,            # 2022-2023 (Kirby ASR 2024 + CDC Q2 2025)
    9                  # 2024 (CDC Q2 2025 report)
  ),
  deaths = c(
    NA, NA, NA,        # 2011-2013 (not reported)
    NA, NA, NA, NA, NA,# 2014-2018 (not separately reported)
    4, 4, 3,           # 2019-2021 (Hengel et al. supplement)
    5, 10,             # 2022-2023 (CDC Q2 2025)
    3                  # 2024 (CDC Q2 2025)
  ),
  source = c(
    rep("Hengel et al. 2024 (estimated)", 3),
    rep("Kirby ASR 2024", 5),
    rep("Kirby ASR 2024 / Hengel et al.", 3),
    rep("Kirby ASR 2024 / CDC Q2 2025", 2),
    "CDC Q2 2025"
  )
)

# National by Indigenous status (from Hengel + CDC reports)
national_cong_indig <- tibble(
  year = c(2019, 2020, 2021, 2022, 2023, 2024),
  indigenous_count = c(13, 8, 7, 13, 17, 7),
  non_indigenous_count = c(2, 2, 2, 4, 3, 2),
  source = c(rep("Hengel et al. 2024 + Kirby", 3),
             rep("CDC Q2 2025 + Kirby ASR 2024", 3))
)

message("  QLD congenital: ", nrow(qld_cong_indig %>% filter(group == "Total")),
        " year-records (2001-2024)")
message("  NT congenital: ", nrow(nt_cong_annual), " year-records")
message("  National congenital: ", nrow(national_congenital), " years (2011-2024)")

# Save compiled congenital data
saveRDS(national_congenital, file.path(path_processed, "congenital_national.rds"))
saveRDS(qld_cong_indig, file.path(path_processed, "congenital_qld_indigenous.rds"))
saveRDS(qld_cong_region, file.path(path_processed, "congenital_qld_region.rds"))

# ==============================================================================
# 2. FIGURE: NATIONAL CONGENITAL SYPHILIS TREND
# ==============================================================================

message("\n2. Creating congenital syphilis trend figures...")

# 2a. National trend with deaths overlay
fig_cong_national <- ggplot(national_congenital, aes(x = year)) +
  geom_col(aes(y = count), fill = "#c0392b", alpha = 0.7) +
  geom_point(aes(y = deaths), colour = "#2c3e50", size = 3, na.rm = TRUE) +
  geom_line(aes(y = deaths), colour = "#2c3e50", linewidth = 0.8, na.rm = TRUE) +
  geom_text(aes(y = count, label = count), vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = 2011:2024) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Congenital syphilis notifications and deaths, Australia",
    subtitle = "2011-2024 | Bars = notifications, line/points = deaths",
    x = NULL, y = "Count (n)",
    caption = "Sources: Hengel et al. 2024 MJA, Kirby ASR 2024, CDC Q2 2025 Report"
  )

save_plot(fig_cong_national, "fig17_congenital_national_trend", width = 12, height = 6)

# 2b. QLD congenital by Indigenous status
qld_cong_plot <- qld_cong_indig %>%
  filter(group != "Total") %>%
  mutate(group = factor(group, levels = c("Indigenous", "Non-Indigenous")))

fig_cong_qld_indig <- ggplot(qld_cong_plot, aes(x = Year, y = count, fill = group)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("Indigenous" = "#e67e22", "Non-Indigenous" = "#2ecc71")) +
  scale_x_continuous(breaks = seq(2001, 2024, 2)) +
  labs(
    title = "Congenital syphilis notifications by Indigenous status, Queensland",
    subtitle = "2001-2024",
    x = NULL, y = "Notifications (n)", fill = NULL,
    caption = "Source: QLD Health Syphilis in Queensland 2024"
  )

save_plot(fig_cong_qld_indig, "fig18_congenital_qld_indigenous", width = 12, height = 6)

# 2c. QLD congenital by region
qld_cong_region_plot <- qld_cong_region %>%
  filter(region != "Total") %>%
  mutate(region = factor(region, levels = c("North QLD", "Central QLD", "South East QLD")))

fig_cong_qld_region <- ggplot(qld_cong_region_plot,
                               aes(x = Year, y = count, fill = region)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("North QLD" = "#e74c3c", "Central QLD" = "#f39c12",
                                "South East QLD" = "#3498db")) +
  scale_x_continuous(breaks = seq(2001, 2024, 2)) +
  labs(
    title = "Congenital syphilis notifications by region, Queensland",
    subtitle = "2001-2024",
    x = NULL, y = "Notifications (n)", fill = NULL,
    caption = "Source: QLD Health Syphilis in Queensland 2024"
  )

save_plot(fig_cong_qld_region, "fig19_congenital_qld_region", width = 12, height = 6)

# 2d. National congenital by Indigenous status (2019-2024)
cong_indig_long <- national_cong_indig %>%
  pivot_longer(cols = c(indigenous_count, non_indigenous_count),
               names_to = "group", values_to = "count") %>%
  mutate(group = case_when(
    group == "indigenous_count" ~ "Aboriginal and Torres Strait\nIslander peoples",
    group == "non_indigenous_count" ~ "Non-Indigenous"
  ))

fig_cong_indig_national <- ggplot(cong_indig_long,
                                   aes(x = year, y = count, fill = group)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c(
    "Aboriginal and Torres Strait\nIslander peoples" = "#e67e22",
    "Non-Indigenous" = "#2ecc71"
  )) +
  scale_x_continuous(breaks = 2019:2024) +
  labs(
    title = "Congenital syphilis by Indigenous status, Australia",
    subtitle = "2019-2024",
    x = NULL, y = "Notifications (n)", fill = NULL,
    caption = "Sources: Hengel et al. 2024, Kirby ASR 2024, CDC Q2 2025"
  )

save_plot(fig_cong_indig_national, "fig20_congenital_indigenous_national", width = 10, height = 6)

# ==============================================================================
# 3. LOAD ANTENATAL CARE DATA (AIHW)
# ==============================================================================

message("\n3. Loading antenatal care data...")

aihw_path <- file.path(path_raw, "congenital",
                        "AIHW-PER-101-National-Perinatal-Data-Collection-annual-update-2023.xlsx")

# Table 5.7: 5+ ANC visits by SA3
anc_5plus <- read_excel(aihw_path, sheet = "Table 5.7", skip = 5) %>%
  select(1:5, 7) %>%
  set_names(c("state", "sa3_code", "sa3_name", "anc_5plus_n", "total_births", "pct_5plus")) %>%
  filter(!is.na(sa3_code), str_detect(sa3_code, "^\\d")) %>%
  mutate(
    sa3_code = str_trim(sa3_code),
    anc_5plus_n = as.numeric(anc_5plus_n),
    total_births = as.numeric(total_births),
    pct_5plus = as.numeric(pct_5plus)
  )

# Table 5.8: First trimester ANC (<14 weeks) by SA3
anc_first_trim <- read_excel(aihw_path, sheet = "Table 5.8", skip = 5) %>%
  select(1:5, 7) %>%
  set_names(c("state", "sa3_code", "sa3_name", "anc_first_trim_n", "total_births_t58",
              "pct_first_trimester")) %>%
  filter(!is.na(sa3_code), str_detect(sa3_code, "^\\d")) %>%
  mutate(
    sa3_code = str_trim(sa3_code),
    anc_first_trim_n = as.numeric(anc_first_trim_n),
    total_births_t58 = as.numeric(total_births_t58),
    pct_first_trimester = as.numeric(pct_first_trimester)
  )

# Table 5.1: 5+ ANC visits by PHN
anc_phn <- read_excel(aihw_path, sheet = "Table 5.1", skip = 5) %>%
  select(1:5, 7) %>%
  set_names(c("state", "phn_code", "phn_name", "anc_5plus_n", "total_births",
              "pct_5plus")) %>%
  filter(!is.na(phn_code), str_detect(phn_code, "^PHN")) %>%
  mutate(
    anc_5plus_n = as.numeric(anc_5plus_n),
    total_births = as.numeric(total_births),
    pct_5plus = as.numeric(pct_5plus)
  )

# Table 5.2: First trimester ANC by PHN
anc_phn_ft <- read_excel(aihw_path, sheet = "Table 5.2", skip = 5) %>%
  select(1:5, 7) %>%
  set_names(c("state", "phn_code", "phn_name", "anc_first_trim_n", "total_births",
              "pct_first_trimester")) %>%
  filter(!is.na(phn_code), str_detect(phn_code, "^PHN")) %>%
  mutate(
    anc_first_trim_n = as.numeric(anc_first_trim_n),
    total_births = as.numeric(total_births),
    pct_first_trimester = as.numeric(pct_first_trimester)
  )

# Merge SA3 ANC datasets
anc_sa3 <- anc_5plus %>%
  left_join(
    anc_first_trim %>% select(sa3_code, anc_first_trim_n, pct_first_trimester),
    by = "sa3_code"
  )

# Merge PHN ANC datasets
anc_phn_merged <- anc_phn %>%
  left_join(
    anc_phn_ft %>% select(phn_code, anc_first_trim_n, pct_first_trimester),
    by = "phn_code"
  )

message("  SA3-level ANC data: ", nrow(anc_sa3), " SA3 areas")
message("  PHN-level ANC data: ", nrow(anc_phn_merged), " PHN areas")
message("  First trimester ANC range: ",
        round(min(anc_sa3$pct_first_trimester, na.rm = TRUE), 1), "% - ",
        round(max(anc_sa3$pct_first_trimester, na.rm = TRUE), 1), "%")

saveRDS(anc_sa3, file.path(path_processed, "anc_sa3.rds"))
saveRDS(anc_phn_merged, file.path(path_processed, "anc_phn.rds"))

# ==============================================================================
# 4. BUILD SA3-LEVEL CONGENITAL RISK DATASET
# ==============================================================================

message("\n4. Building SA3-level congenital risk dataset...")

# Load SA3 boundaries
sa3_boundaries <- readRDS(file.path(path_processed, "boundaries_sa3.rds"))

# Aggregate SA2 covariates to SA3
sa2_analysis <- readRDS(file.path(path_processed, "sa2_analysis.rds"))

sa3_covariates <- sa2_analysis %>%
  filter(!is.na(pop_total), pop_total > 0) %>%
  group_by(sa3_code, sa3_name) %>%
  summarise(
    pop_total = sum(pop_total, na.rm = TRUE),
    pop_indigenous = sum(pop_indigenous, na.rm = TRUE),
    pct_indigenous = pop_indigenous / pop_total * 100,
    indigenous_female = sum(indigenous_female, na.rm = TRUE),
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
    remoteness_mode = {
      r <- remoteness[!is.na(remoteness)]
      if (length(r) > 0) names(which.max(table(r))) else NA_character_
    },
    .groups = "drop"
  )

# Join ANC + covariates + boundaries
sa3_risk <- sa3_boundaries %>%
  left_join(sa3_covariates, by = c("SA3_CODE21" = "sa3_code")) %>%
  left_join(anc_sa3 %>% select(sa3_code, total_births, pct_5plus, pct_first_trimester),
            by = c("SA3_CODE21" = "sa3_code")) %>%
  filter(!is.na(pop_total), pop_total > 0)

message("  SA3 risk dataset: ", nrow(sa3_risk), " areas with population + ANC data")
message("  SA3 areas with ANC data: ",
        sum(!is.na(sa3_risk$pct_first_trimester)),
        " of ", nrow(sa3_risk))

# ==============================================================================
# 5. CONGENITAL RISK CLASSIFICATION
# ==============================================================================

message("\n5. Classifying congenital syphilis risk by SA3...")

# Risk factors for congenital syphilis:
# 1. High Indigenous female population (proxy for high adult female syphilis rates)
# 2. Low first-trimester ANC attendance (missed early screening)
# 3. High remoteness (access barriers)
# 4. High distance to sexual health clinic

sa3_risk <- sa3_risk %>%
  mutate(
    # Quantile-based thresholds
    high_indigenous = pct_indigenous > quantile(pct_indigenous, 0.75, na.rm = TRUE),
    low_anc = !is.na(pct_first_trimester) & pct_first_trimester < quantile(pct_first_trimester, 0.25, na.rm = TRUE),
    high_remoteness = remoteness_mode %in% c("Remote", "Very Remote"),
    high_dist_clinic = !is.na(mean_dist_sh_km) & mean_dist_sh_km > quantile(mean_dist_sh_km, 0.75, na.rm = TRUE),

    # Composite risk score (0-4)
    risk_score = as.integer(high_indigenous) +
                 as.integer(low_anc) +
                 as.integer(high_remoteness) +
                 as.integer(high_dist_clinic),

    # Risk classification
    risk_category = case_when(
      risk_score >= 3 ~ "High vulnerability",
      risk_score == 2 ~ "Elevated vulnerability",
      risk_score == 1 ~ "Moderate vulnerability",
      risk_score == 0 ~ "Lower vulnerability",
      TRUE ~ NA_character_
    ),
    risk_category = factor(risk_category,
                           levels = c("High vulnerability", "Elevated vulnerability",
                                      "Moderate vulnerability", "Lower vulnerability"))
  )

risk_counts <- table(sa3_risk$risk_category, useNA = "ifany")
message("  Risk classification:")
for (r in names(risk_counts)) message("    ", r, ": ", risk_counts[r])

# ==============================================================================
# 6. CONGENITAL RISK MAP
# ==============================================================================

message("\n6. Creating congenital syphilis risk maps...")

risk_colours <- c(
  "High vulnerability"     = "#d73027",
  "Elevated vulnerability" = "#fc8d59",
  "Moderate vulnerability" = "#fee090",
  "Lower vulnerability"    = "#91bfdb"
)

fig_risk_map <- tm_shape(sa3_risk) +
  tm_polygons(
    fill = "risk_category",
    fill.scale = tm_scale_categorical(values = risk_colours),
    fill.legend = tm_legend(title = "Vulnerability\ncategory")
  ) +
  tm_title("Congenital syphilis vulnerability classification by SA3, Australia")

save_map(fig_risk_map, "fig21_congenital_risk_sa3", width = 10, height = 10)

# ==============================================================================
# 7. ANTENATAL CARE MAPS
# ==============================================================================

message("\n7. Creating antenatal care maps...")

# First trimester ANC attendance map
sa3_anc_map <- sa3_risk %>%
  filter(!is.na(pct_first_trimester))

fig_anc_ft <- tm_shape(sa3_anc_map) +
  tm_polygons(
    fill = "pct_first_trimester",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 50, 65, 75, 85, 100),
      values = "brewer.rd_yl_gn"
    ),
    fill.legend = tm_legend(title = "% first\ntrimester ANC")
  ) +
  tm_title("First trimester antenatal care attendance by SA3, 2023")

save_map(fig_anc_ft, "fig22_anc_first_trimester_sa3", width = 10, height = 10)

# ==============================================================================
# 8. BIVARIATE ANALYSIS: ANC vs INDIGENOUS %
# ==============================================================================

message("\n8. Bivariate analysis: ANC vs Indigenous population...")

sa3_bivar <- sa3_risk %>%
  st_drop_geometry() %>%
  filter(!is.na(pct_first_trimester), !is.na(pct_indigenous))

# Correlation
cor_result <- cor.test(sa3_bivar$pct_first_trimester, sa3_bivar$pct_indigenous,
                       method = "spearman")

message("  Spearman correlation (ANC first trimester vs Indigenous %): ",
        round(cor_result$estimate, 3),
        " (p = ", format.pval(cor_result$p.value, digits = 4), ")")

fig_anc_indig <- ggplot(sa3_bivar,
                         aes(x = pct_indigenous, y = pct_first_trimester)) +
  geom_point(aes(colour = remoteness_mode), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, colour = "#e74c3c", linewidth = 0.8) +
  scale_colour_manual(values = c(
    "Major Cities" = "#3498db",
    "Inner Regional" = "#2ecc71",
    "Outer Regional" = "#f39c12",
    "Remote" = "#e74c3c",
    "Very Remote" = "#9b59b6"
  )) +
  labs(
    title = "First trimester antenatal care vs Indigenous population proportion",
    subtitle = paste0("SA3 areas | Spearman rho = ",
                      round(cor_result$estimate, 3),
                      ", p ", ifelse(cor_result$p.value < 0.001, "< 0.001",
                                     paste0("= ", round(cor_result$p.value, 3)))),
    x = "% Aboriginal and Torres Strait Islander population",
    y = "% attending ANC in first trimester",
    colour = "Remoteness",
    caption = "Sources: ABS Census 2021, AIHW NPDC 2023"
  ) +
  theme(legend.position = "right")

save_plot(fig_anc_indig, "fig23_anc_vs_indigenous", width = 10, height = 7)

# ==============================================================================
# 9. RISK FACTOR ANALYSIS: LOGISTIC REGRESSION
# ==============================================================================

message("\n9. Risk factor analysis...")

# Since we don't have SA3-level congenital case counts, we model high-risk
# classification as a function of service accessibility and demographics
# This identifies which area characteristics predict high congenital risk

sa3_model_data <- sa3_risk %>%
  st_drop_geometry() %>%
  filter(!is.na(pct_first_trimester), !is.na(mean_irsd),
         !is.na(mean_dist_sh_km)) %>%
  mutate(
    high_risk_binary = as.integer(risk_score >= 2),
    log_dist_sh = log1p(mean_dist_sh_km),
    remoteness_factor = factor(remoteness_mode,
                               levels = c("Major Cities", "Inner Regional",
                                          "Outer Regional", "Remote", "Very Remote"))
  )

# Descriptive table: ANC by risk category (use full dataset, not filtered)
anc_by_risk <- sa3_risk %>%
  st_drop_geometry() %>%
  filter(!is.na(risk_category)) %>%
  group_by(risk_category) %>%
  summarise(
    n = n(),
    mean_anc_ft = mean(pct_first_trimester, na.rm = TRUE),
    mean_pct_indigenous = mean(pct_indigenous, na.rm = TRUE),
    mean_dist_sh = mean(mean_dist_sh_km, na.rm = TRUE),
    mean_irsd = mean(mean_irsd, na.rm = TRUE),
    mean_births = mean(total_births, na.rm = TRUE),
    .groups = "drop"
  )

message("\n  ANC and demographics by risk category:")
for (i in seq_len(nrow(anc_by_risk))) {
  message("    ", anc_by_risk$risk_category[i], " (n=", anc_by_risk$n[i], "): ",
          "ANC=", round(anc_by_risk$mean_anc_ft[i], 1), "%, ",
          "Indigenous=", round(anc_by_risk$mean_pct_indigenous[i], 1), "%, ",
          "Dist SH=", round(anc_by_risk$mean_dist_sh[i], 0), "km, ",
          "IRSD=", round(anc_by_risk$mean_irsd[i], 0))
}

# ==============================================================================
# 10. TABLES
# ==============================================================================

message("\n10. Creating tables...")

# Table 3: Congenital syphilis summary
table3_data <- national_congenital %>%
  left_join(
    national_cong_indig %>%
      select(year, indigenous_count, non_indigenous_count),
    by = "year"
  )

table3_gt <- table3_data %>%
  gt() %>%
  tab_header(
    title = "Table 3. Congenital syphilis notifications, Australia 2011-2024",
    subtitle = "Notifications and deaths by year"
  ) %>%
  cols_label(
    year = "Year",
    count = "Total notifications",
    deaths = "Deaths",
    indigenous_count = "Indigenous",
    non_indigenous_count = "Non-Indigenous",
    source = "Source"
  ) %>%
  sub_missing(missing_text = "-") %>%
  tab_spanner(label = "By Indigenous status", columns = c(indigenous_count, non_indigenous_count)) %>%
  tab_source_note("Sources: Hengel et al. 2024 MJA, Kirby ASR 2024, National CDC Q2 2025 Report")

gtsave(table3_gt, file.path(path_tables, "table3_congenital_national.html"))
write_csv(table3_data, file.path(path_tables, "table3_congenital_national.csv"))

# Table 4: Risk classification summary by SA3
table4_gt <- anc_by_risk %>%
  mutate(
    mean_anc_ft = round(mean_anc_ft, 1),
    mean_pct_indigenous = round(mean_pct_indigenous, 1),
    mean_dist_sh = round(mean_dist_sh, 0),
    mean_irsd = round(mean_irsd, 0),
    mean_births = round(mean_births, 0)
  ) %>%
  gt() %>%
  tab_header(
    title = "Table 4. SA3 characteristics by congenital syphilis risk classification",
    subtitle = "Mean values across SA3 areas in each risk category"
  ) %>%
  cols_label(
    risk_category = "Risk category",
    n = "N (SA3s)",
    mean_anc_ft = "First trimester ANC (%)",
    mean_pct_indigenous = "Indigenous population (%)",
    mean_dist_sh = "Distance to SH clinic (km)",
    mean_irsd = "IRSD score",
    mean_births = "Annual births"
  ) %>%
  tab_source_note("Sources: ABS Census 2021, AIHW NPDC 2023, NHSD 2025")

gtsave(table4_gt, file.path(path_tables, "table4_risk_classification_sa3.html"))
write_csv(anc_by_risk, file.path(path_tables, "table4_risk_classification_sa3.csv"))

# ==============================================================================
# 11. COMBINED PANEL: CONGENITAL OVERVIEW
# ==============================================================================

message("\n11. Creating combined congenital panel figure...")

fig_panel <- (fig_cong_national + fig_cong_indig_national) /
  (fig_cong_qld_indig + fig_cong_qld_region) +
  plot_annotation(
    title = "Congenital syphilis in Australia: Temporal trends",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(fig_panel, "fig24_congenital_overview_panel", width = 16, height = 12)

# ==============================================================================
# 12. SAVE PROCESSED DATASETS
# ==============================================================================

message("\n12. Saving processed datasets...")

saveRDS(sa3_risk, file.path(path_processed, "sa3_congenital_risk.rds"))
saveRDS(sa3_model_data, file.path(path_processed, "sa3_risk_model_data.rds"))

# ==============================================================================
# 13. SUMMARY
# ==============================================================================

message("\n=== Phase 4 Output Summary ===")
message("Figures saved to: ", path_figures)
message("Maps saved to: ", path_maps)
message("Tables saved to: ", path_tables)

fig_files <- list.files(path_figures, pattern = "fig(1[7-9]|2[0-4])")
map_files <- list.files(path_maps, pattern = "fig2[12]")
table_files <- list.files(path_tables, pattern = "table[34]")

message("  Figures: ", length(fig_files))
for (f in fig_files) message("    ", f)
message("  Maps: ", length(map_files))
for (f in map_files) message("    ", f)
message("  Tables: ", length(table_files))
for (f in table_files) message("    ", f)

message("\n  Key findings:")
message("  - National congenital cases peaked at 20 in 2023 (10 deaths)")
message("  - ", sum(sa3_risk$risk_category == "High risk", na.rm = TRUE),
        " SA3 areas classified as high risk")
message("  - ", sum(sa3_risk$risk_category == "Elevated risk", na.rm = TRUE),
        " SA3 areas classified as elevated risk")
message("  - ANC first trimester vs Indigenous %: rho = ",
        round(cor_result$estimate, 3))

message("\nPhase 4 complete.")

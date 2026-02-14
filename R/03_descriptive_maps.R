# ==============================================================================
# 03_descriptive_maps.R — Phase 2: Descriptive Spatial Analysis
# National trends, state maps, sub-state maps, demographic comparisons
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

message("Loading analysis datasets...")

# Boundaries
ste <- readRDS(file.path(path_processed, "boundaries_ste.rds"))
sa2_sf <- readRDS(file.path(path_processed, "sa2_analysis_sf.rds"))
qld_hhs_boundaries <- readRDS(file.path(path_processed, "boundaries_qld_hhs.rds"))
vic_lga_boundaries <- readRDS(file.path(path_processed, "boundaries_vic_lga.rds"))

# Notifications
state_ts <- readRDS(file.path(path_processed, "state_time_series.rds"))
national_ts <- readRDS(file.path(path_processed, "national_time_series.rds"))
substate <- readRDS(file.path(path_processed, "substate_notifications.rds"))
vic_analysis <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))
qld_hhs_notif <- readRDS(file.path(path_processed, "qld_hhs_notifications.rds"))

# Kirby stratified data
kirby_state_sex <- readRDS(file.path(path_processed, "kirby_state_sex.rds"))
kirby_state_indig <- readRDS(file.path(path_processed, "kirby_state_indigenous.rds"))
kirby_national_sex <- readRDS(file.path(path_processed, "kirby_national_sex.rds"))
kirby_national_age <- readRDS(file.path(path_processed, "kirby_national_age.rds"))
kirby_national_indig <- readRDS(file.path(path_processed, "kirby_national_indigenous.rds"))
kirby_national_remote <- readRDS(file.path(path_processed, "kirby_national_remoteness.rds"))
kirby_remote_indig <- readRDS(file.path(path_processed, "kirby_remoteness_indigenous.rds"))
kirby_remote_sex <- readRDS(file.path(path_processed, "kirby_remoteness_sex.rds"))

# NSW LHD boundaries
nsw_lhd_boundaries <- nswgeo::lhd %>% st_transform(crs_gda2020)

# SA2 analysis (non-spatial)
sa2_analysis <- readRDS(file.path(path_processed, "sa2_analysis.rds"))

# ==============================================================================
# 2. FIGURE 1: NATIONAL SYPHILIS TREND
# ==============================================================================

message("Creating Figure 1: National trend...")

fig1 <- ggplot(national_ts, aes(x = year, y = count)) +
  geom_col(fill = "#c0392b", alpha = 0.8) +
  geom_text(aes(label = scales::comma(count)), vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = 2013:2022) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Infectious syphilis notifications, Australia",
    subtitle = "2013-2022",
    x = NULL, y = "Notifications (n)",
    caption = "Source: Kirby Institute Annual Surveillance Report"
  )

save_plot(fig1, "fig1_national_trend", width = 10, height = 6)

# ==============================================================================
# 3. FIGURE 2a: STATE-LEVEL RATE TRENDS
# ==============================================================================

message("Creating Figure 2a: State rate trends...")

fig2a <- ggplot(state_ts, aes(x = year, y = rate_per_100k, colour = state)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = state_colours) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  labs(
    title = "Infectious syphilis notification rates by state/territory",
    subtitle = "2013-2022",
    x = NULL, y = "Rate per 100,000 population",
    colour = "State/Territory",
    caption = "Source: Kirby Institute"
  ) +
  theme(legend.position = "right")

save_plot(fig2a, "fig2a_state_rate_trends", width = 10, height = 6)

# ==============================================================================
# 4. FIGURE 2b: STATE CHOROPLETH (most recent year)
# ==============================================================================

message("Creating Figure 2b: State choropleth (2022)...")

state_2022 <- state_ts %>%
  filter(year == 2022) %>%
  mutate(ste_name = case_when(
    state == "NSW" ~ "New South Wales",
    state == "VIC" ~ "Victoria",
    state == "QLD" ~ "Queensland",
    state == "SA"  ~ "South Australia",
    state == "WA"  ~ "Western Australia",
    state == "TAS" ~ "Tasmania",
    state == "NT"  ~ "Northern Territory",
    state == "ACT" ~ "Australian Capital Territory"
  ))

# Filter to mainland + Tasmania (exclude Other Territories)
ste_map <- ste %>%
  filter(STE_NAME21 %in% state_2022$ste_name) %>%
  left_join(state_2022, by = c("STE_NAME21" = "ste_name"))

fig2b <- tm_shape(ste_map) +
  tm_polygons(
    fill = "rate_per_100k",
    fill.scale = tm_scale_continuous(
      values = "brewer.yl_or_rd",
      midpoint = NA
    ),
    fill.legend = tm_legend(title = "Rate per\n100,000")
  ) +
  tm_text("state", size = 0.6) +
  tm_title("Infectious syphilis notification rates by state, 2022")

save_map(fig2b, "fig2b_state_choropleth_2022", width = 8, height = 8)

# ==============================================================================
# 5. FIGURE 3: RATES BY SEX AND INDIGENOUS STATUS (national)
# ==============================================================================

message("Creating Figure 3: National rates by sex and Indigenous status...")

# 3a. By sex
national_by_sex <- kirby_national_sex %>%
  filter(sex %in% c("Female", "Male"))

fig3a <- ggplot(national_by_sex, aes(x = year, y = count, fill = sex)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("Female" = "#e74c3c", "Male" = "#3498db")) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Infectious syphilis notifications by sex",
    subtitle = "Australia, 2013-2022",
    x = NULL, y = "Notifications (n)", fill = "Sex",
    caption = "Source: Kirby Institute"
  )

# 3b. By Indigenous status
national_by_indig <- kirby_national_indig %>%
  filter(indigenous_status != "Grand Total",
         indigenous_status != "Not reported") %>%
  mutate(indigenous_status = str_replace(indigenous_status,
    "Aboriginal and/or Torres Strait Islander",
    "Aboriginal and Torres\nStrait Islander peoples"))

fig3b <- ggplot(national_by_indig, aes(x = year, y = count,
                                        fill = indigenous_status)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c(
    "Aboriginal and Torres\nStrait Islander peoples" = "#e67e22",
    "Non-Indigenous" = "#2ecc71"
  )) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Infectious syphilis notifications by Indigenous status",
    subtitle = "Australia, 2013-2022",
    x = NULL, y = "Notifications (n)", fill = NULL,
    caption = "Source: Kirby Institute"
  )

fig3 <- fig3a + fig3b + plot_layout(ncol = 2, widths = c(1, 1))
save_plot(fig3, "fig3_national_sex_indigenous", width = 14, height = 6)

# ==============================================================================
# 6. FIGURE 4: REMOTENESS × INDIGENOUS STATUS RATES
# ==============================================================================

message("Creating Figure 4: Remoteness x Indigenous status...")

# Parse remoteness and Indigenous status from combined population_group
remote_indig <- kirby_remote_indig %>%
  mutate(
    indigenous = ifelse(str_detect(population_group, "Aboriginal"),
                        "Aboriginal and Torres Strait\nIslander peoples",
                        "Non-Indigenous"),
    remoteness_cat = case_when(
      str_detect(population_group, "major cities") ~ "Major Cities",
      str_detect(population_group, "regional") ~ "Regional",
      str_detect(population_group, "remote") ~ "Remote"
    )
  )

fig4 <- ggplot(remote_indig, aes(x = year, y = rate,
                                  colour = remoteness_cat,
                                  linetype = indigenous)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(
    "Major Cities" = "#3498db",
    "Regional" = "#f39c12",
    "Remote" = "#e74c3c"
  )) +
  scale_linetype_manual(values = c(
    "Aboriginal and Torres Strait\nIslander peoples" = "solid",
    "Non-Indigenous" = "dashed"
  )) +
  scale_x_continuous(breaks = seq(2013, 2022, 2)) +
  labs(
    title = "Infectious syphilis rates by remoteness and Indigenous status",
    subtitle = "Australia, 2013-2022",
    x = NULL, y = "Rate per 100,000 population",
    colour = "Remoteness", linetype = "Indigenous status",
    caption = "Source: Kirby Institute"
  ) +
  theme(legend.position = "right")

save_plot(fig4, "fig4_remoteness_indigenous_rates", width = 12, height = 7)

# ==============================================================================
# 7. QLD HHS MAPS — TEMPORAL SMALL MULTIPLES
# ==============================================================================

message("Creating QLD HHS temporal maps...")

# Select key years for small multiples
qld_years <- c(2012, 2015, 2018, 2021, 2024)

qld_map_data <- qld_hhs_notif %>%
  filter(year %in% qld_years, area_name != "Queensland") %>%
  left_join(
    qld_hhs_boundaries %>% rename(area_name = hhs),
    by = "area_name"
  ) %>%
  st_as_sf()

# Fixed breaks for consistent comparison across years
fig5_qld <- tm_shape(qld_map_data) +
  tm_polygons(
    fill = "rate",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 10, 25, 50, 100, 300),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(title = "Rate per\n100,000")
  ) +
  tm_facets(by = "year", ncol = 5, free.coords = FALSE) +
  tm_title("Infectious syphilis rates by HHS, Queensland")

save_map(fig5_qld, "fig5_qld_hhs_small_multiples", width = 16, height = 5)

# ==============================================================================
# 8. VIC LGA MAP — MOST RECENT YEAR + TEMPORAL
# ==============================================================================

message("Creating VIC LGA maps...")

# Clean LGA names for joining: remove "(C)", "(S)", "(RC)" etc.
vic_name_lookup <- vic_lga_boundaries %>%
  st_drop_geometry() %>%
  transmute(
    lga_name_clean = LGA_NAME21,
    lga_name_boundary = LGA_NAME21
  )

# Clean notification names similarly
vic_for_map <- vic_analysis %>%
  mutate(
    lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>%
      str_trim()
  )

# Build VIC spatial data for 2024 (most recent)
vic_2024 <- vic_for_map %>%
  filter(year == 2024) %>%
  group_by(lga_name_clean) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    rate_per_100k = count / population * 100000,
    .groups = "drop"
  )

vic_map_2024 <- vic_lga_boundaries %>%
  left_join(vic_2024, by = c("LGA_NAME21" = "lga_name_clean"))

fig6_vic <- tm_shape(vic_map_2024) +
  tm_polygons(
    fill = "rate_per_100k",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 5, 10, 20, 40, 100),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(title = "Rate per\n100,000")
  ) +
  tm_title("Infectious syphilis rates by LGA, Victoria 2024")

save_map(fig6_vic, "fig6_vic_lga_2024", width = 8, height = 8)

# VIC small multiples (2019-2024)
vic_all_years <- vic_for_map %>%
  group_by(lga_name_clean, year) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    rate_per_100k = count / population * 100000,
    .groups = "drop"
  )

vic_map_all <- vic_lga_boundaries %>%
  cross_join(tibble(year = 2019:2024)) %>%
  left_join(vic_all_years, by = c("LGA_NAME21" = "lga_name_clean", "year"))

fig6b_vic <- tm_shape(vic_map_all) +
  tm_polygons(
    fill = "rate_per_100k",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 5, 10, 20, 40, 100),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(title = "Rate per\n100,000")
  ) +
  tm_facets(by = "year", ncol = 3, free.coords = FALSE) +
  tm_title("Infectious syphilis rates by LGA, Victoria 2019-2024")

save_map(fig6b_vic, "fig6b_vic_lga_small_multiples", width = 12, height = 10)

# ==============================================================================
# 9. NSW LHD MAP — TEMPORAL
# ==============================================================================

message("Creating NSW LHD maps...")

nsw_notif <- substate %>%
  filter(state == "NSW") %>%
  select(area_name, year, count)

# Select key years
nsw_years <- c(2005, 2010, 2015, 2020, 2024)

nsw_map_data <- nsw_lhd_boundaries %>%
  cross_join(tibble(year = nsw_years)) %>%
  left_join(nsw_notif, by = c("lhd_name" = "area_name", "year")) %>%
  # Remove "NSW not otherwise specified"
  filter(!str_detect(lhd_name, "not otherwise"))

fig7_nsw <- tm_shape(nsw_map_data) +
  tm_polygons(
    fill = "count",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 10, 50, 100, 250, 700),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(title = "Notifications\n(n)")
  ) +
  tm_facets(by = "year", ncol = 5, free.coords = FALSE) +
  tm_title("Infectious syphilis notifications by LHD, NSW")

save_map(fig7_nsw, "fig7_nsw_lhd_small_multiples", width = 16, height = 5)

# ==============================================================================
# 10. FIGURE 8: SA2 COVARIATES — SERVICE ACCESSIBILITY
# ==============================================================================

message("Creating Figure 8: Service accessibility overview...")

# Distance to sexual health clinic by remoteness
sa2_remote <- sa2_analysis %>%
  filter(!is.na(remoteness), !is.na(dist_sh_clinic_km))

fig8 <- ggplot(sa2_remote,
               aes(x = factor(remoteness,
                               levels = c("Major Cities", "Inner Regional",
                                          "Outer Regional", "Remote",
                                          "Very Remote")),
                   y = dist_sh_clinic_km)) +
  geom_boxplot(fill = "#3498db", alpha = 0.6, outlier.size = 0.5) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Distance to nearest sexual health clinic by remoteness",
    subtitle = "SA2 centroids, 2021 ASGS boundaries",
    x = NULL, y = "Distance (km, log scale)",
    caption = "Source: NHSD 2025 via AURIN"
  )

save_plot(fig8, "fig8_service_accessibility_remoteness", width = 10, height = 6)

# ==============================================================================
# 11. FIGURE 9: SA2 INDIGENOUS PROPORTION MAP
# ==============================================================================

message("Creating Figure 9: Indigenous population distribution...")

sa2_indig_map <- sa2_sf %>%
  filter(!is.na(pct_indigenous), pop_total > 0)

fig9 <- tm_shape(sa2_indig_map) +
  tm_polygons(
    fill = "pct_indigenous",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 2, 5, 15, 50, 100),
      values = "brewer.oranges"
    ),
    fill.legend = tm_legend(title = "% Aboriginal and\nTorres Strait\nIslander peoples")
  ) +
  tm_title("Aboriginal and Torres Strait Islander population by SA2, 2021 Census")

save_map(fig9, "fig9_indigenous_pct_sa2", width = 10, height = 10)

# ==============================================================================
# 12. TABLE 1: DESCRIPTIVE EPIDEMIOLOGY
# ==============================================================================

message("Creating Table 1: Descriptive epidemiology...")

# State-level summary for most recent year (2022)
table1_state <- state_ts %>%
  filter(year == 2022) %>%
  select(state, rate_per_100k) %>%
  arrange(desc(rate_per_100k))

# By Indigenous status (2022 rates from Kirby state data)
table1_indig <- kirby_state_indig %>%
  filter(year == 2022) %>%
  select(state, indigenous_status, rate_per_100k = rate) %>%
  pivot_wider(names_from = indigenous_status,
              values_from = rate_per_100k,
              names_prefix = "rate_") %>%
  clean_names()

# By remoteness (national 2022)
table1_remote <- kirby_national_remote %>%
  filter(year == 2022, remoteness != "Grand Total") %>%
  select(remoteness, count_2022 = count)

# Combine into formatted table
table1_combined <- table1_state %>%
  left_join(table1_indig, by = "state") %>%
  arrange(desc(rate_per_100k))

table1_gt <- table1_combined %>%
  gt() %>%
  tab_header(
    title = "Table 1. Infectious syphilis notification rates by state/territory and Indigenous status, 2022",
    subtitle = "Rate per 100,000 population"
  ) %>%
  cols_label(
    state = "State",
    rate_per_100k = "Overall rate",
    rate_aboriginal_and_torres_strait_islander_people = "Aboriginal and Torres Strait Islander",
    rate_non_indigenous_people = "Non-Indigenous"
  ) %>%
  fmt_number(columns = where(is.numeric), decimals = 1) %>%
  tab_source_note("Source: Kirby Institute Annual Surveillance Report 2024")

gtsave(table1_gt, file.path(path_tables, "table1_descriptive_epidemiology.html"))

# Also save as CSV
write_csv(table1_combined, file.path(path_tables, "table1_descriptive_epidemiology.csv"))

# Remoteness table
table1_remote_gt <- kirby_remote_indig %>%
  filter(year == 2022) %>%
  mutate(
    indigenous = ifelse(str_detect(population_group, "Aboriginal"),
                        "Aboriginal and Torres Strait Islander",
                        "Non-Indigenous"),
    remoteness_cat = case_when(
      str_detect(population_group, "major cities") ~ "Major Cities",
      str_detect(population_group, "regional") ~ "Regional",
      str_detect(population_group, "remote") ~ "Remote"
    )
  ) %>%
  select(remoteness_cat, indigenous, rate) %>%
  pivot_wider(names_from = indigenous, values_from = rate) %>%
  gt() %>%
  tab_header(
    title = "Table 1b. Infectious syphilis rates by remoteness and Indigenous status, 2022",
    subtitle = "Rate per 100,000 population"
  ) %>%
  cols_label(remoteness_cat = "Remoteness") %>%
  fmt_number(columns = where(is.numeric), decimals = 1) %>%
  tab_source_note("Source: Kirby Institute Annual Surveillance Report 2024")

gtsave(table1_remote_gt, file.path(path_tables, "table1b_remoteness_indigenous.html"))

# ==============================================================================
# 13. SUMMARY
# ==============================================================================

message("\n=== Phase 2 Output Summary ===")
message("Figures saved to: ", path_figures)
message("Maps saved to: ", path_maps)
message("Tables saved to: ", path_tables)

# List outputs
fig_files <- list.files(path_figures, pattern = "fig")
map_files <- list.files(path_maps, pattern = "fig")
table_files <- list.files(path_tables)

message("  Figures: ", length(fig_files))
for (f in fig_files) message("    ", f)
message("  Maps: ", length(map_files))
for (f in map_files) message("    ", f)
message("  Tables: ", length(table_files))
for (f in table_files) message("    ", f)

message("\nPhase 2 complete.")

# ==============================================================================
# 08_manuscript_figures.R — Phase 6 (MSM rewrite): Publication-Ready Figures
# Creates 3 key figures and 2 tables for Victoria-focused MSM paper
# ==============================================================================

setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# Publication theme
theme_pub <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )
theme_set(theme_pub)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("=== Phase 6 (MSM rewrite): Manuscript Figure Preparation ===")
message("Loading data...")

vic_spatial <- readRDS(file.path(path_processed, "vic_lga_spatial_results.rds"))
reg_data <- readRDS(file.path(path_processed, "vic_lga_regression_data.rds"))
vic_analysis <- readRDS(file.path(path_processed, "vic_lga_analysis.rds"))
vic_lga_boundaries <- readRDS(file.path(path_processed, "boundaries_vic_lga.rds"))
model_comparison <- readRDS(file.path(path_processed, "model_comparison.rds"))

# ==============================================================================
# MANUSCRIPT FIGURE 1: Victorian LGA rate map (2024) + temporal small multiples
# ==============================================================================

message("\nFigure 1: Victorian LGA rate maps...")

# Clean LGA names for join
vic_name_lookup <- vic_lga_boundaries %>%
  st_drop_geometry() %>%
  select(LGA_NAME21)

# Panel A: 2024 rate map
vic_2024 <- vic_analysis %>%
  filter(year == 2024) %>%
  mutate(lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>% str_trim())

vic_map_2024 <- vic_lga_boundaries %>%
  left_join(vic_2024 %>% select(lga_name_clean, rate_per_100k),
            by = c("LGA_NAME21" = "lga_name_clean"))

ms_fig1a <- tm_shape(vic_map_2024) +
  tm_polygons(
    fill = "rate_per_100k",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 5, 15, 30, 50, 150),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(title = "Rate per\n100,000")
  ) +
  tm_title("A. Infectious syphilis notification rate by LGA, 2024") +
  tm_layout(frame = FALSE)

save_map(ms_fig1a, "ms_fig1a_vic_lga_2024", width = 7, height = 7, dpi = 600)

# Panel B: Temporal small multiples (2019, 2021, 2024)
vic_years <- c(2019, 2021, 2024)

vic_map_temporal <- vic_analysis %>%
  filter(year %in% vic_years) %>%
  mutate(lga_name_clean = str_remove(area_name, "\\s*\\(.*\\)$") %>% str_trim()) %>%
  left_join(vic_lga_boundaries %>% rename(lga_name_clean = LGA_NAME21),
            by = "lga_name_clean") %>%
  st_as_sf()

ms_fig1b <- tm_shape(vic_map_temporal) +
  tm_polygons(
    fill = "rate_per_100k",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 5, 15, 30, 50, 150),
      values = "brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(
      title = "Rate per\n100,000",
      position = tm_pos_out("right", "center")
    )
  ) +
  tm_facets(by = "year", ncol = 3, free.coords = FALSE) +
  tm_layout(frame = FALSE)

save_map(ms_fig1b, "ms_fig1b_vic_temporal", width = 12, height = 5, dpi = 600)

# ==============================================================================
# MANUSCRIPT FIGURE 2: Spatial analysis — LISA + Gi* (2-panel)
# ==============================================================================

message("Figure 2: Spatial analysis maps...")

# Panel A: LISA clusters
lisa_colours <- c(
  "High-High" = "#d73027",
  "High-Low" = "#fdae61",
  "Low-High" = "#abd9e9",
  "Low-Low" = "#4575b4",
  "Not significant" = "#f0f0f0"
)

ms_fig2a <- tm_shape(vic_spatial) +
  tm_polygons(
    fill = "lisa_cluster",
    fill.scale = tm_scale_categorical(values = lisa_colours),
    fill.legend = tm_legend(title = "LISA cluster")
  ) +
  tm_title("A. Local Moran's I clusters") +
  tm_layout(frame = FALSE)

# Panel B: Gi* hot spots
hotspot_colours <- c(
  "Hot spot (99% CI)" = "#d73027",
  "Hot spot (95% CI)" = "#fc8d59",
  "Hot spot (90% CI)" = "#fee08b",
  "Not significant"   = "#f0f0f0",
  "Cold spot (90% CI)" = "#d9ef8b",
  "Cold spot (95% CI)" = "#91bfdb",
  "Cold spot (99% CI)" = "#4575b4"
)

ms_fig2b <- tm_shape(vic_spatial) +
  tm_polygons(
    fill = "hotspot_class",
    fill.scale = tm_scale_categorical(values = hotspot_colours),
    fill.legend = tm_legend(title = "Gi* class")
  ) +
  tm_title("B. Getis-Ord Gi* hot spots") +
  tm_layout(frame = FALSE)

save_map(ms_fig2a, "ms_fig2a_lisa_clusters", width = 6, height = 6, dpi = 600)
save_map(ms_fig2b, "ms_fig2b_gi_hotspots", width = 6, height = 6, dpi = 600)

# ==============================================================================
# MANUSCRIPT FIGURE 3: GWR local R-squared + distance coefficient (2-panel)
# ==============================================================================

message("Figure 3: GWR spatial variation maps...")

ms_fig3a <- tm_shape(reg_data) +
  tm_polygons(
    fill = "gwr_localR2",
    fill.scale = tm_scale_continuous(
      values = "brewer.yl_gn_bu",
      midpoint = NA
    ),
    fill.legend = tm_legend(title = "Local R\u00b2")
  ) +
  tm_title("A. GWR local R-squared") +
  tm_layout(frame = FALSE)

ms_fig3b <- tm_shape(reg_data) +
  tm_polygons(
    fill = "gwr_coef_dist",
    fill.scale = tm_scale_continuous(
      values = "brewer.rd_yl_bu",
      midpoint = 0
    ),
    fill.legend = tm_legend(title = "Coefficient:\nlog(Distance)")
  ) +
  tm_title("B. GWR local coefficient: distance to sexual health clinic") +
  tm_layout(frame = FALSE)

save_map(ms_fig3a, "ms_fig3a_gwr_local_r2", width = 6, height = 6, dpi = 600)
save_map(ms_fig3b, "ms_fig3b_gwr_coef_distance", width = 6, height = 6, dpi = 600)

# ==============================================================================
# MANUSCRIPT TABLE 1: Victorian LGA descriptive statistics
# ==============================================================================

message("\nCreating publication tables...")
message("Table 1: Victorian LGA descriptive statistics...")

# Summarise VIC LGA characteristics
vic_desc <- vic_spatial %>%
  st_drop_geometry() %>%
  filter(!is.na(mean_annual_rate)) %>%
  summarise(
    n_lgas = n(),
    rate_median = median(mean_annual_rate, na.rm = TRUE),
    rate_iqr_lo = quantile(mean_annual_rate, 0.25, na.rm = TRUE),
    rate_iqr_hi = quantile(mean_annual_rate, 0.75, na.rm = TRUE),
    rate_max = max(mean_annual_rate, na.rm = TRUE),
    pop_median = median(lga_pop, na.rm = TRUE),
    pct_male2044_median = median(pct_male_20_44, na.rm = TRUE),
    density_median = median(pop_density, na.rm = TRUE),
    irsd_median = median(mean_irsd, na.rm = TRUE),
    dist_sh_median = median(mean_dist_sh_km, na.rm = TRUE)
  )

# Summarise by LISA cluster status
cluster_summary <- vic_spatial %>%
  st_drop_geometry() %>%
  filter(!is.na(mean_annual_rate)) %>%
  mutate(cluster_type = ifelse(lisa_cluster == "High-High", "Hot spot cluster", "Other LGAs")) %>%
  group_by(cluster_type) %>%
  summarise(
    n = n(),
    mean_rate = round(mean(mean_annual_rate, na.rm = TRUE), 1),
    median_rate = round(median(mean_annual_rate, na.rm = TRUE), 1),
    mean_pct_male2044 = round(mean(pct_male_20_44, na.rm = TRUE), 1),
    mean_density = round(mean(pop_density, na.rm = TRUE), 0),
    mean_irsd = round(mean(mean_irsd, na.rm = TRUE), 0),
    mean_dist_sh = round(mean(mean_dist_sh_km, na.rm = TRUE), 1),
    .groups = "drop"
  )

ms_table1_gt <- cluster_summary %>%
  gt() %>%
  tab_header(
    title = "Table 1. Characteristics of Victorian LGAs by spatial clustering status",
    subtitle = "Mean annual infectious syphilis rate, 2019-2024"
  ) %>%
  cols_label(
    cluster_type = "",
    n = "N",
    mean_rate = "Mean rate",
    median_rate = "Median rate",
    mean_pct_male2044 = "% Males 20-44",
    mean_density = "Pop. density",
    mean_irsd = "IRSD",
    mean_dist_sh = "Dist. to SH (km)"
  ) %>%
  tab_spanner(label = "Rate per 100,000", columns = c(mean_rate, median_rate)) %>%
  tab_spanner(label = "LGA characteristics (mean)", columns = c(mean_pct_male2044, mean_density, mean_irsd, mean_dist_sh)) %>%
  tab_source_note("Source: VIC Health notifications 2019-2024, ABS Census 2021, NHSD 2025") %>%
  tab_source_note("Hot spot clusters identified by Local Moran's I (p < 0.05)")

gtsave(ms_table1_gt, file.path(path_tables, "ms_table1_vic_descriptive.html"))
write_csv(cluster_summary, file.path(path_tables, "ms_table1_vic_descriptive.csv"))

# --- Manuscript Table 2: Spatial regression ---
message("Table 2: Spatial regression (reading from updated table)...")

table2_data <- read_csv(file.path(path_tables, "table2_spatial_regression.csv"),
                        show_col_types = FALSE)

ms_table2_gt <- table2_data %>%
  gt() %>%
  tab_header(
    title = "Table 2. Spatial regression models: predictors of infectious syphilis rates",
    subtitle = "Victorian LGAs (n=74), mean annual rate 2019-2024, log-transformed"
  ) %>%
  tab_source_note("Coefficients shown as estimate (SE). * p<0.05, ** p<0.01, *** p<0.001") %>%
  tab_source_note("Data sources: VIC Health, ABS Census 2021, NHSD 2025")

gtsave(ms_table2_gt, file.path(path_tables, "ms_table2_regression.html"))

# ==============================================================================
# SUPPLEMENTARY FIGURES LIST (MSM-focused paper)
# ==============================================================================

message("\nCompiling supplementary figure list...")

supp_figures <- tibble(
  figure = c(
    "S1. VIC LGA temporal small multiples (2019-2024, all years)",
    "S2. Moran's I scatterplot",
    "S3. Temporal stability of Moran's I (2019-2024)",
    "S4. GWR local coefficient: % Males aged 20-44",
    "S5. GWR local coefficient: Population density",
    "S6. OLS residual diagnostic plot",
    "S7. GWR local coefficient significance map (distance)",
    "S8. NSW monthly ITS analysis (sensitivity, n=204 months)",
    "S9. NSW seasonal notification patterns by COVID period"
  ),
  source_file = c(
    "fig6b_vic_lga_small_multiples",
    "fig10_moran_scatterplot",
    "temporal_morans_i",
    "fig15_gwr_coef_male2044",
    "fig16_gwr_coef_density",
    "ols_residual_diagnostic",
    "gwr_significance_distance",
    "fig33_nsw_monthly_its",
    "fig34_nsw_seasonal_pattern"
  )
)

write_csv(supp_figures, file.path(path_tables, "supplementary_figure_list.csv"))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== Phase 6 (MSM rewrite) Output Summary ===")

ms_figs <- list.files(path_figures, pattern = "^ms_fig")
ms_maps <- list.files(path_maps, pattern = "^ms_fig")
ms_tables <- list.files(path_tables, pattern = "^ms_table")

message("Manuscript figures: ", length(ms_figs))
for (f in ms_figs) message("  ", f)
message("Manuscript maps: ", length(ms_maps))
for (f in ms_maps) message("  ", f)
message("Manuscript tables: ", length(ms_tables))
for (f in ms_tables) message("  ", f)
message("Supplementary figures: ", nrow(supp_figures))

message("\nFigure mapping (MSM-focused paper):")
message("  Figure 1: ms_fig1a (VIC LGA 2024) + ms_fig1b (temporal 2019/2021/2024)")
message("  Figure 2: ms_fig2a (LISA clusters) + ms_fig2b (Gi* hot spots)")
message("  Figure 3: ms_fig3a (GWR local R2) + ms_fig3b (GWR distance coefficient)")
message("  Table 1: ms_table1_vic_descriptive (LGA characteristics by cluster status)")
message("  Table 2: ms_table2_regression (spatial regression comparison)")

message("\nPhase 6 (MSM rewrite) complete.")

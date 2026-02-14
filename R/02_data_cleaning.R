# ==============================================================================
# 02_data_cleaning.R — Harmonise notifications, build analysis datasets
# Creates unified sub-state notification dataset + SA2 analysis dataset
# ==============================================================================

# Source setup — works whether run from project root or R/ directory
setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. LOAD PROCESSED DATA FROM 01_data_acquisition.R
# ==============================================================================

message("Loading processed data objects...")

# Boundaries
sa2      <- readRDS(file.path(path_processed, "boundaries_sa2.rds"))
sa3      <- readRDS(file.path(path_processed, "boundaries_sa3.rds"))
qld_hhs  <- readRDS(file.path(path_processed, "boundaries_qld_hhs.rds"))
vic_lga  <- readRDS(file.path(path_processed, "boundaries_vic_lga.rds"))

# Demographics
sa2_pop      <- readRDS(file.path(path_processed, "census_sa2_population.rds"))
sa2_indig    <- readRDS(file.path(path_processed, "census_sa2_indigenous.rds"))
seifa        <- readRDS(file.path(path_processed, "seifa_sa2.rds"))
remoteness   <- readRDS(file.path(path_processed, "remoteness_sa2.rds"))
erp_lga      <- readRDS(file.path(path_processed, "erp_lga.rds"))

# Kirby national data
kirby_state_sex       <- readRDS(file.path(path_processed, "kirby_state_sex.rds"))
kirby_state_indig     <- readRDS(file.path(path_processed, "kirby_state_indigenous.rds"))
kirby_national_sex    <- readRDS(file.path(path_processed, "kirby_national_sex.rds"))
kirby_national_age    <- readRDS(file.path(path_processed, "kirby_national_age.rds"))
kirby_national_indig  <- readRDS(file.path(path_processed, "kirby_national_indigenous.rds"))
kirby_national_remote <- readRDS(file.path(path_processed, "kirby_national_remoteness.rds"))

# Sub-state notifications
nsw_lhd  <- readRDS(file.path(path_processed, "nsw_lhd_notifications.rds"))
vic_lga_notif <- readRDS(file.path(path_processed, "vic_lga_notifications.rds"))
qld_hhs_notif <- readRDS(file.path(path_processed, "qld_hhs_notifications.rds"))
nt_district   <- readRDS(file.path(path_processed, "nt_district_notifications.rds"))

# Services
sh_clinics  <- readRDS(file.path(path_processed, "sexual_health_clinics.rds"))
ah_services <- readRDS(file.path(path_processed, "aboriginal_health_services.rds"))

# ==============================================================================
# 2. HARMONISE STATE-LEVEL TIME SERIES (from Kirby)
# ==============================================================================

message("Building state-level time series...")

# Kirby state × sex rates already include Overall, so extract for state totals
state_time_series <- kirby_state_sex %>%
  filter(sex == "Overall") %>%
  select(state, year, rate_per_100k = rate) %>%
  mutate(data_type = "rate")

# National totals from Kirby
national_time_series <- kirby_national_sex %>%
  filter(sex == "Grand Total") %>%
  select(year, count) %>%
  mutate(state = "Australia")

message("  State time series: ", n_distinct(state_time_series$state),
        " states, ", min(state_time_series$year), "-",
        max(state_time_series$year))

# ==============================================================================
# 3. HARMONISE SUB-STATE NOTIFICATIONS
# ==============================================================================

message("Harmonising sub-state notifications...")

# --- 3a. NSW LHD: counts, need LHD populations ---
# NSW LHD populations not in our data — use counts for mapping/clustering
# Can derive approximate populations from SA2 → LHD concordance later

nsw_harmonised <- nsw_lhd %>%
  select(state, geography_level, area_name, year, count) %>%
  mutate(data_type = "count")

# --- 3b. VIC LGA: have counts and population ---

vic_harmonised <- vic_lga_notif %>%
  mutate(
    rate_per_100k = count / population * 100000,
    data_type = "count_and_rate"
  ) %>%
  select(state, geography_level, area_name, year, count, population,
         rate_per_100k, data_type)

# --- 3c. QLD HHS: have rates only (from PDF extraction) ---

qld_harmonised <- qld_hhs_notif %>%
  mutate(
    rate_per_100k = rate,
    data_type = "rate"
  ) %>%
  select(state, geography_level, area_name, year,
         rate_per_100k, data_type)

# --- 3d. NT districts: have counts and rates for 2024 ---

nt_harmonised <- nt_district %>%
  filter(Characteristic %in% c("Darwin_Urban", "Darwin_Rural", "East_Arnhem",
                                "Katherine", "Barkly",
                                "Alice_Springs_Urban", "Alice_Springs_Rural")) %>%
  transmute(
    state = "NT",
    geography_level = "Health District",
    area_name = Characteristic,
    year = 2024L,
    count = Annual_2024_Count,
    data_type = "count"
  )

# --- Combine all sub-state data ---

substate_notifications <- bind_rows(
  nsw_harmonised,
  vic_harmonised %>% select(state, geography_level, area_name, year,
                            count, data_type),
  qld_harmonised %>% select(state, geography_level, area_name, year,
                            rate_per_100k, data_type),
  nt_harmonised
)

message("  Sub-state records: ", nrow(substate_notifications))
message("    NSW LHDs: ", n_distinct(nsw_harmonised$area_name))
message("    VIC LGAs: ", n_distinct(vic_harmonised$area_name))
message("    QLD HHS:  ", n_distinct(qld_harmonised$area_name))
message("    NT districts: ", n_distinct(nt_harmonised$area_name))

# ==============================================================================
# 4. BUILD SA2-LEVEL ANALYSIS DATASET
# ==============================================================================

message("Building SA2 analysis dataset...")

# Join demographics at SA2 level
sa2_analysis <- sa2 %>%
  st_drop_geometry() %>%
  select(sa2_code = SA2_CODE21, sa2_name = SA2_NAME21,
         sa3_code = SA3_CODE21, sa3_name = SA3_NAME21,
         ste_code = STE_CODE21, ste_name = STE_NAME21,
         area_sqkm = AREASQKM21) %>%
  # Add state abbreviation
  mutate(state = case_when(
    ste_code == "1" ~ "NSW",
    ste_code == "2" ~ "VIC",
    ste_code == "3" ~ "QLD",
    ste_code == "4" ~ "SA",
    ste_code == "5" ~ "WA",
    ste_code == "6" ~ "TAS",
    ste_code == "7" ~ "NT",
    ste_code == "8" ~ "ACT",
    ste_code == "9" ~ "OT",  # Other Territories
    TRUE ~ NA_character_
  )) %>%
  # Join Census population
  left_join(sa2_pop, by = "sa2_code") %>%
  # Join Indigenous detail
  left_join(sa2_indig %>% select(sa2_code, indigenous_total, indigenous_female),
            by = "sa2_code") %>%
  # Join SEIFA
  left_join(seifa %>% select(sa2_code, irsd_score, irsd_decile),
            by = "sa2_code") %>%
  # Join remoteness
  left_join(remoteness %>% select(sa2_code, remoteness),
            by = "sa2_code") %>%
  # Calculate percentage Indigenous
  mutate(
    pct_indigenous = ifelse(pop_total > 0,
                            indigenous_total / pop_total * 100, NA_real_)
  )

# Filter to meaningful SA2s (exclude "Migratory - Offshore - Shipping"
# and "No usual address" type areas)
sa2_analysis <- sa2_analysis %>%
  filter(pop_total > 0 | is.na(pop_total),  # Keep areas with population
         !str_detect(sa2_name, "Migratory|No usual address|Offshore"))

message("  SA2 analysis dataset: ", nrow(sa2_analysis), " areas")

# ==============================================================================
# 5. CALCULATE DISTANCE TO NEAREST CLINIC FOR EACH SA2
# ==============================================================================

message("Calculating distance to nearest clinic...")

# Get SA2 centroids
sa2_sf <- sa2 %>%
  filter(SA2_CODE21 %in% sa2_analysis$sa2_code)

sa2_centroids <- st_centroid(sa2_sf)

# Distance to nearest sexual health clinic
if (nrow(sh_clinics) > 0) {
  dist_matrix <- st_distance(sa2_centroids, sh_clinics)
  min_dist_sh <- apply(dist_matrix, 1, min)

  sa2_dist_sh <- tibble(
    sa2_code = sa2_centroids$SA2_CODE21,
    dist_sh_clinic_km = as.numeric(min_dist_sh) / 1000
  )
} else {
  sa2_dist_sh <- tibble(sa2_code = character(), dist_sh_clinic_km = numeric())
}

# Distance to nearest Aboriginal health service
if (nrow(ah_services) > 0) {
  dist_matrix_ah <- st_distance(sa2_centroids, ah_services)
  min_dist_ah <- apply(dist_matrix_ah, 1, min)

  sa2_dist_ah <- tibble(
    sa2_code = sa2_centroids$SA2_CODE21,
    dist_ah_service_km = as.numeric(min_dist_ah) / 1000
  )
} else {
  sa2_dist_ah <- tibble(sa2_code = character(), dist_ah_service_km = numeric())
}

# Distance to nearest of EITHER type
all_services <- readRDS(file.path(path_processed, "health_services.rds"))
if (nrow(all_services) > 0) {
  dist_matrix_any <- st_distance(sa2_centroids, all_services)
  min_dist_any <- apply(dist_matrix_any, 1, min)

  sa2_dist_any <- tibble(
    sa2_code = sa2_centroids$SA2_CODE21,
    dist_any_service_km = as.numeric(min_dist_any) / 1000
  )
} else {
  sa2_dist_any <- tibble(sa2_code = character(), dist_any_service_km = numeric())
}

# Join distances to analysis dataset
sa2_analysis <- sa2_analysis %>%
  left_join(sa2_dist_sh, by = "sa2_code") %>%
  left_join(sa2_dist_ah, by = "sa2_code") %>%
  left_join(sa2_dist_any, by = "sa2_code")

message("  Distance calculated for ", sum(!is.na(sa2_analysis$dist_sh_clinic_km)),
        " SA2 areas")

# ==============================================================================
# 6. BUILD VIC LGA ANALYSIS DATASET (with geometry)
# ==============================================================================

message("Building VIC LGA spatial analysis dataset...")

# Clean VIC LGA names to match boundary file
vic_analysis <- vic_harmonised %>%
  # Get the most recent year with most complete data
  group_by(area_name) %>%
  # Keep all years for time series
  ungroup() %>%
  left_join(
    vic_lga %>%
      st_drop_geometry() %>%
      transmute(
        lga_name_match = LGA_NAME21,
        lga_code = LGA_CODE21
      ),
    by = c("area_name" = "lga_name_match")
  )

message("  VIC LGA analysis: ", n_distinct(vic_analysis$area_name), " LGAs, ",
        nrow(vic_analysis), " records")

# ==============================================================================
# 7. VALIDATION CHECKS
# ==============================================================================

message("\n=== Validation Summary ===")

# Check for missing joins
n_missing_seifa <- sum(is.na(sa2_analysis$irsd_score))
n_missing_remote <- sum(is.na(sa2_analysis$remoteness))
n_missing_pop <- sum(is.na(sa2_analysis$pop_total))
n_missing_dist <- sum(is.na(sa2_analysis$dist_sh_clinic_km))

message("  SA2 areas missing SEIFA: ", n_missing_seifa)
message("  SA2 areas missing remoteness: ", n_missing_remote)
message("  SA2 areas missing population: ", n_missing_pop)
message("  SA2 areas missing clinic distance: ", n_missing_dist)

# Summary by state
state_summary <- sa2_analysis %>%
  group_by(state) %>%
  summarise(
    n_sa2 = n(),
    median_pop = median(pop_total, na.rm = TRUE),
    median_pct_indigenous = median(pct_indigenous, na.rm = TRUE),
    median_dist_sh_km = median(dist_sh_clinic_km, na.rm = TRUE),
    .groups = "drop"
  )

message("\nSA2 summary by state:")
print(as.data.frame(state_summary))

# Sub-state notification summary
substate_summary <- substate_notifications %>%
  group_by(state, geography_level) %>%
  summarise(
    n_areas = n_distinct(area_name),
    year_range = paste(min(year), max(year), sep = "-"),
    .groups = "drop"
  )

message("\nSub-state notification data:")
print(as.data.frame(substate_summary))

# ==============================================================================
# 8. SAVE ANALYSIS DATASETS
# ==============================================================================

message("\nSaving analysis datasets...")

saveRDS(sa2_analysis, file.path(path_processed, "sa2_analysis.rds"))
saveRDS(state_time_series, file.path(path_processed, "state_time_series.rds"))
saveRDS(national_time_series, file.path(path_processed, "national_time_series.rds"))
saveRDS(substate_notifications, file.path(path_processed, "substate_notifications.rds"))
saveRDS(vic_analysis, file.path(path_processed, "vic_lga_analysis.rds"))

# Kirby remoteness × indigenous data already saved by 01_data_acquisition.R

# Save SA2 analysis as spatial dataset (with geometry)
sa2_analysis_sf <- sa2_sf %>%
  rename(sa2_code = SA2_CODE21) %>%
  left_join(sa2_analysis %>% select(-sa2_name, -sa3_code, -sa3_name,
                                     -ste_code, -ste_name, -area_sqkm,
                                     -state),
            by = "sa2_code")

saveRDS(sa2_analysis_sf, file.path(path_processed, "sa2_analysis_sf.rds"))

message("All analysis datasets saved to ", path_processed)
message("Data cleaning complete.")

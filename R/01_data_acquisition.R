# ==============================================================================
# 01_data_acquisition.R — Extract and prepare raw data files
# Reads all raw data from data/raw/ and saves standardised R objects
# ==============================================================================

# Source setup — works whether run from project root or R/ directory
setup_path <- if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R"
source(setup_path)

# ==============================================================================
# 1. GEOGRAPHIC BOUNDARIES
# ==============================================================================

message("Loading boundary files...")

# SA2 boundaries (already GDA2020)
sa2 <- st_read(file.path(path_boundaries, "SA2_2021", "SA2_2021_AUST_GDA2020.shp"),
               quiet = TRUE) %>%
  select(SA2_CODE21, SA2_NAME21, SA3_CODE21, SA3_NAME21,
         SA4_CODE21, SA4_NAME21, STE_CODE21, STE_NAME21, AREASQKM21)

# SA3 boundaries (already GDA2020)
sa3 <- st_read(file.path(path_boundaries, "SA3_2021", "SA3_2021_AUST_GDA2020.shp"),
               quiet = TRUE) %>%
  select(SA3_CODE21, SA3_NAME21, SA4_CODE21, SA4_NAME21,
         STE_CODE21, STE_NAME21, AREASQKM21)

# State boundaries (already GDA2020)
ste <- st_read(file.path(path_boundaries, "STE_2021", "STE_2021_AUST_GDA2020.shp"),
               quiet = TRUE) %>%
  select(STE_CODE21, STE_NAME21, AREASQKM21)

# PHN boundaries — GDA94, reproject to GDA2020
phn <- st_read(file.path(path_boundaries, "PHN_boundaries.geojson"), quiet = TRUE) %>%
  st_transform(crs_gda2020)

# QLD HHS boundaries — Web Mercator, reproject to GDA2020
qld_hhs <- st_read(file.path(path_boundaries, "QLD_HHS_boundaries.geojson"),
                    quiet = TRUE) %>%
  st_transform(crs_gda2020)

# VIC LGA boundaries — filter from national file (already GDA2020)
lga <- st_read(file.path(path_boundaries, "LGA_2021",
                          "LGA_2021_AUST_GDA2020.shp"), quiet = TRUE)
vic_lga <- lga %>%
  filter(STE_CODE21 == "2")  # VIC = state code 2

message("  SA2: ", nrow(sa2), " areas")
message("  SA3: ", nrow(sa3), " areas")
message("  PHN: ", nrow(phn), " areas")
message("  QLD HHS: ", nrow(qld_hhs), " areas")
message("  VIC LGA: ", nrow(vic_lga), " areas")

# ==============================================================================
# 2. DEMOGRAPHIC DATA — Census 2021
# ==============================================================================

message("Extracting Census 2021 data...")

# Extract Census zip if not already done
census_dir <- file.path(path_demographics, "census_2021_sa2")
if (!dir.exists(census_dir)) {
  dir.create(census_dir, recursive = TRUE)
  unzip(file.path(path_demographics, "Census_2021_G07_Indigenous_SA2.zip"),
        exdir = census_dir)
}

census_path <- file.path(census_dir,
                         "2021 Census GCP Statistical Area 2 for AUS")

# G01: Total population by SA2 (includes age groups and Indigenous totals)
g01 <- read_csv(file.path(census_path, "2021Census_G01_AUST_SA2.csv"),
                show_col_types = FALSE)

# Extract key population variables from G01
sa2_population <- g01 %>%
  transmute(
    sa2_code = as.character(SA2_CODE_2021),
    pop_total   = Tot_P_P,
    pop_male    = Tot_P_M,
    pop_female  = Tot_P_F,
    pop_indigenous = Indigenous_P_Tot_P,
    pop_15_29 = Age_15_19_yr_P + Age_20_24_yr_P + Age_25_34_yr_P * 0.5,
    pop_15_44 = Age_15_19_yr_P + Age_20_24_yr_P + Age_25_34_yr_P +
                Age_35_44_yr_P
  ) %>%
  filter(!sa2_code %in% c("AUST", "0"))  # Remove summary rows

# G07: Indigenous population by age group by SA2
g07 <- read_csv(file.path(census_path, "2021Census_G07_AUST_SA2.csv"),
                show_col_types = FALSE)

# Extract Indigenous population for reproductive-age females (15-44)
sa2_indigenous_detail <- g07 %>%
  transmute(
    sa2_code = as.character(SA2_CODE_2021),
    indigenous_total = Tot_Indigenous_P,
    indigenous_female = Tot_Indigenous_F,
    non_indigenous_total = Tot_Non_Indigenous_P
  ) %>%
  filter(!sa2_code %in% c("AUST", "0"))

message("  Census SA2 population: ", nrow(sa2_population), " areas")

# ==============================================================================
# 3. SEIFA (Index of Relative Socio-economic Disadvantage)
# ==============================================================================

message("Loading SEIFA data...")

seifa <- read_excel(file.path(path_demographics, "SEIFA_2021_SA2.xlsx"),
                    sheet = "Table 1", skip = 5) %>%
  clean_names()

# The first two cols are SA2 code and name, followed by IRSD score and decile
seifa_clean <- seifa %>%
  select(1:4) %>%
  set_names(c("sa2_code", "sa2_name", "irsd_score", "irsd_decile")) %>%
  filter(!is.na(sa2_code), !str_detect(sa2_code, "Total|Source|Note")) %>%
  mutate(
    sa2_code = as.character(sa2_code),
    irsd_score = as.numeric(irsd_score),
    irsd_decile = as.integer(irsd_decile)
  )

message("  SEIFA: ", nrow(seifa_clean), " SA2 areas")

# ==============================================================================
# 4. REMOTENESS CLASSIFICATION
# ==============================================================================

message("Loading remoteness classification...")

# SA1-level remoteness (finest available)
ra_sa1 <- read_excel(file.path(path_demographics, "RA_2021_AUST.xlsx"),
                     sheet = 1) %>%
  clean_names()

# Map SA1 codes to SA2 via the first 9 digits
# SA1 codes are 11 digits, SA2 codes are the first 9 digits
# Use RA_NAME for classification (RA_CODE has state-specific sub-codes)
ra_sa2 <- ra_sa1 %>%
  mutate(
    sa2_code = str_sub(sa1_code_2021, 1, 9),
    ra_name = ra_name_2021,
    remoteness = case_when(
      str_detect(ra_name, "Major Cities") ~ "Major Cities",
      str_detect(ra_name, "Inner Regional") ~ "Inner Regional",
      str_detect(ra_name, "Outer Regional") ~ "Outer Regional",
      str_detect(ra_name, "^Remote") ~ "Remote",
      str_detect(ra_name, "Very Remote") ~ "Very Remote",
      TRUE ~ NA_character_
    )
  ) %>%
  # Remove non-geographic SA1s
  filter(!is.na(remoteness), !str_detect(sa2_code, "Z")) %>%
  # For each SA2, assign the most common remoteness category
  group_by(sa2_code) %>%
  summarise(
    remoteness = names(which.max(table(remoteness))),
    .groups = "drop"
  )

message("  Remoteness: ", nrow(ra_sa2), " SA2 areas")

# ==============================================================================
# 5. ERP (LGA level — for VIC analysis)
# ==============================================================================

message("Loading ERP data...")

erp_lga <- read_excel(file.path(path_demographics, "ERP_SA2_2001-2024.xlsx"),
                      sheet = "Table 1", skip = 5) %>%
  clean_names()

# Rename and reshape: first 2 cols are LGA code + name, then years 2001-2024
year_cols <- as.character(2001:2024)
erp_names <- c("lga_code", "lga_name", year_cols)

erp_lga_clean <- erp_lga %>%
  select(1:(2 + length(year_cols))) %>%
  set_names(erp_names) %>%
  filter(!is.na(lga_code), !str_detect(as.character(lga_code), "Total|Source|Note")) %>%
  pivot_longer(cols = all_of(year_cols), names_to = "year",
               values_to = "erp") %>%
  mutate(
    lga_code = as.character(lga_code),
    year = as.integer(year),
    erp = as.numeric(erp)
  )

message("  ERP LGA: ", n_distinct(erp_lga_clean$lga_code), " LGAs x ",
        n_distinct(erp_lga_clean$year), " years")

# ==============================================================================
# 6. NOTIFICATION DATA
# ==============================================================================

message("Loading notification data...")

# --- 6a. Kirby national data (wide format → long) ---

kirby_state_sex <- read_csv(
  file.path(path_notifications, "kirby_rates_by_state_sex.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-c(State, Population_group), names_to = "year",
               values_to = "rate") %>%
  rename(state = State, sex = Population_group) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_state_indigenous <- read_csv(
  file.path(path_notifications, "kirby_rates_by_state_indigenous.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-c(State, Population_group), names_to = "year",
               values_to = "rate") %>%
  rename(state = State, indigenous_status = Population_group) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_national_sex <- read_csv(
  file.path(path_notifications, "kirby_national_counts_by_gender.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-Category, names_to = "year", values_to = "count") %>%
  rename(sex = Category) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_national_age <- read_csv(
  file.path(path_notifications, "kirby_national_counts_by_age.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-Category, names_to = "year", values_to = "count") %>%
  rename(age_group = Category) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_national_indigenous <- read_csv(
  file.path(path_notifications, "kirby_national_counts_by_indigenous.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-Category, names_to = "year", values_to = "count") %>%
  rename(indigenous_status = Category) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_national_remoteness <- read_csv(
  file.path(path_notifications, "kirby_national_counts_by_remoteness.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-Category, names_to = "year", values_to = "count") %>%
  rename(remoteness = Category) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_remoteness_indigenous <- read_csv(
  file.path(path_notifications, "kirby_rates_australia_remoteness_indigenous.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-c(Geography, Population_group), names_to = "year",
               values_to = "rate") %>%
  rename(geography = Geography, population_group = Population_group) %>%
  mutate(year = as.integer(year), source = "kirby")

kirby_remoteness_sex <- read_csv(
  file.path(path_notifications, "kirby_rates_australia_remoteness_sex.csv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-c(Geography, Population_group), names_to = "year",
               values_to = "rate") %>%
  rename(geography = Geography, population_group = Population_group) %>%
  mutate(year = as.integer(year), source = "kirby")


# --- 6b. NSW LHD data ---

nsw_raw <- read_excel(
  file.path(path_notifications, "NSW_syphilis_by_LHD_2000-2025.xls"),
  sheet = 1, skip = 2
)

# Drop unnamed columns (auto-named ...N by readxl) — these are empty spacers
nsw_raw <- nsw_raw %>% select(-starts_with("..."))

# Fill down the Strain column
nsw_raw <- nsw_raw %>% fill(Strain, .direction = "down")

# Pivot LHD columns to long format
lhd_cols <- setdiff(names(nsw_raw), c("Strain", "Year", "Total"))
nsw_lhd <- nsw_raw %>%
  filter(!is.na(Year), Year <= 2024) %>%
  pivot_longer(cols = all_of(lhd_cols), names_to = "lhd",
               values_to = "count") %>%
  mutate(
    year = as.integer(Year),
    strain = str_trim(Strain),
    count = as.integer(count),
    state = "NSW",
    geography_level = "LHD"
  ) %>%
  select(state, geography_level, area_name = lhd, year, strain, count)

# Aggregate: infectious syphilis only (primary analysis)
nsw_infectious <- nsw_lhd %>%
  filter(str_detect(strain, "infectious")) %>%
  group_by(state, geography_level, area_name, year) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

message("  NSW: ", n_distinct(nsw_infectious$area_name), " LHDs, ",
        min(nsw_infectious$year), "-", max(nsw_infectious$year))

# --- 6c. VIC LGA data ---

vic_files <- list.files(path_notifications, pattern = "VIC_syphilis_by_LGA_\\d{4}\\.csv",
                        full.names = TRUE)

vic_lga_data <- map_dfr(vic_files, function(f) {
  yr <- as.integer(str_extract(basename(f), "\\d{4}"))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(year = yr)
})

vic_lga_notifications <- vic_lga_data %>%
  group_by(LGA, year) %>%
  summarise(count = sum(Cases, na.rm = TRUE),
            population = sum(Population_ERP, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    state = "VIC",
    geography_level = "LGA",
    area_name = LGA
  ) %>%
  select(state, geography_level, area_name, year, count, population)

message("  VIC: ", n_distinct(vic_lga_notifications$area_name), " LGAs, ",
        min(vic_lga_notifications$year), "-", max(vic_lga_notifications$year))

# --- 6d. QLD HHS data ---

qld_hhs_rates <- read_csv(
  file.path(path_notifications, "QLD_infectious_syphilis_rate_by_HHS_2012-2024.csv"),
  show_col_types = FALSE
)

qld_hhs_long <- qld_hhs_rates %>%
  filter(!is.na(HHS)) %>%
  pivot_longer(cols = -c(Region, HHS), names_to = "year",
               values_to = "rate") %>%
  mutate(
    year = as.integer(year),
    state = "QLD",
    geography_level = "HHS",
    area_name = HHS
  ) %>%
  select(state, geography_level, area_name, region = Region, year, rate)

message("  QLD: ", n_distinct(qld_hhs_long$area_name), " HHS areas, ",
        min(qld_hhs_long$year), "-", max(qld_hhs_long$year))

# --- 6e. NT health district data ---

nt_district <- read_csv(
  file.path(path_notifications, "NT_infectious_syphilis_2024_by_district.csv"),
  show_col_types = FALSE
)

# NT quarterly time series
nt_quarterly <- read_csv(
  file.path(path_notifications, "NT_infectious_syphilis_quarterly_2019-2024.csv"),
  show_col_types = FALSE
)

message("  NT: district data for 2024, quarterly 2019-2024")

# --- 6f. WA summary ---

wa_summary <- read_csv(
  file.path(path_notifications, "WA_SORG_communique_summary_2024-2025.csv"),
  show_col_types = FALSE
)

message("  WA: statewide summary only (SORG communique)")

# --- 6g. NSW monthly infectious syphilis ---
# Converted from NSW_syphilis_by_month_2009-2026.xls via Python xlrd
# (original XLS unreadable by readxl due to format issues)
nsw_monthly_path <- file.path(path_notifications,
                               "NSW_infectious_syphilis_monthly_2009-2025.csv")
if (file.exists(nsw_monthly_path)) {
  nsw_monthly <- read_csv(nsw_monthly_path, show_col_types = FALSE)
  message("  NSW monthly: ", nrow(nsw_monthly), " observations (",
          min(nsw_monthly$year), "-", max(nsw_monthly$year), ")")
  saveRDS(nsw_monthly, file.path(path_processed, "nsw_monthly_notifications.rds"))
} else {
  message("  NSW monthly: CSV not found — run Python extraction first")
}

# ==============================================================================
# 7. HEALTH SERVICE LOCATIONS
# ==============================================================================

message("Loading health service locations...")

services <- read_csv(
  file.path(path_services, "NHSD_sexual_aboriginal_health_services_2025.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs_gda2020)

sh_clinics <- services %>% filter(service_category == "sexual_health")
ah_services <- services %>% filter(service_category == "aboriginal_health")

message("  Sexual health clinics: ", nrow(sh_clinics))
message("  Aboriginal health services: ", nrow(ah_services))

# ==============================================================================
# 8. SAVE ALL OBJECTS
# ==============================================================================

message("Saving prepared data objects...")

# Boundaries
saveRDS(sa2, file.path(path_processed, "boundaries_sa2.rds"))
saveRDS(sa3, file.path(path_processed, "boundaries_sa3.rds"))
saveRDS(ste, file.path(path_processed, "boundaries_ste.rds"))
saveRDS(phn, file.path(path_processed, "boundaries_phn.rds"))
saveRDS(qld_hhs, file.path(path_processed, "boundaries_qld_hhs.rds"))
saveRDS(vic_lga, file.path(path_processed, "boundaries_vic_lga.rds"))

# Demographics
saveRDS(sa2_population, file.path(path_processed, "census_sa2_population.rds"))
saveRDS(sa2_indigenous_detail, file.path(path_processed, "census_sa2_indigenous.rds"))
saveRDS(seifa_clean, file.path(path_processed, "seifa_sa2.rds"))
saveRDS(ra_sa2, file.path(path_processed, "remoteness_sa2.rds"))
saveRDS(erp_lga_clean, file.path(path_processed, "erp_lga.rds"))

# National/Kirby data
saveRDS(kirby_state_sex, file.path(path_processed, "kirby_state_sex.rds"))
saveRDS(kirby_state_indigenous, file.path(path_processed, "kirby_state_indigenous.rds"))
saveRDS(kirby_national_sex, file.path(path_processed, "kirby_national_sex.rds"))
saveRDS(kirby_national_age, file.path(path_processed, "kirby_national_age.rds"))
saveRDS(kirby_national_indigenous, file.path(path_processed, "kirby_national_indigenous.rds"))
saveRDS(kirby_national_remoteness, file.path(path_processed, "kirby_national_remoteness.rds"))
saveRDS(kirby_remoteness_indigenous, file.path(path_processed, "kirby_remoteness_indigenous.rds"))
saveRDS(kirby_remoteness_sex, file.path(path_processed, "kirby_remoteness_sex.rds"))

# Sub-state notifications
saveRDS(nsw_infectious, file.path(path_processed, "nsw_lhd_notifications.rds"))
saveRDS(vic_lga_notifications, file.path(path_processed, "vic_lga_notifications.rds"))
saveRDS(qld_hhs_long, file.path(path_processed, "qld_hhs_notifications.rds"))
saveRDS(nt_district, file.path(path_processed, "nt_district_notifications.rds"))
saveRDS(nt_quarterly, file.path(path_processed, "nt_quarterly_notifications.rds"))

# Services
saveRDS(services, file.path(path_processed, "health_services.rds"))
saveRDS(sh_clinics, file.path(path_processed, "sexual_health_clinics.rds"))
saveRDS(ah_services, file.path(path_processed, "aboriginal_health_services.rds"))

message("All data saved to ", path_processed)
message("Data acquisition complete.")

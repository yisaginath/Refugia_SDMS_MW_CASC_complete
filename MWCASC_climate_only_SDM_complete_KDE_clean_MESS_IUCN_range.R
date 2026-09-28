# Script: Predicting and projecting biodiversity under climate change
# Author: Yisa Ginath Yuh
# contact: yisaginath80@yahoo.com
##--------------------------------------------------------------------------------------

library(dismo) # for maxent modelling
library(raster) # for raster creation, manipulation, and analysis
library(dplyr) # tidy data
library(magrittr) # for streamlining manipulated data and making codes easily readable with the pipe operator
library(rgdal) # for handling spatial data in different formats and facilitating projection
library(sp) # for handling spatial data, especially vector layers
library(maxnet) # for MaxEnt modeliing
library(ggplot2) # plotting results
library(rnaturalearth)# for handling background map data
library(rnaturalearthdata) # for handling background map data
Sys.setenv("PROJ_NETWORK" = "ON")
library(sf)# for handling spatial data in different formats and facilitating projection
library(terra) #for raster creation, manipulation, and analysis
library(blockCV)# for spatial cross-validation with presence-absence data
library(spatialEco)# provide importnt functions for manipulating, querying, and modeling spatial datasets in ecology
library(pROC)
library(tibble)
library(blockCV)
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)

set.seed(1)

##--------------------------------------------------------------------------------------

# The MaxEnt model fitting requires java, and it is very imperative ti allocate a maximum memory space with java, so as to reduce memory issues when running java operations (e.g., rjava).
# Ensure to set java heap space based on your computer RAM.
# Here, I am setting my Java heap space to 8GB. users can modify to different GB limits
options(java.parameters = '-Xmx8g')


# ============================================================
# 1. SETTINGS AND OUTPUT DIRECTORIES
# ============================================================

taxa_group <- "mammals_IUCN_range"

base_output_dir <- file.path(
  "D:/yisa/UI/uillinois_sdm",
  taxa_group
)

dirs <- list(
  kde_raster_dir = file.path(base_output_dir, "kde_rasters"),
  kde_group_plot_dir = file.path(base_output_dir, "kde_group_plots"),
  kde_species_plot_dir = file.path(base_output_dir, "kde_species_plots"),
  current_fold_dir = file.path(base_output_dir, "current_fold_predictions"),
  current_avg_dir = file.path(base_output_dir, "current_average_predictions"),
  current_binary_dir = file.path(base_output_dir, "current_binary_predictions"),
  future_avg_dir = file.path(base_output_dir, "future_average_predictions_full_extent"),
  future_binary_dir = file.path(base_output_dir, "future_binary_predictions_full_extent"),
  mess_dir = file.path(base_output_dir, "MESS_future_novelty"),
  variable_importance_dir = file.path(base_output_dir, "variable_importance"),
  performance_dir = file.path(base_output_dir, "model_performance"),
  
  # ----------------------------------------------------------
  # NEW:
  # Directory for saving fitted MaxEnt model objects
  # ----------------------------------------------------------
  model_object_dir = file.path(base_output_dir, "model_objects")
)

lapply(
  dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

safe_name <- function(x) {
  x %>%
    trimws() %>%
    gsub("[^A-Za-z0-9]+", "_", .) %>%
    gsub("_+$", "", .) %>%
    gsub("^_+", "", .)
}


# Ensure MESS directory exists
if (is.null(dirs$mess_dir) || is.na(dirs$mess_dir)) {
  stop("ERROR: dirs$mess_dir is not defined properly.")
}

dir.create(
  dirs$mess_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Ensure model-object directory exists
if (
  is.null(dirs$model_object_dir) ||
  is.na(dirs$model_object_dir)
) {
  
  stop(
    "ERROR: dirs$model_object_dir is not defined properly."
  )
}

dir.create(
  dirs$model_object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. READ STUDY AREA AND GREAT LAKES
# ============================================================

studyarea <- st_read(
  "D:/yisa/UI/study_area/study_area_projected.shp",
  quiet = TRUE
)

great_lake <- st_read(
  "D:/yisa/UI/great_lakes/great_lakes_CASC.shp",
  quiet = TRUE
)


# ============================================================
# 3. READ CURRENT CLIMATE
# ============================================================

climate_current_dir <- "D:/yisa/UI/climate_data/historical"

current_files <- list.files(
  climate_current_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

climate_current <- rast(
  current_files
)

names(climate_current) <- gsub(
  "masked_",
  "",
  tools::file_path_sans_ext(
    basename(current_files)
  )
)

template <- climate_current[[1]]

studyarea_clim <- st_transform(
  studyarea,
  crs(template)
)

great_lake_clim <- st_transform(
  great_lake,
  crs(template)
)


# ------------------------------------------------------------
# Create one dissolved study-area boundary
# ------------------------------------------------------------

studyarea_union <- st_union(
  st_geometry(
    st_make_valid(studyarea_clim)
  )
)

studyarea_union_sf <- st_sf(
  study_area = "CASC",
  geometry = studyarea_union
)


# ============================================================
# 4. READ FUTURE CLIMATE STACKS
# ============================================================

base_climate_dir <- "D:/yisa/UI/climate_data"

future_dirs <- list(
  "2050_ssp126" = file.path(
    base_climate_dir,
    "2041_2070_ssp126"
  ),
  "2050_ssp245" = file.path(
    base_climate_dir,
    "2041_2070_ssp245"
  ),
  "2050_ssp585" = file.path(
    base_climate_dir,
    "2041_2070_ssp585"
  ),
  "2100_ssp126" = file.path(
    base_climate_dir,
    "2071_2100_ssp126"
  ),
  "2100_ssp245" = file.path(
    base_climate_dir,
    "2071_2100_ssp245"
  ),
  "2100_ssp585" = file.path(
    base_climate_dir,
    "2071_2100_ssp585"
  )
)

read_future_stack <- function(folder) {
  
  f <- list.files(
    folder,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  r <- rast(
    f
  )
  
  names(r) <- gsub(
    "masked_",
    "",
    tools::file_path_sans_ext(
      basename(f)
    )
  )
  
  r <- r[[names(climate_current)]]
  
  return(r)
}

future_climate_stacks <- lapply(
  future_dirs,
  read_future_stack
)


# ============================================================
# 5. READ DEM, ALIGN TO CLIMATE, COMPUTE SLOPE
# ============================================================

dem_raw <- rast(
  "D:/yisa/UI/study_area/DEM/DEM_cropped_to_study_area_30s.tif"
)

dem_1km <- project(
  dem_raw,
  template,
  method = "bilinear"
)

names(dem_1km) <- "elevation"

slope_1km <- terrain(
  dem_1km,
  v = "slope",
  unit = "degrees",
  neighbors = 8
)

names(slope_1km) <- "slope"


# ============================================================
# 6. GREAT LAKES MASK FUNCTION
# ============================================================

mask_great_lakes <- function(
    x,
    great_lake_sf
) {
  
  lake_vect <- vect(
    great_lake_sf
  )
  
  lake_raster <- rasterize(
    lake_vect,
    x[[1]],
    field = 1,
    background = NA
  )
  
  x_masked <- mask(
    x,
    lake_raster,
    maskvalues = 1,
    updatevalue = NA
  )
  
  return(
    x_masked
  )
}


# ============================================================
# 7. BUILD CURRENT AND FUTURE ENVIRONMENTAL STACKS
# ============================================================

current_env_full <- c(
  climate_current,
  dem_1km,
  slope_1km
)

current_env_full <- crop(
  current_env_full,
  vect(studyarea_clim)
)

current_env_full <- mask(
  current_env_full,
  vect(studyarea_clim)
)

current_env_full <- mask_great_lakes(
  current_env_full,
  great_lake_clim
)

future_env_stacks <- lapply(
  future_climate_stacks,
  function(x) {
    
    fx <- c(
      x,
      dem_1km,
      slope_1km
    )
    
    fx <- crop(
      fx,
      vect(studyarea_clim)
    )
    
    fx <- mask(
      fx,
      vect(studyarea_clim)
    )
    
    fx <- mask_great_lakes(
      fx,
      great_lake_clim
    )
    
    names(fx) <- names(
      current_env_full
    )
    
    return(fx)
  }
)

cat(
  "\nCurrent environmental variables:\n"
)

print(
  names(current_env_full)
)


# ============================================================
# 8. READ MAMMAL OBSERVATIONS AND SPECIES SUMMARY
# ============================================================

species_observations <- read.csv(
  "D:/yisa/UI/recent_cleaned_species_observation/mammal_observations_cleaned.csv"
)

species_summary <- read.csv(
  "D:/yisa/UI/recent_cleaned_species_observation/mammal_summary.csv"
)


species_observations <- species_observations %>%
  mutate(
    scientific_name = paste(
      genus,
      species
    ),
    scientific_name_clean = trimws(
      tolower(scientific_name)
    ),
    group = taxa_group
  )

species_summary <- species_summary %>%
  mutate(
    scientific_name = paste(
      genus,
      species
    ),
    scientific_name_clean = trimws(
      tolower(scientific_name)
    ),
    species_safe = safe_name(
      scientific_name
    ),
    mean_auc = NA_real_,
    mean_threshold = NA_real_,
    n_final_records = NA_integer_,
    n_background = NA_integer_
  ) %>%
  filter(
    n_MWCASC >= 20
  )

modelPerformanceAllSpecies <- species_summary

for (v in names(current_env_full)) {
  
  modelPerformanceAllSpecies[[v]] <- NA_real_
  
}


# ============================================================
# 8A. READ MATCHED IUCN SPECIES RANGE
# ============================================================

# ------------------------------------------------------------
# IMPORTANT:
#
# Current calibration, background sampling and current
# predictions will now be restricted to these matched IUCN
# species ranges.
#
# The future predictions will NOT be restricted to these
# current ranges. Future predictions remain across the full
# CASC study extent.
# ------------------------------------------------------------

matched_range_file <- file.path(
  "D:/yisa/UI/study_area/IUCN_species_CASC",
  "mammals_IUCN_CASC_names_matched.shp"
)


matched_iucn_range <- st_read(
  matched_range_file,
  quiet = TRUE
)


# ------------------------------------------------------------
# Verify scientific-name field
# ------------------------------------------------------------

if (
  !"sci_name" %in% names(
    matched_iucn_range
  )
) {
  
  stop(
    "ERROR: Matched IUCN range does not contain sci_name."
  )
}


# ------------------------------------------------------------
# Standardize range names for matching
# ------------------------------------------------------------

matched_iucn_range <- matched_iucn_range %>%
  mutate(
    sci_name = trimws(
      as.character(sci_name)
    ),
    scientific_name_clean = trimws(
      tolower(sci_name)
    )
  )


# ------------------------------------------------------------
# Repair geometries
# ------------------------------------------------------------

matched_iucn_range <- st_make_valid(
  matched_iucn_range
)


# ------------------------------------------------------------
# Remove empty geometries
# ------------------------------------------------------------

matched_iucn_range <- matched_iucn_range[
  !st_is_empty(matched_iucn_range),
]


# ------------------------------------------------------------
# Transform matched ranges to environmental-raster CRS
# ------------------------------------------------------------

matched_iucn_range <- st_transform(
  matched_iucn_range,
  crs(template)
)


# ------------------------------------------------------------
# Print matching information
# ------------------------------------------------------------

cat(
  "\n====================================================\n",
  "MATCHED IUCN RANGE DATA\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Range shapefile:\n",
  matched_range_file,
  "\n\n",
  sep = ""
)

cat(
  "Number of polygon features: ",
  nrow(matched_iucn_range),
  "\n",
  sep = ""
)

cat(
  "Number of unique IUCN species: ",
  length(
    unique(
      matched_iucn_range$scientific_name_clean
    )
  ),
  "\n",
  sep = ""
)


# ============================================================
# 9. CREATE GROUP-LEVEL KDE BIAS SURFACE
# ============================================================

environmental_mask <- current_env_full[[1]]

obs_sf_all <- species_observations %>%
  filter(
    !is.na(longitude),
    !is.na(latitude)
  ) %>%
  st_as_sf(
    coords = c(
      "longitude",
      "latitude"
    ),
    crs = 4326
  ) %>%
  st_transform(
    crs(template)
  )

inside_study <- st_intersects(
  obs_sf_all,
  studyarea_clim,
  sparse = FALSE
)

obs_sf_all <- obs_sf_all[
  which(
    rowSums(inside_study) > 0
  ),
]

obs_sf_sample <- obs_sf_all %>%
  slice_sample(
    n = min(
      100000,
      nrow(.)
    )
  )

coords <- st_coordinates(
  obs_sf_sample
)

coords <- unique(
  coords
)

kde_window <- as.owin(
  st_union(
    st_geometry(
      studyarea_clim
    )
  )
)

pp <- ppp(
  x = coords[, 1],
  y = coords[, 2],
  window = kde_window,
  check = TRUE
)

kde <- density.ppp(
  pp,
  edge = TRUE,
  at = "pixels",
  eps = terra::res(
    environmental_mask
  )
)

kde_vals <- kde$v[
  nrow(kde$v):1,
]

kde_raster <- rast(
  kde_vals,
  extent = ext(
    kde$xrange[1],
    kde$xrange[2],
    kde$yrange[1],
    kde$yrange[2]
  ),
  crs = crs(template)
)

kde_raster <- project(
  kde_raster,
  environmental_mask
)

kde_raster <- crop(
  kde_raster,
  environmental_mask
)

kde_raster <- mask(
  kde_raster,
  environmental_mask
)

kde_minmax <- minmax(
  kde_raster
)

kde_bias <- (
  kde_raster - kde_minmax[1, 1]
) /
  (
    kde_minmax[2, 1] -
      kde_minmax[1, 1]
  )

names(kde_bias) <- "kde_bias"

writeRaster(
  kde_bias,
  file.path(
    dirs$kde_raster_dir,
    "kde_bias_mammals_0_1_scaled.tif"
  ),
  overwrite = TRUE
)

png(
  file.path(
    dirs$kde_group_plot_dir,
    "kde_group_mammals.png"
  ),
  width = 1800,
  height = 1600,
  res = 200
)

plot(
  kde_bias,
  main = "Group-level KDE bias surface: mammals"
)

plot(
  st_geometry(
    obs_sf_sample
  ),
  add = TRUE,
  pch = 20,
  cex = 0.25
)

dev.off()


# ============================================================
# 10. MODEL EACH MAMMAL SPECIES
# ============================================================

for (
  i in seq_len(
    nrow(
      modelPerformanceAllSpecies
    )
  )
) {
  
  scientific_name_i <-
    modelPerformanceAllSpecies$scientific_name[i]
  
  scientific_name_clean_i <-
    modelPerformanceAllSpecies$scientific_name_clean[i]
  
  species_safe_i <-
    modelPerformanceAllSpecies$species_safe[i]
  
  message(
    "\n================================================="
  )
  
  message(
    "Processing: ",
    scientific_name_i
  )
  
  message(
    "================================================="
  )
  
  
  # ------------------------------------------------------------
  # Extract observations for current species
  # ------------------------------------------------------------
  
  species_i_obs <- species_observations %>%
    filter(
      scientific_name_clean ==
        scientific_name_clean_i
    ) %>%
    filter(
      !is.na(longitude),
      !is.na(latitude)
    )
  
  
  if (
    "square_id" %in% names(
      species_i_obs
    )
  ) {
    
    species_i_obs <- species_i_obs %>%
      filter(
        !is.na(square_id)
      ) %>%
      group_by(
        square_id
      ) %>%
      slice(1) %>%
      ungroup()
    
  } else {
    
    species_i_obs <- species_i_obs %>%
      distinct(
        longitude,
        latitude,
        .keep_all = TRUE
      )
  }
  
  
  species_i_sf <- species_i_obs %>%
    st_as_sf(
      coords = c(
        "longitude",
        "latitude"
      ),
      crs = 4326
    ) %>%
    st_transform(
      crs(template)
    )
  
  
  # ------------------------------------------------------------
  # Retain observations within CASC study area
  # ------------------------------------------------------------
  
  inside_species <- st_intersects(
    species_i_sf,
    studyarea_clim,
    sparse = FALSE
  )
  
  species_i_sf <- species_i_sf[
    which(
      rowSums(
        inside_species
      ) > 0
    ),
  ]
  
  
  # ============================================================
  # NEW:
  # EXTRACT MATCHED IUCN RANGE FOR CURRENT SPECIES
  # ============================================================
  
  species_range_i <- matched_iucn_range %>%
    filter(
      scientific_name_clean ==
        scientific_name_clean_i
    )
  
  
  # ------------------------------------------------------------
  # Skip species with no matched IUCN range
  # ------------------------------------------------------------
  
  if (
    nrow(species_range_i) == 0
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": no matched IUCN range found."
    )
    
    next
  }
  
  
  # ------------------------------------------------------------
  # Repair range geometry
  # ------------------------------------------------------------
  
  species_range_i <- st_make_valid(
    species_range_i
  )
  
  
  species_range_i <- species_range_i[
    !st_is_empty(species_range_i),
  ]
  
  
  if (
    nrow(species_range_i) == 0
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": IUCN range geometry is empty."
    )
    
    next
  }
  
  
  # ------------------------------------------------------------
  # Dissolve all IUCN polygons belonging to this species
  # ------------------------------------------------------------
  
  species_range_geom <- st_union(
    st_geometry(
      species_range_i
    )
  )
  
  
  species_range <- st_sf(
    scientific_name = scientific_name_i,
    geometry = species_range_geom
  )
  
  
  # ------------------------------------------------------------
  # Ensure range is constrained to CASC study area
  # ------------------------------------------------------------
  
  species_range <- suppressWarnings(
    st_intersection(
      species_range,
      studyarea_union_sf
    )
  )
  
  
  species_range <- st_make_valid(
    species_range
  )
  
  
  species_range <- species_range[
    !st_is_empty(species_range),
  ]
  
  
  if (
    nrow(species_range) == 0
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": IUCN range does not overlap study area."
    )
    
    next
  }
  
  
  # ============================================================
  # NEW:
  # RETAIN CURRENT OCCURRENCES WITHIN THE MATCHED IUCN RANGE
  # ============================================================
  
  inside_iucn_range <- lengths(
    st_intersects(
      species_i_sf,
      species_range
    )
  ) > 0
  
  
  species_i_sf <- species_i_sf[
    inside_iucn_range,
  ]
  
  
  # ------------------------------------------------------------
  # Require at least 20 occurrence records within IUCN range
  # ------------------------------------------------------------
  
  if (
    nrow(species_i_sf) < 20
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": fewer than 20 records within matched IUCN range."
    )
    
    next
  }
  
  
  modelPerformanceAllSpecies$n_final_records[i] <-
    nrow(
      species_i_sf
    )
  
  
  # ------------------------------------------------------------
  # Species modelling extent = MATCHED IUCN RANGE
  #
  # PREVIOUSLY:
  # Species range was a 50-km buffer around occurrences.
  #
  # NOW:
  # Current calibration and prediction are restricted to the
  # matched IUCN range.
  # ------------------------------------------------------------
  
  current_env_i <- crop(
    current_env_full,
    vect(species_range)
  )
  
  current_env_i <- mask(
    current_env_i,
    vect(species_range)
  )
  
  
  # ------------------------------------------------------------
  # Check that current environmental data exist in range
  # ------------------------------------------------------------
  
  current_valid_cells <- global(
    !is.na(
      current_env_i[[1]]
    ),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  
  if (
    is.na(current_valid_cells) ||
    current_valid_cells == 0
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": no environmental cells within IUCN range."
    )
    
    next
  }
  
  
  # ------------------------------------------------------------
  # Extract mean environment at occurrence points
  # ------------------------------------------------------------
  
  env_vals_presence <- terra::extract(
    current_env_i,
    vect(species_i_sf)
  ) %>%
    dplyr::select(
      -ID
    )
  
  
  env_means <- colMeans(
    env_vals_presence,
    na.rm = TRUE
  )
  
  
  modelPerformanceAllSpecies[
    i,
    names(current_env_full)
  ] <- env_means
  
  
  # ------------------------------------------------------------
  # Species-specific KDE
  #
  # KDE is now cropped and masked to matched IUCN range.
  # ------------------------------------------------------------
  
  kde_i <- crop(
    kde_bias,
    current_env_i
  )
  
  kde_i <- mask(
    kde_i,
    current_env_i[[1]]
  )
  
  
  if (
    all(
      is.na(
        values(kde_i)
      )
    )
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": empty KDE raster."
    )
    
    next
  }
  
  
  png(
    file.path(
      dirs$kde_species_plot_dir,
      paste0(
        "kde_species_",
        species_safe_i,
        ".png"
      )
    ),
    width = 1800,
    height = 1600,
    res = 200
  )
  
  
  plot(
    kde_i,
    main = paste(
      "KDE bias + observations:",
      scientific_name_i
    )
  )
  
  
  plot(
    st_geometry(
      species_range
    ),
    add = TRUE,
    border = "black",
    lwd = 2
  )
  
  
  plot(
    st_geometry(
      species_i_sf
    ),
    add = TRUE,
    pch = 20,
    col = "red",
    cex = 0.5
  )
  
  
  dev.off()
  
  
  # ------------------------------------------------------------
  # Background sampling using KDE bias
  #
  # Since kde_i is masked to the matched IUCN range,
  # all background points are sampled within the IUCN range.
  # ------------------------------------------------------------
  
  range_pixels <- global(
    (
      (kde_i * 0) + 1
    ),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  
  backg_n <- if (
    range_pixels < 100000
  ) {
    
    round(
      range_pixels * 0.20
    )
    
  } else if (
    range_pixels < 250000
  ) {
    
    round(
      range_pixels * 0.10
    )
    
  } else {
    
    round(
      range_pixels * 0.05
    )
  }
  
  
  #backg_n <- max(backg_n, 1000)
  
  
  # ------------------------------------------------------------
  # Ensure background sample is valid
  # ------------------------------------------------------------
  
  backg_n <- min(
    backg_n,
    floor(
      range_pixels
    )
  )
  
  
  if (
    is.na(backg_n) ||
    backg_n < 1
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": invalid number of background points."
    )
    
    next
  }
  
  
  modelPerformanceAllSpecies$n_background[i] <-
    backg_n
  
  
  backg <- dismo::randomPoints(
    mask = raster(kde_i),
    n = backg_n,
    prob = TRUE
  )
  
  
  backg_df <- as.data.frame(
    backg
  )
  
  
  colnames(
    backg_df
  ) <- c(
    "x",
    "y"
  )
  
  
  backg_sf <- st_as_sf(
    backg_df,
    coords = c(
      "x",
      "y"
    ),
    crs = crs(template)
  ) %>%
    mutate(
      p = 0
    )
  
  
  pres_sf <- species_i_sf %>%
    mutate(
      p = 1
    )
  
  
  pres_backg_sf <- rbind(
    pres_sf[, "p"],
    backg_sf[, "p"]
  )
  
  
  pb_extract <- terra::extract(
    current_env_i,
    vect(pres_backg_sf)
  ) %>%
    dplyr::select(
      -ID
    ) %>%
    as.data.frame()
  
  
  good_rows <- complete.cases(
    pb_extract
  )
  
  
  pres_backg_sf <- pres_backg_sf[
    good_rows,
  ]
  
  
  pb_extract <- pb_extract[
    good_rows,
  ]
  
  
  if (
    length(
      unique(
        pres_backg_sf$p
      )
    ) < 2
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": invalid presence/background data."
    )
    
    next
  }
  
  
  # Reference environment used for MESS
  # This represents the environmental space available to
  # the current model within the matched IUCN range.
  
  mess_reference <- pb_extract
  
  
  # ------------------------------------------------------------
  # Spatial block cross-validation
  # ------------------------------------------------------------
  
  pres_backg_sf$X <-
    st_coordinates(
      pres_backg_sf
    )[, 1]
  
  
  pres_backg_sf$Y <-
    st_coordinates(
      pres_backg_sf
    )[, 2]
  
  
  sb <- spatialBlock(
    speciesData = pres_backg_sf,
    species = "p",
    rasterLayer = raster::stack(
      current_env_i
    ),
    selection = "systematic",
    rows = 10,
    cols = 10,
    k = 5,
    biomod2Format = TRUE,
    verbose = FALSE
  )
  
  
  fold_id <- sb$foldID
  
  
  AUCs_i <- rep(
    NA_real_,
    5
  )
  
  
  thresholds_i <- rep(
    NA_real_,
    5
  )
  
  
  models_i <- vector(
    "list",
    5
  )
  
  
  preds_i <- vector(
    "list",
    5
  )
  
  
  # ------------------------------------------------------------
  # NEW:
  # Create species-specific model-object folder
  # ------------------------------------------------------------
  
  species_model_dir <- file.path(
    dirs$model_object_dir,
    species_safe_i
  )
  
  
  dir.create(
    species_model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ------------------------------------------------------------
  # Fit MaxEnt models
  # ------------------------------------------------------------
  
  all_var_imp <- list()
  
  
  for (
    k in 1:5
  ) {
    
    trainSet <- which(
      fold_id != k
    )
    
    
    testSet <- which(
      fold_id == k
    )
    
    
    if (
      length(trainSet) == 0 ||
      length(testSet) == 0
    ) {
      
      next
    }
    
    
    trainEnv <- pb_extract[
      trainSet,
    ]
    
    
    testEnv <- pb_extract[
      testSet,
    ]
    
    
    train_p <- pres_backg_sf$p[
      trainSet
    ]
    
    
    test_p <- pres_backg_sf$p[
      testSet
    ]
    
    
    if (
      length(
        unique(train_p)
      ) < 2 ||
      length(
        unique(test_p)
      ) < 2
    ) {
      
      message(
        "Fold ",
        k,
        " skipped: lacks both presence and background."
      )
      
      next
    }
    
    
    maxent_model <- dismo::maxent(
      x = trainEnv,
      p = train_p,
      args = c(
        "outputformat=cloglog",
        "betamultiplier=1",
        "maximumiterations=1000"
      )
    )
    
    
    models_i[[k]] <- maxent_model
    
    
    # ========================================================
    # NEW:
    # SAVE INDIVIDUAL FOLD MODEL OBJECT
    # ========================================================
    
    saveRDS(
      maxent_model,
      file = file.path(
        species_model_dir,
        paste0(
          "maxent_model_",
          species_safe_i,
          "_fold_",
          k,
          ".rds"
        )
      )
    )
    
    
    var_imp <- maxent_model@results
    
    
    var_imp_rows <- grep(
      "contribution$|permutation.importance$",
      rownames(var_imp),
      value = TRUE
    )
    
    
    # extract fold importance
    
    fold_imp <- var_imp[
      var_imp_rows,
      ,
      drop = FALSE
    ]
    
    
    # save fold result
    
    write.csv(
      fold_imp,
      file.path(
        dirs$variable_importance_dir,
        paste0(
          "variable_importance_",
          species_safe_i,
          "_fold_",
          k,
          ".csv"
        )
      )
    )
    
    
    # store for aggregation
    
    fold_imp_df <- data.frame(
      variable = rownames(
        fold_imp
      ),
      importance = fold_imp[, 1],
      fold = k
    )
    
    
    all_var_imp[[k]] <- fold_imp_df
    
    
    # --------------------------------------------------------
    # Current prediction is ONLY within matched IUCN range
    # because current_env_i is range-masked.
    # --------------------------------------------------------
    
    pred_k <- terra::predict(
      current_env_i,
      maxent_model,
      na.rm = TRUE
    )
    
    
    preds_i[[k]] <- pred_k
    
    
    writeRaster(
      pred_k,
      file.path(
        dirs$current_fold_dir,
        paste0(
          "current_habitat_suitability_",
          species_safe_i,
          "_fold_",
          k,
          ".tif"
        )
      ),
      overwrite = TRUE
    )
    
    
    eval_k <- dismo::evaluate(
      p = testEnv[
        test_p == 1,
      ],
      a = testEnv[
        test_p == 0,
      ],
      model = maxent_model
    )
    
    
    AUCs_i[k] <- eval_k@auc
    
    
    thresholds_i[k] <- tryCatch(
      dismo::threshold(
        eval_k,
        "spec_sens"
      ),
      error = function(e) {
        NA_real_
      }
    )
    
    
    message(
      "Fold ",
      k,
      " AUC: ",
      round(
        AUCs_i[k],
        3
      )
    )
  }
  
  
  valid_folds <- which(
    !sapply(
      preds_i,
      is.null
    ) &
      !is.na(
        AUCs_i
      )
  )
  
  
  if (
    length(valid_folds) == 0
  ) {
    
    message(
      "Skipping ",
      scientific_name_i,
      ": no valid folds."
    )
    
    next
  }
  
  
  # ------------------------------------------------------------
  # AUC-weighted current average prediction
  #
  # Prediction remains restricted to matched IUCN range.
  # ------------------------------------------------------------
  
  auc_weights <- (
    AUCs_i[valid_folds] - 0.5
  )^2
  
  
  if (
    sum(
      auc_weights,
      na.rm = TRUE
    ) == 0
  ) {
    
    auc_weights <- rep(
      1,
      length(valid_folds)
    )
  }
  
  
  current_pred_stack <- rast(
    preds_i[
      valid_folds
    ]
  )
  
  
  current_pred_avg <- app(
    current_pred_stack,
    fun = function(x) {
      
      if (
        all(
          is.na(x)
        )
      ) {
        
        return(
          NA
        )
        
      } else {
        
        return(
          weighted.mean(
            x,
            w = auc_weights,
            na.rm = TRUE
          )
        )
      }
    }
  )
  
  
  names(
    current_pred_avg
  ) <- "current_weighted_avg_suitability"
  
  
  writeRaster(
    current_pred_avg,
    file.path(
      dirs$current_avg_dir,
      paste0(
        "current_weighted_avg_habitat_suitability_",
        species_safe_i,
        ".tif"
      )
    ),
    overwrite = TRUE
  )
  
  
  avg_auc <- mean(
    AUCs_i[
      valid_folds
    ],
    na.rm = TRUE
  )
  
  
  modelPerformanceAllSpecies$mean_auc[i] <-
    avg_auc
  
  
  valid_thresholds <- which(
    !is.na(
      thresholds_i
    ) &
      !is.na(
        AUCs_i
      )
  )
  
  
  if (
    length(valid_thresholds) > 0
  ) {
    
    threshold_weights <- (
      AUCs_i[
        valid_thresholds
      ] - 0.5
    )^2
    
    
    if (
      sum(
        threshold_weights,
        na.rm = TRUE
      ) == 0
    ) {
      
      species_threshold <- mean(
        thresholds_i[
          valid_thresholds
        ],
        na.rm = TRUE
      )
      
    } else {
      
      species_threshold <- weighted.mean(
        thresholds_i[
          valid_thresholds
        ],
        w = threshold_weights,
        na.rm = TRUE
      )
    }
    
  } else {
    
    species_threshold <- NA_real_
  }
  
  
  modelPerformanceAllSpecies$mean_threshold[i] <-
    species_threshold
  
  
  if (
    !is.na(
      species_threshold
    )
  ) {
    
    current_binary <-
      current_pred_avg >
      species_threshold
    
    
    writeRaster(
      current_binary,
      file.path(
        dirs$current_binary_dir,
        paste0(
          "current_binary_habitat_suitability_",
          species_safe_i,
          ".tif"
        )
      ),
      overwrite = TRUE
    )
  }
  
  
  # ============================================================
  # NEW:
  # SAVE COMPLETE SPECIES MODEL OBJECT
  # ============================================================
  
  species_model_object <- list(
    
    scientific_name =
      scientific_name_i,
    
    scientific_name_clean =
      scientific_name_clean_i,
    
    species_safe =
      species_safe_i,
    
    taxonomic_group =
      taxa_group,
    
    predictor_names =
      names(current_env_full),
    
    valid_folds =
      valid_folds,
    
    models =
      models_i,
    
    fold_auc =
      AUCs_i,
    
    fold_threshold =
      thresholds_i,
    
    auc_weights =
      auc_weights,
    
    mean_auc =
      avg_auc,
    
    final_threshold =
      species_threshold,
    
    n_occurrences =
      nrow(species_i_sf),
    
    n_background =
      backg_n,
    
    calibration_extent =
      "matched IUCN species range",
    
    future_projection_extent =
      "full CASC study area",
    
    matched_iucn_name =
      unique(
        species_range_i$sci_name
      ),
    
    model_arguments = c(
      "outputformat=cloglog",
      "betamultiplier=1",
      "maximumiterations=1000"
    )
  )
  
  
  saveRDS(
    species_model_object,
    file = file.path(
      species_model_dir,
      paste0(
        "maxent_model_object_",
        species_safe_i,
        "_ALL_FOLDS.rds"
      )
    )
  )
  
  
  message(
    "Model object saved: ",
    file.path(
      species_model_dir,
      paste0(
        "maxent_model_object_",
        species_safe_i,
        "_ALL_FOLDS.rds"
      )
    )
  )
  
  
  # ------------------------------------------------------------
  # Future projections across full study extent
  #
  # IMPORTANT:
  #
  # THIS PART REMAINS UNCHANGED IN SPATIAL EXTENT.
  #
  # Future environmental stacks are NOT cropped or masked to
  # the current IUCN range.
  #
  # Future predictions therefore identify suitable habitat
  # anywhere across the full CASC study area.
  #
  # Future binary uses current threshold.
  # ------------------------------------------------------------
  
  for (
    scenario_name in names(
      future_env_stacks
    )
  ) {
    
    message(
      "Future projection: ",
      scientific_name_i,
      " | ",
      scenario_name
    )
    
    
    # ------------------------------------------------------------
    # MESS analysis for future environmental novelty
    # Negative MESS values = novel climate conditions
    # ------------------------------------------------------------
    
    scenario_mess_dir <- file.path(
      dirs$mess_dir,
      scenario_name
    )
    
    
    dir.create(
      scenario_mess_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    future_mess <- dismo::mess(
      x = raster::stack(
        future_env_stacks[[scenario_name]]
      ),
      v = mess_reference
    )
    
    
    raster::writeRaster(
      future_mess,
      filename = file.path(
        scenario_mess_dir,
        paste0(
          "MESS_",
          species_safe_i,
          "_",
          scenario_name,
          ".tif"
        )
      ),
      overwrite = TRUE
    )
    
    
    future_novel <- future_mess < 0
    
    
    raster::writeRaster(
      future_novel,
      filename = file.path(
        scenario_mess_dir,
        paste0(
          "novel_climate_binary_",
          species_safe_i,
          "_",
          scenario_name,
          ".tif"
        )
      ),
      overwrite = TRUE
    )
    
    
    future_preds <- list()
    
    
    for (
      k in valid_folds
    ) {
      
      # --------------------------------------------------------
      # IMPORTANT:
      # Predict on full future environmental stack.
      # No IUCN range mask is applied here.
      # --------------------------------------------------------
      
      future_preds[[as.character(k)]] <- terra::predict(
        future_env_stacks[[scenario_name]],
        models_i[[k]],
        na.rm = TRUE
      )
    }
    
    
    future_pred_stack <- rast(
      future_preds
    )
    
    
    future_pred_avg <- app(
      future_pred_stack,
      fun = function(x) {
        
        if (
          all(
            is.na(x)
          )
        ) {
          
          return(
            NA
          )
          
        } else {
          
          return(
            weighted.mean(
              x,
              w = auc_weights,
              na.rm = TRUE
            )
          )
        }
      }
    )
    
    
    names(
      future_pred_avg
    ) <- paste0(
      "future_",
      scenario_name
    )
    
    
    scenario_avg_dir <- file.path(
      dirs$future_avg_dir,
      scenario_name
    )
    
    
    scenario_binary_dir <- file.path(
      dirs$future_binary_dir,
      scenario_name
    )
    
    
    dir.create(
      scenario_avg_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    dir.create(
      scenario_binary_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    writeRaster(
      future_pred_avg,
      file.path(
        scenario_avg_dir,
        paste0(
          "future_weighted_avg_habitat_suitability_",
          species_safe_i,
          "_",
          scenario_name,
          ".tif"
        )
      ),
      overwrite = TRUE
    )
    
    
    if (
      !is.na(
        species_threshold
      )
    ) {
      
      future_binary <-
        future_pred_avg >
        species_threshold
      
      
      writeRaster(
        future_binary,
        file.path(
          scenario_binary_dir,
          paste0(
            "future_binary_habitat_suitability_",
            species_safe_i,
            "_",
            scenario_name,
            ".tif"
          )
        ),
        overwrite = TRUE
      )
    }
    
    
    # ----------------------------------------------------------
    # Remove large scenario objects before next future scenario
    # ----------------------------------------------------------
    
    rm(
      future_mess,
      future_novel,
      future_preds,
      future_pred_stack,
      future_pred_avg
    )
    
    
    if (
      exists(
        "future_binary"
      )
    ) {
      
      rm(
        future_binary
      )
    }
    
    
    gc()
  }
  
  
  # ============================================================
  # 11. SAVE MODEL PERFORMANCE TABLE
  # ============================================================
  
  performance_path <- file.path(
    dirs$performance_dir,
    paste0(
      "model_performance_summary_",
      taxa_group,
      ".csv"
    )
  )
  
  
  write.csv(
    modelPerformanceAllSpecies,
    performance_path,
    row.names = FALSE
  )
  
  
  message(
    "Model performance table saved to: ",
    performance_path
  )
  
  
  # ------------------------------------------------------------
  # Combine all successfully generated variable importance
  # results
  # ------------------------------------------------------------
  
  all_var_imp_valid <- all_var_imp[
    !sapply(
      all_var_imp,
      is.null
    )
  ]
  
  
  if (
    length(
      all_var_imp_valid
    ) > 0
  ) {
    
    # combine all folds
    
    all_var_imp_df <- do.call(
      rbind,
      all_var_imp_valid
    )
    
    
    # average across folds
    
    final_var_imp <- aggregate(
      importance ~ variable,
      data = all_var_imp_df,
      FUN = mean
    )
    
    
    # save final species-level importance
    
    write.csv(
      final_var_imp,
      file.path(
        dirs$variable_importance_dir,
        paste0(
          "variable_importance_",
          species_safe_i,
          "_ALL_FOLDS.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # Clean large species-specific objects before next species
  # ------------------------------------------------------------
  
  rm(
    species_i_obs,
    species_i_sf,
    species_range_i,
    species_range_geom,
    species_range,
    current_env_i,
    env_vals_presence,
    kde_i,
    backg,
    backg_df,
    backg_sf,
    pres_sf,
    pres_backg_sf,
    pb_extract,
    mess_reference,
    sb,
    fold_id,
    preds_i,
    current_pred_stack,
    current_pred_avg,
    species_model_object
  )
  
  
  if (
    exists(
      "current_binary"
    )
  ) {
    
    rm(
      current_binary
    )
  }
  
  
  gc()
}

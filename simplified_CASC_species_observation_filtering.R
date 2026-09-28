
install.packages("rgdal")
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
library(randomForest)
library(PresenceAbsence)
library(corrplot)
library(tidyr)
library(readr)
library(RStoolbox) # has rasterPCA
library(data.table)


# refilter species observation data

# amphibians

# species summary
#forest <- rast("F:/Uottawa_data/new_predictors_occurrence/updated2/forest_proportion_1km_with_SF.tif")
#forest
#ncell(forest)

# use climate data here
temp <- rast("F:/Uillinois_data/analysis/habitat_suitability/climate_data/historical/masked_MAT.tif")
temp
ncell(temp)

species_observations <- read.csv("F:/Uillinois_data/analysis/habitat_suitability/species_observ_new_preparation/amphibian_observations.csv")
head(species_observations)
species_observations


# Ensure CRS consistency
crs_observations <- "EPSG:4326"  # Assuming the amphibian data is in WGS84
if (crs(temp) != crs_observations) {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
  species_points <- project(species_points, crs(temp))  # Reproject to match raster
} else {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
}

species_points
temp

# Extract raster cell index (square_id)
species_observations$square_id <- cells(temp, species_points)[, "cell"]
head(species_observations)


write.csv(species_observations, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/amphibian_observations_cleaned.csv", row.names = FALSE)


# Step 1: Group and summarize
species_summary <- species_observations %>%
  group_by(genus, species,  scientific_name, group) %>%
  summarise(n_MWCASC = n(), .groups = 'drop')

# Step 2: Reorder columns
species_summary <- species_summary %>%
  dplyr::select(genus, species, group,  scientific_name, n_MWCASC)
species_summary

#View(species_summary)

# Step 3: Export to CSV
write.csv(species_summary, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/amphibian_summary.csv", row.names = FALSE)




# reptiles

# use climate data here
temp <- rast("F:/Uillinois_data/analysis/habitat_suitability/climate_data/historical/masked_MAT.tif")
temp
ncell(temp)

species_observations <- read.csv("F:/Uillinois_data/analysis/habitat_suitability/species_observ_new_preparation/reptile_observations.csv")
head(species_observations)
species_observations


# Ensure CRS consistency
crs_observations <- "EPSG:4326"  # Assuming the amphibian data is in WGS84
if (crs(temp) != crs_observations) {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
  species_points <- project(species_points, crs(temp))  # Reproject to match raster
} else {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
}

species_points
temp

# Extract raster cell index (square_id)
species_observations$square_id <- cells(temp, species_points)[, "cell"]
head(species_observations)


write.csv(species_observations, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/reptile_observations_cleaned.csv", row.names = FALSE)


# Step 1: Group and summarize
species_summary <- species_observations %>%
  group_by(genus, species,  scientific_name, group) %>%
  summarise(n_MWCASC = n(), .groups = 'drop')

# Step 2: Reorder columns
species_summary <- species_summary %>%
  dplyr::select(genus, species, group,  scientific_name, n_MWCASC)
species_summary

#View(species_summary)

# Step 3: Export to CSV
write.csv(species_summary, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/reptile_summary.csv", row.names = FALSE)



# mammals

# use climate data here
temp <- rast("F:/Uillinois_data/analysis/habitat_suitability/climate_data/historical/masked_MAT.tif")
temp
ncell(temp)

species_observations <- read.csv("F:/Uillinois_data/analysis/habitat_suitability/species_observ_new_preparation/mammal_observations.csv")
head(species_observations)
species_observations


# Ensure CRS consistency
crs_observations <- "EPSG:4326"  # Assuming the amphibian data is in WGS84
if (crs(temp) != crs_observations) {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
  species_points <- project(species_points, crs(temp))  # Reproject to match raster
} else {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
}

species_points
temp

# Extract raster cell index (square_id)
species_observations$square_id <- cells(temp, species_points)[, "cell"]
head(species_observations)


write.csv(species_observations, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/mammal_observations_cleaned.csv", row.names = FALSE)


# Step 1: Group and summarize
species_summary <- species_observations %>%
  group_by(genus, species,  scientific_name, group) %>%
  summarise(n_MWCASC = n(), .groups = 'drop')

# Step 2: Reorder columns
species_summary <- species_summary %>%
  dplyr::select(genus, species, group,  scientific_name, n_MWCASC)
species_summary

#View(species_summary)

# Step 3: Export to CSV
write.csv(species_summary, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/mammal_summary.csv", row.names = FALSE)




# birds
# use climate data here
temp <- rast("F:/Uillinois_data/analysis/habitat_suitability/climate_data/historical/masked_MAT.tif")
temp
ncell(temp)

species_observations <- read.csv("F:/Uillinois_data/GBIF/birds_breeding_April_July_combined.csv")
head(species_observations)
species_observations
#species_observations <- test

# Ensure CRS consistency
crs_observations <- "EPSG:4326"  # Assuming the amphibian data is in WGS84
if (crs(temp) != crs_observations) {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
  species_points <- project(species_points, crs(temp))  # Reproject to match raster
} else {
  species_points <- vect(species_observations, geom = c("longitude", "latitude"), crs = crs_observations)
}

species_points
temp

# Extract raster cell index (square_id)
species_observations$square_id <- cells(temp, species_points)[, "cell"]
head(species_observations)


write.csv(species_observations, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/bird_observations_cleaned.csv", row.names = FALSE)


# alternative processing in chunks to avoid memory failures

library(terra)
library(dplyr)

# -----------------------------
# External disk temp folder
# -----------------------------
terraOptions(
  tempdir = "F:/R_temp_terra",
  memfrac = 0.3,
  progress = 1
)

dir.create("F:/R_temp_terra", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Raster
# -----------------------------
#temp <- rast("F:/Uillinois_data/path_to_your_raster/temp.tif")

# -----------------------------
# Observation data
# -----------------------------
#species_observations <- read.csv(
# "F:/Uillinois_data/GBIF/birds_breeding_April_July_combined.csv"
#)

# -----------------------------
# Chunk settings
# -----------------------------
chunk_size <- 50000
n <- nrow(species_observations)

species_observations$square_id <- NA_integer_

# -----------------------------
# Process points in chunks
# -----------------------------
for (start_i in seq(1, n, by = chunk_size)) {
  
  end_i <- min(start_i + chunk_size - 1, n)
  
  cat("\nProcessing rows:", start_i, "to", end_i, "\n")
  
  chunk_df <- species_observations[start_i:end_i, ]
  
  # Convert only this chunk to SpatVector
  chunk_points <- vect(
    chunk_df,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  # Project chunk to raster CRS
  chunk_points <- project(chunk_points, crs(temp))
  
  # Extract cell ID
  chunk_cells <- cells(temp, chunk_points)
  
  species_observations$square_id[start_i:end_i] <- chunk_cells[, "cell"]
  
  rm(chunk_df, chunk_points, chunk_cells)
  gc()
}

# -----------------------------
# Save final output
# -----------------------------
write.csv(
  species_observations,
  "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/bird_observations_cleaned.csv", row.names = FALSE
)



# Step 1: Group and summarize
species_summary <- species_observations %>%
  group_by(genus, species,  scientific_name, group) %>%
  summarise(n_MWCASC = n(), .groups = 'drop')

# Step 2: Reorder columns
species_summary <- species_summary %>%
  dplyr::select(genus, species, group,  scientific_name, n_MWCASC)
species_summary

#View(species_summary)

# Step 3: Export to CSV
write.csv(species_summary, "F:/Uillinois_data/analysis/habitat_suitability/recent_cleaned_species_observation/bird_summary.csv", row.names = FALSE)


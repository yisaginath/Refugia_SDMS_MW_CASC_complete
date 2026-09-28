# produce climate only SDMs

studyarea <- st_read(
  "D:/yisa/UI/study_area/study_area_projected.shp",
  quiet = TRUE
)
studyarea

IUCN_amphibian <- st_read(
  "D:/global_mammal_range/North_america/north_america_outputs/amphibians_north_america_terrestrial_only.shp",
  quiet = TRUE
)

IUCN_amphibian


IUCN_reptile <- st_read(
  "D:/global_mammal_range/North_america/north_america_outputs/reptiles/north_america_outputs/reptiles_north_america_terrestrial_only.shp",
  quiet = TRUE
)

IUCN_reptile


IUCN_mammal <- st_read(
  "D:/global_mammal_range/North_america/north_america_outputs/mammals/north_america_outputs/mammals_north_america_terrestrial_only.shp",
  quiet = TRUE
)

IUCN_mammal



IUCN_bird <- st_read(
  "D:/global_mammal_range/North_america/north_america_outputs/birds/north_america_outputs/birds_north_america_terrestrial_terrestrial_only_breeding.shp",
  quiet = TRUE
)

IUCN_bird



# ============================================================
# CLIP IUCN SPECIES RANGES TO CASC STUDY AREA
#
# Taxonomic groups:
#   1. Amphibians
#   2. Reptiles
#   3. Mammals
#   4. Birds
#
# Outputs:
#   D:/yisa/UI/study_area/IUCN_species_CASC
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGE
# ============================================================

library(sf)


# ============================================================
# 2. DEFINE INPUT AND OUTPUT PATHS
# ============================================================

study_area_file <- paste0(
  "D:/yisa/UI/study_area/",
  "study_area_projected.shp"
)


amphibian_file <- paste0(
  "D:/global_mammal_range/North_america/",
  "north_america_outputs/",
  "amphibians_north_america_terrestrial_only.shp"
)


reptile_file <- paste0(
  "D:/global_mammal_range/North_america/",
  "north_america_outputs/reptiles/",
  "north_america_outputs/",
  "reptiles_north_america_terrestrial_only.shp"
)


mammal_file <- paste0(
  "D:/global_mammal_range/North_america/",
  "north_america_outputs/mammals/",
  "north_america_outputs/",
  "mammals_north_america_terrestrial_only.shp"
)


bird_file <- paste0(
  "D:/global_mammal_range/North_america/",
  "north_america_outputs/birds/",
  "north_america_outputs/",
  "birds_north_america_terrestrial_terrestrial_only_breeding.shp"
)


# ------------------------------------------------------------
# New output folder
# ------------------------------------------------------------

output_dir <- file.path(
  dirname(study_area_file),
  "IUCN_species_CASC"
)


if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
  
}


cat(
  "\nOutput folder:\n",
  output_dir,
  "\n",
  sep = ""
)


# ============================================================
# 3. READ STUDY AREA
# ============================================================

studyarea <- st_read(
  study_area_file,
  quiet = TRUE
)


cat(
  "\n====================================================\n",
  "STUDY AREA\n",
  "====================================================\n",
  sep = ""
)


print(studyarea)

print(
  st_crs(studyarea)
)


# ============================================================
# 4. PREPARE STUDY AREA
# ============================================================

# ------------------------------------------------------------
# Repair possible geometry problems
# ------------------------------------------------------------

studyarea <- st_make_valid(
  studyarea
)


# ------------------------------------------------------------
# Dissolve the two country polygons into one study-area
# boundary.
#
# This prevents the same species polygon from being split
# unnecessarily according to the Canada-USA boundary.
# ------------------------------------------------------------

studyarea_boundary <- st_union(
  st_geometry(studyarea)
)


studyarea_boundary <- st_sf(
  study_area = "CASC",
  geometry = studyarea_boundary
)


# ------------------------------------------------------------
# Ensure the dissolved boundary retains the original CRS
# ------------------------------------------------------------

st_crs(studyarea_boundary) <- st_crs(studyarea)


cat(
  "\n====================================================\n",
  "DISSOLVED CASC STUDY-AREA BOUNDARY\n",
  "====================================================\n",
  sep = ""
)


print(studyarea_boundary)


# ============================================================
# 5. FUNCTION TO CLIP ONE IUCN TAXON
# ============================================================

clip_iucn_range <- function(
    input_file,
    taxon_name,
    output_file,
    study_boundary
) {
  
  
  cat(
    "\n====================================================\n",
    "PROCESSING: ",
    toupper(taxon_name),
    "\n====================================================\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Read IUCN range data
  # ----------------------------------------------------------
  
  species_range <- st_read(
    input_file,
    quiet = TRUE
  )
  
  
  cat(
    "\nOriginal number of range features: ",
    nrow(species_range),
    "\n",
    sep = ""
  )
  
  
  if ("sci_name" %in% names(species_range)) {
    
    cat(
      "Original number of species: ",
      length(
        unique(
          species_range$sci_name
        )
      ),
      "\n",
      sep = ""
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Remove completely empty geometries
  # ----------------------------------------------------------
  
  species_range <- species_range[
    !st_is_empty(species_range),
  ]
  
  
  # ----------------------------------------------------------
  # Repair invalid geometries
  # ----------------------------------------------------------
  
  species_range <- st_make_valid(
    species_range
  )
  
  
  # ----------------------------------------------------------
  # Transform IUCN ranges to study-area CRS
  #
  # The CASC study area uses a projected Lambert Azimuthal
  # Equal Area CRS, whereas the original IUCN ranges use
  # geographic WGS 84.
  # ----------------------------------------------------------
  
  species_range <- st_transform(
    species_range,
    st_crs(study_boundary)
  )
  
  
  # ----------------------------------------------------------
  # Remove any geometries that became empty during repair
  # ----------------------------------------------------------
  
  species_range <- species_range[
    !st_is_empty(species_range),
  ]
  
  
  # ----------------------------------------------------------
  # First perform spatial subsetting.
  #
  # This avoids running st_intersection() on features located
  # far outside the CASC region.
  # ----------------------------------------------------------
  
  overlap_index <- lengths(
    st_intersects(
      species_range,
      study_boundary
    )
  ) > 0
  
  
  species_overlap <- species_range[
    overlap_index,
  ]
  
  
  cat(
    "Features overlapping study area before clipping: ",
    nrow(species_overlap),
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Stop cleanly if no species overlap
  # ----------------------------------------------------------
  
  if (nrow(species_overlap) == 0) {
    
    warning(
      paste(
        "No",
        taxon_name,
        "ranges overlap the study area."
      )
    )
    
    return(NULL)
    
  }
  
  
  # ----------------------------------------------------------
  # Clip species ranges exactly to study-area boundary
  # ----------------------------------------------------------
  
  clipped <- suppressWarnings(
    
    st_intersection(
      species_overlap,
      study_boundary
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Remove empty geometries
  # ----------------------------------------------------------
  
  clipped <- clipped[
    !st_is_empty(clipped),
  ]
  
  
  # ----------------------------------------------------------
  # Retain polygon geometries only
  #
  # Intersection can occasionally produce geometry
  # collections or line fragments along boundaries.
  # ----------------------------------------------------------
  
  clipped <- suppressWarnings(
    
    st_collection_extract(
      clipped,
      "POLYGON"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Repair final geometries
  # ----------------------------------------------------------
  
  clipped <- st_make_valid(
    clipped
  )
  
  
  clipped <- clipped[
    !st_is_empty(clipped),
  ]
  
  
  # ----------------------------------------------------------
  # Remove the temporary study-area attribute
  # ----------------------------------------------------------
  
  if ("study_area" %in% names(clipped)) {
    
    clipped$study_area <- NULL
    
  }
  
  
  # ----------------------------------------------------------
  # Print final information
  # ----------------------------------------------------------
  
  cat(
    "\nClipped number of range features: ",
    nrow(clipped),
    "\n",
    sep = ""
  )
  
  
  if ("sci_name" %in% names(clipped)) {
    
    cat(
      "Number of species retained: ",
      length(
        unique(
          clipped$sci_name
        )
      ),
      "\n",
      sep = ""
    )
    
  }
  
  
  cat(
    "Output CRS:\n"
  )
  
  
  print(
    st_crs(clipped)
  )
  
  
  # ----------------------------------------------------------
  # Write shapefile
  # ----------------------------------------------------------
  
  st_write(
    clipped,
    output_file,
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  
  cat(
    "\nSaved:\n",
    output_file,
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Remove large temporary objects from memory
  # ----------------------------------------------------------
  
  rm(
    species_range,
    species_overlap
  )
  
  
  gc()
  
  
  return(clipped)
  
}


# ============================================================
# 6. CLIP AMPHIBIAN IUCN RANGES
# ============================================================

amphibian_CASC <- clip_iucn_range(
  
  input_file = amphibian_file,
  
  taxon_name = "amphibians",
  
  output_file = file.path(
    output_dir,
    "amphibians_IUCN_CASC.shp"
  ),
  
  study_boundary = studyarea_boundary
  
)


# ============================================================
# 7. VERIFY AMPHIBIAN OUTPUT
# ============================================================

if (!is.null(amphibian_CASC)) {
  
  print(
    amphibian_CASC
  )
  
}


# ============================================================
# 8. CLIP REPTILE IUCN RANGES
# ============================================================

reptile_CASC <- clip_iucn_range(
  
  input_file = reptile_file,
  
  taxon_name = "reptiles",
  
  output_file = file.path(
    output_dir,
    "reptiles_IUCN_CASC.shp"
  ),
  
  study_boundary = studyarea_boundary
  
)


# ============================================================
# 9. VERIFY REPTILE OUTPUT
# ============================================================

if (!is.null(reptile_CASC)) {
  
  print(
    reptile_CASC
  )
  
}


# ============================================================
# 10. CLIP MAMMAL IUCN RANGES
# ============================================================

mammal_CASC <- clip_iucn_range(
  
  input_file = mammal_file,
  
  taxon_name = "mammals",
  
  output_file = file.path(
    output_dir,
    "mammals_IUCN_CASC.shp"
  ),
  
  study_boundary = studyarea_boundary
  
)


# ============================================================
# 11. VERIFY MAMMAL OUTPUT
# ============================================================

if (!is.null(mammal_CASC)) {
  
  print(
    mammal_CASC
  )
  
}


# ============================================================
# 12. CLIP BIRD IUCN RANGES
# ============================================================

bird_CASC <- clip_iucn_range(
  
  input_file = bird_file,
  
  taxon_name = "birds",
  
  output_file = file.path(
    output_dir,
    "birds_IUCN_CASC.shp"
  ),
  
  study_boundary = studyarea_boundary
  
)


# ============================================================
# 13. VERIFY BIRD OUTPUT
# ============================================================

if (!is.null(bird_CASC)) {
  
  print(
    bird_CASC
  )
  
}


# ============================================================
# 14. SUMMARY OF SAVED OUTPUTS
# ============================================================

cat(
  "\n\n====================================================\n",
  "ALL IUCN RANGE CLIPPING COMPLETED\n",
  "====================================================\n",
  "Output folder:\n",
  output_dir,
  "\n\n",
  "Expected shapefiles:\n",
  "  amphibians_IUCN_CASC.shp\n",
  "  reptiles_IUCN_CASC.shp\n",
  "  mammals_IUCN_CASC.shp\n",
  "  birds_IUCN_CASC.shp\n",
  "====================================================\n",
  sep = ""
)


# ------------------------------------------------------------
# List shapefiles written to output folder
# ------------------------------------------------------------

saved_shapefiles <- list.files(
  output_dir,
  pattern = "\\.shp$",
  full.names = TRUE
)


print(
  saved_shapefiles
)


# ============================================================
# 15. FINAL SPECIES COUNTS
# ============================================================

cat(
  "\n====================================================\n",
  "NUMBER OF UNIQUE SPECIES WITH RANGES IN CASC AREA\n",
  "====================================================\n",
  sep = ""
)


if (!is.null(amphibian_CASC)) {
  
  cat(
    "Amphibians: ",
    length(unique(amphibian_CASC$sci_name)),
    "\n",
    sep = ""
  )
  
}


if (!is.null(reptile_CASC)) {
  
  cat(
    "Reptiles: ",
    length(unique(reptile_CASC$sci_name)),
    "\n",
    sep = ""
  )
  
}


if (!is.null(mammal_CASC)) {
  
  cat(
    "Mammals: ",
    length(unique(mammal_CASC$sci_name)),
    "\n",
    sep = ""
  )
  
}


if (!is.null(bird_CASC)) {
  
  cat(
    "Birds: ",
    length(unique(bird_CASC$sci_name)),
    "\n",
    sep = ""
  )
  
}


cat(
  "====================================================\n"
)




# ============================================================
# PLOT CLIPPED AMPHIBIAN IUCN RANGES ACROSS CASC STUDY AREA
# ============================================================

library(sf)


# ============================================================
# 1. DEFINE FILE PATHS
# ============================================================

study_area_file <- paste0(
  "D:/yisa/UI/study_area/",
  "study_area_projected.shp"
)


amphibian_clip_file <- paste0(
  "D:/yisa/UI/study_area/",
  "IUCN_species_CASC/",
  "amphibians_IUCN_CASC.shp"
)


# ============================================================
# 2. READ STUDY AREA
# ============================================================

studyarea <- st_read(
  study_area_file,
  quiet = TRUE
)


cat(
  "\n====================================================\n",
  "STUDY AREA\n",
  "====================================================\n",
  sep = ""
)


print(
  studyarea
)


# ============================================================
# 3. READ CLIPPED AMPHIBIAN RANGES
# ============================================================

amphibian_CASC <- st_read(
  amphibian_clip_file,
  quiet = TRUE
)


cat(
  "\n====================================================\n",
  "CLIPPED AMPHIBIAN IUCN RANGES\n",
  "====================================================\n",
  sep = ""
)


print(
  amphibian_CASC
)


# ============================================================
# 4. VERIFY CRS
# ============================================================

cat(
  "\nStudy-area CRS:\n"
)

print(
  st_crs(studyarea)
)


cat(
  "\nAmphibian-range CRS:\n"
)

print(
  st_crs(amphibian_CASC)
)


cat(
  "\nDo CRS match? ",
  st_crs(studyarea) == st_crs(amphibian_CASC),
  "\n",
  sep = ""
)


# ============================================================
# 5. VERIFY NUMBER OF SPECIES
# ============================================================

cat(
  "\n====================================================\n",
  "AMPHIBIAN RANGE SUMMARY\n",
  "====================================================\n",
  sep = ""
)


cat(
  "Number of clipped polygon features: ",
  nrow(amphibian_CASC),
  "\n",
  sep = ""
)


cat(
  "Number of unique amphibian species: ",
  length(
    unique(
      amphibian_CASC$sci_name
    )
  ),
  "\n",
  sep = ""
)


# ============================================================
# 6. DISSOLVE STUDY AREA FOR OUTER BOUNDARY
# ============================================================

studyarea_boundary <- st_union(
  st_geometry(studyarea)
)


# ============================================================
# 7. PLOT CLIPPED AMPHIBIAN RANGES
# ============================================================

plot(
  st_geometry(studyarea_boundary),
  col = "grey95",
  border = "black",
  lwd = 2,
  main = paste0(
    "IUCN Amphibian Ranges within CASC Study Area\n",
    "n = ",
    length(unique(amphibian_CASC$sci_name)),
    " species"
  ),
  axes = TRUE
)


# ------------------------------------------------------------
# Add all clipped amphibian range polygons
# ------------------------------------------------------------

plot(
  st_geometry(amphibian_CASC),
  add = TRUE,
  col = adjustcolor(
    "forestgreen",
    alpha.f = 0.25
  ),
  border = adjustcolor(
    "darkgreen",
    alpha.f = 0.35
  ),
  lwd = 0.3
)


# ------------------------------------------------------------
# Redraw study-area boundary on top
# ------------------------------------------------------------

plot(
  st_geometry(studyarea_boundary),
  add = TRUE,
  border = "black",
  lwd = 2
)


# ============================================================
# 8. ADD LEGEND
# ============================================================

legend(
  "bottomleft",
  legend = c(
    "CASC study area",
    "IUCN amphibian ranges"
  ),
  fill = c(
    "grey95",
    adjustcolor(
      "forestgreen",
      alpha.f = 0.25
    )
  ),
  border = c(
    "black",
    "darkgreen"
  ),
  bty = "n",
  cex = 0.9
)



mammal_observations <- read.csv(
"D:/yisa/UI/recent_cleaned_species_observation/mammal_observations_cleaned.csv"
)

mammal_observations


amphibian_observations <- read.csv(
  "D:/yisa/UI/recent_cleaned_species_observation/amphibian_observations_cleaned.csv"
)

amphibian_observations

reptile_observations <- read.csv(
  "D:/yisa/UI/recent_cleaned_species_observation/reptile_observations_cleaned.csv"
)

reptile_observations




# ============================================================
# MATCH IUCN RANGE SCIENTIFIC NAMES TO OBSERVATION DATA
#
# Taxonomic groups:
#   1. Mammals
#   2. Amphibians
#   3. Reptiles
#
# Matching basis:
#   observation genus + species
#
# Existing clipped IUCN range folder:
#   D:/yisa/UI/study_area/IUCN_species_CASC
#
# New outputs:
#   mammals_IUCN_CASC_names_matched.shp
#   amphibians_IUCN_CASC_names_matched.shp
#   reptiles_IUCN_CASC_names_matched.shp
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

library(sf)
library(dplyr)
library(stringr)


# ============================================================
# 2. DEFINE FOLDERS
# ============================================================

range_dir <- paste0(
  "D:/yisa/UI/study_area/",
  "IUCN_species_CASC"
)


observation_dir <- paste0(
  "D:/yisa/UI/",
  "recent_cleaned_species_observation"
)


# ============================================================
# 3. DEFINE IUCN RANGE FILES
# ============================================================

mammal_range_file <- file.path(
  range_dir,
  "mammals_IUCN_CASC.shp"
)


amphibian_range_file <- file.path(
  range_dir,
  "amphibians_IUCN_CASC.shp"
)


reptile_range_file <- file.path(
  range_dir,
  "reptiles_IUCN_CASC.shp"
)


# ============================================================
# 4. DEFINE OBSERVATION FILES
# ============================================================

mammal_obs_file <- file.path(
  observation_dir,
  "mammal_observations_cleaned.csv"
)


amphibian_obs_file <- file.path(
  observation_dir,
  "amphibian_observations_cleaned.csv"
)


reptile_obs_file <- file.path(
  observation_dir,
  "reptile_observations_cleaned.csv"
)


# ============================================================
# 5. READ OBSERVATION DATA
# ============================================================

mammal_observations <- read.csv(
  mammal_obs_file,
  stringsAsFactors = FALSE
)


amphibian_observations <- read.csv(
  amphibian_obs_file,
  stringsAsFactors = FALSE
)


reptile_observations <- read.csv(
  reptile_obs_file,
  stringsAsFactors = FALSE
)


# ============================================================
# 6. READ CLIPPED IUCN RANGE DATA
# ============================================================

mammal_range <- st_read(
  mammal_range_file,
  quiet = TRUE
)


amphibian_range <- st_read(
  amphibian_range_file,
  quiet = TRUE
)


reptile_range <- st_read(
  reptile_range_file,
  quiet = TRUE
)


# ============================================================
# 7. FUNCTION TO STANDARDIZE SCIENTIFIC NAMES
# ============================================================

standardize_scientific_name <- function(x) {
  
  x <- as.character(x)
  
  # Remove leading/trailing spaces
  x <- str_trim(x)
  
  # Replace multiple spaces with one space
  x <- str_squish(x)
  
  # Convert to lowercase for matching only
  x <- str_to_lower(x)
  
  return(x)
}


# ============================================================
# 8. FUNCTION TO PREPARE OBSERVATION SPECIES NAMES
# ============================================================

prepare_observation_names <- function(obs_data) {
  
  
  # ----------------------------------------------------------
  # Check required fields
  # ----------------------------------------------------------
  
  required_columns <- c(
    "genus",
    "species"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    names(obs_data)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Observation data are missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Construct scientific name from genus + species
  # ----------------------------------------------------------
  
  obs_names <- obs_data %>%
    
    mutate(
      
      genus = str_squish(
        as.character(genus)
      ),
      
      species = str_squish(
        as.character(species)
      ),
      
      observation_sci_name = paste(
        genus,
        species
      ),
      
      match_name = standardize_scientific_name(
        observation_sci_name
      )
      
    ) %>%
    
    # Remove incomplete names
    filter(
      !is.na(genus),
      !is.na(species),
      genus != "",
      species != "",
      genus != "NA",
      species != "NA"
    ) %>%
    
    # One record per unique species name
    distinct(
      match_name,
      .keep_all = TRUE
    ) %>%
    
    select(
      match_name,
      observation_sci_name
    )
  
  
  return(
    obs_names
  )
}


# ============================================================
# 9. PREPARE UNIQUE OBSERVATION SPECIES LISTS
# ============================================================

mammal_obs_names <- prepare_observation_names(
  mammal_observations
)


amphibian_obs_names <- prepare_observation_names(
  amphibian_observations
)


reptile_obs_names <- prepare_observation_names(
  reptile_observations
)


# ============================================================
# 10. PRINT OBSERVATION SPECIES COUNTS
# ============================================================

cat(
  "\n====================================================\n",
  "UNIQUE SPECIES IN OBSERVATION DATA\n",
  "====================================================\n",
  sep = ""
)


cat(
  "Mammals: ",
  nrow(mammal_obs_names),
  "\n",
  sep = ""
)


cat(
  "Amphibians: ",
  nrow(amphibian_obs_names),
  "\n",
  sep = ""
)


cat(
  "Reptiles: ",
  nrow(reptile_obs_names),
  "\n",
  sep = ""
)


# ============================================================
# 11. FUNCTION TO MATCH AND CORRECT IUCN RANGE NAMES
# ============================================================

match_range_names <- function(
    range_data,
    observation_names,
    taxon_name
) {
  
  
  cat(
    "\n\n====================================================\n",
    "MATCHING ",
    toupper(taxon_name),
    "\n",
    "====================================================\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Ensure sci_name exists
  # ----------------------------------------------------------
  
  if (!"sci_name" %in% names(range_data)) {
    
    stop(
      paste0(
        taxon_name,
        " range data do not contain a sci_name column."
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Store original IUCN scientific name
  # ----------------------------------------------------------
  
  range_data <- range_data %>%
    
    mutate(
      
      iucn_name = as.character(
        sci_name
      ),
      
      match_name = standardize_scientific_name(
        sci_name
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Match IUCN names against observation genus + species
  # ----------------------------------------------------------
  
  range_data <- range_data %>%
    
    left_join(
      observation_names,
      by = "match_name"
    )
  
  
  # ----------------------------------------------------------
  # Identify whether each polygon's species matched
  # ----------------------------------------------------------
  
  range_data <- range_data %>%
    
    mutate(
      
      name_match = ifelse(
        !is.na(observation_sci_name),
        "matched",
        "unmatched"
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Replace sci_name with observation scientific name
  # where a match exists.
  #
  # Otherwise retain original IUCN name.
  # ----------------------------------------------------------
  
  range_data <- range_data %>%
    
    mutate(
      
      sci_name = ifelse(
        name_match == "matched",
        observation_sci_name,
        iucn_name
      )
      
    )
  
  
  # ==========================================================
  # SPECIES-LEVEL MATCH SUMMARY
  # ==========================================================
  
  species_summary <- range_data %>%
    
    st_drop_geometry() %>%
    
    select(
      iucn_name,
      observation_sci_name,
      sci_name,
      name_match
    ) %>%
    
    distinct()
  
  
  matched_species <- species_summary %>%
    
    filter(
      name_match == "matched"
    )
  
  
  unmatched_species <- species_summary %>%
    
    filter(
      name_match == "unmatched"
    )
  
  
  cat(
    "\nIUCN species in clipped range data: ",
    nrow(species_summary),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Species matched to observations: ",
    nrow(matched_species),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Species not matched: ",
    nrow(unmatched_species),
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Print matched species
  # ----------------------------------------------------------
  
  cat(
    "\n----------------------------------------------------\n",
    "MATCHED SPECIES\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  
  print(
    matched_species
  )
  
  
  # ----------------------------------------------------------
  # Print unmatched IUCN species
  # ----------------------------------------------------------
  
  cat(
    "\n----------------------------------------------------\n",
    "UNMATCHED IUCN RANGE SPECIES\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  
  print(
    unmatched_species
  )
  
  
  # ----------------------------------------------------------
  # Determine observation species absent from IUCN range data
  # ----------------------------------------------------------
  
  matched_observation_names <- unique(
    matched_species$observation_sci_name
  )
  
  
  observations_without_range <- observation_names %>%
    
    filter(
      !observation_sci_name %in%
        matched_observation_names
    )
  
  
  cat(
    "\n----------------------------------------------------\n",
    "OBSERVATION SPECIES WITHOUT AN EXACT IUCN NAME MATCH\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  
  print(
    observations_without_range
  )
  
  
  # ----------------------------------------------------------
  # Remove temporary matching field.
  #
  # Keep:
  #   sci_name = final matched/corrected name
  #   iucn_name = original IUCN name
  #   name_match = match status
  #
  # Shapefiles have a 10-character field-name limitation,
  # so long fields would otherwise be abbreviated.
  # ----------------------------------------------------------
  
  range_data$match_name <- NULL
  
  
  # Shorten fields before shapefile export
  names(range_data)[
    names(range_data) == "observation_sci_name"
  ] <- "obs_name"
  
  
  names(range_data)[
    names(range_data) == "name_match"
  ] <- "matchstat"
  
  
  names(range_data)[
    names(range_data) == "iucn_name"
  ] <- "iucn_name"
  
  
  # ----------------------------------------------------------
  # Return everything
  # ----------------------------------------------------------
  
  return(
    list(
      
      range = range_data,
      
      summary = species_summary,
      
      matched = matched_species,
      
      unmatched = unmatched_species,
      
      observations_without_range =
        observations_without_range
      
    )
  )
}


# ============================================================
# 12. MATCH MAMMAL NAMES
# ============================================================

mammal_match <- match_range_names(
  
  range_data = mammal_range,
  
  observation_names = mammal_obs_names,
  
  taxon_name = "mammals"
  
)


mammal_range_matched <- mammal_match$range


# ============================================================
# 13. MATCH AMPHIBIAN NAMES
# ============================================================

amphibian_match <- match_range_names(
  
  range_data = amphibian_range,
  
  observation_names = amphibian_obs_names,
  
  taxon_name = "amphibians"
  
)


amphibian_range_matched <- amphibian_match$range


# ============================================================
# 14. MATCH REPTILE NAMES
# ============================================================

reptile_match <- match_range_names(
  
  range_data = reptile_range,
  
  observation_names = reptile_obs_names,
  
  taxon_name = "reptiles"
  
)


reptile_range_matched <- reptile_match$range


# ============================================================
# 15. DEFINE OUTPUT SHAPEFILES
# ============================================================

mammal_output <- file.path(
  range_dir,
  "mammals_IUCN_CASC_names_matched.shp"
)


amphibian_output <- file.path(
  range_dir,
  "amphibians_IUCN_CASC_names_matched.shp"
)


reptile_output <- file.path(
  range_dir,
  "reptiles_IUCN_CASC_names_matched.shp"
)


# ============================================================
# 16. SAVE MATCHED MAMMAL RANGE
# ============================================================

st_write(
  mammal_range_matched,
  mammal_output,
  delete_layer = TRUE,
  quiet = TRUE
)


cat(
  "\nSaved mammal range:\n",
  mammal_output,
  "\n",
  sep = ""
)


# ============================================================
# 17. SAVE MATCHED AMPHIBIAN RANGE
# ============================================================

st_write(
  amphibian_range_matched,
  amphibian_output,
  delete_layer = TRUE,
  quiet = TRUE
)


cat(
  "\nSaved amphibian range:\n",
  amphibian_output,
  "\n",
  sep = ""
)


# ============================================================
# 18. SAVE MATCHED REPTILE RANGE
# ============================================================

st_write(
  reptile_range_matched,
  reptile_output,
  delete_layer = TRUE,
  quiet = TRUE
)


cat(
  "\nSaved reptile range:\n",
  reptile_output,
  "\n",
  sep = ""
)


# ============================================================
# 19. SAVE MATCHING REPORTS AS CSV
# ============================================================

write.csv(
  mammal_match$summary,
  file.path(
    range_dir,
    "mammals_IUCN_observation_name_matching.csv"
  ),
  row.names = FALSE
)


write.csv(
  amphibian_match$summary,
  file.path(
    range_dir,
    "amphibians_IUCN_observation_name_matching.csv"
  ),
  row.names = FALSE
)


write.csv(
  reptile_match$summary,
  file.path(
    range_dir,
    "reptiles_IUCN_observation_name_matching.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SAVE UNMATCHED SPECIES LISTS
# ============================================================

write.csv(
  mammal_match$unmatched,
  file.path(
    range_dir,
    "mammals_IUCN_unmatched_species.csv"
  ),
  row.names = FALSE
)


write.csv(
  amphibian_match$unmatched,
  file.path(
    range_dir,
    "amphibians_IUCN_unmatched_species.csv"
  ),
  row.names = FALSE
)


write.csv(
  reptile_match$unmatched,
  file.path(
    range_dir,
    "reptiles_IUCN_unmatched_species.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAVE OBSERVATION SPECIES WITHOUT IUCN MATCH
# ============================================================

write.csv(
  mammal_match$observations_without_range,
  file.path(
    range_dir,
    "mammal_observation_species_without_IUCN_match.csv"
  ),
  row.names = FALSE
)


write.csv(
  amphibian_match$observations_without_range,
  file.path(
    range_dir,
    "amphibian_observation_species_without_IUCN_match.csv"
  ),
  row.names = FALSE
)


write.csv(
  reptile_match$observations_without_range,
  file.path(
    range_dir,
    "reptile_observation_species_without_IUCN_match.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. READ SAVED SHAPEFILES BACK FOR VERIFICATION
# ============================================================

mammal_saved <- st_read(
  mammal_output,
  quiet = TRUE
)


amphibian_saved <- st_read(
  amphibian_output,
  quiet = TRUE
)


reptile_saved <- st_read(
  reptile_output,
  quiet = TRUE
)


# ============================================================
# 23. VERIFY FINAL SCIENTIFIC NAMES
# ============================================================

cat(
  "\n\n====================================================\n",
  "FINAL MAMMAL SCIENTIFIC NAMES\n",
  "====================================================\n",
  sep = ""
)


print(
  sort(
    unique(
      mammal_saved$sci_name
    )
  )
)


cat(
  "\n\n====================================================\n",
  "FINAL AMPHIBIAN SCIENTIFIC NAMES\n",
  "====================================================\n",
  sep = ""
)


print(
  sort(
    unique(
      amphibian_saved$sci_name
    )
  )
)


cat(
  "\n\n====================================================\n",
  "FINAL REPTILE SCIENTIFIC NAMES\n",
  "====================================================\n",
  sep = ""
)


print(
  sort(
    unique(
      reptile_saved$sci_name
    )
  )
)


# ============================================================
# 24. FINAL SUMMARY
# ============================================================

cat(
  "\n\n====================================================\n",
  "SCIENTIFIC-NAME MATCHING COMPLETED\n",
  "====================================================\n",
  "\nNew range shapefiles:\n",
  "\n",
  basename(mammal_output),
  "\n",
  basename(amphibian_output),
  "\n",
  basename(reptile_output),
  "\n",
  "\nOutput folder:\n",
  range_dir,
  "\n",
  "====================================================\n",
  sep = ""
)





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

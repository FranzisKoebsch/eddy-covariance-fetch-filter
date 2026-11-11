#################################################################################
## Eddy Covariance Fetch Filter
# Flags half hour intervals whose footprint extends beyond a user-specified area of interest (AOI)
# based on the 1D (crosswind-integrated) footprint estimates in the EddyPro FullOutput file
# Adds the fetch flag as additional column to the EddyPro FullOutput file

# Author:     Franziska Koebsch
# Date:       2026-11-10
# Repository: https://github.com/FranzisKoebsch/eddy-covariance-fetch-filter
# License:    MIT

## 1. Input
# EddyPro FullOutput file containing wind direction and the 90% footprint extent (x_90%)
# from the footprint model selected during EddyPro flux processing 
# AOI as kmz polygon file
# EC tower position (lat/lon)

## 2. Processing:
# This script creates an azimuthal fetch table: distance between EC tower position and the AOI boundaries
# for each 5° wind direction sector. Table is stored as .csv.
# Fetch filter: Compares the along-wind footprint extent (x_90%) with the fetch of the corresponding wind direction sector.
# Outcome of the fetch filters are appended as additional column 'fetch_flag' to the EddyPro FullOutput File 

# Flagging policy (following the 0-1-2 QC scheme of Mauder & Foken, 2004):
#   0 = along-wind footprint extent (x_90%) lies within the AOI
#   2 = along-wind footprint extent (x_90%) extends beyond the AOI
#   NA = missing footprint estimate

## 3. Output
# - Fetch table (.csv)
# - modified EddyPro fulloutput .csv with additional column 'fetch_flag', original file structure is pertained
# - optional: plot of tower position and AOI on downloaded ESRI Satellite base map (.jpg)

# ==============================================================================
# --- 0. Dependency Check ---
required_packages <- c("sf", "sp", "maptiles", "data.table")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Package", pkg, "not found. Installing..."))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# --- 1. USER SETTINGS ---
# ==============================================================================
# EddyPro FullOutput file
eddypro_input <- 'path/to/eddypro_fulloutput.csv'

# Area of interest (AOI) polygon as .kmz file
# Expected Coordinate System: Geographic Coordinates (lon/lat, decimal degrees), WGS84 datum (= EPSG:4326)
aoi_kmz <- 'path/to/aoi.kmz'

# Output directory
output_dir <- 'path/to/output/directory'

# Tower coordinates as Lon/Lat in decimal degrees (WGS84, EPSG:4326)
# Example coordinates
tower_lon <- 12.752202 
tower_lat <- 53.825566 

# Plot tower position and AOI on a basemap of satellite map tiles (Esri World Imagery)
plot_aoi <- T # TRUE/FALSE

# ==============================================================================
# --- 2. Check and Create Directories  ---

if (eddypro_input == "" || !exists("eddypro_input")) {
  stop("'eddypro_input' is not defined! /nPlease provide the path and file name of your eddypro fulloutput file in the User Settings")
}

if (aoi_kmz == "" || !exists("aoi_kmz")) {
  stop("'AOI File' is not defined! /nPlease provide the path and file name of your AOI polygon file in the User Settings. File must be a .kmz")
}

if(!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# --- 3. Input Data: Reading and Coordinate Transformation  ---

# Create tower point (WGS84, lon/lat) and transform to projected coordinates (ETRS89 / UTM, metres) for distance calculations 
tower_coord <- st_sfc(st_point(c(tower_lon, tower_lat)), crs = 4326)
# Determine UTM Zone and EPSG code from tower coordinates (WGS84 / UTM, Europe)
utm_zone <- floor((tower_lon + 180) / 6) + 1
epsg_utm <- 25800 + utm_zone   # ETRS89 / UTM north, e.g. zone 32 -> EPSG:25832
tower_coord_etrs <- st_transform(tower_coord, epsg_utm)
tower_coord_etrs_x <- st_coordinates(tower_coord_etrs)[1]
tower_coord_etrs_y <- st_coordinates(tower_coord_etrs)[2]
# Transform to webmercator (EPSG:3857) for plotting on web basemaps
tower_webmercator <- st_transform(st_as_sf(tower_coord), 3857)


# Unzip and read AOI polygon 
aoi_file <- unzip(aoi_kmz, exdir = tempdir()) 
aoi <- st_read(aoi_file)
# Transform to projected coordinates (ETRS89 / UTM zone, metres) for distance calculations
aoi_etrs <- st_transform(aoi, epsg_utm)
# Transform to webmercator (EPSG:3857) for plotting on web basemaps
aoi_webmercator <- st_transform(aoi_etrs, 3857)

# ==============================================================================
# --- 4. Fetch Filter Parameters  ---

bin_size <- 5 # size of wind sectors (default: 5°)
angles <- seq(0, 355, by = bin_size) # Angles from 0 (North) to 355 degrees (clockwise)
ray_length <- 1000 # A sufficiently large length in meters to ensure intersection with the AOI boundaries (default: 1000m)

# Create an sf object containing all rays 
# Each ray is a LineString from the tower to a distant point in a specific direction.
all_rays_sfc <- st_sfc(lapply(angles, function(angle_deg) {
  angle_rad <- angle_deg * pi / 180 # Convert to radians
  # Calculate end point of ray based on angle from North (meteorological convention: 0=N, 90=E)
  end_x <- tower_coord_etrs_x + ray_length * sin(angle_rad)
  end_y <- tower_coord_etrs_y + ray_length * cos(angle_rad)
  st_linestring(rbind(c(tower_coord_etrs_x, tower_coord_etrs_y), c(end_x, end_y)))
}), crs=epsg_utm)

# Transform to webmercator (EPSG:3857) for plotting on web basemaps
all_rays_sfc_webmercator <- st_transform(all_rays_sfc, 3857) 

# ==============================================================================
# --- 5. Creates Azimuthal Fetch Table  ---
# Calculate distance between EC tower position and the AOI boundaries

# Initialize a data frame to store results
distances_to_border <- data.frame(angle_deg = angles, distance_m = NA)

# Loop through each ray 
for (i in c(1:length(all_rays_sfc))) {
  current_ray <- all_rays_sfc[i]
  
  # Find geometric intersection between the current ray and the AOI polygon 
  intersection_geom <- st_intersection(current_ray, aoi_etrs)
  
  # Check if an intersection occurred and if the resulting geometry is not empty
  if (length(intersection_geom) > 0 && !st_is_empty(intersection_geom)) {
    # Cast the intersection to POINTs. This handles cases where intersection_geom might be
    # MULTIPOINT (multiple intersections) or LINESTRING (ray overlaps the AOI).
    intersection_points <- st_cast(intersection_geom, "POINT")
    
    # Filter out any empty points that might have resulted from the cast
    intersection_points <- intersection_points[!st_is_empty(intersection_points)]
    
    # Transform to webmercator (EPSG:3857) for plotting on web basemaps
    intersection_points_webmercator <- st_transform(intersection_points, 3857) 
    
    if (length(intersection_points) > 0) {
      # Calculate distance from the tower to each intersection point found
      distances_from_tower <- st_distance(tower_coord_etrs, intersection_points)
      
      # The "closest" intersection for a ray refers to the EC tower,
      # So we take the maximum distance found along this ray.
      max_dist <- max(distances_from_tower)
      
      # Store the maximum distance (in meters) for the current angle
      distances_to_border$distance_m[i] <- as.numeric(max_dist)
    }
  }
} # end of loop for (i in c(1:length(all_rays_sfc))) {

# Write Fetch table (.csv)
write.csv(distances_to_border, file.path(output_dir, "FetchTable.csv"), row.names = F)

# ==============================================================================
# --- 6. Fetch Filter  ---
# Read original Eddypro FullOutput File
# Header lines are kept as raw text to write them back unchanged.
header <- fread(eddypro_input, sep = ",", header = FALSE, nrows = 3,
                fill = TRUE, colClasses = "character")
dat <- fread(eddypro_input, skip = 3, header = FALSE, sep = ",", na.strings = "-9999")   
names(dat) <- as.character(header[1])

# Set max fetch based on wind direction increments of 5%
fetch.max <- rep(NA, length = nrow(dat))
for(i in c(1:length(fetch.max))){
  if(!is.na(dat$wind_dir[i])){
    fetch.max[i] <- distances_to_border$distance_m[which.min(abs(dat$wind_dir[i]-distances_to_border$angle_deg))]
    }
  }

# Apply fetch flag: Half-hour intervals are flagged (2) if the along-wind footprint extent (x_90%) exceeds the fetch of the corresponding wind direction sector.
dat$fetch_flag <- ifelse(dat$`x_90%` != -9999 & fetch.max > dat$`x_90%`, 0, 2)

# Add fetch flag to header 
header[, fetch_flag := c("fetch_filter", "fetch_flag", "[#]")]

# Save modified EddyPro Output file: write header, then append data block
# Create new file name: Original EddyPro FullOutput file + suffix '_fetchflag' + '.csv'
file_name <- strsplit(basename(eddypro_input), ".csv")[[1]]
file_name <- paste0(file_name, "_fetchflag.csv")
fwrite(header, file.path(output_dir, file_name), col.names = FALSE, sep = ",", quote = FALSE, na = "")
fwrite(dat, file.path(output_dir, file_name), append = TRUE, col.names = FALSE,
       sep = ",", quote = FALSE, na = "-9999")

# ==============================================================================
# --- 7. Plot Tower Position and AOI  ---
# Optional, specified by plot_aoi == T/F

if(plot_aoi == T) {

# Bounding box for tile download in Web Mercator ---
bbox_webmercator <- st_bbox(aoi_webmercator)
# add 10% padding
padding_factor <- 0.10 # 10% padding
x_range_wm <- diff(bbox_webmercator[c("xmin", "xmax")])
y_range_wm <- diff(bbox_webmercator[c("ymin", "ymax")])
bbox_webmercator_padded <- c(
  xmin = as.numeric(bbox_webmercator["xmin"]) - padding_factor * x_range_wm,
  ymin = as.numeric(bbox_webmercator["ymin"]) - padding_factor * y_range_wm,
  xmax = as.numeric(bbox_webmercator["xmax"]) + padding_factor * x_range_wm,
  ymax = as.numeric(bbox_webmercator["ymax"]) + padding_factor * y_range_wm
)
names(bbox_webmercator_padded) <- c("xmin", "ymin", "xmax", "ymax")
bbox_webmercator_padded <- st_bbox(bbox_webmercator_padded, crs = 3857)

# Download satellite map tiles (Esri World Imagery) as basemap for plotting;
# tiles are cached in 'maptiles_cache' (relative to working directory) for re-use
tiles <- get_tiles(bbox_webmercator_padded,
                   provider = "Esri.WorldImagery",
                   zoom = 17, 
                   crop = FALSE,
                   cachedir =  file.path(output_dir, "maptiles_cache")) # Optional: cache tiles for faster re-use

# Plot EC Tower position and AOI on basemap 
# Save map as JPEG
jpeg(file.path(output_dir, "aoi_map.jpg"),
     width = 20, height = 20, units = "cm", res = 300, quality = 95)

plot_tiles(tiles)
plot(tower_webmercator,
     pch = 19, # Solid circle (highly visible)
     col = "red", # Bright red for visibility
     cex = 1,    # Size
     add = TRUE)

plot(aoi_webmercator, add=T, col=rgb(1,1,1, alpha=0), border = "white", lwd=2) 
dev.off()
} # end of if(plot_aoi == T) 


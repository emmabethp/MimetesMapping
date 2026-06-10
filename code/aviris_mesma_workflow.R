# AVIRIS-NG imaging spectroscopy workflow:
# 1) read AVIRIS-NG data (.tif)
# 2) crop to area of interest (.geojson)
# 3) calculate vegetation indices from band combinations
# 4) extract index values to points
# 5) run MESMA using luna

library(terra)
library(luna)

# ---- User inputs -------------------------------------------------------------
aviris_tif <- "/absolute/path/to/aviris_ng_cube.tif"
aoi_geojson <- "/absolute/path/to/aoi.geojson"
points_file <- "/absolute/path/to/sample_points.geojson" # points for extraction

# Example AVIRIS-NG band choices (update to your wavelength-to-band mapping):
band_nir <- 50L
band_red <- 30L
band_green <- 20L
band_swir <- 100L

# Endmember names and spectra should be updated for your application.
# Here we extract spectra from sample points and use them as a simple example.

# ---- Read and crop raster ----------------------------------------------------
aviris <- rast(aviris_tif)
aoi <- vect(aoi_geojson)

# Crop first for speed, then mask to AOI boundary.
aviris_crop <- crop(aviris, aoi)
aviris_mask <- mask(aviris_crop, aoi)

# ---- Vegetation indices ------------------------------------------------------
nir <- aviris_mask[[band_nir]]
red <- aviris_mask[[band_red]]
green <- aviris_mask[[band_green]]
swir <- aviris_mask[[band_swir]]
names(nir) <- "nir"
names(red) <- "red"
names(green) <- "green"
names(swir) <- "swir"

# NDVI = (NIR - RED) / (NIR + RED)
ndvi <- (nir - red) / (nir + red)
names(ndvi) <- "NDVI"

# EVI = 2.5 * (NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1)
# If you have a blue band, substitute it below. We use green as a placeholder.
evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * green + 1)
names(evi) <- "EVI"

# NDWI (Gao-style using NIR and SWIR) = (NIR - SWIR) / (NIR + SWIR)
ndwi <- (nir - swir) / (nir + swir)
names(ndwi) <- "NDWI"

veg_indices <- c(ndvi, evi, ndwi)

# ---- Extract to points -------------------------------------------------------
pts <- vect(points_file)

# Extract vegetation index values and selected AVIRIS bands to points.
index_at_points <- extract(veg_indices, pts)
bands_at_points <- extract(c(nir, red, green, swir), pts)

# Combine outputs (ID column comes from terra::extract)
point_data <- merge(index_at_points, bands_at_points, by = "ID", suffixes = c("_idx", "_bands"))

# ---- MESMA analysis with luna -----------------------------------------------
# Build a simple library matrix from extracted band values
# (rows = samples, columns = bands). In practice, replace with known endmembers.
spectral_library <- as.matrix(na.omit(point_data[, c("nir", "red", "green", "swir")]))

# MESMA requires image spectra and candidate endmember spectra.
# Convert the cropped raster to a matrix where rows are pixels and columns are bands.
img_values <- values(aviris_mask, mat = TRUE)
img_values <- img_values[complete.cases(img_values), , drop = FALSE]

# Example MESMA call. Adjust arguments to your luna version/data:
# - image: matrix of pixel spectra
# - emlib: matrix of endmember spectra
# - n: number of endmembers in each model
mesma_result <- mesma(
  image = img_values,
  emlib = spectral_library,
  n = 2
)

# ---- Optional outputs --------------------------------------------------------
# Write vegetation index stack to disk.
writeRaster(
  veg_indices,
  filename = "/absolute/path/to/output/vegetation_indices.tif",
  overwrite = TRUE
)

# Save point-level extracted values and MESMA outputs.
write.csv(point_data, "/absolute/path/to/output/point_extractions.csv", row.names = FALSE)
saveRDS(mesma_result, "/absolute/path/to/output/mesma_result.rds")

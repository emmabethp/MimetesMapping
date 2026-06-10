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
out_indices_tif <- "/absolute/path/to/output/vegetation_indices.tif"
out_points_csv <- "/absolute/path/to/output/point_extractions.csv"
out_mesma_rds <- "/absolute/path/to/output/mesma_result.rds"

# Example AVIRIS-NG band choices. Update these to your wavelength-to-band mapping
# from the AVIRIS-NG metadata/header for your specific scene
# (e.g., NIR ~860 nm, RED ~660 nm, GREEN ~560 nm, SWIR ~1600 nm).
band_nir <- 50L
band_red <- 30L
band_green <- 20L
band_swir <- 100L
band_blue <- NA_integer_ # set to a valid band index to calculate standard EVI

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

# NDWI (Gao-style using NIR and SWIR) = (NIR - SWIR) / (NIR + SWIR)
ndwi <- (nir - swir) / (nir + swir)
names(ndwi) <- "NDWI"

# EVI can only be computed with a valid BLUE band.
veg_indices <- c(ndvi, ndwi)
if (!is.na(band_blue)) {
  blue <- aviris_mask[[band_blue]]
  names(blue) <- "blue"
  evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  names(evi) <- "EVI"
  veg_indices <- c(veg_indices, evi)
}

# ---- Extract to points -------------------------------------------------------
pts <- vect(points_file)

# Extract vegetation index values and selected AVIRIS bands to points.
index_at_points <- extract(veg_indices, pts)
bands_at_points <- extract(c(nir, red, green, swir), pts)
spectra_at_points <- extract(aviris_mask, pts)

# Combine outputs (ID column comes from terra::extract)
point_data <- merge(index_at_points, bands_at_points, by = "ID")

# ---- MESMA analysis with luna -----------------------------------------------
# Build a spectral library matrix from extracted full spectra at points
# (rows = samples, columns = AVIRIS bands). Replace with known endmembers as needed.
spectral_library <- as.matrix(
  na.omit(spectra_at_points[, !names(spectra_at_points) %in% "ID", drop = FALSE])
)
if (nrow(spectral_library) < 1) {
  stop("No valid spectra extracted for spectral library; provide endmember samples.")
}

# MESMA requires image spectra and candidate endmember spectra.
# Convert the cropped raster to a matrix where rows are pixels and columns are bands.
# For very large scenes this can be memory-intensive; consider tiled/chunked workflows.
# Threshold is total raster values (ncell x nlyr), not bytes directly.
# Approximate memory for numeric storage is values x 8 bytes:
# 5e7 values x 8 bytes = ~400 MB (decimal); tune to your available RAM.
max_values_in_memory <- 5e7
if ((ncell(aviris_mask) * nlyr(aviris_mask)) > max_values_in_memory) {
  stop("Scene is too large for full in-memory MESMA example; use a chunked workflow.")
}
img_values <- values(aviris_mask, mat = TRUE)
img_values <- img_values[complete.cases(img_values), , drop = FALSE]

# Example MESMA call. Adjust arguments to your luna version/data:
# - image: matrix of pixel spectra
# - emlib: matrix of endmember spectra
# - n: number of endmembers in each model (set based on expected material mixing)
#   (2 is a simple default; increase when pixels are expected to mix more materials)
mesma_result <- mesma(
  image = img_values,
  emlib = spectral_library,
  n = 2
)
# mesma_result typically contains model fit information (e.g., abundances/error),
# but structure can vary by luna version; inspect with str(mesma_result).

# ---- Optional outputs --------------------------------------------------------
# Write vegetation index stack to disk.
writeRaster(
  veg_indices,
  filename = out_indices_tif,
  overwrite = TRUE
)

# Save point-level extracted values and MESMA outputs.
write.csv(point_data, out_points_csv, row.names = FALSE)
saveRDS(mesma_result, out_mesma_rds)

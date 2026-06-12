# AVIRIS-NG imaging spectroscopy workflow:
# 1) read AVIRIS-NG data (.tif)
# 2) crop to area of interest (.geojson)
# 3) calculate vegetation indices from band combinations
# 4) extract index values to points
# 5) run MESMA using luna

# install.packages('luna', repos='https://rspatial.r-universe.dev')

library(tidyverse)
library(terra)
library(luna)

# ---- User inputs -------------------------------------------------------------
aviris_tif <- "/Users/jasper/Documents/Datasets/BioSCape/ANG/CapePeninsula/SilvermineEast/ang20231126t085417_012_L2A_OE_0b4f48b4_RFL_ORT_QL/ang20231126t085417_012_L2A_OE_0b4f48b4_RFL_ORT.nc"
aviris_rgb_tif <- "/Users/jasper/Documents/Datasets/BioSCape/ANG/CapePeninsula/SilvermineEast/ang20231126t085417_012_L2A_OE_0b4f48b4_RFL_ORT_QL/ang20231126t085417_012_L2A_OE_0b4f48b4_RFL_ORT_QL.tif"
aoi_geojson <- "data/Silvermine_Hons_2026.geojson"
points_file <- "data/SilvermineEM.kml" # points for extraction
out_indices_tif <- "/absolute/path/to/output/vegetation_indices.tif"
out_points_csv <- "/absolute/path/to/output/point_extractions.csv"
out_mesma_rds <- "/absolute/path/to/output/mesma_result.rds"

# Example AVIRIS-NG band choices. Update these to your wavelength-to-band mapping
# from the AVIRIS-NG metadata/header for your specific scene
# (e.g., NIR ~860 nm, RED ~660 nm, GREEN ~560 nm, SWIR ~1600 nm).
band_nir <- 98L
band_red <- 59L
band_green <- 36L
band_swir <- 174L
band_blue <- 15L # set to a valid band index to calculate standard EVI

# Endmember names and spectra should be updated for your application.
# Here we extract spectra from sample points and use them as a simple example.

# ---- Read and crop raster ----------------------------------------------------
aviris <- rast(aviris_tif)
aoi <- project(vect(aoi_geojson), crs(aviris))

# Crop first for speed, then mask to AOI boundary.
aviris_mask <- mask(crop(aviris, aoi), aoi, filename = "data/aviris_cropped.tif", overwrite = TRUE)

# ---- Vegetation indices ------------------------------------------------------

# Extract RGB bands to plot a true colour image
red <- aviris_mask[[band_red]]
green <- aviris_mask[[band_green]]
blue <- aviris_mask[[band_blue]]
names(red) <- "red"
names(green) <- "green"
names(blue) <- "blue"

plotRGB(c(red, green, blue), stretch = "lin", main = "AVIRIS True Color Composite")

# Extract NIR and SWIR bands for index calculations
nir <- aviris_mask[[band_nir]]
swir <- aviris_mask[[band_swir]]

names(nir) <- "nir"
names(swir) <- "swir"

# NDVI = (NIR - RED) / (NIR + RED)
ndvi <- (nir - red) / (nir + red)
names(ndvi) <- "NDVI"

# NDWI (Gao-style using NIR and SWIR) = (NIR - SWIR) / (NIR + SWIR)
ndwi <- (nir - swir) / (nir + swir)
names(ndwi) <- "NDWI"

# EVI can only be computed with a valid BLUE band.
evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
names(evi) <- "EVI"

veg_indices <- c(ndvi, evi, ndwi)
plot(veg_indices, nc = 3)

plot(evi, main = "EVI from AVIRIS")

# ---- Extract to points -------------------------------------------------------
pts <- project(vect(points_file), crs(aviris))

# Extract vegetation index values and selected AVIRIS bands to points.
index_at_points <- terra::extract(veg_indices, pts) |>
  mutate(EndMember = pts$Name) |>
  na.omit()
bands_at_points <- terra::extract(c(nir, red, green, swir), pts) |>
  mutate(EndMember = pts$Name) |>
  na.omit()
spectra_at_points <- terra::extract(aviris_mask, pts) |>
  mutate(EndMember = pts$Name) |>
  na.omit()

# Plot spectral library
spectra_at_points |>
  group_by(EndMember, ID) |>
  pivot_longer(cols = starts_with("reflectance"), names_to = "Band", values_to = "Reflectance") |>
  mutate(Band = parse_number(str_split_i(Band, "=", 2))) |>
  ggplot(aes(x = Band, y = Reflectance, group = ID, color = EndMember)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = EndMember), fun = mean, geom = "line", linewidth = 1) +
  annotate(geom = "rect", xmin = 1350, xmax =  1420 , ymin = -Inf, ymax = Inf, 
            fill = "grey", alpha = 0.3) + #1263-1562 nm 
  annotate(geom = "rect", xmin = 1785, xmax =  1975 , ymin = -Inf, ymax = Inf, 
            fill = "grey", alpha = 0.3) + #1761-1958 nm
  theme_minimal() +
  ylim(0,0.5) +
  facet_wrap(~EndMember) + # facet_wrap(~Forest) +
  labs(title = "Spectral library of Silvermine Endmembers") +
  xlab("Wavelength (nm)")

# Calculate the mean spectrum for each endmember
spectra_at_points |>
  group_by(EndMember) |>
  pivot_longer(cols = starts_with("reflectance"), names_to = "Band", values_to = "Reflectance") |>
  group_by(EndMember, Band) |>
  summarise(mean_reflectance = mean(Reflectance, na.rm = TRUE)) |>
  pivot_wider(names_from = Band, values_from = mean_reflectance)

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



matplot(t(spectral_library), type = "l", lty = 1, col = rainbow(nrow(spectral_library)),
        xlab = "Band Index", ylab = "Reflectance", main = "Spectral Library from Sample Points")

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

# AVIRIS-NG imaging spectroscopy workflow:
# 1) read AVIRIS-NG data (.tif)
# 2) crop to area of interest (.geojson)
# 3) calculate vegetation indices from band combinations
# 4) extract index values to points
# 5) run MESMA using luna

# install.packages('luna', repos='https://rspatial.r-universe.dev')
library(tidyverse) # For general data handling
library(broom) # For tidying model outputs
library(terra) # For handling raster and vector data
library(luna) # A new companion package for terra, but specifically designed for remote sensing analyses
library(mapview) # A useful package for interactive spatial data visualization
library(spectrolab) # For handling and analyzing leaf-level spectral data
library(patchwork) # For combining multiple ggplots into one figure


# ---- User inputs -------------------------------------------------------------
setwd("~/Desktop/GIT/MimetesMapping")
aviris_nc <- "~/Desktop/GIT/BigData/ang20231126t085417_012_L2A_OE_0b4f48b4_RFL_ORT.nc"
aoi_geojson <- "~/Desktop/GIT/MimetesMapping/data/Silvermine_Hons_2026.geojson"
points_file <- "~/Desktop/GIT/MimetesMapping/data/SilvermineEM.kml" # points for extraction
mimetes_points <- "~/Desktop/GIT/MimetesMapping/data/mimetes_canopy.geojson"
leucospermum_points <- "~/Desktop/GIT/MimetesMapping/data/leucospermum_canopy.geojson"

# Example AVIRIS-NG band choices. Update these to your wavelength-to-band mapping
# from the AVIRIS-NG metadata/header for your specific scene
# (e.g., NIR ~860 nm, RED ~660 nm, GREEN ~560 nm, SWIR ~1600 nm).
band_nir <- 98L # see names(aviris)[98]
band_red <- 59L
band_green <- 36L
band_swir <- 174L
band_blue <- 15L 

# Endmember names and spectra should be updated for your application.
# Here we extract spectra from sample points and use them as a simple example.

# ---- Read and crop raster ----------------------------------------------------
aviris_big <- rast(aviris_nc)
aoi <- project(vect(aoi_geojson), crs(aviris_big))

# Crop first for speed, then mask to AOI boundary.
aviris <- mask( # Mask raster by area of interest (AOI) boundary
  crop(aviris_big, aoi) # Crop to AOI before masking (more efficient apparently)
  [[grep("reflectance",names(aviris_big))]], # keep only reflectance bands for analysis (remove water vapour and aerosol bands); adjust if needed
  aoi, # AOI - read in as a separate vector file
  filename = "data/aviris_cropped.tif", # save masked raster to disk so R doesn't have to keep it in memory - more efficient
  overwrite = TRUE) # allow overwriting existing file if one already exists

# ---- Vegetation indices ------------------------------------------------------

# Extract RGB bands to plot a true colour image
red <- aviris[[band_red]]
green <- aviris[[band_green]]
blue <- aviris[[band_blue]]
names(red) <- "red"
names(green) <- "green"
names(blue) <- "blue"

#plot 3 bands
plotRGB(c(red, green, blue), stretch = "lin", main = "AVIRIS True Color Composite")
viewRGB(raster::stack(c(red, green, blue)), r = 1, g = 2, b = 3) #interactive

# Extract NIR and SWIR bands and calculate indices using raster maths
nir <- aviris[[band_nir]]
swir <- aviris[[band_swir]]
names(nir) <- "nir"
names(swir) <- "swir"

# NDVI = (NIR - RED) / (NIR + RED)
ndvi <- (nir - red) / (nir + red)
names(ndvi) <- "NDVI"

# NDWI (Gao-style using NIR and SWIR) = (NIR - SWIR) / (NIR + SWIR)
ndwi <- (nir - swir) / (nir + swir)
names(ndwi) <- "NDWI"

# EVI requires a BLUE band, in addition to red and nir.
evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
names(evi) <- "EVI"

veg_indices <- c(ndvi, evi, ndwi)
plot(veg_indices, nc = 3) #, range = c(-1,1)) # Plot the indices in a 3-panel plot. Could add a common range of -1 to 1 for all indices, because cover types like the sea and city areas have very high values that would make it hard to see the variation in the vegetation areas if we plotted them with their own ranges.





plot(evi, main = "EVI from AVIRIS")

#evi <- clamp(evi, lower = 0, upper = 0.7) # clamp values to 0-0.7 for better visualization if needed

mapview(evi, alpha.regions = 0.4, map.types = "Esri.WorldImagery")
hist(evi, main = "Histogram of EVI values", xlab = "EVI", breaks = 50)

mapview(app(evi>0.3, as.integer), alpha.regions = 0.4, map.types = "Esri.WorldImagery")


# ---- Extract to points -------------------------------------------------------
pts <- project(vect(points_file), crs(aviris))

mapview(pts, label = pts$Name, cex = 5, map.types = "Esri.WorldImagery")

# Extract vegetation index values and selected AVIRIS bands to points.
terra::extract(veg_indices, pts) |>
  mutate(EndMember = pts$Name) |>
  na.omit() |>
  pivot_longer(cols = c("NDVI", "EVI", "NDWI"), names_to = "Index", values_to = "Value") |>
  ggplot(aes(x = EndMember, y = Value, fill = EndMember)) +
  geom_boxplot() +
  facet_wrap(~ Index, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Define a matrix of ranges and numericclass: [from, to, becomes_integer_ID]
reclass_matrix <- matrix(c(0.7,  1,  1, # 1 = shrub
                           0.4, 0.7,  2, # 2 = low vegetation
                           0, 0.4, 3, # 3 = rock/path
                           -1, 0, 4), # 4 = water/sea
                         ncol=3, byrow=TRUE)

# Apply the classification
# include.lowest=TRUE ensures the absolute minimum value is counted
rc <- classify(ndvi, reclass_matrix, include.lowest=TRUE)

plot(rc, main = "Classified by NDVI")
mapview(rc, col.regions = rainbow(4), alpha.regions = 0.4, map.types = "Esri.WorldImagery")

# Extract reflectance spectra at point locations to make a "spectral library" of the different cover classes
spectra_at_points <- terra::extract(aviris, pts) |>
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
           fill = "black", alpha = 0.7) + #1263-1562 nm 
  annotate(geom = "rect", xmin = 1785, xmax =  1975 , ymin = -Inf, ymax = Inf, 
           fill = "black", alpha = 0.7) + #1761-1958 nm
  theme_minimal() +
  ylim(0,0.5) +
  facet_wrap(~EndMember) + # facet_wrap(~Forest) +
  labs(title = "Spectral library of Silvermine Endmembers") +
  xlab("Wavelength (nm)")



#0rrr

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
spectra_at_points |> filter(EndMember == "shrub") |>
  select(starts_with("reflectance")) |> 
  t() |>
  matplot(type = "l", col = rainbow(20), lty = 1,
        xlab = "Band", ylab = "Reflectance", main = "Spectral library of Silvermine Endmembers")

spectra_at_points |>
  group_by(EndMember, ID) |>
  pivot_longer(cols = starts_with("reflectance"), names_to = "Band", values_to = "Reflectance") |>
  mutate(Band = parse_number(str_split_i(Band, "=", 2))) |>
  ggplot(aes(x = Band, y = Reflectance, group = ID, color = ID)) +
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
  pivot_wider(names_from = Band, values_from = mean_reflectance) |>
  column_to_rownames(var = "EndMember") ->
  spectral_library
  
spectral_library |>
  t() |>
  matplot(type = "l", col = rainbow(20), lty = 1,
          xlab = "Band", ylab = "Reflectance", main = "Spectral library of Silvermine Endmembers")

spectra_at_points[c(1,6,11,16),] |>
  remove_rownames() |>
  column_to_rownames(var = "EndMember") |>
  select(starts_with("reflectance")) -> spectral_library

# # Combine outputs (ID column comes from terra::extract)
# point_data <- merge(index_at_points, bands_at_points, by = "ID")

# ---- MESMA analysis with luna -----------------------------------------------
# Build a spectral library matrix from extracted full spectra at points
# (rows = samples, columns = AVIRIS bands). Replace with known endmembers as needed.
# spectral_library <- as.matrix(
#   na.omit(spectra_at_points[, !names(spectra_at_points) %in% "ID", drop = FALSE])
# )
# if (nrow(spectral_library) < 1) {
#   stop("No valid spectra extracted for spectral library; provide endmember samples.")
# }

# Select the first spectrum from each class for our spectral library
spectra_at_points[c(1,6,11,16),] |>
  remove_rownames() |>
  column_to_rownames(var = "EndMember") |>
  select(starts_with("reflectance")) -> spectral_library

# Run MESMA
mesma_result <- mesma(aviris, em = spectral_library[,], iterate = 400)
plot(mesma_result)
mapview::mapview(mesma_result$shrub, alpha.regions = 0.4, map.types = "Esri.WorldImagery")

shrub <- app(mesma_result$shrub>0.7, as.integer)

mapview(shrub, alpha.regions = 0.4, map.types = "Esri.WorldImagery")

# Stack rasters and plot pairwise relationships
pairs(c(mesma_result$shrub, evi, ndvi))

# Get points for each species
mpts <- vect(mimetes_points)
lpts <- vect(leucospermum_points)
cpts <- rbind(mpts, lpts) # combine into one vector object

# Extract spectra for each species
spectra_species <- terra::extract(aviris_big, cpts) |> # Note that I extract from the full AVIRIS-NG data here, rather than the cropped version, because some points are outside of the cropped area. This is a good example of why it is important to check your data and make sure that your points are within the area of interest before cropping or masking your raster data.
  mutate(Species = cpts$Name) |>
  na.omit()

# Plot spectral signatures for each species
spectra_species |>
  group_by(Species, ID) |>
  pivot_longer(cols = starts_with("reflectance"), names_to = "Band", values_to = "Reflectance") |>
  mutate(Band = parse_number(str_split_i(Band, "=", 2))) |>
  ggplot(aes(x = Band, y = Reflectance, group = ID, color = Species)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = Species), fun = mean, geom = "line", linewidth = 1) +
  annotate(geom = "rect", xmin = 1350, xmax =  1420 , ymin = -Inf, ymax = Inf, 
           fill = "black", alpha = 0.7) + #1263-1562 nm 
  annotate(geom = "rect", xmin = 1785, xmax =  1975 , ymin = -Inf, ymax = Inf, 
           fill = "black", alpha = 0.7) + #1761-1958 nm
  theme_minimal() +
  # ylim(0,0.5) +
  # facet_wrap(~Species) + # facet_wrap(~Forest) +
  labs(title = "Mimetes vs Leucospermum spectra") +
  xlab("Wavelength (nm)")

# Do PCA
hyp.pca <- prcomp(spectra_species |> select(starts_with("reflectance")), center = TRUE, scale. = TRUE)
#screeplot(hyp.pca)

#autoplot(hyp.pca, data = mlib, colour = "Forest")

# Augment data with PCA results and original groupings
pca_data <- augment(hyp.pca, spectra_species$Species)

# Plot with ellipses
pca_data |> ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = data)) +
  geom_point() +
  stat_ellipse() + # Adds 95% confidence ellipses
  theme_minimal() +
  labs(title = "PCA grouping by species",
       x = "PC1",
       y = "PC2")

# Select the first spectrum from each class for our spectral library
spectra_species[c(1,6),] |>
  remove_rownames() |>
  column_to_rownames(var = "Species") |>
  select(starts_with("reflectance")) -> spp_spectral_library

# Mask AVIRIS imagery to only include pixels that are classified as "shrub" from the previous MESMA run. This will allow us to focus on the areas where we know there are shrubs and see if we can differentiate between the two species.
aviris_shrub <- mask(aviris, shrub, maskvalue = 0)


# Run MESMA
shrub_mesma_result <- mesma(aviris_shrub, em = spp_spectral_library[,], iterate = 400)
plot(shrub_mesma_result)
hist(shrub_mesma_result[[c(1,2)]])
plot(shrub_mesma_result$Mimetes > shrub_mesma_result$Leucospermum, legend = FALSE, main = "Mimetes (yellow) vs Leucospermum (purple)")
plot(shrub_mesma_result$Mimetes < shrub_mesma_result$Leucospermum, legend = FALSE, main = "Mimetes (purple) vs Leucospermum (yellow)")

#Orrrrr

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
#img_values <- values(aviris_mask, mat = TRUE)
#img_values <- img_values[complete.cases(img_values), , drop = FALSE]

# Example MESMA call. Adjust arguments to your luna version/data:
# - image: matrix of pixel spectra
# - emlib: matrix of endmember spectra
# - n: number of endmembers in each model (set based on expected material mixing)
#   (2 is a simple default; increase when pixels are expected to mix more materials)
mesma_result <- mesma(aviris_mask, em = spectral_library[,], iterate = 400)

mesma_result

plot(mesma_result)

mapview::mapview(mesma_result$shrub, alpha.regions = 0.35, map.types = "Esri.WorldImagery")

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



#Leaf spectroscopy

## Get data
sdat = read_spectra(path = "data/spectra", extract_metadata = TRUE)
sdat

# Normalise spectra
sdat <- normalize(sdat)

# Extract metadata from filename
mdat <- separate_wider_delim(as.data.frame(names(sdat)), 
                             cols = "names(sdat)", 
                             delim = "_", 
                             names = c("Species", "Replicate")) 

## PCA

# Prep data
lib <- data.frame(as.matrix(sdat$value)) |>
  bind_cols(mdat) # Note that binds by row number, so the order of the metadata must match the order of the spectra in the spectrolab object

# Do PCA
hyp.pca <- prcomp(lib |> select(starts_with("X")), center = TRUE, scale. = TRUE)

# Augment data with PCA results and original groupings
pca_data <- augment(hyp.pca, lib)

# Plot PCA axes 1 and 2 with ellipses
pca_data |>
  #select(Species, starts_with(".fitted")) |>
  #pivot_longer(cols = starts_with(".fitted"), names_to = "PC", values_to = "Value") |>
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = Species)) +
  geom_point() +
  stat_ellipse() + # Adds 95% confidence ellipses
  theme_minimal() +
  labs(title = "PCA of leaf spectra by species",
       x = "PC1",
       y = "PC2")

#It looks like there is some degree of differentiation between the three species along the first two principal components. This suggests that it should be possible to differentiate between the two focal species using hyperspectral data, especially when the variation in the higher PC axes is also considered, but further analyses and validation would be needed to confirm this.

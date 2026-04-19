# =============================================================================
# Regenerar figura 4b (incertidumbre) con nueva paleta magma
# Carga el modelo ya ajustado (results_ipp_2m.rds) y regenera solo la figura.
# =============================================================================
library(INLA); library(inlabru); library(sf); library(terra)
library(ggplot2); library(viridis); library(scales); library(ggspatial)
library(ggnewscale)

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("--file=", args)]
if (length(file_arg) > 0) {
  script_path <- sub("--file=", "", file_arg[1])
  BASE <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  BASE <- normalizePath("..", mustWork = FALSE)
}

DATA <- file.path(BASE, "DATA")
FIG  <- file.path(BASE, "FIG")
cat("BASE:", BASE, "\n")

# ── 1. Cargar modelo ─────────────────────────────────────────────────────────
cat("Cargando modelo...\n")
fit_2m <- readRDS(file.path(BASE, "results_ipp_2m.rds"))

# ── 2. Cargar covariables y reconstruir rasters procesados ──────────────────
cat("Cargando covariables...\n")
pts_wgs <- st_read(file.path(DATA, "MenM_VdeA_dem.gpkg"), quiet = TRUE)
pts     <- st_transform(pts_wgs, crs = 9377)
pts_coords <- st_coordinates(pts)

slope_2m  <- rast(file.path(DATA, "cov2m_pendiente.tif"))
aspect_2m <- rast(file.path(DATA, "cov2m_aspecto.tif"))
cov_cob   <- rast(file.path(DATA, "cobertura_2024.tif"))
cov_geo   <- rast(file.path(DATA, "zoni_geotecnica.tif"))

domain_buf <- st_convex_hull(st_union(pts)) |> st_buffer(500) |> st_as_sf()
st_crs(domain_buf) <- 9377
dem_2m <- rast(file.path(DATA, "dem_2m.tif"))
dem_mask <- as.polygons(ifel(!is.na(dem_2m), 1, NA)) |>
  st_as_sf() |> st_union() |> st_as_sf()
st_crs(dem_mask) <- 9377
domain_final <- st_intersection(domain_buf, dem_mask)
domain_final <- st_sf(geometry = st_union(st_geometry(domain_final)), crs = 9377)

# Reescalar covariables
sc <- function(r) { m <- global(r,"mean",na.rm=TRUE)[[1]]; s <- global(r,"sd",na.rm=TRUE)[[1]]; (r-m)/s }
slope_sc_r  <- sc(slope_2m);  names(slope_sc_r)  <- "slope_sc"
aspect_sc_r <- sc(aspect_2m); names(aspect_sc_r) <- "aspect_sc"
cov_sc <- c(slope_sc_r, aspect_sc_r)

cov_cob_r <- project(cov_cob, slope_2m, method="near")
cov_geo_r <- project(cov_geo, slope_2m, method="near")
cob_antropico <- ifel(cov_cob_r == 1, 1L, 0L); names(cob_antropico) <- "cob_antropico"
geo_finos_bin <- ifel(cov_geo_r %in% c(4,5,17,11,12,13,14), 1L, 0L); names(geo_finos_bin) <- "geo_finos_bin"
cov_cat <- c(cob_antropico, geo_finos_bin)

# Reconstruir malla para obtener extent
bnd  <- fmesher::fm_as_segm(domain_final)
mesh <- fmesher::fm_mesh_2d(
  boundary = bnd, max.edge = c(100, 500), offset = c(200, 2000),
  cutoff = 20, crs = sf::st_crs(9377))
cat(sprintf("  Malla: %d nodos\n", mesh$n))

mesh_ext <- terra::ext(
  min(mesh$loc[,1]) - 200, max(mesh$loc[,1]) + 200,
  min(mesh$loc[,2]) - 200, max(mesh$loc[,2]) + 200)

fill_na_focal <- function(r) {
  f <- terra::focal(r, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  f <- terra::focal(f, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  f <- terra::focal(f, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  terra::ifel(is.na(f), 0, f)
}
cov_sc_ext  <- terra::extend(cov_sc, mesh_ext)
tmp_sc <- tempfile(fileext=".tif"); writeRaster(cov_sc_ext, tmp_sc, overwrite=TRUE, datatype="FLT4S")
cov_sc_ext <- rast(tmp_sc)
cov_sc_filled <- terra::rast(lapply(names(cov_sc_ext), function(nm) fill_na_focal(cov_sc_ext[[nm]])))
names(cov_sc_filled) <- names(cov_sc)

cov_cat_ext <- terra::extend(cov_cat, mesh_ext)
tmp_cat <- tempfile(fileext=".tif"); writeRaster(cov_cat_ext, tmp_cat, overwrite=TRUE, datatype="INT1U")
cov_cat_ext <- rast(tmp_cat)
cov_cat_filled <- terra::rast(lapply(names(cov_cat_ext), function(nm) fill_na_focal(cov_cat_ext[[nm]])))
names(cov_cat_filled) <- names(cov_cat)

# ── 3. Grilla de predicción ──────────────────────────────────────────────────
cat("Creando grilla de predicción...\n")
grid_full <- st_make_grid(domain_final, cellsize=100, what="centers") |>
  st_as_sf() |> st_filter(domain_final)
grid_full$geometry <- grid_full$x
grid_full <- st_set_geometry(grid_full, "geometry"); st_crs(grid_full) <- 9377
coords_g  <- st_coordinates(grid_full)

grid_full$slope_sc      <- terra::extract(cov_sc_filled[["slope_sc"]],  coords_g)[,1]
grid_full$aspect_sc     <- terra::extract(cov_sc_filled[["aspect_sc"]], coords_g)[,1]
grid_full$geo_finos_bin <- terra::extract(cov_cat_filled[["geo_finos_bin"]], coords_g)[,1]
grid_full$cob_antropico <- terra::extract(cov_cat_filled[["cob_antropico"]], coords_g)[,1]
for (nm in c("slope_sc","aspect_sc","geo_finos_bin","cob_antropico"))
  grid_full[[nm]][is.na(grid_full[[nm]])] <- 0

cat(sprintf("  Grilla: %d celdas\n", nrow(grid_full)))

# ── 4. Predicción ────────────────────────────────────────────────────────────
cat("Prediciendo (n.samples=500)... esto puede tardar varios minutos\n")
pred_formula <- ~ exp(Intercept + slope_sc + aspect_sc + geo_finos_bin + cob_antropico)
pi_full <- predict(fit_2m, newdata=grid_full, formula=pred_formula, n.samples=500)
grid_full$lambda_mean <- pi_full$mean
grid_full$lambda_sd   <- pi_full$sd
cat(sprintf("  λ medio = %.2e  |  sd medio = %.2e\n",
            mean(grid_full$lambda_mean, na.rm=TRUE),
            mean(grid_full$lambda_sd,   na.rm=TRUE)))

# ── 5. Hillshade ─────────────────────────────────────────────────────────────
cat("Construyendo hillshade...\n")
dem_hs   <- aggregate(dem_2m, fact = 5, fun = "mean")
slope_r  <- terrain(dem_hs, "slope",  unit = "radians")
aspect_r <- terrain(dem_hs, "aspect", unit = "radians")
hs_nw <- shade(slope_r, aspect_r, angle=45, direction=315, normalize=TRUE)
hs_ne <- shade(slope_r, aspect_r, angle=30, direction= 45, normalize=TRUE)
hs_sw <- shade(slope_r, aspect_r, angle=30, direction=225, normalize=TRUE)
hs_se <- shade(slope_r, aspect_r, angle=30, direction=135, normalize=TRUE)
hs_r  <- (hs_nw*2 + hs_ne + hs_sw + hs_se) / 5
hs_wgs <- project(hs_r, "EPSG:4326")
hs_df  <- as.data.frame(hs_wgs, xy=TRUE); names(hs_df)[3] <- "hillshade"

hs_lyr   <- function(alpha=0.38)
  geom_raster(data=hs_df, aes(x=x, y=y, fill=hillshade),
              inherit.aes=FALSE, alpha=alpha, interpolate=TRUE)
hs_scale <- function()
  scale_fill_gradientn(colors=grey.colors(256,start=0.0,end=0.85), guide="none")
map_crs  <- function() coord_sf(crs=sf::st_crs(4326), expand=FALSE)
map_cart <- function(loc_s="br", loc_n="tr") list(
  annotation_scale(location=loc_s, width_hint=0.25,
                   style="ticks", line_width=0.8, text_cex=0.65),
  annotation_north_arrow(location=loc_n, which_north="true",
    style=north_arrow_fancy_orienteering(text_size=7),
    height=unit(1.2,"cm"), width=unit(1.2,"cm")))
map_theme <- function()
  theme_bw(base_size=9) +
  theme(legend.key.height=unit(0.8,"cm"),
        axis.text=element_text(size=7), axis.title=element_text(size=8))

# ── 6. Figura 4b — Incertidumbre con paleta magma ────────────────────────────
cat("Generando figura 4b (incertidumbre, paleta magma)...\n")
sd_rdf <- {
  sp_v   <- vect(grid_full)
  r_tmpl <- rast(ext(sp_v), resolution=100, crs="EPSG:9377")
  r      <- rasterize(sp_v, r_tmpl, field="lambda_sd", fun=mean)
  r_wgs  <- project(r, "EPSG:4326", method="bilinear")
  out    <- as.data.frame(r_wgs, xy=TRUE); names(out)[3] <- "value"
  out[!is.na(out$value) & out$value > 0, ]
}

p_sd <- ggplot() +
  hs_lyr(1.0) + hs_scale() + ggnewscale::new_scale_fill() +
  geom_raster(data=sd_rdf, aes(x=x, y=y, fill=value),
              alpha=0.35, interpolate=TRUE) +
  scale_fill_viridis_c(
    option = "magma",
    name   = "sd[λ|y]",
    trans  = "log10", labels = scales::scientific) +
  map_crs() + map_cart() +
  labs(x="Longitud", y="Latitud") + map_theme()

png(file.path(FIG,"2m_05_incertidumbre.png"), width=2400, height=2000, res=300)
print(p_sd); dev.off()
cat("  ✓ FIG/2m_05_incertidumbre.png regenerada con paleta magma\n")

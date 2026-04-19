# =============================================================================
# IPP — DESLIZAMIENTOS VALLE DE ABURRÁ — COVARIABLES A 2 m
# =============================================================================
# Pasos:
#   1. Cargar puntos y definir dominio
#   2. Recortar DEM al dominio (evita procesar 641M píxeles)
#   3. Remuestrear DEM recortado 1m → 2m
#   4. Derivar pendiente, aspecto y TWI solo sobre el área recortada
#   5. Extraer covariables en los puntos de deslizamiento (verificación)
#   6. Construir malla triangular (max.edge = 100m)
#   7. Ajustar modelo IPP con inlabru/INLA
#
# NOTA TWI: el TWI requiere acumulación de flujo; se usa un buffer de 5 km
# alrededor del dominio para capturar la cuenca aguas arriba antes del recorte
# final. Esto es mucho menor que el DEM completo (27k×23k píxeles).
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── Dependencias ──────────────────────────────────────────────────────────────
for (pkg in c("INLA","inlabru","sf","terra","ggplot2","gridExtra",
              "viridis","scales","ggspatial","ggnewscale","pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "INLA") {
      install.packages("INLA",
        repos = c(INLA = "https://inla.r-inla-download.org/R/stable",
                  CRAN = "https://cloud.r-project.org"), dep = TRUE)
    } else install.packages(pkg)
  }
}

library(INLA); library(inlabru); library(sf); library(terra)
library(ggplot2); library(gridExtra); library(grid)
library(viridis); library(scales); library(ggspatial)
library(ggnewscale); library(pROC)

# ── Rutas ─────────────────────────────────────────────────────────────────────
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  BASE <- normalizePath(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."),
                        mustWork = FALSE)
} else {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("--file=", args)]
  BASE     <- normalizePath(file.path(dirname(sub("--file=","",file_arg[1])), ".."),
                             mustWork = FALSE)
}
DATA <- file.path(BASE, "DATA")
FIG  <- file.path(BASE, "FIG")
dir.create(FIG, showWarnings = FALSE)
cat("BASE:", BASE, "\n")


# =============================================================================
# 1. PUNTOS Y DOMINIO (primero, para recortar el DEM)
# =============================================================================
cat("\n=== 1. Cargando puntos y definiendo dominio ===\n")
pts_wgs <- st_read(file.path(DATA, "MenM_VdeA_dem.gpkg"), quiet = TRUE)
pts     <- st_transform(pts_wgs, crs = 9377)
if (attr(pts, "sf_column") != "geometry")
  pts <- st_sf(st_drop_geometry(pts), geometry = st_geometry(pts), crs = 9377)
pts_coords <- st_coordinates(pts)
cat("  Puntos:", nrow(pts), "\n")

# Dominio = convex hull de los puntos + 500m
domain_buf <- st_convex_hull(st_union(pts)) |> st_buffer(500) |> st_as_sf()
st_crs(domain_buf) <- 9377

# Buffer extra de 5km para capturar cuenca aguas arriba en el cálculo del TWI
domain_twi_buf <- st_convex_hull(st_union(pts)) |> st_buffer(5000) |> st_as_sf()
st_crs(domain_twi_buf) <- 9377


# =============================================================================
# 2 & 3. DEM 2m Y MORFOMÉTRICAS — omitir si los archivos ya existen
# =============================================================================
f_dem   <- file.path(DATA, "dem_2m.tif")
f_slope <- file.path(DATA, "cov2m_pendiente.tif")
f_asp   <- file.path(DATA, "cov2m_aspecto.tif")
f_twi   <- file.path(DATA, "cov2m_twi.tif")

if (all(file.exists(f_dem, f_slope, f_asp, f_twi))) {
  cat("\n=== 2-3. Cargando rasters 2m ya existentes ===\n")
  dem_2m    <- rast(f_dem)
  slope_2m  <- rast(f_slope)
  aspect_2m <- rast(f_asp)
  twi_2m    <- rast(f_twi)
  cat("  dem_2m:    ", dim(dem_2m)[1:2],   "\n")
  cat("  slope_2m:  ", dim(slope_2m)[1:2], "\n")
  cat("  aspect_2m: ", dim(aspect_2m)[1:2],"\n")
  cat("  twi_2m:    ", dim(twi_2m)[1:2],   "\n")
} else {

cat("\n=== 2. Recortando DEM al dominio + buffer TWI (5 km) ===\n")
t0 <- proc.time()

dem_1m <- rast(file.path(DATA, "dem_medellin_1m.tif"))
cat("  DEM completo:", dim(dem_1m)[1:2], "=",
    format(ncell(dem_1m), big.mark=","), "píxeles\n")

# Recortar rectangular (sin mask): el DEM ya tiene NAs fuera del área válida
dem_crop_1m <- crop(dem_1m, vect(domain_twi_buf))
rm(dem_1m); gc()
cat("  DEM recortado:", dim(dem_crop_1m)[1:2], "=",
    format(ncell(dem_crop_1m), big.mark=","), "píxeles\n")

dem_2m <- aggregate(dem_crop_1m, fact = 2, fun = "mean")
names(dem_2m) <- "dem"
rm(dem_crop_1m); gc()
cat("  DEM 2m:", dim(dem_2m)[1:2], "=",
    format(ncell(dem_2m), big.mark=","), "píxeles\n")
writeRaster(dem_2m, f_dem, overwrite=TRUE, datatype="FLT4S",
            gdal=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=6"))
cat(sprintf("  Tiempo: %.1f s\n", (proc.time()-t0)["elapsed"]))

cat("\n=== 3a. Pendiente 2m ===\n"); t0 <- proc.time()
slope_2m <- terrain(dem_2m, v="slope",  unit="degrees", neighbors=8)
names(slope_2m) <- "slope"
writeRaster(slope_2m, f_slope, overwrite=TRUE, datatype="FLT4S",
            gdal=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=6"))
cat(sprintf("  Rango: %.2f – %.2f °  |  Tiempo: %.1f s\n",
            minmax(slope_2m)[1], minmax(slope_2m)[2], (proc.time()-t0)["elapsed"]))

cat("\n=== 3b. Aspecto 2m ===\n"); t0 <- proc.time()
aspect_2m <- terrain(dem_2m, v="aspect", unit="degrees", neighbors=8)
names(aspect_2m) <- "aspect"
writeRaster(aspect_2m, f_asp, overwrite=TRUE, datatype="FLT4S",
            gdal=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=6"))
cat(sprintf("  Tiempo: %.1f s\n", (proc.time()-t0)["elapsed"]))

cat("\n=== 3c. TWI 2m ===\n"); t0 <- proc.time()
cat("  Dirección de flujo D8...\n")
flowdir_2m <- terrain(dem_2m, v="flowdir", neighbors=8)
cat(sprintf("  Dirección: %.1f s\n", (proc.time()-t0)["elapsed"]))
t1 <- proc.time()
cat("  Acumulación de flujo...\n")
flowacc_2m <- flowAccumulation(flowdir_2m, progress=TRUE)
cat(sprintf("  Acumulación: %.1f s\n", (proc.time()-t1)["elapsed"]))
res_m     <- res(dem_2m)[1]
A_esp     <- flowacc_2m * res_m
slope_rad <- slope_2m * (pi/180)
tan_slope <- tan(slope_rad)
tan_slope <- ifel(tan_slope < 0.001, 0.001, tan_slope)
twi_2m    <- log(A_esp / tan_slope)
names(twi_2m) <- "twi"
rm(flowdir_2m, flowacc_2m, A_esp, slope_rad, tan_slope); gc()
writeRaster(twi_2m, f_twi, overwrite=TRUE, datatype="FLT4S",
            gdal=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=6"))
cat(sprintf("  Rango TWI: %.3f – %.3f  |  Tiempo total: %.1f s\n",
            minmax(twi_2m)[1], minmax(twi_2m)[2], (proc.time()-t0)["elapsed"]))

} # fin else rasters existentes


# =============================================================================
# 4. EXTRACCIÓN EN PUNTOS DE DESLIZAMIENTO (verificación)
# =============================================================================
cat("\n=== 4. Covariables 2m en puntos de deslizamiento ===\n")
pts$slope_2m  <- terra::extract(slope_2m,  pts_coords)[, 1]
pts$aspect_2m <- terra::extract(aspect_2m, pts_coords)[, 1]
pts$twi_2m    <- terra::extract(twi_2m,    pts_coords)[, 1]

cat("  Estadísticas en puntos de deslizamiento:\n")
print(summary(data.frame(slope = pts$slope_2m,
                          aspect = pts$aspect_2m,
                          twi    = pts$twi_2m)))
cat("  NAs: slope =", sum(is.na(pts$slope_2m)),
    " aspecto =", sum(is.na(pts$aspect_2m)),
    " TWI =", sum(is.na(pts$twi_2m)), "\n")


# =============================================================================
# 5. DOMINIO FINAL (intersección con máscara del DEM 2m recortado)
# =============================================================================
cat("\n=== 5. Dominio final ===\n")
dem_mask <- as.polygons(ifel(!is.na(dem_2m), 1, NA)) |>
  st_as_sf() |> st_union() |> st_as_sf()
st_crs(dem_mask) <- 9377

domain_final <- st_intersection(domain_buf, dem_mask)
domain_final <- st_sf(geometry = st_union(st_geometry(domain_final)), crs = 9377)
cat("  Área dominio:", round(as.numeric(st_area(domain_final))/1e6, 2), "km²\n")


# =============================================================================
# 6. HILLSHADE PARA MAPAS
# =============================================================================
# Agregamos el DEM 2m a 10m (factor 5) para el hillshade
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
rm(slope_r, aspect_r, hs_nw, hs_ne, hs_sw, hs_se, hs_r); gc()

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
    height=unit(1.2,"cm"), width=unit(1.2,"cm"))
)
map_theme <- function()
  theme_bw(base_size=9) +
  theme(legend.key.height=unit(0.8,"cm"),
        axis.text=element_text(size=7), axis.title=element_text(size=8))


# =============================================================================
# 7. MALLA TRIANGULAR (max.edge = 100m)
# =============================================================================
# Con celdas de 2m, max.edge=100m captura bien la variación espacial de λ(s)
# sin el costo computacional de 50m (que generaría ~4× más nodos).
cat("\n=== 7. Construyendo malla triangular (max.edge = 100 m) ===\n")
t0 <- proc.time()

bnd  <- fmesher::fm_as_segm(domain_final)
mesh <- fmesher::fm_mesh_2d(
  boundary = bnd,
  max.edge = c(100, 500),    # interior ≤100m, buffer exterior ≤500m
  offset   = c(200, 2000),   # margen interior 200m, exterior 2km
  cutoff   = 20,             # mínima distancia entre nodos = 20m
  crs      = sf::st_crs(9377)
)
cat(sprintf("  Nodos: %d  |  Triángulos: %d  |  Tiempo: %.1f s\n",
            mesh$n, nrow(mesh$graph$tv), (proc.time()-t0)["elapsed"]))

# Figura: malla
mesh_loc  <- mesh$loc[,1:2]; mesh_tv <- mesh$graph$tv
mesh_tris <- st_sf(geometry=st_sfc(
  lapply(seq_len(nrow(mesh_tv)), function(i) {
    idx <- c(mesh_tv[i,], mesh_tv[i,1])
    st_polygon(list(mesh_loc[idx,]))
  }), crs=9377))

png(file.path(FIG,"2m_01_mesh.png"), width=2400, height=2000, res=300)
ggplot() +
  hs_lyr(1.0) + hs_scale() +
  geom_sf(data=mesh_tris, fill=NA, color="steelblue", lwd=0.2, alpha=0.4) +
  geom_sf(data=pts, size=0.7, color="#d62728", alpha=0.6) +
  geom_sf(data=domain_final, fill=NA, color="black", lwd=0.6) +
  map_crs() + map_cart() +
  labs(x="Longitud", y="Latitud",
       caption=paste0("Malla 2m | Nodos: ",mesh$n,
                      "  Triángulos: ",nrow(mesh$graph$tv))) +
  theme_bw(base_size=10) +
  theme(axis.text=element_text(size=7), plot.caption=element_text(size=7))
dev.off()
cat("  ✓ Fig 2m_01: malla\n")


# =============================================================================
# 8. COVARIABLES — RECORTE FINAL, ESCALA Y RELLENO DE NAs
# =============================================================================
cat("\n=== 8. Escalando covariables y rellenando NAs ===\n")

# Recortar morfométricas al dominio final (sin el buffer TWI de 5km)
slope_2m  <- crop(slope_2m,  vect(domain_final), extend=TRUE)
aspect_2m <- crop(aspect_2m, vect(domain_final), extend=TRUE)
twi_2m    <- crop(twi_2m,    vect(domain_final), extend=TRUE)

# Covariables categóricas (cobertura, geotecnia)
cov_cob   <- rast(file.path(DATA, "cobertura_2024.tif"))
cov_geo   <- rast(file.path(DATA, "zoni_geotecnica.tif"))
# Remuestrear al template de slope_2m (misma extensión 2m, mismo CRS)
# Esto reduce de 641M píxeles (1m) a ~84M (2m dominio) antes de operar
cov_cob_r <- project(cov_cob, slope_2m, method="near")
cov_geo_r <- project(cov_geo, slope_2m, method="near")

cob_antropico <- ifel(cov_cob_r == 1, 1L, 0L)
names(cob_antropico) <- "cob_antropico"

geo_finos_bin <- ifel(cov_geo_r %in% c(4,5,17,11,12,13,14), 1L, 0L)
names(geo_finos_bin) <- "geo_finos_bin"

cat("  cob_antropico prop =", round(mean(values(cob_antropico),na.rm=TRUE),3), "\n")
cat("  geo_finos_bin prop =", round(mean(values(geo_finos_bin),na.rm=TRUE),3), "\n")

# Escalar covariables continuas (media=0, sd=1)
sc <- function(r) { m <- global(r,"mean",na.rm=TRUE)[[1]]; s <- global(r,"sd",na.rm=TRUE)[[1]]; (r-m)/s }
slope_sc_r  <- sc(slope_2m);  names(slope_sc_r)  <- "slope_sc"
aspect_sc_r <- sc(aspect_2m); names(aspect_sc_r) <- "aspect_sc"
cov_sc      <- c(slope_sc_r, aspect_sc_r)
rm(slope_sc_r, aspect_sc_r); gc()

# Extender extent para cubrir los nodos del buffer exterior de la malla
mesh_ext <- terra::ext(
  min(mesh$loc[,1]) - 200, max(mesh$loc[,1]) + 200,
  min(mesh$loc[,2]) - 200, max(mesh$loc[,2]) + 200
)

fill_na_focal <- function(r) {
  # Dos pasadas de w=3 (equivale a rango ~5 celdas) más estable que w=11
  f <- terra::focal(r, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  f <- terra::focal(f, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  f <- terra::focal(f, w=3, fun=mean, na.policy="only", na.rm=TRUE)
  terra::ifel(is.na(f), 0, f)
}

cat("  Extendiendo y rellenando NAs (continuas)...\n")
cov_sc_ext  <- terra::extend(cov_sc, mesh_ext)
# Escribir a disco para que focal opere sobre raster respaldado en archivo
tmp_sc      <- tempfile(fileext=".tif")
writeRaster(cov_sc_ext, tmp_sc, overwrite=TRUE, datatype="FLT4S")
cov_sc_ext  <- rast(tmp_sc)
cov_sc_filled <- terra::rast(lapply(names(cov_sc_ext),
                                    function(nm) fill_na_focal(cov_sc_ext[[nm]])))
names(cov_sc_filled) <- names(cov_sc)
rm(cov_sc_ext); gc()

cat("  Extendiendo y rellenando NAs (categóricas)...\n")
cov_cat     <- c(cob_antropico, geo_finos_bin)
cov_cat_ext <- terra::extend(cov_cat, mesh_ext)
tmp_cat     <- tempfile(fileext=".tif")
writeRaster(cov_cat_ext, tmp_cat, overwrite=TRUE, datatype="INT1U")
cov_cat_ext <- rast(tmp_cat)
cov_cat_filled <- terra::rast(lapply(names(cov_cat_ext),
                                     function(nm) fill_na_focal(cov_cat_ext[[nm]])))
names(cov_cat_filled) <- names(cov_cat)
rm(cov_cat_ext); gc()

# Verificar sin NAs en nodos de la malla
vals_chk <- terra::extract(cov_sc_filled[["slope_sc"]], mesh$loc[,1:2])[,1]
cat("  NAs en nodos tras relleno:", sum(is.na(vals_chk)), "de", length(vals_chk), "\n")

# Pre-extraer covariables escaladas en los puntos de deslizamiento
pts$slope_sc      <- terra::extract(cov_sc_filled[["slope_sc"]],  pts_coords)[,1]
pts$aspect_sc     <- terra::extract(cov_sc_filled[["aspect_sc"]], pts_coords)[,1]
pts$cob_antropico <- terra::extract(cov_cat_filled[["cob_antropico"]], pts_coords)[,1]
pts$geo_finos_bin <- terra::extract(cov_cat_filled[["geo_finos_bin"]], pts_coords)[,1]
cat("  Covariables extraídas en", nrow(pts), "puntos. NAs slope:",
    sum(is.na(pts$slope_sc)), "\n")


# =============================================================================
# 9. MODELO IPP
# =============================================================================
# Liberar objetos grandes antes de INLA para reducir presión de memoria
rm(dem_2m, slope_2m, aspect_2m, twi_2m, cov_cob, cov_cob_r, cov_geo, cov_geo_r,
   cob_antropico, geo_finos_bin, cov_sc, cov_cat)
gc()
cat(sprintf("  RAM antes de INLA: %.0f MB usados\n",
            sum(gc()[,2]) ))

# log λ(s) = β₀ + β₁·pendiente + β₂·aspecto + β₃·geo_finos_bin + β₄·cob_antropico
# Covariables morfométricas derivadas del DEM remuestreado a 2m. TWI omitido (no significativo).
cat("\n=== 9. Ajustando modelo IPP (DEM 2m, sin TWI) ===\n")
t0 <- proc.time()

components <- ~ Intercept(1) +
  slope_sc     (cov_sc_filled[["slope_sc"]],       model="linear") +
  aspect_sc    (cov_sc_filled[["aspect_sc"]],      model="linear") +
  cob_antropico(cov_cat_filled[["cob_antropico"]], model="linear") +
  geo_finos_bin(cov_cat_filled[["geo_finos_bin"]], model="linear")

fit_2m <- tryCatch(
  lgcp(
    components = components,
    data       = pts,
    formula    = geometry ~ .,
    samplers   = domain_final,
    domain     = list(geometry = mesh),
    options    = list(
      control.fixed = list(mean=0, prec=1.0),
      control.inla  = list(int.strategy="ccd", strategy="laplace"),
      verbose = TRUE
    )
  ),
  error = function(e) { cat("ERROR en lgcp():", conditionMessage(e), "\n"); stop(e) }
)

cat(sprintf("  ✓ Modelo ajustado  |  Tiempo: %.1f min\n",
            (proc.time()-t0)["elapsed"]/60))

saveRDS(fit_2m, file=file.path(BASE, "results_ipp_2m.rds"))
cat("  Guardado en results_ipp_2m.rds\n")

# Resumen de coeficientes
fixed_2m <- fit_2m$summary.fixed
cat("\n=== EFECTOS FIJOS — IPP 2m ===\n")
print(round(fixed_2m[, c("mean","sd","0.025quant","0.975quant","mode")], 4))

var_map <- c(slope_sc="Pendiente", aspect_sc="Aspecto",
             cob_antropico="Urbano", geo_finos_bin="Suelos finos")
cat("\n--- Resumen ---\n")
for (v in c("slope_sc","aspect_sc","geo_finos_bin","cob_antropico"))
  if (v %in% rownames(fixed_2m))
    cat(sprintf("  %-18s %8.3f  IC95=[%7.3f, %7.3f]\n",
                var_map[v], fixed_2m[v,"mean"],
                fixed_2m[v,"0.025quant"], fixed_2m[v,"0.975quant"]))


# =============================================================================
# 10. PREDICCIÓN EN GRILLA (100m)
# =============================================================================
cat("\n=== 10. Grilla de predicción (100 m) ===\n")
t0 <- proc.time()

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

cat(sprintf("  Grilla: %d celdas  |  Tiempo: %.1f s\n",
            nrow(grid_full), (proc.time()-t0)["elapsed"]))

pred_formula <- ~ exp(Intercept + slope_sc + aspect_sc +
                        geo_finos_bin + cob_antropico)
cat("  Prediciendo intensidad...\n")
pi_full <- predict(fit_2m, newdata=grid_full, formula=pred_formula, n.samples=500)
grid_full$lambda_mean <- pi_full$mean
grid_full$lambda_sd   <- pi_full$sd
cat(sprintf("  ✓ λ medio = %.2e\n", mean(grid_full$lambda_mean, na.rm=TRUE)))


# =============================================================================
# 11. EVALUACIÓN (hexbins + AUC)
# =============================================================================
cat("\n=== 11. Evaluación ===\n")
hex_grid <- st_make_grid(domain_final, cellsize=500, square=FALSE) |>
  st_as_sf() |> st_filter(domain_final)
hex_grid$cell_id <- seq_len(nrow(hex_grid))
hex_grid$area_m2 <- as.numeric(st_area(hex_grid))

pts_in_hex  <- st_join(pts, hex_grid, join=st_within)
obs_count   <- as.data.frame(table(pts_in_hex$cell_id))
names(obs_count) <- c("cell_id","n_obs")
obs_count$cell_id <- as.integer(as.character(obs_count$cell_id))

grid_in_hex <- st_join(grid_full, hex_grid, join=st_within)
pred_by_hex <- aggregate(lambda_mean ~ cell_id, data=grid_in_hex, FUN=mean)

hex_compare <- merge(hex_grid, obs_count,    by="cell_id", all.x=TRUE)
hex_compare <- merge(hex_compare, pred_by_hex, by="cell_id", all.x=TRUE)
hex_compare$n_obs[is.na(hex_compare$n_obs)] <- 0
hex_compare$n_exp <- hex_compare$lambda_mean * hex_compare$area_m2
hex_compare$n_exp[is.na(hex_compare$n_exp)] <- 0
hex_compare$pearson_res <- with(hex_compare,
  ifelse(n_exp > 0 & is.finite(n_exp), (n_obs - n_exp)/sqrt(n_exp), NA))

roc_data <- hex_compare[!is.na(hex_compare$n_exp) & is.finite(hex_compare$n_exp), ]
roc_obj  <- roc(as.numeric(roc_data$n_obs > 0), roc_data$n_exp, quiet=TRUE)
auc_val  <- auc(roc_obj)
max_n    <- quantile(c(hex_compare$n_obs, hex_compare$n_exp), 0.99, na.rm=TRUE)
hex_sc   <- hex_compare[hex_compare$n_exp <= max_n & is.finite(hex_compare$n_exp), ]
r2_val   <- round(cor(hex_sc$n_obs, hex_sc$n_exp, use="complete.obs")^2, 3)
cat(sprintf("  AUC = %.3f  |  r² = %.3f\n", auc_val, r2_val))


# =============================================================================
# 12. FIGURAS
# =============================================================================
cat("\n=== 12. Figuras ===\n")

lambda_rdf <- {
  sp_v   <- vect(grid_full)
  r_tmpl <- rast(ext(sp_v), resolution=100, crs="EPSG:9377")
  r      <- rasterize(sp_v, r_tmpl, field="lambda_mean", fun=mean)
  r_wgs  <- project(r, "EPSG:4326", method="bilinear")
  out    <- as.data.frame(r_wgs, xy=TRUE); names(out)[3] <- "value"
  out[!is.na(out$value) & out$value > 0, ]
}
sd_rdf <- {
  sp_v   <- vect(grid_full)
  r_tmpl <- rast(ext(sp_v), resolution=100, crs="EPSG:9377")
  r      <- rasterize(sp_v, r_tmpl, field="lambda_sd", fun=mean)
  r_wgs  <- project(r, "EPSG:4326", method="bilinear")
  out    <- as.data.frame(r_wgs, xy=TRUE); names(out)[3] <- "value"
  out[!is.na(out$value) & out$value > 0, ]
}

# Efectos fijos
coef_df <- data.frame(variable=rownames(fixed_2m), mean=fixed_2m$mean,
                      lo=fixed_2m$`0.025quant`, hi=fixed_2m$`0.975quant`)
coef_df <- coef_df[coef_df$variable != "Intercept", ]
coef_df$label <- factor(var_map[coef_df$variable],
  levels=c("Pendiente","Aspecto","Suelos finos","Urbano"))
coef_df <- coef_df[!is.na(coef_df$label), ]

png(file.path(FIG,"2m_02_efectos_fijos.png"), width=2400, height=1600, res=300)
ggplot(coef_df, aes(x=label, y=mean, ymin=lo, ymax=hi)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray60") +
  geom_pointrange(size=0.9, linewidth=1.3, color="#e6550d") +
  labs(x=NULL, y=expression(hat(beta)~"(log-intensidad)")) +
  theme_bw(base_size=10)
dev.off(); cat("  ✓ 2m_02: efectos fijos\n")

# Intensidad predicha — paleta verde (bajo) → rojo (alto)
p_lambda <- ggplot() +
  hs_lyr(1.0) + hs_scale() + ggnewscale::new_scale_fill() +
  geom_raster(data=lambda_rdf, aes(x=x, y=y, fill=value),
              alpha=0.35, interpolate=TRUE) +
  scale_fill_gradientn(
    colors = c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c"),
    name   = "λ(s)\n(eventos/m²)",
    trans  = "log10", labels = scales::scientific) +
  geom_sf(data=pts, size=0.4, color="white", alpha=0.4) +
  map_crs() + map_cart() +
  labs(x="Longitud", y="Latitud") + map_theme()
png(file.path(FIG,"2m_04_intensidad.png"), width=2400, height=2000, res=300)
print(p_lambda); dev.off(); cat("  ✓ 2m_04: intensidad\n")

# Incertidumbre — paleta magma (violeta → amarillo), diferente de la intensidad
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
print(p_sd); dev.off(); cat("  ✓ 2m_05: incertidumbre\n")

# Obs vs predicho
png(file.path(FIG,"2m_08_obs_vs_predicho.png"), width=2000, height=2000, res=300)
ggplot(hex_sc, aes(x=n_exp, y=n_obs)) +
  geom_abline(slope=1, intercept=0, color="red", lty=2, lwd=0.8) +
  geom_point(alpha=0.5, size=1.2, color="#3182bd") +
  geom_smooth(method="lm", color="#e6550d", se=TRUE, lwd=0.8) +
  labs(x="Conteo esperado E[N(A)|y]", y="Conteo observado N(A)") +
  annotate("text", x=max_n*0.7, y=max_n*0.1,
           label=paste0("r² = ",r2_val), size=3.5, color="#e6550d") +
  theme_bw(base_size=10)
dev.off(); cat("  ✓ 2m_08: obs vs predicho\n")

# Residuos
lim_res <- quantile(abs(hex_compare$pearson_res), 0.99, na.rm=TRUE)
p_res <- ggplot() +
  hs_lyr(1.0) + hs_scale() + ggnewscale::new_scale_fill() +
  geom_sf(data=hex_compare, aes(fill=pearson_res), color=NA, alpha=0.25) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#d73027",
                       midpoint=0, limits=c(-lim_res,lim_res),
                       name="Residuo\nPearson", na.value="gray90") +
  geom_sf(data=pts, size=0.4, color="black", alpha=0.3) +
  map_crs() + map_cart() +
  labs(x="Longitud", y="Latitud") + map_theme()
png(file.path(FIG,"2m_09_residuos.png"), width=2400, height=2000, res=300)
print(p_res); dev.off(); cat("  ✓ 2m_09: residuos\n")

# Mapas obs vs esperado (dos paneles)
p_obs <- ggplot() +
  hs_lyr(1.0) + hs_scale() + ggnewscale::new_scale_fill() +
  geom_sf(data=hex_compare, aes(fill=n_obs), color=NA, alpha=0.35) +
  scale_fill_gradientn(
    colors=c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c"),
    name="N observado") +
  geom_sf(data=pts, size=0.3, color="white", alpha=0.4) +
  map_crs() + map_cart(loc_s="bl", loc_n="tl") +
  labs(x="Longitud", y="Latitud") + map_theme()

p_exp <- ggplot() +
  hs_lyr(1.0) + hs_scale() + ggnewscale::new_scale_fill() +
  geom_sf(data=hex_compare, aes(fill=n_exp), color=NA, alpha=0.35) +
  scale_fill_gradientn(
    colors=c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c"),
    name="N esperado", trans="log1p") +
  map_crs() + map_cart() +
  labs(x="Longitud", y="Latitud") + map_theme()

png(file.path(FIG,"2m_10_obs_vs_predicho_mapas.png"), width=4000, height=2000, res=300)
gridExtra::grid.arrange(p_obs, p_exp, ncol=2); dev.off()
cat("  ✓ 2m_10: mapas obs vs predicho\n")

# ROC
png(file.path(FIG,"2m_12_roc.png"), width=1800, height=1800, res=300)
plot(roc_obj, col="#3182bd", lwd=2)
text(0.3, 0.1, paste0("AUC = ",round(auc_val,3)), cex=1.2, col="#3182bd")
abline(0, 1, lty=2, col="gray")
dev.off(); cat("  ✓ 2m_12: ROC  AUC =",round(auc_val,3),"\n")


# =============================================================================
# 13. RESUMEN FINAL
# =============================================================================
cat("\n============================================================\n")
cat("IPP DEM 2m — Valle de Aburrá\n")
cat("============================================================\n")
cat("DEM:               2m (recortado al dominio + 5km buffer, agregado desde 1m)\n")
cat("Malla max.edge:    100m interior / 500m exterior\n")
cat(sprintf("Malla nodos:       %d\n", mesh$n))
cat(sprintf("Malla triángulos:  %d\n", nrow(mesh$graph$tv)))
cat(sprintf("Eventos totales:   %d\n", nrow(pts)))
cat(sprintf("Área dominio:      %.2f km²\n", as.numeric(st_area(domain_final))/1e6))
cat("\n--- Efectos fijos ---\n")
for (v in c("slope_sc","aspect_sc","geo_finos_bin","cob_antropico"))
  if (v %in% rownames(fixed_2m))
    cat(sprintf("  %-18s %8.3f  IC95=[%7.3f, %7.3f]\n",
                var_map[v], fixed_2m[v,"mean"],
                fixed_2m[v,"0.025quant"], fixed_2m[v,"0.975quant"]))
cat(sprintf("\nAUC (ROC): %.3f\n", auc_val))
cat(sprintf("r²  (hex):  %.3f\n", r2_val))
cat("Resultados: results_ipp_2m.rds\n")
cat("Figuras:    FIG/2m_*.png\n")
cat("============================================================\n")

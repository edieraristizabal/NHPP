# =============================================================================
# LOG-GAUSSIAN COX PROCESS (LGCP) — DESLIZAMIENTOS VALLE DE ABURRÁ
# Método: SPDE-INLA (Lindgren et al. 2011) + inlabru
# =============================================================================
#
# FUNDAMENTO ESTADÍSTICO
# ─────────────────────────────────────────────────────────────────────────────
# Un LGCP modela la ubicación de eventos (deslizamientos) como una realización
# de un proceso puntual de Poisson no homogéneo cuya intensidad λ(s) varía en
# el espacio de forma log-gaussiana:
#
#   log λ(s) = β₀ + β₁·elev(s) + β₂·pend(s) + β₃·asp(s) + β₄·log_acc(s)
#              + ξ(s)
#
# donde ξ(s) es un Campo Aleatorio Gaussiano (GRF) estacionario con
# covarianza de Matérn que captura la autocorrelación espacial residual
# (clustering no explicado por las covariables).
#
# La aproximación SPDE (Stochastic PDE) de Lindgren et al. (2011) convierte
# el GRF continuo en un proceso Markoviano Gaussiano (GMRF) definido sobre
# una malla triangular, lo que hace el cálculo tratable con INLA.
#
# ESTRUCTURA DE CARPETAS
#   CODE/  ← este script
#   DATA/  ← GeoPackage de puntos + rasters de covariables
#   FIG/   ← todas las figuras de salida
# =============================================================================


# ─── 0. LIBRERÍAS ─────────────────────────────────────────────────────────────
# INLA: motor de inferencia variacional (Integrated Nested Laplace Approximation)
# inlabru: interfaz de alto nivel para modelos de proceso puntual con INLA
# sf: manejo de datos vectoriales (puntos, polígonos)
# terra: manejo de rasters (covariables topográficas)
# ggplot2 + patchwork + viridis: visualización científica

options(repos = c(CRAN = "https://cloud.r-project.org"))

# inlabru 2.10+ requiere INLA >= 23.1.31; forzar actualización si la versión es antigua
if (!requireNamespace("INLA", quietly = TRUE) ||
    packageVersion("INLA") < "23.1.31") {
  install.packages("INLA",
    repos = c(INLA = "https://inla.r-inla-download.org/R/stable", CRAN = "https://cloud.r-project.org"),
    dep = TRUE
  )
}
if (!requireNamespace("inlabru",  quietly = TRUE)) install.packages("inlabru")
if (!requireNamespace("sf",       quietly = TRUE)) install.packages("sf")
if (!requireNamespace("terra",    quietly = TRUE)) install.packages("terra")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("viridis",  quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("scales",   quietly = TRUE)) install.packages("scales")
if (!requireNamespace("pROC",     quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("ggspatial",  quietly = TRUE)) install.packages("ggspatial")
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")

library(INLA)
library(inlabru)
library(sf)
library(terra)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(scales)
library(ggspatial)
library(ggnewscale)

# ─── 1. RUTAS ──────────────────────────────────────────────────────────────────
# El script se ejecuta desde CODE/; las rutas son relativas al directorio padre
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  BASE <- normalizePath(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."),
                        mustWork = FALSE)
} else {
  # Ruta cuando se ejecuta con Rscript fuera de RStudio
  args      <- commandArgs(trailingOnly = FALSE)
  file_arg  <- args[grep("--file=", args)]
  script_path <- sub("--file=", "", file_arg[1])
  BASE <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
}

DATA  <- file.path(BASE, "DATA")
FIG   <- file.path(BASE, "FIG")
dir.create(FIG, showWarnings = FALSE)

cat("BASE:", BASE, "\nDATA:", DATA, "\nFIG:", FIG, "\n")


# ─── 2. CARGA DE DATOS ────────────────────────────────────────────────────────

# 2.1 Puntos de deslizamientos con covariables extraídas (EPSG:4326)
pts_wgs <- st_read(file.path(DATA, "MenM_VdeA_dem.gpkg"), quiet = TRUE)
cat("Puntos cargados:", nrow(pts_wgs), "\n")
cat("Columnas:", paste(names(pts_wgs), collapse = " | "), "\n")

# Reproyectar a MAGNA-SIRGAS (metros) — necesario para la malla SPDE
pts <- st_transform(pts_wgs, crs = 9377)
# inlabru busca la columna "geometry"; renombrar si viene con otro nombre del gpkg
if (attr(pts, "sf_column") != "geometry") {
  pts <- st_sf(st_drop_geometry(pts), geometry = st_geometry(pts), crs = 9377)
}

# 2.2 Rasters de covariables a 10m (MAGNA-SIRGAS)
cov_elev    <- rast(file.path(DATA, "cov_elevacion.tif"))
cov_slope   <- rast(file.path(DATA, "cov_pendiente.tif"))
cov_aspect  <- rast(file.path(DATA, "cov_aspecto.tif"))
cov_logacc  <- rast(file.path(DATA, "cov_log_area_acum.tif"))

# cov_elev se usa solo para definir el dominio (dem_mask); no entra al modelo
cat("Rasters cargados. Resolución:", res(cov_elev), "m\n")
cat("Extensión:", as.character(ext(cov_elev)), "\n")


# ─── 3. DOMINIO DE ESTUDIO ────────────────────────────────────────────────────
#
# El dominio Ω es el área dentro de la cual se modelan los deslizamientos.
# Se define como la envolvente convexa de los puntos ampliada 500m,
# intersectada con la extensión del DEM (zona con datos válidos).
#
# Usar un dominio bien definido es crítico para el LGCP porque la intensidad
# integrada ∫Ω λ(s)ds aparece explícitamente en la verosimilitud de Poisson.

# Convex hull de los puntos + buffer 500m
domain <- st_convex_hull(st_union(pts)) |>
  st_buffer(500) |>
  st_as_sf()
st_crs(domain) <- 9377

# Máscara del DEM (pixels con dato válido) como polígono de dominio real
dem_mask <- as.polygons(ifel(!is.na(cov_elev), 1, NA)) |>
  st_as_sf() |>
  st_union() |>
  st_as_sf()
st_crs(dem_mask) <- 9377

# Intersección: dominio = convex hull ∩ área del DEM
domain_final <- st_intersection(domain, dem_mask)
# inlabru 2.10+ requiere que el sampler tenga la geometría nombrada "geometry"
domain_final <- st_sf(geometry = st_union(st_geometry(domain_final)), crs = 9377)

cat("Área del dominio:", round(as.numeric(st_area(domain_final)) / 1e6, 2), "km²\n")

# ─── HILLSHADE Y HELPERS PARA MAPAS ──────────────────────────────────────────
# Hillshade multidireccional: combina 4 ángulos de sol para un relieve más realista.
# Un solo ángulo produce sombras planas; la combinación enfatiza crestas/valles
# desde múltiples direcciones, similar a la representación de terreno en cartografía.
slope_r  <- terrain(cov_elev, "slope",  unit = "radians")
aspect_r <- terrain(cov_elev, "aspect", unit = "radians")

hs_nw  <- shade(slope_r, aspect_r, angle = 45, direction = 315, normalize = TRUE)  # NO
hs_ne  <- shade(slope_r, aspect_r, angle = 30, direction =  45, normalize = TRUE)  # NE
hs_sw  <- shade(slope_r, aspect_r, angle = 30, direction = 225, normalize = TRUE)  # SO
hs_se  <- shade(slope_r, aspect_r, angle = 30, direction = 135, normalize = TRUE)  # SE

# Combinación ponderada: dirección principal NO (peso 2) + tres secundarias
hs_r   <- (hs_nw * 2 + hs_ne + hs_sw + hs_se) / 5

hs_wgs <- project(hs_r, "EPSG:4326")
hs_df  <- as.data.frame(hs_wgs, xy = TRUE)
names(hs_df)[3] <- "hillshade"

# Capa de hillshade reutilizable
hs_lyr <- function(alpha = 0.38) {
  geom_raster(data = hs_df, aes(x = x, y = y, fill = hillshade),
              inherit.aes = FALSE, alpha = alpha, interpolate = TRUE)
}
hs_scale <- function() {
  scale_fill_gradientn(colors = grey.colors(256, start = 0.0, end = 1.0),
                       guide = "none")
}
map_crs <- function() coord_sf(crs = sf::st_crs(4326), expand = FALSE)
map_cart <- function(loc_s = "br", loc_n = "tr") {
  list(
    annotation_scale(location = loc_s, width_hint = 0.25,
                     style = "ticks", line_width = 0.8, text_cex = 0.65),
    annotation_north_arrow(location = loc_n, which_north = "true",
      style = north_arrow_fancy_orienteering(text_size = 7),
      height = unit(1.2, "cm"), width = unit(1.2, "cm"))
  )
}


# ─── 4. MALLA SPDE ───────────────────────────────────────────────────────────
#
# La malla triangular es el soporte del GRF aproximado por SPDE.
# Parámetros clave:
#   max.edge : tamaño máximo de las aristas [interior, extensión exterior]
#              Regla práctica: ≤ 1/5 del rango espacial esperado
#   offset   : extensión de la malla más allá del dominio (evita efecto borde)
#   cutoff   : distancia mínima entre nodos (evita triángulos degenerados)
#
# Con datos en metros y extensión ~24km, usamos max.edge = c(1000, 3000) m

# Figura 0: Mapa del inventario de deslizamientos de masa
png(file.path(FIG, "00_inventario_deslizamientos.png"), width = 2400, height = 2400, res = 300)
ggplot() +
  hs_lyr(0.28) + hs_scale() +
  geom_sf(data = domain_final, fill = NA, color = "black", lwd = 0.5) +
  geom_sf(data = pts, size = 0.7, color = "#d62728", alpha = 0.7, shape = 21,
          fill = "#d62728", stroke = 0.1) +
  map_crs() +
  map_cart() +
  labs(x = "Longitud", y = "Latitud") +
  theme_bw(base_size = 10) +
  theme(axis.text  = element_text(size = 7),
        axis.title = element_text(size = 8))
dev.off()
cat("✓ Fig 00: inventario\n")

bnd <- fmesher::fm_as_segm(domain_final)

mesh <- fmesher::fm_mesh_2d(
  boundary = bnd,
  max.edge = c(2000, 5000),   # aristas interiores ≤2km, exteriores ≤5km
  offset   = c(2000, 6000),   # margen 2km interior, 6km exterior
  cutoff   = 1000,            # mínima distancia entre nodos = 1km
  crs      = sf::st_crs(9377) # CRS requerido por inlabru 2.10+ para lgcp()
)

cat("Malla SPDE:\n")
cat("  Nodos:", mesh$n, "\n")
cat("  Triángulos:", nrow(mesh$graph$tv), "\n")

# Figura 1: Malla + puntos de observación
# Convertir triangulación a sf para plotear con geom_sf (compatible con todas las versiones)
mesh_loc <- mesh$loc[, 1:2]
mesh_tv  <- mesh$graph$tv
mesh_tris <- st_sf(
  geometry = st_sfc(
    lapply(seq_len(nrow(mesh_tv)), function(i) {
      idx <- c(mesh_tv[i, ], mesh_tv[i, 1])
      st_polygon(list(mesh_loc[idx, ]))
    }),
    crs = 9377
  )
)

png(file.path(FIG, "01_mesh_spde.png"), width = 2400, height = 2000, res = 300)
ggplot() +
  hs_lyr(0.22) + hs_scale() +
  geom_sf(data = mesh_tris, fill = NA, color = "steelblue", lwd = 0.3, alpha = 0.4) +
  geom_sf(data = pts, size = 0.7, color = "#d62728", alpha = 0.6) +
  geom_sf(data = domain_final, fill = NA, color = "black", lwd = 0.6) +
  map_crs() +
  map_cart() +
  labs(x = "Longitud", y = "Latitud",
       caption = paste0("Rojo: deslizamientos observados  |  Azul: triangulación SPDE  |  Nodos: ",
                        mesh$n, "  Triángulos: ", nrow(mesh$graph$tv))) +
  theme_bw(base_size = 10) +
  theme(axis.text  = element_text(size = 7),
        plot.caption = element_text(size = 7))
dev.off()
cat("✓ Fig 01: malla SPDE\n")


# ─── 5. DEFINICIÓN DEL PROCESO SPDE ──────────────────────────────────────────
#
# El GRF ξ(s) se parametriza con la covarianza de Matérn:
#   Cov(ξ(s), ξ(s')) = (σ²/Γ(ν)2^(ν-1)) · (κ||s-s'||)^ν · K_ν(κ||s-s'||)
#
# donde:
#   κ = √(8ν) / ρ   (ρ = rango práctico: distancia donde la correlación ≈ 0.13)
#   σ² = varianza marginal
#   ν = 1 (suavidad del campo, α = ν + d/2 = 2 en 2D → SPDE de orden 2)
#
# Priors PC (Penalized Complexity) de Simpson et al. (2017):
#   P(ρ < ρ₀) = p_ρ   → prior que penaliza rangos pequeños (sobreajuste)
#   P(σ > σ₀) = p_σ   → prior que penaliza varianzas grandes
#
# Priors ajustados según diagnóstico de confusión espacial:
#   prior.range = c(2000, 0.5) → P(ρ < 2km) = 0.5  (prior centrado en 2km)
#   prior.sigma = c(2.0,  0.1) → P(σ > 2)   = 0.1  (penaliza varianzas grandes)
# Con priors permisivos (ρ₀=5km, p=0.01; σ₀=1.5, p=0.01) el GRF convergía
# a σ²≈933, absorbiendo todos los efectos de las covariables. Priors más
# informativos devuelven σ²≈0.83 y ρ≈0.6km, valores físicamente razonables.

spde <- inla.spde2.pcmatern(
  mesh  = mesh,
  alpha = 2,
  prior.range = c(5000, 0.01),  # P(ρ < 5km) = 0.01
  prior.sigma = c(2.0,  0.05)   # P(σ > 2)   = 0.05
)

cat("SPDE definido (Matérn, α=2, PC priors)\n")


# ─── 6. MUESTREO DE COVARIABLES EN PUNTOS Y GRILLA ───────────────────────────
#
# Para el LGCP con inlabru, las covariables se proveen como objetos SpatRaster.
# inlabru las muestrea automáticamente en los puntos de integración y en los
# puntos de observación.

# ── Calcular Índice Topográfico de Humedad (TWI) ──────────────────────────────
# TWI = ln(A / tan(β))
#   A  = área acumulada (m²), ya disponible como exp(cov_logacc)
#   β  = pendiente en radianes
# Pendiente mínima = 0.01° para evitar tan(0) → división por cero en valles planos
pts_coords  <- st_coordinates(pts)
slope_rad   <- cov_slope * (pi / 180)
slope_safe  <- ifel(slope_rad < (0.01 * pi / 180), 0.01 * pi / 180, slope_rad)
cov_twi     <- cov_logacc - log(tan(slope_safe))
cat("TWI calculado. Rango:",
    round(global(cov_twi, "min", na.rm=TRUE)[[1]], 2), "–",
    round(global(cov_twi, "max", na.rm=TRUE)[[1]], 2), "\n")

# Escalar covariables del modelo (media=0, sd=1)
# Covariables del modelo: pendiente, TWI, aspecto
sc <- function(r) {
  m <- global(r, "mean", na.rm = TRUE)[[1]]
  s <- global(r, "sd",   na.rm = TRUE)[[1]]
  (r - m) / s
}
cov_sc <- c(sc(cov_slope), sc(cov_twi), sc(cov_aspect))
names(cov_sc) <- c("slope_sc", "twi_sc", "aspect_sc")

# Verificar covariables en los puntos de deslizamiento
cov_at_pts <- as.data.frame(terra::extract(cov_sc, pts_coords))
names(cov_at_pts) <- names(cov_sc)
cat("Estadísticas de covariables escaladas en los puntos:\n")
print(summary(cov_at_pts))

cat("Covariables escaladas (media=0, sd=1)\n")

# ── Rellenar NAs en los rasters escalados ─────────────────────────────────────
# Los nodos del buffer exterior de la malla caen fuera del DEM y obtienen NA.
# Si inlabru evalúa el raster en esos nodos y recibe NA, los trata como 0,
# sesgando los coeficientes. Se propagan los valores de los vecinos más cercanos
# mediante focal iterativo hasta eliminar todos los NAs.
cat("Extendiendo y rellenando NAs en covariables...\n")
fill_na_focal <- function(r) {
  filled <- r
  iter   <- 0
  while (any(is.na(values(filled, na.rm = FALSE))) && iter < 100) {
    filled <- terra::focal(filled, w = 3, fun = mean,
                           na.policy = "only", na.rm = TRUE)
    iter <- iter + 1
  }
  filled
}

# Extender el extent del raster para cubrir TODOS los nodos de la malla,
# incluidos los del buffer exterior que caen fuera del DEM.
# Sin esto, terra::extract devuelve NA en esos nodos y el modelo colapsa.
mesh_ext <- terra::ext(
  min(mesh$loc[, 1]) - 500,
  max(mesh$loc[, 1]) + 500,
  min(mesh$loc[, 2]) - 500,
  max(mesh$loc[, 2]) + 500
)
cov_sc_ext <- terra::extend(cov_sc, mesh_ext)  # añade celdas NA en el margen

cov_sc_filled <- terra::rast(lapply(names(cov_sc_ext), function(nm) {
  fill_na_focal(cov_sc_ext[[nm]])
}))
names(cov_sc_filled) <- names(cov_sc)

# Verificar: ningún nodo debe quedar con NA
vals_check <- terra::extract(cov_sc_filled[["slope_sc"]], mesh$loc[, 1:2])[, 1]
cat("  NAs en nodos de la malla tras relleno:", sum(is.na(vals_check)),
    "de", length(vals_check), "\n")

# ── Pre-extracción de covariables en los puntos de deslizamiento ──────────────
pts$slope_sc  <- terra::extract(cov_sc_filled[["slope_sc"]],  pts_coords)[, 1]
pts$twi_sc    <- terra::extract(cov_sc_filled[["twi_sc"]],    pts_coords)[, 1]
pts$aspect_sc <- terra::extract(cov_sc_filled[["aspect_sc"]], pts_coords)[, 1]

cat("Covariables pre-extraídas en", nrow(pts), "puntos de deslizamiento\n")
cat("  NAs en puntos:", sum(is.na(pts$slope_sc)), "(slope),",
    sum(is.na(pts$twi_sc)), "(TWI),",
    sum(is.na(pts$aspect_sc)), "(aspecto)\n")


# ─── 7. ESPECIFICACIÓN DEL MODELO LGCP ───────────────────────────────────────
#
# Modelo:
#   log λ(s) = intercepto
#              + β_slope  · slope_sc(s)
#              + β_twi    · twi_sc(s)
#              + β_aspect · aspect_sc(s)
#              + ξ(s)                    ← campo espacial latente (GRF)
#
# Covariables:
#   slope_sc  : pendiente estandarizada (°, media=0, sd=1)
#   twi_sc    : Índice Topográfico de Humedad estandarizado
#               TWI = ln(area_acum / tan(pendiente_rad))
#               Valores altos → zonas de convergencia hídrica (fondos de quebrada)
#               Valores bajos → laderas convexas con drenaje rápido
#   aspect_sc : aspecto estandarizado (°, 0=N, 90=E, 180=S, 270=O)
#
# Se usa cov_sc_filled (sin NAs) para que inlabru evalúe correctamente las
# covariables en los nodos del buffer exterior (fuera del DEM).

components <- ~ Intercept(1) +
  slope_sc (cov_sc_filled[["slope_sc"]],  model = "linear") +
  twi_sc   (cov_sc_filled[["twi_sc"]],    model = "linear") +
  aspect_sc(cov_sc_filled[["aspect_sc"]], model = "linear") +
  field(geometry, model = spde)   # campo espacial latente ξ(s)

cat("Modelo LGCP especificado\n")
cat("  Efectos fijos: pendiente, TWI, aspecto\n")
cat("  Efecto aleatorio espacial: GRF Matérn vía SPDE\n")


# ─── 8. AJUSTE DEL MODELO ─────────────────────────────────────────────────────
#
# lgcp() en inlabru:
#   - data      : puntos de observación (proceso puntual)
#   - samplers  : dominio Ω de integración
#   - domain    : especifica la malla para la integración del proceso puntual
#
# Internamente, inlabru discretiza el dominio en puntos de integración usando
# la estrategia de cuadratura de Berman & Turner (1992), donde la intensidad
# integrada ∫Ω λ(s)ds ≈ Σᵢ wᵢ·λ(sᵢ) con pesos wᵢ proporcionales al área
# de cada celda de la teselación de Voronoi de la malla.
#
# INLA aproxima la posterior conjunta p(θ|y) mediante:
#   p(θ|y) ≈ Σₖ p(θ|xₖ,y) · Δxₖ   (suma sobre hiperparámetros)
# usando la aproximación de Laplace para cada término.

cat("Ajustando modelo LGCP... (puede tomar varios minutos)\n")

fit <- lgcp(
  components = components,
  data       = pts,
  formula    = geometry ~ .,   # inlabru 2.10+: especificar explícitamente la col. de geometría
  samplers   = domain_final,
  domain     = list(geometry = mesh),
  options    = list(
    control.fixed = list(
      mean = 0, prec = 0.01  # prior Normal(0, σ=10) en cada β
    ),
    control.mode = list(
      # Valores iniciales físicamente razonables para los hiperparámetros del GRF.
      # Evita que el optimizador converja al mínimo degenerado (ρ→0, σ²→∞).
      # theta[1] = log(rango_interno) ≈ log(5000) = 8.52
      # theta[2] = log(sigma_interno) ≈ log(1.5)  = 0.41
      theta   = c(8.52, 0.41),
      restart = FALSE
    ),
    control.inla = list(
      int.strategy = "ccd",    # integra sobre grilla de hiperparámetros (evita mínimos locales)
      strategy     = "laplace"
    ),
    verbose = TRUE
  )
)

cat("✓ Modelo ajustado\n")
print(summary(fit))


# ─── 9. RESULTADOS: EFECTOS FIJOS ─────────────────────────────────────────────
#
# Los coeficientes β son los efectos marginales en la escala log-intensidad.
# Un β_slope > 0 significa que a mayor pendiente, mayor intensidad de deslizamientos.
# Los intervalos de credibilidad al 95% provienen directamente de la posterior.

fixed_summary <- fit$summary.fixed
cat("\n=== EFECTOS FIJOS (escala log-intensidad) ===\n")
print(round(fixed_summary[, c("mean","sd","0.025quant","0.975quant","mode")], 4))

# Figura 2: Coeficientes con intervalos de credibilidad al 95%
coef_df <- data.frame(
  variable = rownames(fixed_summary),
  mean     = fixed_summary$mean,
  lo       = fixed_summary$`0.025quant`,
  hi       = fixed_summary$`0.975quant`
)
coef_df$variable <- factor(coef_df$variable,
  levels = c("Intercept","slope_sc","twi_sc","aspect_sc"),
  labels = c("Intercepto","Pendiente\n(estand.)","TWI\n(estand.)","Aspecto\n(estand.)")
)
# Color: positivo=azul, negativo=rojo, no significativo=gris
coef_df$sig <- ifelse(coef_df$lo > 0, "Positivo (sig.)",
               ifelse(coef_df$hi < 0, "Negativo (sig.)", "No significativo"))

png(file.path(FIG, "02_efectos_fijos.png"), width = 2400, height = 1600, res = 300)
ggplot(coef_df[-1,], aes(x = variable, y = mean, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.8, linewidth = 1.2) +
  scale_color_manual(values = c("Positivo (sig.)" = "#1f78b4",
                                "Negativo (sig.)" = "#e31a1c",
                                "No significativo" = "#636363"),
                     name = "Efecto") +
  labs(x = NULL, y = expression(hat(beta) ~ "(log-intensidad)")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()
cat("✓ Fig 02: efectos fijos\n")


# ─── 10. RESULTADOS: HIPERPARÁMETROS SPDE ─────────────────────────────────────
#
# El GRF tiene dos hiperparámetros:
#   rango ρ: distancia a la que la correlación espacial decae a ~0.13
#            Un rango pequeño → variación espacial muy local (rugoso)
#            Un rango grande  → variación suave y extendida
#
#   σ (desv. estándar marginal): magnitud de la variabilidad espacial residual
#            σ grande → el campo espacial domina sobre las covariables

spde_result <- inla.spde2.result(fit, "field", spde)  # fit hereda clase "inla" en inlabru 2.12

# Distribuciones marginales de rango y σ
post_range <- data.frame(
  x = inla.emarginal(function(x) x, spde_result$marginals.range.nominal[[1]]),
  marginal = spde_result$marginals.range.nominal[[1]]
)

range_stats <- inla.zmarginal(spde_result$marginals.range.nominal[[1]])
sigma_stats <- inla.zmarginal(spde_result$marginals.variance.nominal[[1]])

cat("\n=== HIPERPARÁMETROS SPDE ===\n")
cat("Rango (ρ):\n")
cat("  Media:", round(range_stats$mean / 1000, 2), "km\n")
cat("  IC 95%: [",
    round(range_stats$quant0.025 / 1000, 2), ",",
    round(range_stats$quant0.975 / 1000, 2), "] km\n")
cat("Varianza (σ²):\n")
cat("  Media:", round(sigma_stats$mean, 4), "\n")
cat("  IC 95%: [",
    round(sigma_stats$quant0.025, 4), ",",
    round(sigma_stats$quant0.975, 4), "]\n")

# Figura 3: Distribuciones posteriores de hiperparámetros
pdf_range <- as.data.frame(spde_result$marginals.range.nominal[[1]])
pdf_sigma <- as.data.frame(spde_result$marginals.variance.nominal[[1]])
pdf_range$x <- pdf_range$x / 1000  # m → km

p_range <- ggplot(pdf_range, aes(x, y)) +
  geom_area(fill = "#3182bd", alpha = 0.4) +
  geom_line(color = "#3182bd", lwd = 1) +
  geom_vline(xintercept = range_stats$mean/1000, lty=2, color="red") +
  labs(x = "ρ (km)", y = "Densidad posterior") +
  theme_bw(base_size = 10)

p_sigma <- ggplot(pdf_sigma, aes(x, y)) +
  geom_area(fill = "#e6550d", alpha = 0.4) +
  geom_line(color = "#e6550d", lwd = 1) +
  geom_vline(xintercept = sigma_stats$mean, lty=2, color="red") +
  labs(x = expression(sigma^2), y = "Densidad posterior") +
  theme_bw(base_size = 10)

png(file.path(FIG, "03_hiperparametros_spde.png"), width=2400, height=1200, res=300)
gridExtra::grid.arrange(p_range, p_sigma, ncol = 2)
dev.off()
cat("✓ Fig 03: hiperparámetros SPDE\n")


# ─── 11. PREDICCIÓN: INTENSIDAD ESPACIAL ─────────────────────────────────────
#
# Predecimos λ(s) = exp(η(s)) sobre una grilla regular que cubre el dominio.
# La predicción posterior tiene media E[λ(s)|y] y desviación estándar sd[λ(s)|y].
#
# Interpretación física:
#   λ(s) = número esperado de deslizamientos por unidad de área (m⁻²)
#   λ(s) × A_celda = número esperado de eventos en una celda de área A
#
# Usamos la función predict() de inlabru para obtener la distribución posterior
# completa en cada celda de la grilla.

# Crear grilla de predicción dentro del dominio (resolución 200m)
grid_pred <- st_make_grid(domain_final, cellsize = 200, what = "centers") |>
  st_as_sf() |>
  st_filter(domain_final)
grid_pred$geometry <- grid_pred$x
grid_pred <- st_set_geometry(grid_pred, "geometry")
st_crs(grid_pred) <- 9377

# Muestrear covariables escaladas en la grilla (desde raster sin NAs)
grid_coords <- st_coordinates(grid_pred)
cov_grid <- as.data.frame(terra::extract(cov_sc_filled, grid_coords))
names(cov_grid) <- names(cov_sc_filled)
grid_pred <- cbind(grid_pred, cov_grid)

cat("Grilla de predicción:", nrow(grid_pred), "celdas (200m)\n")
cat("Prediciendo intensidad posterior...\n")

# predict() estima la distribución posterior de log λ(s) en cada punto de la grilla
pred_intensity <- predict(
  fit,
  newdata = grid_pred,
  formula = ~ exp(Intercept + slope_sc + twi_sc + aspect_sc + field),
  n.samples = 1000   # muestras de la posterior para integración de Monte Carlo
)

grid_pred$lambda_mean <- pred_intensity$mean   # media posterior de λ(s)
grid_pred$lambda_sd   <- pred_intensity$sd     # incertidumbre en λ(s)
grid_pred$lambda_lo   <- pred_intensity$q0.025 # IC 95% inferior
grid_pred$lambda_hi   <- pred_intensity$q0.975 # IC 95% superior

cat("✓ Predicción completada\n")


# ─── 12. CAMPO ESPACIAL LATENTE ───────────────────────────────────────────────
#
# El campo ξ(s) captura la autocorrelación espacial residual:
# valores positivos → zonas con más deslizamientos de lo que predicen las covariables
# valores negativos → zonas con menos deslizamientos de lo esperado
# Este "exceso" puede reflejar factores no observados (litología, vegetación, etc.)

pred_field <- predict(
  fit,
  newdata = grid_pred,
  formula = ~ field,
  n.samples = 1000
)
grid_pred$field_mean <- pred_field$mean
grid_pred$field_sd   <- pred_field$sd


# ─── 13. MAPAS DE RESULTADOS ──────────────────────────────────────────────────

# Función auxiliar para mapas ggplot uniformes
map_theme <- function() {
  theme_bw(base_size = 9) +
  theme(
    legend.key.height = unit(0.8, "cm"),
    axis.text  = element_text(size = 7),
    axis.title = element_text(size = 8)
  )
}

# Figura 4: Intensidad predicha (media posterior)
p_lambda <- ggplot() +
  hs_lyr(0.22) + hs_scale() +
  geom_sf(data = grid_pred, aes(color = lambda_mean), size = 0.4, shape = 15) +
  scale_color_viridis_c(option = "plasma", name = "λ(s)\n(eventos/m²)",
                        trans = "log10", labels = scientific) +
  geom_sf(data = pts, size = 0.4, color = "white", alpha = 0.4) +
  map_crs() + map_cart() +
  labs(x = "Longitud", y = "Latitud") +
  map_theme()

# Figura 5: Incertidumbre (desviación estándar posterior)
p_sd <- ggplot() +
  hs_lyr(0.22) + hs_scale() +
  geom_sf(data = grid_pred, aes(color = lambda_sd), size = 0.4, shape = 15) +
  scale_color_viridis_c(option = "cividis", name = "sd[λ(s)|y]",
                        trans = "log10", labels = scientific) +
  map_crs() + map_cart() +
  labs(x = "Longitud", y = "Latitud") +
  map_theme()

# Figura 6: Campo espacial latente (ξ(s))
lim_field <- max(abs(grid_pred$field_mean), na.rm = TRUE)
p_field <- ggplot() +
  hs_lyr(0.22) + hs_scale() +
  geom_sf(data = grid_pred, aes(color = field_mean), size = 0.4, shape = 15) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#d73027",
                        midpoint = 0, limits = c(-lim_field, lim_field),
                        name = "ξ(s)") +
  geom_sf(data = pts, size = 0.4, color = "black", alpha = 0.3) +
  map_crs() + map_cart() +
  labs(x = "Longitud", y = "Latitud") +
  map_theme()

png(file.path(FIG, "04_mapa_intensidad_predicha.png"), width=2400, height=2000, res=300)
print(p_lambda)
dev.off()
cat("✓ Fig 04: intensidad predicha\n")

png(file.path(FIG, "05_mapa_incertidumbre.png"), width=2400, height=2000, res=300)
print(p_sd)
dev.off()
cat("✓ Fig 05: incertidumbre\n")

png(file.path(FIG, "06_mapa_campo_latente.png"), width=2400, height=2000, res=300)
print(p_field)
dev.off()
cat("✓ Fig 06: campo latente ξ(s)\n")

# Figura 7: Panel 3×1 — comparación espacial completa
png(file.path(FIG, "07_panel_mapas.png"), width=4800, height=2000, res=300)
gridExtra::grid.arrange(p_lambda, p_field, p_sd, ncol = 3)
dev.off()
cat("✓ Fig 07: panel de mapas\n")


# ─── 14. COMPARACIÓN REAL vs SIMULADO ────────────────────────────────────────
#
# Estrategia: dividir el dominio en celdas (hexbin) y comparar el conteo
# observado de deslizamientos con el conteo esperado bajo el modelo.
#
# Conteo esperado en celda i = ∫_{Aᵢ} λ(s)ds ≈ λ̄ᵢ × Aᵢ
# donde λ̄ᵢ es la intensidad media en la celda.
#
# Un buen ajuste → puntos cercanos a la diagonal en el scatter observado vs esperado.

# Crear hexgrid de ~500m sobre el dominio
hex_grid <- st_make_grid(domain_final, cellsize = 500, square = FALSE) |>
  st_as_sf() |>
  st_filter(domain_final)
hex_grid$cell_id <- seq_len(nrow(hex_grid))
hex_grid$area_m2 <- as.numeric(st_area(hex_grid))

# Conteo observado por celda
pts_in_hex <- st_join(pts, hex_grid, join = st_within)
obs_count  <- as.data.frame(table(pts_in_hex$cell_id))
names(obs_count) <- c("cell_id", "n_obs")
obs_count$cell_id <- as.integer(as.character(obs_count$cell_id))

# Intensidad predicha media por celda
grid_in_hex <- st_join(grid_pred, hex_grid, join = st_within)
pred_by_hex <- aggregate(lambda_mean ~ cell_id, data = grid_in_hex, FUN = mean)

hex_compare <- merge(hex_grid, obs_count,  by = "cell_id", all.x = TRUE)
hex_compare <- merge(hex_compare, pred_by_hex, by = "cell_id", all.x = TRUE)
hex_compare$n_obs[is.na(hex_compare$n_obs)] <- 0

# Conteo esperado = λ̄ × área de la celda
hex_compare$n_exp <- hex_compare$lambda_mean * hex_compare$area_m2
hex_compare$n_exp[is.na(hex_compare$n_exp)] <- 0

# Residuos de Pearson: (obs - esp) / sqrt(esp)
hex_compare$pearson_res <- with(hex_compare,
  ifelse(n_exp > 0, (n_obs - n_exp) / sqrt(n_exp), NA)
)

# Figura 8: Scatter obs vs predicho
max_n <- max(c(hex_compare$n_obs, hex_compare$n_exp), na.rm = TRUE)
png(file.path(FIG, "08_obs_vs_predicho.png"), width=2000, height=2000, res=300)
ggplot(hex_compare, aes(x = n_exp, y = n_obs)) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2, lwd = 0.8) +
  geom_point(alpha = 0.5, size = 1.2, color = "#3182bd") +
  geom_smooth(method = "lm", color = "#e6550d", se = TRUE, lwd = 0.8) +
  scale_x_continuous(limits = c(0, max_n)) +
  scale_y_continuous(limits = c(0, max_n)) +
  labs(
    x = "Conteo esperado bajo el modelo E[N(A)|y]",
    y = "Conteo observado N(A)"
  ) +
  annotate("text", x = max_n*0.7, y = max_n*0.1,
           label = paste0("r² = ", round(cor(hex_compare$n_obs,
                                              hex_compare$n_exp, use="complete.obs")^2, 3)),
           size = 3.5, color = "#e6550d") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold"))
dev.off()
cat("✓ Fig 08: observado vs predicho\n")

# Figura 9: Mapa de residuos de Pearson
png(file.path(FIG, "09_mapa_residuos.png"), width=2400, height=2000, res=300)
lim_res <- quantile(abs(hex_compare$pearson_res), 0.99, na.rm=TRUE)
ggplot() +
  hs_lyr(0.22) + hs_scale() +
  ggnewscale::new_scale_fill() +
  geom_sf(data = hex_compare, aes(fill = pearson_res), color = NA, alpha = 0.80) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#d73027",
                       midpoint=0, limits=c(-lim_res, lim_res),
                       name="Residuo\nPearson", na.value="gray90") +
  geom_sf(data=pts, size=0.4, color="black", alpha=0.3) +
  map_crs() + map_cart() +
  labs(x = "Longitud", y = "Latitud") +
  map_theme()
dev.off()
cat("✓ Fig 09: mapa de residuos\n")

# Figura 10: Mapa de conteo observado vs esperado (doble panel)
p_obs <- ggplot() +
  hs_lyr(0.22) + hs_scale() +
  ggnewscale::new_scale_fill() +
  geom_sf(data = hex_compare, aes(fill = n_obs), color = NA, alpha = 0.80) +
  scale_fill_viridis_c(option="plasma", name="N observado") +
  geom_sf(data=pts, size=0.3, color="white", alpha=0.4) +
  map_crs() + map_cart(loc_s="bl", loc_n="tl") +
  labs(x = "Longitud", y = "Latitud") + map_theme()

p_exp <- ggplot() +
  hs_lyr(0.22) + hs_scale() +
  ggnewscale::new_scale_fill() +
  geom_sf(data = hex_compare, aes(fill = n_exp), color = NA, alpha = 0.80) +
  scale_fill_viridis_c(option="plasma", name="N esperado") +
  map_crs() + map_cart() +
  labs(x = "Longitud", y = "Latitud") + map_theme()

png(file.path(FIG, "10_obs_vs_predicho_mapas.png"), width=4000, height=2000, res=300)
gridExtra::grid.arrange(p_obs, p_exp, ncol = 2)
dev.off()
cat("✓ Fig 10: mapas comparativos\n")


# ─── 15. DIAGNÓSTICO ESTADÍSTICO ──────────────────────────────────────────────
#
# DIC (Deviance Information Criterion): análogo bayesiano del AIC.
#   DIC = -2·log p(y|θ̄) + 2·p_D   (p_D = nº efectivo de parámetros)
#   Menor DIC → mejor ajuste relativo.
#
# WAIC (Widely Applicable Information Criterion): alternativa más robusta al DIC.
#
# CPO (Conditional Predictive Ordinate): probabilidad de observar el dato i
#   dado el resto del conjunto de datos. CPO bajo → observación influyente.
#   LCPO = -Σlog(CPO_i): análogo a la log-verosimilitud Leave-One-Out.

cat("\n=== DIAGNÓSTICO ESTADÍSTICO ===\n")
cat("DIC: ", round(fit$dic$dic,  2), "\n")
cat("p_D: ", round(fit$dic$p.eff, 2), "(nº efectivo de parámetros)\n")
cat("WAIC:", round(fit$waic$waic, 2), "\n")

# Log CPO (Leave-One-Out)
lcpo <- -mean(log(fit$cpo$cpo), na.rm = TRUE)
cat("LCPO (negativo → mejor):", round(lcpo, 4), "\n")
cat("Fallos CPO:", sum(fit$cpo$failure > 0, na.rm=TRUE), "de", nrow(pts), "\n")

# Figura 11: Distribución del CPO
png(file.path(FIG, "11_diagnostico_cpo.png"), width=2400, height=1600, res=300)
cpo_df <- data.frame(cpo = fit$cpo$cpo, pit = fit$cpo$pit)
p_cpo <- ggplot(cpo_df, aes(x=log(cpo))) +
  geom_histogram(fill="#3182bd", color="white", bins=40) +
  labs(x="log(CPO)", y="Frecuencia") +
  theme_bw(base_size=10)

p_pit <- ggplot(cpo_df, aes(x=pit)) +
  geom_histogram(fill="#e6550d", color="white", bins=20, boundary=0) +
  geom_hline(yintercept=nrow(pts)/20, lty=2, color="gray50") +
  labs(x="PIT", y="Frecuencia") +
  theme_bw(base_size=10)
gridExtra::grid.arrange(p_cpo, p_pit, ncol = 2)
dev.off()
cat("✓ Fig 11: diagnóstico CPO/PIT\n")

# Figura 12: Curva ROC espacial (ROC basada en intensidad predicha)
# Convertir el problema a clasificación binaria: celda con al menos 1 deslizamiento
library(pROC)
roc_data <- hex_compare[!is.na(hex_compare$n_exp), ]
roc_obj  <- roc(as.numeric(roc_data$n_obs > 0), roc_data$n_exp)
auc_val  <- auc(roc_obj)

png(file.path(FIG, "12_curva_roc.png"), width=1800, height=1800, res=300)
plot(roc_obj, col="#3182bd", lwd=2)
text(0.3, 0.1, paste0("AUC = ", round(auc_val, 3)), cex=1.2, col="#3182bd")
abline(0, 1, lty=2, col="gray")
dev.off()
cat("✓ Fig 12: curva ROC  AUC =", round(auc_val, 3), "\n")


# ─── 16. RESUMEN FINAL ────────────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("RESUMEN DEL MODELO LGCP — DESLIZAMIENTOS VALLE DE ABURRÁ\n")
cat(strrep("=", 60), "\n\n")

cat("Eventos observados:", nrow(pts), "\n")
cat("Área del dominio:", round(as.numeric(st_area(domain_final))/1e6, 2), "km²\n")
cat("Intensidad media: ~",
    round(nrow(pts) / as.numeric(st_area(domain_final)) * 1e6, 4),
    "eventos/km²\n\n")

cat("--- Efectos fijos (β en escala log-intensidad) ---\n")
for (nm in rownames(fixed_summary)) {
  cat(sprintf("  %-20s  media=%6.3f  IC95=[%6.3f, %6.3f]\n",
              nm,
              fixed_summary[nm, "mean"],
              fixed_summary[nm, "0.025quant"],
              fixed_summary[nm, "0.975quant"]))
}

cat("\n--- Hiperparámetros SPDE ---\n")
cat(sprintf("  Rango ρ:    %.1f km  IC95=[%.1f, %.1f] km\n",
            range_stats$mean/1000,
            range_stats$quant0.025/1000,
            range_stats$quant0.975/1000))
cat(sprintf("  Varianza σ²: %.4f   IC95=[%.4f, %.4f]\n",
            sigma_stats$mean,
            sigma_stats$quant0.025,
            sigma_stats$quant0.975))

cat("\n--- Ajuste del modelo ---\n")
cat("  DIC:", round(fit$dic$dic, 2), "\n")
cat("  WAIC:", round(fit$waic$waic, 2), "\n")
cat("  LCPO:", round(lcpo, 4), "\n")
cat("  AUC (ROC):", round(auc_val, 3), "\n")

cat("\nFiguras generadas en:", FIG, "\n")
cat(strrep("=", 60), "\n")

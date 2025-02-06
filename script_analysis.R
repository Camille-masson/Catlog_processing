### 0. LIBRARIES AND CONSTANTS ###
#--------------------------------#

# Nettoyage de la mémoire
gc()

# Chargement de la configuration
source("config.R")

# Définition de l'année d'analyse
YEAR <- 2024
ALPAGES_TOTAL <- list(
  "2022" = c("Ane-et-Buyant", "Cayolle", "Combe-Madame", "Grande-Fesse", "Jas-des-Lievres", "Lanchatra", "Pelvas", "Sanguiniere", "Viso"),
  "2023" = c("Cayolle", "Crouzet", "Grande-Cabane", "Lanchatra", "Rouanette", "Sanguiniere", "Vacherie-de-Roubion", "Viso"),
  "2024" = c("Viso", "Cayolle", "Sanguiniere")
)
ALPAGES <- ALPAGES_TOTAL[[as.character(YEAR)]]

### 1. GPS TRAJECTORIES FILTERING ###
#-----------------------------------#

if (FALSE) {  # Mettre TRUE pour exécuter
  library(terra)
  source(file.path(functions_dir, "Functions_filtering.R"))
  
  raw_data_dir <- file.path(data_dir, paste0("Colliers_", YEAR, "_brutes"))
  alpages <- c("Sanguiniere", "Cayolle")
  output_file <- file.path(data_dir, paste0("Donnees_brutes_", YEAR, "_simplifiees.gpkg"))
  
  lapply(alpages, function(alpage) {
    collar_dir <- file.path(raw_data_dir, alpage)
    collar_files <- list.files(collar_dir, full.names = TRUE)
    lapply(collar_files, function(collar_f) {
      collar_ID <- strsplit(basename(collar_f), split = "_")[[1]][1]
      load_catlog_data(collar_f) %>%
        slice(which(row_number() %% 30 == 10)) %>%
        mutate(ID = collar_ID, date = lubridate::format_ISO8601(date)) %>%
        vect(geom = c("lon", "lat"), crs = CRS_WSG84)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>%
    writeVector(filename = output_file, overwrite = TRUE)
}

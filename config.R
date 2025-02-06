### CONFIGURATION DU PROJET ###

# Définition des chemins 
root_dir <- getwd()  # Récupère le chemin du projet automatiquement
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "outputs")
raster_dir <- file.path(root_dir, "raster")
functions_dir <- file.path(root_dir, "Functions")

# Chargement des librairies
library(tidyverse) # inclut ggplot2, dplyr, etc.
library(lubridate)
library(sf)
library(sp)
library(terra)
library(viridisLite)

# Chargement automatique des fonctions
# Charger automatiquement les fichiers de fonction
function_files <- list.files(functions_dir, pattern = "\\.R$", full.names = TRUE)

# Vérifier si des fichiers existent avant de les sourcer
if (length(function_files) > 0) {
  for (file in function_files) {
    print(paste("Sourcing:", file))  # Debugging: Affiche le fichier en cours de chargement
    tryCatch(
      source(file),
      error = function(e) {
        message(paste("⚠️ Erreur lors du chargement de :", file))
        message(e)
      }
    )
  }
} else {
  warning("⚠ Aucun fichier de fonction trouvé dans Functions/")
}
sapply(function_files, source)

# Paramètres globaux
YEAR <- 2024
ncores <- parallel::detectCores() /3 # Utilisation optimale des cœurs CPU


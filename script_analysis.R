### 0. LIBRARIES AND CONSTANTS ###
#--------------------------------#

# Chargement de la configuration
source("config.R")

# Définition de l'année d'analyse
YEAR <- 2022
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
  alpages <- c("Viso")
  output_file <- file.path(output_dir, paste0("Donnees_brutes_", YEAR, "_simplifiees.gpkg"))
  
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




### 1.2 BJONERAAS FILTER CALIBRATION ###
#--------------------------------------#

if (TRUE) {  # Remplace (F) par (TRUE) si tu veux exécuter directement
  
  # Chargement des fonctions et librairies nécessaires
  source(file.path(functions_dir, "Functions_filtering.R"))
  source(file.path(functions_dir, "Functions_map_plot.R"))
  library(dplyr)
  library(terra)
  
  # Définition des chemins
  raw_data_dir <- file.path(data_dir, paste0("Colliers_", YEAR, "_brutes"))
  AIF <- file.path(data_dir, paste0(YEAR, "_infos_alpages.csv"))
  
  # Chargement des informations sur les alpages
  AIF_data <- read.csv(AIF, sep = ",", header = TRUE, encoding = "UTF-8")
  
  # L’alpage à traiter
  alpage <- "Viso"
  
  # Création du fichier PDF pour enregistrer les graphiques
  output_pdf <- file.path(output_dir, "Filtering_calibration.pdf")
  pdf(output_pdf, width = 9, height = 9)
  
  # Liste des fichiers de l'alpage concerné
  collar_dir <- file.path(raw_data_dir, alpage)
  files <- list.files(collar_dir, full.names = TRUE)
  
  # Sélection des 3 premiers fichiers pour le test (évite erreur si < 3 fichiers)
  files <- files[1:min(3, length(files))]
  
  # Vérification qu'il y a bien des fichiers à traiter
  if (length(files) == 0) {
    stop("Aucun fichier GPS trouvé pour l'alpage: ", alpage)
  }
  
  # Chargement des données GPS
  data <- do.call(rbind, lapply(files, function(file) {
    data <- load_catlog_data(file)
    data$ID <- file
    return(data)
  }))
  
  # Récupération des dates de suivi
  beg_date <- as.POSIXct(get_alpage_info(alpage, AIF_data, "date_pose"), tz = "GMT", format = "%d/%m/%Y %H:%M")
  end_date <- as.POSIXct(get_alpage_info(alpage, AIF_data, "date_retrait"), tz = "GMT", format = "%d/%m/%Y %H:%M")
  
  # Filtrage des données
  data <- date_filter(data, beg_date, end_date)
  
  # Vérification qu'il reste des données après filtrage
  if (nrow(data) == 0) {
    stop("Erreur : Aucune donnée après filtrage pour l'alpage ", alpage)
  }
  
  # Projection en Lambert 93
  CRS_L93 <- "EPSG:2154"
  data_xy <- data %>%
    terra::vect(crs = "EPSG:4326") %>%
    terra::project(CRS_L93) %>%
    as.data.frame(geom = "XY")
  
  # Vérification si des données existent après transformation
  if (nrow(data_xy) == 0) {
    stop("Erreur : Aucune donnée après projection pour l'alpage ", alpage)
  }
  
  # Histogramme des périodes d’échantillonnage
  temps <- diff(data_xy$date)
  temps <- as.numeric(temps, units = "mins")
  hist(temps, nclass = 30, main = "Histogramme des périodes d’échantillonnage", xlab = "Temps (min)")
  
  # Histogramme des distances parcourues
  dist <- sqrt(diff(data_xy$x)^2 + diff(data_xy$y)^2)
  hist(dist, nclass = 30, main = "Histogramme des distances parcourues", xlab = "Distance (m)", xaxt = "n")
  
  ### Application du filtre avec plusieurs paramètres ###
  medcrits <- c(750, 500, 750)
  meancrits <- c(500, 500, 350)
  spikesps <- c(1500, 1500, 1500)
  spikecoss <- c(-0.95, -0.95, -0.95)
  
  for (i in seq_along(medcrits)) {
    trajectories <- position_filter(
      data,
      medcrit = medcrits[i],
      meancrit = meancrits[i],
      spikesp = spikesps[i],
      spikecos = spikecoss[i]
    )
    
    minmax_xy <- get_minmax_L93(trajectories[!(trajectories$R1error | trajectories$R2error), ], buffer = 100)
    
    trajectories$errors <- 1
    trajectories$errors[trajectories$R1error] <- 2
    trajectories$errors[trajectories$R2error] <- 3
    
    pal <- c("#56B4E9", "red", "black")
    
    print(
      ggplot(trajectories, aes(x, y, col = errors)) +
        geom_path(size = 0.2) +
        geom_point(size = 0.3) +
        coord_equal() +
        xlim(minmax_xy$x_min, minmax_xy$x_max) +
        ylim(minmax_xy$y_min, minmax_xy$y_max) +
        ggtitle(paste0(
          "medcrit = ", medcrits[i], ", meancrit = ", meancrits[i], 
          ", spikesp = ", spikesps[i], ", spikecos = ", spikecoss[i]
        )) +
        scale_colour_gradientn(colors = pal, guide = "legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error"))
    )
    
    # Zoom sur les 5 premiers jours
    trajectories <- trajectories %>% filter(date < beg_date + 3600 * 24 * 5)
    
    print(
      ggplot(trajectories, aes(x, y, col = errors)) +
        geom_path(size = 0.2) +
        geom_point(size = 0.3) +
        coord_equal() +
        ggtitle(paste0(
          "5 jours seulement, medcrit = ", medcrits[i], 
          ", meancrit = ", meancrits[i], ", spikesp = ", spikesps[i], 
          ", spikecos = ", spikecoss[i]
        )) +
        scale_colour_gradientn(colors = pal, guide = "legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error"))
    )
  }
  
  # Fermeture du fichier PDF
  dev.off()
  
  print("Calibration du filtre Bjorneraas terminée.")
}


### 1.3 FILTERING CATLOG DATA ###
#--------------------------------------#


# Chargement des fonctions nécessaires
source(file.path(functions_dir, "Functions_filtering.R"))

# ENTREES
# Un dossier contenant les trajectoires brutes (fichiers CSV issus des colliers Catlog), rangées dans des sous-dossiers au nom des alpages
raw_data_dir <- file.path(data_dir, paste0("Colliers_", YEAR, "_brutes"))

# Un fichier contenant les informations sur chaque individu équipé (dates de pose et de retrait, paramètres des colliers)
IIF <- file.path(raw_data_dir, paste0(YEAR, "_colliers_poses.csv"))


IFF <- read.csv(IIF)

# Les alpages à traiter
alpages <- c("Viso")

# SORTIES
# Un fichier .RDS contenant les trajectoires filtrées, avec coordonnées en Lambert93
output_rds_file <- file.path(output_dir, paste0("Catlog_", YEAR, "_filtered.rds"))

# Un fichier .csv contenant les performances des colliers (pourcentage de points éliminés, colliers défectueux...)
indicator_file <- file.path(output_dir, paste0("catlog_", YEAR, "_filtering.csv"))

# Boucle sur chaque alpage sélectionné
for (alpage in alpages) {
  print(paste("WORKING ON ALPAGE", alpage))
  
  # Définition du dossier contenant les données brutes de l'alpage
  collar_dir <- file.path(raw_data_dir, alpage)
  
  # Liste des fichiers CSV des trajectoires
  collar_files <- list.files(collar_dir, pattern = ".csv", full.names = TRUE)
  
  # Récupération des paramètres de filtrage depuis le fichier AIF
  medcrit <- get_alpage_info(alpage, AIF, "medcrit")
  meancrit <- get_alpage_info(alpage, AIF, "meancrit")
  spikesp <- get_alpage_info(alpage, AIF, "spikesp")
  spikecos <- as.numeric(gsub(",", ".", get_alpage_info(alpage, AIF, "spikecos")))
  
  print(paste0("Bjorneraas filter parameters: medcrit=", medcrit, 
               ", meancrit=", meancrit, ", spikesp=", spikesp, 
               ", spikecos=", spikecos))
  
  # Filtrage des trajectoires et calcul des indicateurs de performance
  indicators <- lapply(collar_files, function(collar) {
    filter_one_collar(load_catlog_data(collar),
                      collar, output_rds_file, alpage, beg_date, end_date, IIF,
                      bjoneraas.medcrit = medcrit,
                      bjoneraas.meancrit = meancrit,
                      bjoneraas.spikesp = spikesp,
                      bjoneraas.spikecos = spikecos)
  }) %>%
    do.call(rbind, .)
  
  # Agrégation des performances des colliers ayant fonctionné correctement
  indicators_tot <- indicators %>%
    filter(worked_until_end == 1) %>% # On exclut les colliers défectueux pour le calcul global
    add_row(name = paste("TOTAL", alpage), 
            worked_until_end = sum(.$worked_until_end), 
            nloc = NA,
            R1error = NA, 
            R2error = NA,
            error_perc = sum(.$nloc * .$error_perc) / sum(.$nloc), 
            localisation_rate = mean(.$localisation_rate))
  
  # Ajout des statistiques globales aux indicateurs
  indicators <- rbind(indicators, indicators_tot[nrow(indicators_tot),])
  
  # Écriture des résultats dans le fichier CSV des performances des colliers
  write.table(indicators, file = indicator_file, append = TRUE, sep = ',', row.names = FALSE, col.names = FALSE)
}

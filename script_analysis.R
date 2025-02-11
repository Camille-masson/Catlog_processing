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


### 1.2 BJONERAAS FILTER CALIBRATION
if (F) {
  # Chargement des fonctions nécessaires
  source(file.path(functions_dir, "Functions_filtering.R"))
  source(file.path(functions_dir, "Functions_map_plot.R"))
  # ENTREES
  # Un dossier contenant les trajectoires brutes, au format csv issu des colliers catlog, rangées dans des sous-dossiers au nom de leurs alpages
  raw_data_dir = file.path(data_dir,paste0("Colliers_",YEAR,"_brutes"))
  # Un data.frame contenant les dates de pose et de retrait des colliers, Doit contenir les colonnes  "alpage", "date_pose" et "date_retrait"
    AIF <- file.path(data_dir, paste0(YEAR,"_infos_alpages.csv"))
  #Vérification du format du fichier (souvent mal formaté, attention csv UTF8)
  AIF_data <- read.csv(AIF, sep = ",", header = TRUE, row.names = NULL, check.names = FALSE, encoding = "UTF-8")
  
  # L’alpage devant être traité
  alpage = "Viso"
  
  #Fichier pdf de sortie
  pdf(file.path(output_dir,paste0("Filtering_calibration_",YEAR,"_",alpage,".pdf")), width = 9, height = 9)
  
  # Liste des fichiers de données brutes pour l'alpage
  files <- list.files(file.path(raw_data_dir, alpage), full.names = TRUE)
  files <- files[1:3]  # Sélection des trois premiers fichiers
  
  # Chargement et concaténation des données
  data <- do.call(rbind, lapply(files, function(file) { 
    data <- load_catlog_data(file)
    data$ID <- file
    return(data) 
  }))
  
  # Récupération des dates de pose et de retrait du collier
  beg_date = as.POSIXct(get_alpage_info(alpage, AIF, "date_pose"), tz="GMT", format="%d/%m/%Y %H:%M")
  end_date = as.POSIXct(get_alpage_info(alpage, AIF, "date_retrait"), tz="GMT", format="%d/%m/%Y %H:%M")
  data = date_filter(data, beg_date, end_date)
  
  # Projection des données en Lambert93
  data_xy <- data %>%
    terra::vect(crs="EPSG:4326") %>%
    terra::project("EPSG:2154") %>%
    as.data.frame(geom = "XY") %>%
    head()
  
  
  # Histogramme des périodes d’échantillonnage
  temps <- diff(data_xy$date)
  temps <- as.numeric(temps, units = "mins")
  hist(temps, nclass = 30)
  
  # Histogramme des distances parcourues
  dist <- sqrt(diff(data_xy$x)^2+diff(data_xy$y)^2)
  h <- hist(dist, nclass = 30, xlab='Distance (m)', xaxt="n")
  
  ### Test de différents filtres
  # Parameters sets to be tested
  medcrits = c(750, 500, 750)
  meancrits = c(500, 500, 350)
  spikesps = c(1500, 1500, 1500)
  spikecoss = c(-0.95, -0.95, -0.95)
  
  for (i in 1:length(medcrits)) {
    trajectories <- position_filter(data, medcrit=medcrits[i], meancrit=meancrits[i], spikesp=spikesps[i], spikecos=spikecoss[i])
    
    minmax_xy = get_minmax_L93(trajectories[!(trajectories$R1error | trajectories$R2error ),], buffer = 100)
    
    trajectories$errors = 1
    trajectories$errors[trajectories$R1error] = 2
    trajectories$errors[trajectories$R2error] = 3
    
    pal <- c("#56B4E9", "red", "black")
    print(ggplot(trajectories, aes(x, y, col = errors)) +
            geom_path(size = 0.2) +
            geom_point(size = 0.3) +
            coord_equal() +
            xlim(minmax_xy$x_min, minmax_xy$x_max) + ylim(minmax_xy$y_min, minmax_xy$y_max) +
            ggtitle(paste0("medcrit = ", medcrits[i], ", meancrit = ", meancrits[i], ", spikesp = ", spikesps[i], ", spikecos = ", spikecoss[i])) +
            scale_colour_gradientn(colors=pal, guide="legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error")))
    
    trajectories <- trajectories %>%
      filter(date < beg_date + 3600*24*5)
    print(ggplot(trajectories, aes(x, y, col = errors)) +
            geom_path(size = 0.2) +
            geom_point(size = 0.3) +
            coord_equal() +
            ggtitle(paste0("5 days only, medcrit = ", medcrits[i], ", meancrit = ", meancrits[i], ", spikesp = ", spikesps[i], ", spikecos = ", spikecoss[i])) +
            scale_colour_gradientn(colors=pal, guide="legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error")))
  }
  dev.off()
}

##OK NORMALEMENT


#En travaux 

### 1.3 FILTERING CATLOG DATA
if (F) {
  source(file.path(functions_dir, "Functions_filtering.R"))
  # ENTREES
  # Un dossier contenant les trajectoires brutes, au format csv issu des colliers catlog, rangées dans des sous-dossiers au nom de leurs alpages. Coordonnées en WSG84. Le nom des fichiers, sous la forme "ID_quelquechose.csv", sera utilisé pour déterminer l’ID du collier qui doit comporter 3 caractères.
  raw_data_dir = file.path(data_dir,paste0("Colliers_",YEAR,"_brutes"))
  # Un fichier contenant les informations sur chaque individu équipé, les dates de pose et de retrait des colliers, ainsi que la proportion de temps pour laquelle les colliers sont programmés pour être allumés (18h par jour = 0.75). Doit contenir les colonnes "Collier", "Alpage", "Espece", "Race", "date_pose", "date_retrait" et "proportion_jour_allume"
  IIF = file.path(raw_data_dir, paste0(YEAR,"_colliers_poses.csv"))
  
  #Load and vérife data collier pose (format)
  IFF_data <- read.csv(IIF, stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Les alpages à traiter
  alpages = c("Viso")
  # SORTIES
  # Un .RDS contenant les trajectoires filtrées (les nouvelles trajectoires sont ajoutées à la suite des trajectoires traitées précédemment). Coordonnées en Lambert93.
  output_rds_file = file.path(output_dir, paste0("Catlog_",YEAR,"_filtered_",alpages,".rds"))
  # Un .csv contenant les performances des colliers (pourcentages de points éliminés à chaque étape, colliers défectueux...)
  indicator_file = file.path(output_dir, paste0(YEAR,"_filtering_",alpages,".csv"))
  
  
  for (alpage in alpages) {
    print(paste("WORKING ON ALPAGE", alpage))
    collar_dir <- file.path(raw_data_dir, alpage)
    collar_files <- list.files(collar_dir, pattern = ".csv", full.names = TRUE)
    
    medcrit = get_alpage_info(alpage, AIF, "medcrit")
    meancrit = get_alpage_info(alpage, AIF, "meancrit")
    spikesp = get_alpage_info(alpage, AIF, "spikesp")
    spikecos = as.numeric(gsub(",", ".", get_alpage_info(alpage, AIF, "spikecos")))
    print(paste0("Bjorneraas filter parameters: medcrit=",medcrit,", meancrit=", meancrit, ", spikesp=", spikesp, ", spikecos=", spikecos))
    print(collar_files)
    
    
    # Filtrage des trajectoires et calcul des indicateurs
    indicators <- lapply(collar_files, function(collar) {
      filter_one_collar(
        load_catlog_data(collar),  
        basename(collar),  # On passe uniquement le nom du fichier
        output_rds_file, alpage, beg_date, end_date, IIF,
        bjoneraas.medcrit = medcrit,
        bjoneraas.meancrit = meancrit,
        bjoneraas.spikesp = spikesp,
        bjoneraas.spikecos = spikecos
      )
    }) %>%
      do.call(rbind, .)
    
    indicators_tot = indicators %>%
      filter(worked_until_end == 1) %>% # to compute performance indicators at the alpage level, we remove defective collars
      add_row(name = paste("TOTAL", alpage), worked_until_end = sum(.$worked_until_end), nloc = NA,
              R1error = NA, R2error = NA,
              error_perc = sum(.$nloc*.$error_perc)/sum(.$nloc), localisation_rate = mean(.$localisation_rate))
    indicators = rbind(indicators, indicators_tot[nrow(indicators_tot),])
    
    write.table(indicators, file=indicator_file, append = T, sep=',', row.names=F, col.names=F)
  }
}

str(IFF_data)

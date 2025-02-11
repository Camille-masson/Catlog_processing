### 0. LIBRARIES AND CONSTANTS ###
#--------------------------------#

# Chargement de la configuration
source("config.R")

# Définition de l'année d'analyse
YEAR <- 0000
ALPAGES_TOTAL <- list(
  "0000" = c("Alpage_demo"),
  "2022" = c("Ane-et-Buyant", "Cayolle", "Combe-Madame", "Grande-Fesse", "Jas-des-Lievres", "Lanchatra", "Pelvas", "Sanguiniere", "Viso"),
  "2023" = c("Cayolle", "Crouzet", "Grande-Cabane", "Lanchatra", "Rouanette", "Sanguiniere", "Vacherie-de-Roubion", "Viso"),
  "2024" = c("Viso", "Cayolle", "Sanguiniere")
)
ALPAGES <- ALPAGES_TOTAL[[as.character(YEAR)]]

### 1. GPS TRAJECTORIES FILTERING ###
#-----------------------------------#


### 1.1 Simplification en GPKG


if (FALSE) {  # Mettre TRUE pour exécuter
  library(terra)
  source(file.path(functions_dir, "Functions_filtering.R"))
  
  raw_data_dir <- file.path(data_dir, paste0("Colliers_", YEAR, "_brutes"))
  alpages <- c("Alpage_demo")
  output_file <- file.path(output_dir, paste0("Donnees_brutes_", YEAR,"_",alpages,"_simplifiees.gpkg"))
  
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


### 3. HMM FITTING ### 
#--------------------#
if (F) {
  library(snow)
  library(stats)
  # Movement modelling packages
  library(momentuHMM)
  # library(foieGras)
  library(adehabitatLT)
  library(adehabitatHR)
  # Libraries RMarkdown
  library(knitr)
  library(rmarkdown)
  source(file.path(functions_dir, "Functions_HMM_fitting.R"))
  
  # ENTREES
  # Un .RDS contenant les trajectoires filtrées
  input_rds_file = file.path(output_dir, paste0("Catlog_",YEAR,"_filtered_",alpage,".rds"))
  
  # Un data.frame contenant la correspondance entre colliers et alpages. Doit contenir les colonnes  "ID", "Alpage" et "Periode d’echantillonnage"
  individual_info_file <- file.path(data_dir, paste0("Colliers_", YEAR, "_brutes"), paste0(YEAR, "_colliers_poses.csv"))
  
  # Les alpages à traiter
  alpages = ALPAGES
  
  # SORTIES
  # Un .RDS contenant les trajectoires catégorisées par comportement (les nouvelles trajectoires sont ajoutées à la suite des trajectoires traitées précédemment)
  output_rds_file = file.path(output_dir, paste0("Catlog_",YEAR,"_",alpage,"_viterbi.rds"))
  
  ### LOADING DATA FOR ANALYSES
  data = readRDS(input_rds_file)
  data = data[data$species == "brebis",]
  data = data[data$alpage %in% alpages,]
  
  ### HMM FIT
  run_parameters = list(
    # Model
    model = "HMM",
    
    # Resampling
    resampling_ratio = 5,
    resampling_first_index = 0,
    rollavg = FALSE,
    rollavg_convolution = c(0.15, 0.7, 0.15),
    knownRestingStates = FALSE,
    
    # Observation distributions (step lengths and turning angles)
    dist = list(step = "gamma", angle = "vm"),
    # Design matrices to be used for the probability distribution parameters of each data stream
    DM = list(angle=list(mean = ~1, concentration = ~1)),
    # Covariants formula
    covariants = ~cos(hour*3.141593/12), # ~1 if no covariants used
    
    # 3-state HMM
    Par0 = list(step = c(10, 25, 50, 10, 15, 40), angle = c(tan(pi/2), tan(0/2), tan(0/2), log(0.5), log(0.5), log(3))),
    fixPar = list(angle = c(tan(pi/2), tan(0/2), tan(0/2), NA, NA, NA))
  )
  run_parameters = scale_step_parameters_to_resampling_ratio(run_parameters)
  
  startTime = Sys.time()
  results = par_HMM_fit(data, run_parameters, ncores = ncores, individual_info_file, sampling_period = 120, output_dir)
  endTime = Sys.time()
  
  ### SUMMARIZE MODEL FITTING BY ALPAGE in rmarkdown PDFs (A revoir)
  parameters_df <- parameters_to_data.frame(run_parameters)
  
  individual_IDs <- sapply(results, function(hmm) hmm$data$ID[1])
  individual_alpages <- get_individual_alpage(individual_IDs, individual_info_file)
  for (alpage in unique(individual_alpages)) {
    res_index  <- individual_alpages == alpage
    data_alpage <- do.call("rbind", lapply(results[res_index], function(result) result$data))
    results_df <- do.call("rbind", lapply(results[res_index], hmm_result_to_data.frame))
    
    runs_to_pdf(alpage, parameters_df, results_df, data_alpage , paste(round(difftime(endTime, startTime, units='mins'),2), "min"), output_dir, paste0(output_dir,"1 HMM Fit/HMM_fitting_",alpage,".pdf"), show_performance_indicators = FALSE)
  }
  
  ### SAVE RESULTING TRAJECTORIES
  data_hmm <- do.call("rbind", lapply(results, function(result) result$data))
  viterbi_trajectory_to_rds(data_hmm, output_rds_file, individual_info_file)
}







### 4.1. FLOCK STOCKING RATE (charge) BY DAY AND BY STATE ###
#-------------------------------------------------------------#
if (F) {
  library(adehabitatHR)
  library(data.table)
  library(snow)
  source(file.path(functions_dir, "Functions_map_plot.R"))
  source(file.path(functions_dir, "Functions_flock_density.R"))
  
  # ENTREES
  # Un .RDS contenant les trajectoires catégorisées par comportement
  input_rds_file = file.path(output_dir, paste0("Catlog_",YEAR,"_",alpage,"_viterbi.rds"))
  
  # Un data.frame contenant les tailles de troupeaux, le pourcentage de temps d’allumage des colliers sur une journée et le chemin de la carte de phénologie phenOTB (dont la grille sert de base à l’analyse raster).
  # Un data.frame contenant les évolutions des tailles de troupeaux en fonction de la date (une ligne par évolution).
  # Doit contenir les colonnes  "alpage", "date_debut_periode", "taille_totale_troupeau"
  flock_size_file <- file.path(data_dir, paste0(YEAR,"_tailles_troupeaux.csv"))
  #A reprendre mais vérif les inputs
  Tailles_troupeaux <- read.csv(flock_size_file, stringsAsFactors = FALSE, encoding = "UTF-8")
  # Les alpages à traiter
  alpages = ALPAGES
  
  
  
  
  
  
  
  
  
  # Définition du dossier principal pour stocker les résultats des chargements
  save_dir <- file.path(data_dir, "Chargements_calcules")
  
  # Création du dossier principal s'il n'existe pas
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Fonction pour créer un chemin spécifique pour chaque alpage et année
  get_alpage_save_dir <- function(alpage, year) {
    alpage_dir <- file.path(save_dir, paste0(alpage, "_", year))
    
    # Création du sous-dossier pour l'alpage s'il n'existe pas
    if (!dir.exists(alpage_dir)) dir.create(alpage_dir, recursive = TRUE)
    
    return(alpage_dir)
  }
  
  # Préfixes des fichiers de sortie
  state_daily_rds_prefix <- paste0("by_day_and_state_", YEAR, "_")
  daily_rds_prefix <- paste0("by_day_", YEAR, "_")
  state_rds_prefix <- paste0("by_state_", YEAR, "_")
  total_rds_prefix <- paste0("total_", YEAR, "_")
  
  
  h <- 25 # Distance caractéristique pour calculer le chargement, écart-type de la gaussienne 2D sur laquelle chaque point est "dilué"
  
  # Boucle sur les alpages
  for (alpage in alpages) {
    
    # Définition du dossier de stockage spécifique à l'alpage
    alpage_save_dir <- get_alpage_save_dir(alpage, YEAR)  # Créé automatiquement le dossier "Viso_2022/" par exemple
    
    # Récupération de la taille du troupeau et du temps de fonctionnement du collier
    flock_sizes <- get_flock_size_through_time(alpage, flock_size_file)
    prop_time_collar_on <- get_alpage_info(alpage, AIF, "proportion_jour_allume")
    
    # Chargement des données
    data <- readRDS(input_rds_file)
    data <- data[data$alpage == alpage,]
    
    # Chargement du raster de phénologie, recadré selon les données
    pheno_t0 <- get_raster_cropped_L93(
      get_alpage_info(alpage, AIF, "chemin_carte_phenologie"), 
      get_minmax_L93(data, 100), 
      reproject = TRUE, 
      band = 2, 
      as = "SpatialPixelDataFrame"
    )
    
    # Sauvegarde des charges journalières et par état dans le dossier dédié à l'alpage
    flock_load_by_day_and_state_to_rds_kernelbb(
      data, 
      pheno_t0, 
      alpage_save_dir,  # ⬅️ Utilisation du bon dossier
      file.path(alpage_save_dir, paste0(state_daily_rds_prefix, alpage, ".rds")), 
      flock_sizes, 
      prop_time_collar_on
    ) 
    
    # Libération de la mémoire
    rm(data)
    
    # Vérification : Chargement des résultats
    charge <- readRDS(file.path(alpage_save_dir, paste0(state_daily_rds_prefix, alpage, ".rds")))
  }
  
    # By state
    charge_state <- charge %>%
      group_by(x,y,state) %>%
      summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
      as.data.frame()
    save_file <- paste0(save_dir,state_rds_prefix,alpage,".rds")
    saveRDS(charge_state, file = save_file)
    rm(charge_state)
    
    # By day
    # charge_day <- charge %>%
    #     group_by(x,y,day) %>%
    #     summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
    #     as.data.frame()
    # save_file <- paste0(save_dir,daily_rds_prefix,alpage,".rds")
    # saveRDS(charge_day, file = save_file)
    # rm(charge_day)
    # Alternative if the data.frame is too large :
    charge_day = lapply(unique(charge$day), function(d) {charge  %>%
        filter(day == d) %>%
        group_by(x,y,day) %>%
        summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') } )
    charge_day = as.data.frame(rbindlist(charge_day, use.names=TRUE))
    save_file <- paste0(save_dir,daily_rds_prefix,alpage,".rds")
    saveRDS(charge_day, file = save_file)
    rm(charge_day)
    
    # Total
    charge_tot <- charge %>%
      group_by(x,y) %>%
      summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
      as.data.frame()
    save_file <- paste0(save_dir,total_rds_prefix,alpage,".rds")
    saveRDS(charge_tot, file = save_file)
    rm(charge_tot)
    
    rm(charge)
  }
}




### 4.2. IF NEEDED: READAPT FLOCK STOCKING RATE TO NEW FLOCK SIZES ###
#--------------------------------------------------------------------#
if (F) {
  # Permet de recalculer un chargement déjà calculer pour l’adapter à une nouvelle estimation des évolutions de la taille
  # du troupeau, sans pour autant relancer le lourd calcul du home-range.
  source("Functions/Functions_flock_density.R")
  # ENTREES
  # Un data.frame contenant les tailles de troupeaux, le pourcentage de temps d’allumage des colliers sur une journée
  # et le chemin de la carte de phénologie phenOTB (dont la grille sert de base à l’analyse raster).
  # Un data.frame contenant les évolutions des tailles de troupeaux en fonction de la date (une ligne par).
  # Doit contenir les colonnes  "alpage", "date_debut_periode", "taille_totale_troupeau"
  flock_size_file <- paste0(data_dir,YEAR,"_tailles_troupeaux.csv")
  # Les alpages à traiter
  alpages = ALPAGES
  alpages = c("Viso")
  
  # SORTIES
  # Dossier de sortie
  save_dir = paste0(data_dir,"Chargements_calcules/")
  # Un .RDS par alpage contenant les charges journalières par comportement
  state_daily_rds_prefix = paste0("by_day_and_state_",YEAR,"_")
  # Un .RDS par alpage contenant les charges journalières
  daily_rds_prefix = paste0("by_day_",YEAR,"_")
  # Un .RDS par alpage contenant les charges par comportement
  state_rds_prefix = paste0("by_state_",YEAR,"_")
  # Un .RDS par alpage contenant la charge totale sur toute la saison
  total_rds_prefix = paste0("total_",YEAR,"_")
  
  
  h <- 25 # Distance caractéristique pour calculer le chargement, écart-type de la gaussienne 2D sur laquelle chaque point est "dilué"
  
  for (alpage in alpages) {
    flock_sizes <- get_flock_size_through_time(alpage, flock_size_file)
    prop_time_collar_on <- get_alpage_info(alpage, AIF, "proportion_jour_allume")
    
    charge <- readRDS(paste0(save_dir,state_daily_rds_prefix,alpage,"_NEW.rds"))
    s2 = sum(charge$Charge)/100
    s2
    charge <- readRDS(paste0(save_dir,state_daily_rds_prefix,alpage,".rds"))
    s1 = sum(charge$Charge)/100
    s1
    
    sum(flock_sizes[1:284])*0.75
    
    # Recompute stocking rates by state and day
    charge <- charge %>%
      group_by(day) %>%
      group_modify(function(charge_init, d) {
        recompute_daily_flock_load_by_state(charge_init, flock_sizes[as.numeric(d)], prop_time_collar_on)},
        .keep = FALSE ) %>%
      ungroup() %>%
      as.data.frame()
    save_file <- paste0(save_dir,state_daily_rds_prefix,alpage,"_NEW.rds")
    saveRDS(charge, file = save_file)
    
    # By state
    charge_state <- charge %>%
      group_by(x,y,state) %>%
      summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
      as.data.frame()
    save_file <- paste0(save_dir,state_rds_prefix,alpage,".rds")
    saveRDS(charge_state, file = save_file)
    rm(charge_state)
    
    # By day
    charge_day <- charge %>%
      group_by(x,y,day) %>%
      summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
      as.data.frame()
    save_file <- paste0(save_dir,daily_rds_prefix,alpage,".rds")
    saveRDS(charge_day, file = save_file)
    rm(charge_day)
    # Alternative if the data.frame is too large :
    # charge_day = lapply(unique(charge$day), function(d) {charge  %>%
    #                                                         filter(day == d) %>%
    #                                                         group_by(x,y,day) %>%
    #                                                         summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') } )
    # charge_day = as.data.frame(rbindlist(charge_day, use.names=TRUE))
    # save_file <- paste0(save_dir,daily_rds_prefix,alpage,".rds")
    # saveRDS(charge_day, file = save_file)
    # rm(charge_day)
    
    # Total
    charge_tot <- charge %>%
      group_by(x,y) %>%
      summarise(Charge = sum(Charge, na.rm = T), .groups = 'drop') %>%
      as.data.frame()
    save_file <- paste0(save_dir,total_rds_prefix,alpage,".rds")
    saveRDS(charge_tot, file = save_file)
    rm(charge_tot)
    
    rm(charge)
  }
}





 
### FLOCK LOAD COMPUTING ###
#--------------------------#

# OLD: compute UD using kernel around points (bivariate normal or Epanechnikov kernel)
flock_load_from_daily_and_state_UD <- function(UD, n_points_state, n_points_total, flock_size, prop_time_collar_on) {
    # INPUTS
    #   UD : the kernel computed by adehabitatHR for this day and this state
    #   n_points_state : number of relocation attributed to this state this day
    #   n_points_total : total number of relocations this day
    #   flock_size : number of animals in the flock
    # OUTPUT
    #   a data.frame with Charge and x and y coordinates

    UD <- as.data.frame.estUD(UD)
    colnames(UD) <- c("Charge", "x", "y")
    UD$Charge <- (UD$Charge
        / sum(UD$Charge)
        * n_points_state
        / n_points_total
        * flock_size
        * prop_time_collar_on
        * 10000 / abs(min(diff(unique(UD$x)))) / abs(min(diff(unique(UD$y))))) # To get a result in animals/ha for each pixel
    return(UD)
}

flock_load_by_day_and_state_to_rds_kernelUD <- function(data, grid, save_dir, save_rds_name, flock_size, prop_time_collar_on) {
    # Compute flock load by day and state
    # INPUTS
    #    data : a dataframe of animal relocations with x, y and time columns
    #    grid : a spatialPixelDataFrame object, the load will be computed on the same grid
    #    save_dir : path to the directory the results will be saved to.
    #    save_rds_name : name of the rds file the results will be saved to.
    #    flock_size : the number of animals in the flock
    # OUPUTS
    #    a rds by state file containing a data.frame with "x", "y", "day" and "Charge" columns 

    data$day <- yday(data$time)
    data <- SpatialPointsDataFrame(data[, c("x", "y")], data)
    proj4string(data) <- CRS_L93

    
    for (state in unique(data$state)) {
        save_rds_prefix = paste0(save_dir, save_rds_name,"_",state,"_")
        if (file.exists(paste0(save_dir,save_rds_name))) {
            file.remove(paste0(save_dir,save_rds_name)) #Delete file if it exists
        }

        days = unique(data$day)

        removed_days =  numeric(0) # kernelUD needs at least 5 relocations, days with less relocations are removed
        for (d in days) {
            if(nrow(data[data$state == state & data$day == d,]) < 5) {
                data = data[!(data$state == state & data$day == d),]
                removed_days = append(removed_days, d)
            }
        }

        hr <- kernelUD(data[data$state == state, "day"], grid = grid, h = h, kern = "epa") # kern :a character string. If "bivnorm", a bivariate normal kernel is used. If "epa", an Epanechnikov kernel is used.
        charge_jour <- flock_load_from_daily_and_state_UD(hr[[as.character(days[!(days %in% removed_days)][1])]], 1, 1, flock_size, prop_time_collar_on)
        ordre = order(charge_jour$x, charge_jour$y) # To get the right x and y orders

        for (d in days) {
            print(paste0("Day ",d,"/",max(unique(data$day))," state ",state))
            if (d %in% removed_days) { # If the day was removed, we consider the load is equal to 0 everywhere
                charge_jour <- pheno_t0 %>%
                                as.data.frame(xy = T) %>%
                                dplyr::select("x", "y") %>%
                                mutate(Charge = 0)
                charge_jour <- charge_jour[order(charge_jour$x, charge_jour$y), ]
            } else {
                charge_jour <- flock_load_from_daily_and_state_UD(hr[[as.character(d)]], 
                                    n_points_state = sum(data$day == d & data$state == state),
                                    n_points_total = sum(data$day == d), flock_size, prop_time_collar_on)
                charge_jour <- charge_jour[ordre, ]
            }
            charge_jour$day <- d
            charge_jour$state <- state
            saveRDS(charge_jour, file = paste0(save_rds_prefix,d,".rds"))
        }
        rm(charge_jour)
        rm(hr)
    }

    # Merge the .rds files in one
    files = list.files(save_dir, pattern=paste0("^",save_rds_name), full.names=T)
    charge = as.data.frame(rbindlist(lapply(files, readRDS), use.names=TRUE))
    saveRDS(charge, paste0(save_dir,save_rds_name))
    lapply(files, file.remove)
}


# NEW: compute UD using a Brownian Bridge kernel around trajectories
flock_load_from_daily_and_state_UD_kernelbb <- function(UD, n_points_state, n_points_total, flock_size, prop_time_collar_on) {
    # INPUTS
    #   UD : the kernel computed by adehabitatHR for this day and this state
    #   n_points_state : number of relocation attributed to this state this day
    #   n_points_total : total number of relocations this day
    #   flock_size : number of animals in the flock
    # OUTPUT
    #   a data.frame with Charge and x and y coordinates

    if (is.list(UD)) {
        UD <- lapply(UD, as.data.frame.estUD) %>%
                     do.call(rbind, .)
    }
    else {
        UD <- as.data.frame.estUD(UD)
    }

    UD %>% rename(Charge=ud, x=Var2, y=Var1) %>%
            group_by(x,y) %>%
            summarise(Charge=sum(Charge, na.rm = T), .groups="drop") %>% # it is necessary to remove NA, because if there is only one point in the burst, the corresponding UD is NaN everywhere
            mutate(Charge = (Charge
                            / sum(Charge)
                            * n_points_state
                            / n_points_total
                            * flock_size
                            * prop_time_collar_on
                            * 10000 / abs(min(diff(unique(x)))) / abs(min(diff(unique(y)))))) %>% # To get a result in animals/ha for each pixel
            return()
}

flock_load_by_day_and_state_to_rds_kernelbb <- function(data, grid, save_dir, save_rds_name, flock_sizes, prop_time_collar_on) {
    # Compute flock load by day and state
    # INPUTS
    #    data : a dataframe of animal relocations with x, y and time columns
    #    grid : a spatialPixelDataFrame object, the load will be computed on the same grid
    #    save_dir : path to the directory the results will be saved to.
    #    save_rds_name : name of the rds file the results will be saved to.
    #    flock_sizes : a 365-long vector containing the number of animals in the flock for each year_day (values for days without GPS relocations are not used)
    # OUPUTS
    #    a rds by state file containing a data.frame with "x", "y", "day" and "Charge" columns 

    data$day <- yday(data$time)
    # data <- SpatialPointsDataFrame(data[, c("x", "y")], data)
    # proj4string(data) <- CRS_L93
    if (file.exists(paste0(save_dir,save_rds_name))) {
        file.remove(paste0(save_dir,save_rds_name)) #Delete file if it exists
    }
    Tmax = 42 #in minutes for split_at_gap
    hmin = 15
    # Ds <- list(Repos = 0.7, #for BRB
    #            Paturage = 3,
    #            Deplacement = 8)
    Ds <- list(Repos = 1.25,
                Paturage = 3,
                Deplacement = 4.5)
    
    for (state in unique(data$state)) {
        save_rds_prefix = paste0(save_dir, save_rds_name,"_",state,"_")

        days = unique(data$day)

        data_state = data %>%
                filter(state == !!state)
        ltr = as.ltraj(xy = data_state[c("x", "y")], date = as.POSIXct(data_state$time), id = data_state$ID)

        # print(paste("DIFFUSION COEFFICIENTS COMPUTED FOR STATE",state))
        # print(liker(ltr, rangesig1 = c(0,10), sig2 = 10))
        # print(BRB.D(ltr, Tmax = Tmax, Lmin = 0, habitat = NULL, activity = NULL))

        # Prepare cluster for parallelised computation of daily utilisation distribution
        clus <- makeCluster(ncores, outfile='') # outfile='' is verbose option
        clusterExport(clus, list("days", "state", "save_rds_prefix", "data", "data_state", "pheno_t0", "flock_sizes", "prop_time_collar_on"), envir = environment())
        clusterExport(clus, as.list(lsf.str(.GlobalEnv))) # export all manually loaded functions
        clusterCall(clus, function() {
            # For the pipes
            suppressPackageStartupMessages(library(tidyverse)) # includes ggplot2 and dplyr among others
            suppressPackageStartupMessages(library(adehabitatHR))
            options(warn=0)
        })

        parLapply(clus, days,
                    function(d) {
                        print(paste0("Day ",d,"/",max(days)," state ",state))

                        ltr = data_state %>%
                            filter(day == d) %>%
                            mutate(time = as.POSIXct(time))

                        if( nrow(ltr) > 0 ) {
                            ltr = split_at_gap(ltr, max_gap = Tmax) #split the trajectory of each individual into bursts of the specified state (otherwise, every bursts are bridged together when using kernelbb)
                            ltr = as.ltraj(xy = ltr[c("x", "y")], date = ltr$time, id = ltr$ID)

                            hr <- kernelbb(ltr, sig1 = Ds[[state]], sig2=10, grid = grid, same4all = FALSE, byburst = TRUE,
                                    extent = 0.5, nalpha = 25)

                            # hr <- BRB(ltr, Ds[[state]], Tmax = Tmax, Lmin = 1, hmin = hmin, type="UD",
                            # filtershort=TRUE, grid = 20, b=FALSE, same4all=TRUE, extent=1, tau = NULL,
                            # boundary=NULL)

                            charge_jour <- flock_load_from_daily_and_state_UD_kernelbb(hr,
                                                n_points_state = sum(data$day == d & data$state == state),
                                                n_points_total = sum(data$day == d), flock_sizes[d], prop_time_collar_on)

                            charge_jour$day <- d
                            charge_jour$state <- state
                            saveRDS(charge_jour, file = paste0(save_rds_prefix,d,".rds"))
                        }
                } )
        stopCluster(clus)
    }

    # Merge the .rds files in one
    files = list.files(save_dir, pattern=paste0("^",save_rds_name), full.names=T)
    charge = as.data.frame(rbindlist(lapply(files, readRDS), use.names=TRUE))
    saveRDS(charge, paste0(save_dir,save_rds_name))
    lapply(files, file.remove)
}


split_at_gap <- function(data, max_gap = 60, shortest_track = 0) {
    #' Courtesy Théo Michelot
    #' Split track at gaps
    #'
    #' @param data Data frame with (at least) columns for "ID" and "time"
    #' @param max_gap Longest allowed gap, in minutes (track will be split at longer gaps)
    #' @param shortest_track Shortest track to keep after splitting, in minutes. Shorter
    #' tracks will be removed from the output data set.
    #'
    #' @return Data frame with identical structure as input, where ID column
    #' has been replaced by new ID for split tracks. Old ID still accessible as
    #' ID_old column

    # Number of tracks
    n_tracks <- length(unique(data$ID))

    # Save old ID and reinitialise ID column
    data$ID_old <- data$ID
    data$ID <- character(nrow(data))

    # Loop over tracks (i.e., over IDs)
    for(i_track in 1:n_tracks) {
        # Indices for this track
        ind_this_track <- which(data$ID_old == unique(data$ID_old)[i_track])
        track_length <- length(ind_this_track)

        # Time intervals in min
        dtimes <- difftime(data$time[ind_this_track[-1]],
                           data$time[ind_this_track[-track_length]],
                           units = "mins")

        # Indices of gaps longer than max_gap
        ind_gap <- c(0, which(dtimes > max_gap), track_length)

        # Create new ID based on split track
        subtrack_ID <- rep(1:(length(ind_gap) - 1), diff(ind_gap)) #la diff permet d’avoir la longueur de chaque subtrack, puis le rep va permettre de répéter l’indice ordinal de cette subtrack autant de fois que cette longueur
        data$ID[ind_this_track] <- paste0(data$ID_old[ind_this_track], "-", subtrack_ID)
    }

    # Only keep sub-tracks longer than some duration
    track_lengths <- sapply(unique(data$ID), function(id) {
        ind <- which(data$ID == id)
        difftime(data$time[ind[length(ind)]], data$time[ind[1]], units = "min")
    })
    ID_keep <- names(track_lengths)[which(track_lengths >= shortest_track)]
    data <- subset(data, ID %in% ID_keep)

    return(data)
}


get_flock_size_through_time <- function(alpage, flock_size_file) {
    # INPUTS
    #   alpage: the name of the pasture of interest
    #   flock_size_file: path to a CSV file containing the evolutions of flock sizes (one line is the new size of the flock from the specified date)
    #                    Must contain les colonnes  "alpage", "date_debut_periode", "taille_totale_troupeau"
    # OUTPUT
    #   a 365-long vector containing the number of animals in the flock for each year_day (values for days without GPS relocations are not used)
    dates_sizes <- read.csv(flock_size_file, header=TRUE, sep=",") %>%
        filter(alpage == !!alpage) %>%
        mutate(yday = yday(as.POSIXct(date_debut_periode, tz="GMT", format="%d/%m/%Y"))) %>%
        arrange(yday)

    sizes = rep(0, 1, 365)
    for (i in 1:nrow(dates_sizes)) {
        sizes[dates_sizes$yday[i]:365] = dates_sizes$taille_totale_troupeau[i]
    }
    
    return(sizes)
}


recompute_daily_flock_load_by_state <- function(charge_d, flock_size_d, prop_time_collar_on) {
    charge_tot_init = sum(charge_d$Charge)
    charge_tot_fin = flock_size_d * 
                        prop_time_collar_on *
                        10000 / abs(min(diff(unique(charge_d$x)))) / abs(min(diff(unique(charge_d$y)))) # To get a result in animals/ha for each pixel
    n_states = length(unique(charge_d$state))

    for(s in unique(charge$state)) {
        I_s = (charge_d$state == s)
        charge_d$Charge[I_s] = charge_d$Charge[I_s] / charge_tot_init * charge_tot_fin
    }
    return(charge_d)
}

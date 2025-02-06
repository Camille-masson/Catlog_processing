### FUNCTIONS ###
#---------------#

generate_one_frame <- function(data, initial_time, time_per_frame, remanance_time, output_file, hmm_scale, ombrage, charge) {
    print(initial_time)
    data1 = data[(data$time >= initial_time) & (data$time < initial_time + time_per_frame), ]
    title = format(max(data1$time)+2*60*60, format = "%d/%m/%Y %H:%M")

    if(nrow(data1)>0) {
        data2 = data[(data$time >= initial_time - remanance_time) & (data$time < initial_time + time_per_frame), ]
        instant <- data1 %>%
                    group_by(ID) %>%
                    filter(time == max(time)) %>%
                    as.data.frame()
        
        charge = charge[charge$day == max(data1$day)-1, c("x","y","cumulative_charge")]
        
        hour_colour_scale = scale_colour_gradient2(low = "blue", mid = "#ff0000", high = "blue", midpoint = 12, breaks = c(0, 6, 12, 18, 24),
                                                    limits = c(5, 23),
                                                    name = "Dernières 24 heures", guide="colorbar")
        
        pdf(file = NULL) # do not show plot
        ggplot(data1, aes(x, y, group=ID)) + # l’argument group permet de lier les différents points indépendament de leur couleurs (sinon il y a un chemin par couleur). Pour tout lier, on peut mettre group="all" ou n’importe quelle autre chaîne de caractères
            # Map background
            geom_raster(data=ombrage, aes(x=x, y=y, fill=Ombrage), alpha=0.5, inherit.aes = FALSE) +
            scale_fill_gradientn(colours=c("black","white"), guide="none") +
            coord_equal() +
            # Cumulative charge up to the day before
            new_scale("fill") +
            scale_fill_viridis(option = "E", direction = -1,
                                      oob = scales::squish, # squish anything that exceeds the limits (out-of-bounds/oob values) to the nearest extreme
                                      limits = c(1, 500),
                                      guide = triangle_colourbar_up(), name = "Charge depuis le début de la saison (brebis.jour/ha)") +
            scale_alpha(range=c(0, 0.7), limits = c(0, 20), oob = scales::squish, guide="none") +
            geom_tile(data=charge, aes(x=x, y=y, fill = cumulative_charge, alpha = cumulative_charge), inherit.aes = FALSE) +
            # Trajectory over the last 24 hours, coloured according to the daytime
            new_scale("alpha") +
            geom_path(data=data2, linewidth = 0.5, aes(alpha=time, col=hour)) +
            scale_alpha(range=c(0, 0.8), guide="none") +
            hour_colour_scale +
            scale_alpha_continuous(range=c(0.07 ,0.7), guide="none") +
            # Trajectory over the last 30 minutes, coloured according to the behavioural state, and last animal position
            new_scale("colour") +
            geom_path(linewidth = 0.7, aes(col = state)) + 
            # geom_point(size = 0.7, aes(col = state), show.legend=FALSE) +
            geom_point(data = instant, size = 2.2, col = "black", aes(shape = ID), show.legend=FALSE) +
            geom_point(data = instant, size = 1.5, aes(shape = ID, col = state), show.legend=FALSE) +
            hmm_scale +
            scale_shape_manual(values=c(0, 1, 2, 5, 6, 11, 7, 9, 10, 12, 13, 14, 0)) +
            # Plot elements
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5)) + # legend.key.size = unit(1, 'cm')
            labs(colour = "") +
            xlab("") +
            ylab("")
        ggsave(output_file, width = 14, height = 7)
        dev.off() # close the current graphical device to remove plot from memory
    }
}

generate_all_frames_parallelized <- function(data, cumulative_charge, ombrage, time_per_frame, remanance_time, output_folder) {
    # Crop charge to Ombrage extent
    xy_minmax = get_minmax_L93(ombrage)
    charge_month <- cumulative_charge %>%
                    filter(day %in% (unique(data$day)-1) ) %>% # -1 to get the charge of the last day of the month before
                    filter((x>xy_minmax$x_min) & (x<xy_minmax$x_max) & (y>xy_minmax$y_min) & (y<xy_minmax$y_max))
                    
    # Generating images one by one
    data$state = as.factor(data$state)
    frames_info = data.frame(initial_time = format(seq(min(data$time), max(data$time), time_per_frame), format="%Y-%m-%d %H:%M:%S %Z"))
    frames_info$index = 1:nrow(frames_info)

    startTime <- Sys.time()
    clus <- makeCluster(ncores) # outfile='' is verbose option
    clusterExport(clus, append(as.list(lsf.str(.GlobalEnv)), list("data", "output_dir", "alpage", "time_per_frame", "remanance_time",
                    "hmm_scale", "ombrage", "charge_month", "output_folder")), envir = environment())
    # Load libraries
    clusterCall(clus, function() {
        library(tidyverse) # includes ggplot2 and dplyr among others
        theme_set(theme_bw()) # theme for ggplot2
        library(ggnewscale)
        library(viridis)
        library(grid)
        library(gtable)
        # GIS packages
        library(terra)
        library(tidyterra)
        library(maptiles)
    })

    parRapply(clus, frames_info, function(fi) {
        output_file = paste0(output_folder,as.numeric(fi[2]),".png")
        generate_one_frame(data, as.POSIXct(fi[1], tz="GMT"), time_per_frame, remanance_time, output_file,
                            hmm_scale, ombrage, charge_month)
    })

    stopCluster(clus)
    endTime <- Sys.time()
    print(paste("+++ Cluster total excecution time :", round(endTime - startTime,2), "min +++"))
}


pngs_to_mp4 <- function(pngs_folder, output_video_file) {
    # Generates and saves video, removes temporary files

    imgs <- list.files(path=pngs_folder, pattern="*.png")
    numbers = as.numeric(regmatches(imgs, regexpr("[0-9]+", imgs)))
    imgs = imgs[order(numbers)]
    imgs <- paste0(output_dir,alpage,"/animation/", imgs)

    av::av_encode_video(imgs, framerate = 6,
                        output = output_video_file)

    unlink(pngs_folder, recursive=TRUE)
}
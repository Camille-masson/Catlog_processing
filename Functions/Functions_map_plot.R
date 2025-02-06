library(viridis)
library(terra)
source("Functions/Functions_ggplot_custom.R")
theme_set(theme_bw()) # theme for ggplot2

if (!exists("BACKGROUND_TYPE")) { # If no user-defined background-type was chosen, then use the default one
    # BACKGROUND_TYPE = "BDALTI" # using a hillshade raster that was made with IGN’s BDALTI
    BACKGROUND_TYPE = "ESRI" # based on maptiles library
}


### MAP BACKGROUND FUNCTIONS ###
#******************************#
if(BACKGROUND_TYPE == "BDALTI") {
    get_ombrage  <- function(minmax_xy_L93) {
        ombrage <- get_raster_cropped_L93(paste0(raster_dir,"BDALTI/Fond_topo.tif"), minmax_xy_L93, as = "data.frame")
        colnames(ombrage)[3] = "Ombrage"
        return(ombrage)
    }

    add_ombrage_layer <- function(ombrage, UP) {
        minmax_xy_L93 = get_minmax_L93(ombrage)
        return( list(geom_raster(data = ombrage, aes(x = x, y = y, fill = Ombrage), alpha = 0.7, inherit.aes = FALSE, show.legend = FALSE),
                     scale_fill_gradientn(colours = c("black", "white")),
                     coord_equal(),
                     new_scale("fill"),
                     if(is.null(UP)) {NULL} else {geom_spatvector(data = UP, linewidth=0.7, linetype="dotdash", alpha=0.6, fill=NA)}, # if UP is unspecified or NULL, then it is not plotted
                     xlab(NULL), 
                     ylab(NULL),
                     scale_x_continuous(expand = c(0, 0), limits=c(minmax_xy_L93$x_min, minmax_xy_L93$x_max), guide=guide_axis(angle=45)),
                     scale_y_continuous(expand = c(0, 0), limits=c(minmax_xy_L93$y_min, minmax_xy_L93$y_max)),
                     # scale_y_continuous(labels = function(x) paste(round(x/1000), "km")),
                     theme(axis.text.x = element_text(size=7),
                           axis.text.y = element_text(size=7),
                           legend.title = element_text(size=8),
                           plot.title = element_text(hjust=0) ) ))
    }
} else if (BACKGROUND_TYPE == "ESRI") {
    # library(OpenStreetMap)
    library(maptiles)
    library(tidyterra)
    
    get_ombrage  <- function(minmax_xy_L93) {
        e1 <- ext(minmax_xy_L93$x_min, minmax_xy_L93$x_max, minmax_xy_L93$y_min, minmax_xy_L93$y_max) %>%
                vect(CRS_L93)
        e <- e1 %>% terra::project(CRS_WSG84)


        ombrage <- get_tiles(e, crop = TRUE, provider = "Esri.WorldTopoMap", zoom=14) %>% # see https://wiki.openstreetmap.org/wiki/Zoom_levels for zoom
                terra::project(CRS_L93) %>%
                crop(e1)
        
        return(ombrage)
    }

    add_ombrage_layer <- function(ombrage, UP) {
        minmax_xy_L93 = get_minmax_L93(ombrage)
        return( list(coord_sf(crs = CRS_L93),
                     geom_spatraster_rgb(data = ombrage, alpha = 0.8),
                     if(is.null(UP)) {NULL} else {geom_spatvector(data = UP, inherit.aes=F, linewidth=0.7, linetype="dotdash", alpha=0.6, fill=NA)}, # if UP is unspecified or NULL, then it is not plotted
                     xlab(NULL), 
                     ylab(NULL),
                     scale_x_continuous(expand = c(0, 0), limits=c(minmax_xy_L93$x_min, minmax_xy_L93$x_max), guide=guide_axis(angle=45)),
                     scale_y_continuous(expand = c(0, 0), limits=c(minmax_xy_L93$y_min, minmax_xy_L93$y_max)),
                     # scale_y_continuous(labels = function(x) paste(round(x/1000), "km")),
                     theme(axis.text.x = element_text(size=7),
                           axis.text.y = element_text(size=7),
                           legend.title = element_text(size=8),
                           plot.title = element_text(hjust=0) ),
                     annotate(geom = 'text', label = "Tiles © Esri", x = Inf, y = -Inf, hjust = 1.1, vjust = -1, color="#505050", size=3.5) )) # Credit the map
    }
}

get_UP_polygon <- function(alpage, alpage_info_file, UP_file) {
    UP_name1 = get_alpage_info(alpage, alpage_info_file, "nom1_UP")

    UP = vect(UP_file) %>%
        filter(nom1 == UP_name1) %>%
    return()
}

get_explored_space_polygon <- function(charge_tot, threshold = 0.1) {
    charge_tot <- mutate(Charge = (Charge > 0.1)) %>%
        rast(type="xyz", crs=CRS_L93) %>% 
        as.polygons(e) %>% 
        filter(Charge == 1) %>%
        return()
}


### RASTER FUNCTIONS ###
#**********************#

get_raster_cropped_L93 <- function(raster_file, minmax_xy_L93, reproject = FALSE, reproject_method = "bilinear", band = 1, as) {
    # reproject_method : "bilinear" for interpolation, 'ngb' for nearest neighbourgh
    # INPUTS :
    #   minmax_xy_l93 = a list of x_min, x_max, y_min, y_max in L93 to which the raster should be cropped
    #   as = the desired type of object, either "SpatialPixelDataFrame", "data.frame" or "spatRaster"
    
    raster <- terra::rast(raster_file)
    raster <- terra::subset(raster, band)

    e <- ext(minmax_xy_L93$x_min, minmax_xy_L93$x_max, minmax_xy_L93$y_min, minmax_xy_L93$y_max)

    if (crs(raster) == "") {
        crs(raster) <- CRS_L93
    } else if (reproject) {
        raster <- project(raster, y = CRS_L93, method = reproject_method)
    }

    raster <- terra::crop(raster, e)
    
    if(as == "SpatialPixelDataFrame") {
        raster <- as.data.frame(raster, xy=TRUE, na.rm = FALSE)
        raster <- SpatialPixelsDataFrame(raster[c(1,2)], raster[3])
        proj4string(raster) <- CRS_L93
    } else if (as == "data.frame") {
        raster <- as.data.frame(raster, xy=TRUE, na.rm = FALSE)
    }

    return(raster)
}

plot_raster_over_ombrage <- function(data, col_name, scale, alpha = 0.7, ombrage, title, label, UP=NULL) {
    col_name <- as.symbol(col_name)
    print(ggplot(data) +
        add_ombrage_layer(ombrage, UP) +
        geom_raster(aes(x, y, fill = !!col_name), alpha = alpha) +
        scale +
        guides(alpha = FALSE) +
        ggtitle(title) +
        labs(fill = label) +
        colourbar_right())
}

plot_raster_over_ombrage_varying_alpha <- function(data, col_name, scale, ombrage, title, label, UP=NULL) {
    col_name <- as.symbol(col_name)
    print(ggplot(data) +
        add_ombrage_layer(ombrage, UP) +
        geom_tile(aes(x, y, fill = !!col_name, alpha = !!col_name)) +
        scale +
        scale_alpha(trans = "log10") +
        guides(alpha = FALSE) +
        ggtitle(title) +
        labs(fill = label))
}

charge_labels <- function(x) {
  l = length(x)
  x = as.character(round(x))
  if(l == 2) { x[l] = ""
  } else { x[l] = paste0(">", x[l]) }
  return(x)
}


plot_charge <- function(data, ombrage, title = "", charge_max = 1000, UP=NULL, col_name = "Charge", lower_threshold = 10, charge_label = "brebis.jours/ha") {
    scale = scale_fill_viridis_b(option = "D",
                                    n.breaks = ifelse(charge_max<800, 7, 9),
                                    nice.breaks = FALSE,
                                    labels = charge_labels,
                                    oob = scales::squish, # squish anything that exceeds the limits (out-of-bounds/oob values) to the nearest extreme
                                    limits = c(0, charge_max),
                                    show.limits = TRUE)
    plot_raster_over_ombrage(data[data[col_name] >= lower_threshold, ], col_name,
                    scale,
                    alpha = 0.7, ombrage, title, charge_label, UP)
}

plot_perc <- function(data, col_name, ombrage, title, UP=NULL) {
  plot_raster_over_ombrage(data, col_name,
                    scale = scale_fill_viridis_b(option = "plasma", direction = -1,
                    limits = c(0, 100),
                    breaks = c(0,20,40,60,80,100), labels = c("0 %","20 %","40 %","60 %","80 %","100 %") ),
                    alpha = 0.7, ombrage, title, NULL, UP)
}

set_scale_max_level <- function(data, quantile = 0.99, round = 100) {
    return(max(round(quantile(data, quantile)/round)*round, round))
}

crop_data.frame <- function(data, cropping_polygon, names_col = NULL, values_col = "Charge") {
    if(is.null(names_col)) {
        rast(data, type="xyz", crs=CRS_L93) %>%
            crop(cropping_polygon, mask=TRUE) %>%
            as.data.frame(xy=TRUE, na.rm=FALSE) %>%
            return()
    } else {
        names = unique(data[,names_col])
        data %>%
            pivot_wider(names_from = names_col, values_from = values_col) %>%
            rast(type="xyz", crs=CRS_L93) %>%
            crop(cropping_polygon, mask=TRUE) %>%
            as.data.frame(xy=TRUE, na.rm=FALSE) %>%
            pivot_longer(cols=names, names_to=names_col, values_to=values_col) %>%
            as.data.frame() %>%
            return()
    }
}

### VECTOR FUNCTIONS ###
#**********************#

plot_trajectory_over_ombrage <- function(data, col_name, scale, ombrage, title, label, UP=NULL) {
    col_name <- as.symbol(col_name)
    print(ggplot(data, aes(x, y, col = !!col_name)) +
            add_ombrage_layer(ombrage, UP) +
            geom_path(linewidth = 0.3) +
            geom_point(size = 0.3) +
            scale +
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5)) +
            labs(colour = label))
}

plot_polygons_over_ombrage <- function(data, col_name, scale, alpha = 0.7, ombrage, title, label, UP=NULL) {
    col_name <- as.symbol(col_name)
    print(ggplot(data) +
            add_ombrage_layer(ombrage, UP) +
            geom_spatvector(aes(fill = !!col_name), alpha = alpha, colour = "white", size = 0.2) +
            coord_sf(crs = CRS_L93, datum = CRS_L93) +
            scale +
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5)) +
            labs(fill = label) )
}

plot_charge_by_polygon <- function(data, col_name = "Charge", charge_max = 1000, ombrage, title = "", label = "brebis.jours/ha", UP) {
    plot_polygons_over_ombrage(data[as.data.frame(data[,col_name]) >= 0.1, ], col_name,
                    scale_fill_viridis(option = "H",
                                      oob = scales::squish, # squish anything that exceeds the limits (out-of-bounds/oob values) to the nearest extreme
                                      limits = c(0, charge_max)),
                    alpha = 0.7, ombrage, title, label, UP)
}

get_minmax_L93 <- function(data, buffer = 0) {
    if("SpatRaster" %in% class(data)) {
        ext = as.list(ext(data))
        return(list(x_min = ext$xmin - buffer,
                x_max = ext$xmax + buffer,
                y_min = ext$ymin - buffer,
                y_max = ext$ymax + buffer))
    }
    return(list(x_min = min(data$x, na.rm=T) - buffer,
                x_max = max(data$x, na.rm=T) + buffer,
                y_min = min(data$y, na.rm=T) - buffer,
                y_max = max(data$y, na.rm=T) + buffer))
}

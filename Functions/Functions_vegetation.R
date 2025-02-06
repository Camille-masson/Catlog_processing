# Typologie STrouMPH v1
typology_STrouMPH <- data.frame(id=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 0),
                       name = c("Pelouses de mode nival", "Formations mixtes nivales/thermiques",
                                "Pelouses intermédiaires de l’alpin", "Pelouses intermédiaires du subalpin",
                                "Nardaies denses du subalpin", "Queyrellins", "Pelouses productives",
                                "Pelouses en bombement de l’alpin", "Pelouses thermiques écorchées",
                                "Pelouses thermiques enherbées", "Pelouses thermiques à Brachypode pénné",
                                "Pelouses thermiques méditerranéo-montagnardes", "Pelouses nitrophiles",
                                "Pelouses humides", "Éboulis à ressource pastorale", "Sous-bois pastoraux",
                                "Landes", "Formations minérales", "Forêts non pastorales", "Megaphorbiaies et Aulnaies",
                                "Autres"))
typology_STrouMPH <- data.frame(id=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 0, 910, 1112, 1518, 1619, 16921),
                       name = c("Pelouses nivales", "Formations mixtes nivales/thermiques",
                                "P. intermédiaires de l’alpin", "P. intermédiaires du subalpin",
                                "Nardaies denses du subalpin", "Queyrellins", "Pelouses productives",
                                "P. en bombement de l’alpin", "P. thermiques écorchées",
                                "P. thermiques enherbées", "P. th. à Brachypode pénné",
                                "P. th. méditerranéo-montagnardes", "Pelouses nitrophiles",
                                "Pelouses humides", "Éboulis à ressource pastorale", "Sous-bois pastoraux",
                                "Landes", "Formations minérales", "Forêts non pastorales", "Megaphorbiaies et Aulnaies",
                                "Autres", "P. thermiques fus.", "P. th. montagnardes", "Minéral fus.", "Forêts fus.", "Ligneux hauts fus."))
typology_STrouMPH_abbreviations <- data.frame(id=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 0, 910, 1112, 1518, 1619, 16921),
                       name = c("P. nivales", "F. mixtes",
                                "P. int. alpin", "P. int. subalpin",
                                "Nardaies", "Queyrellins", "P. prod.",
                                "P. en bomb.", "P. th. éc.",
                                "P. th. en.", "P. th. Brachypode",
                                "P. th. méd-mont.", "P. nitrophiles",
                                "P. humides", "Éboulis", "Sous-bois",
                                "Landes", "F. minérales", "Forêts", "Megaphorb. et Aul.",
                                "Autres", "P. th. fus.", "P. th. mont.", "Minéral fus.", "Forêts fus.", "Ligneux hauts fus."))
habitat_scale_STrouMPH <- scale_fill_manual(values = c("Pelouses nivales" = "#2eacff",
                                "Formations mixtes nivales/thermiques" = "#62c0ff",
                                "P. intermédiaires de l’alpin" = "orange",
                                "P. intermédiaires du subalpin" =  "#ffd17c",
                                "Nardaies denses du subalpin" = "#2efd62",
                                "Queyrellins" = "#7bff00",
                                "Pelouses productives" = "#c47f00",
                                "P. en bombement de l’alpin" = "#ff0000",
                                "P. thermiques écorchées" = "#ff3c00",
                                "P. thermiques enherbées" = "#ff784f",
                                "P. th. à Brachypode pénné" = "#fdcebf",
                                "P. th. méditerranéo-montagnardes" = "#fff1ec",
                                "Pelouses nitrophiles" = "#ffe600",
                                "Pelouses humides" = "#576dee",
                                "Éboulis à ressource pastorale" = "#aaaaaa",
                                "Sous-bois pastoraux" = "#138630",
                                "Landes" = "#9456ff",
                                "Formations minérales" = "#414141",
                                "Forêts non pastorales" = "#0a6357",
                                "Megaphorbiaies et Aulnaies" = "#0124ee",
                                "Autres" = "#cccccc",
                                "P. thermiques fus." = "#ff3c00",
                                "P. th. montagnardes" = "#fdcebf",
                                "Minéral fus." = "#414141",
                                "Forêts fus." = "#138630",
                                "Ligneux hauts fus." = "#138630"))
# Typologie Télédec v1
typology_teledetection <- data.frame(id = c(1, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29),
                       name = c("Pelouses nivales", "Eboulis végétalisés", "Forets", "Landes ligneuses",
                       "Surfaces rocheuses", "Montagnard très productif", "Montagnard moyennement productif",
                       "Montagnard peu productif", "Alpin très productif", "Alpin moyennement productif",
                       "Alpin peu productif", "Subalpin très productif", "Subalpin moyennement productif",
                       "Subalpin peu productif"))
habitat_scale_teledetection <- scale_fill_manual(values = c("Pelouses nivales" = "#0099ff",
                        "Eboulis végétalisés" = "#a16a6a",
                        "Forets" = "#138630",
                        "Landes ligneuses" = "#2efd62",
                        "Surfaces rocheuses" = "#5a3838",
                        "Montagnard très productif" = "orange",
                        "Montagnard moyennement productif" = "#ffd17c",
                        "Montagnard peu productif" = "#c47f00",
                        "Alpin très productif" = "#ff0000",
                        "Alpin moyennement productif" = "#ff784f",
                        "Alpin peu productif" = "#ff3c00",
                        "Subalpin très productif" = "#ffe600",
                        "Subalpin moyennement productif" = "#eee69f",
                        "Subalpin peu productif" = "#b19f00"))


###  VEGETATION DATA ###
#**********************#

get_file_extension = function(file_name) {
    file_ext = strsplit(file_name, split = ".", fixed = T)[[1]]
    return(file_ext[length(file_ext)])
}

get_vegetation_rasterized <- function(vegetation_file, vegetation_typology_name, grid, cropping_polygon=NULL){
    # Get the vegetation map as a raster over a specified grid
    # INPUTS
    #    vegetation_file : the path to the .shp or .tif vegetation_map
    #    vegetation_typology_name : Either "typology_STrouMPH" or "typologie_teledetection"
    #    grid : a data.frame with "x" and "y" columns containing the Lambert93 coordinates grid the final raster should be aligned on (can contain duplicated and some missing coordinates)
    #    cropping_polygon : a terra::spatvect containing a polygon to crop the vegetation map
    # OUPUTS
    #    a raster of the vegetation map aligned on grid
    
    file_ext = get_file_extension(vegetation_file)
    grid = terra::rast(grid)
    crs(grid) = CRS_L93

    if(file_ext == "shp") { # If we get a vector file, it needs to be rasterized
        habitats <- terra::vect(vegetation_file)
        habitats <- terra::rasterize(habitats, grid, "Typ_pst", background=0)
        habitats[which(is.na(habitats[]))] = 0
    } else if (file_ext == "tif") { # If we get a raster file, it needs to be extrapolated on the right grid
        habitats <- terra::rast(vegetation_file)
        habitats <- terra::project(x = habitats, y = grid, method="near")
    }
    if (!is.null(cropping_polygon)) {habitats <- terra::crop(habitats, cropping_polygon, mask=TRUE)} # If a cropping polygon is provided, then the vegetation map is cropped
    habitats <- as.data.frame(habitats, xy = TRUE, na.rm = F)
    colnames(habitats) <- c("x", "y", "vegetation_type")
    habitats <- habitats[order(habitats$x, habitats$y), ]

    # Replace habitat types numbers by habitat names from the right typology
    typology = get_typology(vegetation_typology_name)
    habitats$vegetation_type <- as.character(sapply(habitats$vegetation_type, function(id) typology$name[which(typology$id==id)] ))
    habitats$vegetation_type[habitats$vegetation_type == "character(0)"] <- NA

    return(habitats)
}

get_vegetation_polygons <- function(vegetation_file, vegetation_typology_name, grid, cropping_polygon=NULL) {
    # Get the vegetation map as polygons
    # The polygons will be cropped to the grid extent
    # If the imput map is a raster, it will be reprojected on the grid before beeing vectorized
    # INPUTS
    #    vegetation_file : the path to the .shp or .tif vegetation_map
    #    vegetation_typology_name : Either "typology_STrouMPH" or "typologie_teledetection"
    #    grid : a data.frame with "x" and "y" columns containing a Lambert93 coordinates grid (can contain duplicated and some missing coordinates), which extend is used to crop the map
    #    cropping_polygon : a terra::spatvect containing a polygon to crop the vegetation map
    # OUPUTS
    #    a raster of the vegetation map aligned on grid
    
    file_ext <- get_file_extension(vegetation_file)
    grid <- terra::rast(grid)
    crs(grid) <- CRS_L93

    if(file_ext == "shp") { # If we get a vector file, no need to transform to polygons
        vegetation_units <- terra::vect(vegetation_file) %>%
            terra::project(y = CRS_L93) %>%
            terra::crop(y = grid)
        vegetation_units <- vegetation_units[,length(names(vegetation_units))] # only keep the vegetation type attribute
    } else if (file_ext == "tif") { # If we get a raster file, it needs to be vectorized to polygons
        vegetation_units <- terra::rast(vegetation_file) %>%
            terra::project(y = grid, method="near") %>%
            as.polygons(na.rm=T, dissolve=T) %>% # convert raster to polygons
            terra::disagg() # disolve polygons (one by habitat type) into individual vegetation units
    }
    if (!is.null(cropping_polygon)) {vegetation_units <- terra::crop(vegetation_units, cropping_polygon)} # If a cropping polygon is provided, then the vegetation map is cropped
    names(vegetation_units) <- "vegetation_type"

    # Replace by habitat names from the correct typology
    typology = get_typology(vegetation_typology_name)
    vegetation_units$vegetation_type <- as.character(sapply(vegetation_units$vegetation_type, function(id) typology$name[which(typology$id==id)] ))
    vegetation_units$vegetation_type[vegetation_units$vegetation_type == "character(0)"] <- NA

    return(vegetation_units)
}

get_typology <- function(vegetation_typology_name){
    # INPUTS
    #    vegetation_typology_name : Either "typology_STrouMPH" or "typologie_teledetection"
    # OUPUTS
    #    a data.frame containing the name of each habitat in correspondance to their id number, columns "id" and "name"


    if (vegetation_typology_name == "typology_STrouMPH") {
        return (typology_STrouMPH)
    } else if (vegetation_typology_name == "typology_STrouMPH_abbreviations") {
        return (typology_STrouMPH_abbreviations)
    } else if (vegetation_typology_name == "typology_teledetection") {
        return (typology_teledetection)
    }

}

get_habitat_scale <- function(vegetation_typology_name){
    # INPUTS
    #    vegetation_typology_name : Either "typology_STrouMPH" or "typologie_teledetection"
    # OUPUTS
    #    a ggplot2 fill scale with one colour per habitat name
    if (vegetation_typology_name == "typology_STrouMPH") {
        return (habitat_scale_STrouMPH)
    } else if (vegetation_typology_name == "typology_teledetection") {
        return (habitat_scale_teledetection)
    }
}


###  FIGURES ###
#**************#

barplot_load_by_vegetation <- function(charge_tot, habitats, ombrage) {
    habitats <- habitats[order(habitats$x, habitats$y), ]
    charge_tot <- charge_tot[order(charge_tot$x, charge_tot$y), ] 

    pixel_surface = get_pixel_surface(charge_tot)
    charge_tot$Charge <- charge_tot$Charge * pixel_surface / 10000 # going from to animals.days/ha to animals.days/pixel

    charge_tot$vegetation_type = habitats$vegetation_type
    charge_tot = charge_tot[!is.na(charge_tot$vegetation_type), ]
    habitats = habitats[!is.na(habitats$vegetation_type), ]

    habitat_labs <- get_habitat_labs(habitats)

    charge_tot <- charge_tot %>% 
        mutate(vegetation_type = factor(vegetation_type)) %>%
        group_by(vegetation_type) %>%
        summarize(Charge = sum(Charge)) %>%
        ungroup()

    ggplot(charge_tot, aes(y=vegetation_type, x=Charge, fill=vegetation_type)) +
        geom_col() +
        xlab("Présence du troupeau (brebis.jours)") +
        # xlab("Flock presence (sheep.days)") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        habitat_scale +
        guides(fill = FALSE)
}

barplot_load_by_vegetation_state <- function(charge_state, habitats, ombrage) {
    habitats <- habitats[order(habitats$x, habitats$y), ]
    charge_state <- charge_state[order(charge_state$state, charge_state$x, charge_state$y), ]

    pixel_surface = get_pixel_surface(charge_state)
    charge_state$Charge <- charge_state$Charge * pixel_surface / 10000 # going from to animals.days/ha to animals.days/pixel

    charge_state$vegetation_type = habitats$vegetation_type
    charge_state = charge_state[!is.na(charge_state$vegetation_type), ]
    habitats = habitats[!is.na(habitats$vegetation_type), ]

    habitat_labs <- get_habitat_labs(habitats)

    charge_state %>% 
        mutate(vegetation_type = factor(vegetation_type), state = factor(state)) %>%
        group_by(vegetation_type, state) %>%
        summarize(Charge = sum(Charge)) %>%
        ggplot(aes(y=vegetation_type, x=Charge, fill=state)) +
        geom_col() +
        xlab("Flock presence (sheeps.days)") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        hmm_scale_fill +
        theme(legend.position = c(0.8, .2))
}

boxplot_load_by_vegetation <- function(charge_tot, habitats, ombrage) {
    habitats <- habitats[order(habitats$x, habitats$y), ]
    charge_tot <- charge_tot[order(charge_tot$x, charge_tot$y), ]

    pixel_surface = get_pixel_surface(charge_state)

    charge_tot = charge_tot[!is.na(habitats$vegetation_type), ]
    habitats = habitats[!is.na(habitats$vegetation_type), ]

    habitat_labs <- get_habitat_labs(habitats)

    ggplot(charge_tot, aes(y=habitats$vegetation_type, x=Charge, fill=habitats$vegetation_type)) +
        geom_boxplot(outlier.shape = NA, varwidth = T) +
        xlab("Flock load (sheep.days/ha)") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        habitat_scale +
        guides(fill = F, size = F, color = guide_legend("Mean flock load", override.aes = list(linetype = 0, size = 0.8))) +
        coord_cartesian(xlim = c(0, quantile(charge_tot$Charge, 0.95))) + # contrairement à scale_x_continuous(xlim =...), coord_cartesian ne supprime pas les données hors graphe avant de calculer les statistiques
        stat_summary(fun.y = mean, mapping=aes(geom="point", color="darkred"), size=0.4,
                    position = position_dodge2(width = 0.75), show.legend = T) +
        scale_color_manual(values="darkred", labels="") +
        theme(legend.position = c(0.8, .1))
}

get_habitat_labs <- function(habitats) {
    # Creates labels for plotting bar_charts of data by habitat, the label being the habitat name and its surface in ha.
    # INPUTS
    #    habitats : a data.frame representing a raster (with Lambert93 x and y columns) OR a terra SpatVector object, containing a "vegetation_type" column
    # OUPUTS
    #    a data.frame with a vegetation_type column, and a lab column containing the habitat names and surface
    if(class(habitats) == "SpatVector") {
        habitat_labs <- habitats %>%
                    as.data.frame() %>%
                    mutate(area = terra::expanse(habitats, unit="ha")) %>%
                    group_by(vegetation_type) %>%
                    summarise(area=sum(area)) %>%
                    rename("lab"="area") %>%
                    mutate(lab = paste0(vegetation_type,"\n", round(lab)," ha")) %>%
                    as.data.frame()
    } else if (class(habitats) == "data.frame") {
        pixel_surface = get_pixel_surface(habitats)
        habitat_labs <- as.data.frame(table(habitats$vegetation_type))
        colnames(habitat_labs) <- c("vegetation_type", "lab")
        habitat_labs$lab <- paste0(habitat_labs$vegetation_type,"\n", round(habitat_labs$lab*pixel_surface/10000), " ha")
    }
    return(habitat_labs)
}

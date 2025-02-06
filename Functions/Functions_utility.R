get_main_ID <- function(data) {
    # INPUTS:
    # data : GPS data, with ID field identifying individual trajectories being a combination in the form "animal_ID-subtrajectory_ID"
    # OUTPUT:
    # a vector containing only the animal ID of each row

    return(sapply(as.character(data$ID), function(x) { x = strsplit(x, split = "-")
                          return(x[[1]][1])} ))
}

save_append_replace_IDs <- function(data_save, file_path) {
    # If file_path does not exist yet, save_append_replace will create it and write data_save inside.
    # If file_path exists, all data_save lines will be appended at the end of the file_path data.frame
    # Any line of file_path with an ID present in data_save will be removed from file_path prior to the append of data_save
    # INPUTS :
    #   data_save : a data.frame containing an ID column, and the same columns as the file_path rds object
    #   file_path : path to the rds file the data.frame will be saved to, must contain an ID column

    if (!file.exists(file_path)) {
        saveRDS(data_save, file = file_path)
        print(paste("Data written to newly created", file_path, "file."))
    } else {
        data_file = readRDS(file_path)
        same_IDs_indices = data_file$ID %in% unique(data_save$ID)
        same_IDs = unique(data_file$ID[same_IDs_indices])
        data_file = data_file[!same_IDs_indices, ]
        data_file = rbind(data_file, data_save)
        saveRDS(data_file, file = file_path)
        print(paste("Data appended to", file_path, "file."))
        if (length(same_IDs > 0)) {
            same_IDs = paste(same_IDs, collapse = ', ')
            print(paste0("*** IDs ", same_IDs," were replaced in ", file_path))
        }
    }
}

save_append <- function(data_save, file_path) {
    # If file_path does not exist yet, save_append_replace will create it and write data_save inside.
    # If file_path exists, all data_save lines will be appended at the end of the file_path data.frame
    # INPUTS :
    #   data_save : a data.frame containing the same columns as the file_path rds object
    #   file_path : path to the rds file the data.frame will be saved to

    if (!file.exists(file_path)) {
        saveRDS(data_save, file = file_path)
        print(paste("Data written to newly created", file_path, "file."))
    } else {
        data_file = readRDS(file_path)
        data_file = rbind(data_file, data_save)
        saveRDS(data_file, file = file_path)
        print(paste("Data appended to", file_path, "file."))
    }
}

get_individual_info <- function(ID, individual_info_file, col_name) {
    individual_info <- read.csv(individual_info_file, header=TRUE, sep=",")
    info <- individual_info[as.numeric(sapply(ID, function(id) which(individual_info$Collier==id))), col_name]

    return(info)
}

get_individual_alpage <- function(ID, individual_info_file) {
    return(get_individual_info(ID, individual_info_file, "Alpage"))
}

get_alpage_info  <- function(alpage, alpage_info_file, col_name) {
    alpage_infos <- read.csv(alpage_info_file, header=TRUE, sep=",")
    info<- alpage_infos[alpage_infos$alpage==alpage, col_name]

    return(info)
}

get_pixel_surface   <- function(data) {
    # Compute the pixel surface of a raster contained in a data.frame
    # INPUTS
    #    data : a data.frame containing a raster, with "x" and "y" columns
    # OUPUTS
    #    the pixel surface, in the same unit as x and y

    y = sort(data$x)
    x = sort(data$x)
    pixel_surface = min(abs(diff(unique(x))))*min(abs(diff(unique(y))))

    return(pixel_surface)
}

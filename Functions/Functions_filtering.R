source("Functions/Functions_clean_Bjorneraas.R")


### FUNCTIONS ###
#---------------#

load_catlog_data  <- function(input_file) {
    # INPUT :
    #     input_file : path to a csv file containing a Catlog trajectory
    # OUTPUT :
    #     a data.frame containing trajectory, colnames being c("lat", "lon", "Altitude", "actX", "actY", "date")

    traject <- read.csv(input_file, skip=6, header=TRUE)
    colnames(traject)[c(3, 4)] = c("lat", "lon")
    traject$date = as.POSIXct(paste(traject$Date, traject$Time), tz="GMT", format="%m/%d/%Y %H:%M:%S")
    traject <- traject %>% dplyr::select(-c(Date, Time, Satellites, HDOP, PDOP, Temperature..C., Speed..km.h., TTFF, SNR.AVG, SNR.MAX))
    return(traject)
}

load_followit_data  <- function(input_file) {
    # INPUT :
    #     input_file : path to a csv file containing a Folowit trajectory
    # OUTPUT :
    #     a data.frame containing trajectory, colnames being c("lat", "lon", "Altitude", "actX", "actY", "date")

    traject <- read.csv(input_file, header=TRUE, sep="\t")
    traject <- traject[-c(1),]
    traject$time = as.POSIXct(paste(traject$Date, traject$Time), format = "%Y %m %d %H:%M:%S", tz = "UTC")
    traject <- subset(traject, select = -c(X.1, X.2, X2D3D, TTF, DOP, SVs, FOM, Date, Time))
    colnames(traject) <- c("lat", "lon", "Altitude", "actX", "actY", "date")
    traject$lat = as.numeric(traject$lat)
    traject$lon = as.numeric(traject$lon)
    traject <- traject[!is.na(traject$lat) & !is.na(traject$lon),]
    return(traject)
}

followit_correct_date  <- function(traject, correct_beg_date) {
    date_diff <- date(correct_beg_date) - date(traject$date[1])
    date(traject$date) <- date(traject$date) + date_diff
    return(traject)
}

date_filter <- function(traject, beg_date, end_date) {
    previous_length = length(traject$date)
    traject <- traject %>% filter(between(date, beg_date, end_date))
    new_length = length(traject$date)
    print(paste(previous_length-new_length, "points were removed by date filter.", new_length,"points remain."))
    return(traject)
}

WSG84_speed <- function(i,x) {
    if (i==1) {
        return(NA)
    }
    else {
        distance = 6371008 * acos( sin(pi/180*x$lat[i])*sin(pi/180*x$lat[i-1]) + cos(pi/180*x$lat[i])*cos(pi/180*x$lat[i-1])*cos(pi/180*(x$lon[i-1]-x$lon[i])) )
    
        time = as.numeric(difftime(time1 = x$date[i], time2 = x$date[i-1], units = "secs"))
        return(distance/time)
    }
}

position_filter <- function(traject, medcrit=750, meancrit=500, spikesp=1500, spikecos=(-0.97)) {
    # Using Bjoneraas2010’s script
    
    WSG84_CRS = CRS("+proj=longlat +datum=WGS84")
    L93_CRS = CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs")
    traject <- SpatialPointsDataFrame(traject[,c("lon","lat")],
                                      traject,    #the R object to convert
                                      proj4string = WSG84_CRS)
    
    traject <- spTransform(traject, L93_CRS)
    traject@data$x = traject@coords[,1]
    traject@data$y = traject@coords[,2]
    
    # medcrit et meancrit sont en m
    # spikesp est en m/h
    scr = GPS.screening.wrp(id=traject$ID, x=traject@data$x, y=traject@data$y, 
                            da=traject$date, medcrit, meancrit, spikesp, spikecos)
    
    # GPS.screening.wrp mixes up the order of the individual trajectories (but not the points within each trajectory)
    # therefore the need for the loop
    traject$R1error <- TRUE
    traject$R2error <- TRUE
    for (collar in unique(traject$ID)) {
        print(paste("Nombre de R1errors pour", collar, ":", sum(scr$R1error[scr$id == collar], na.rm=TRUE)))
        print(paste("Nombre de R2errors pour", collar, ":", sum(scr$R2error[scr$id == collar], na.rm=TRUE)))
        traject$R1error[traject$ID==collar] <- scr$R1error[scr$id==collar]
        traject$R2error[traject$ID==collar] <- scr$R2error[scr$id==collar]
    }
    
    
    return(traject@data)
}

filter_one_collar <- function(traject, collar_file, output_rds_file, alpage_name, beg_date, end_date, individual_info_file,
                              bjoneraas.medcrit, bjoneraas.meancrit, bjoneraas.spikesp, bjoneraas.spikecos, sampling_period = 120) {
    # Filters the relocation from one collar based on date and speed/position (Bjoneraas2010).
    # The filtered trajectory is appended to output_rds_file
    # INPUTS
    #    traject : a data.frame containing the raw relocations of the collar, with columns date, lat and lon
    #    collar_ID : the ID of the collar
    #    output_rds_file : the rds file the filtered trajectory should be appended to
    #    alpage_name : name of the alpage
    #    beg_date : the date from which relocations are to be kept
    #    end_date : the date until which relocations are to be kept
    #    bjoneraas.medcrit, bjoneraas.meancrit, bjoneraas.spikesp, bjoneraas.spikecos : parameter of the Bjorneraas relocation errors filter
    #    sampling_period : theoretical time between two consecutive relocations, in seconds
    # OUPUTS
    #    a one-row data.frame of performance indicators of the collar: name (collar ID), worked_until_end (1 if the collar didn’t stop working until 24 hours before beg_date), nloc (number of relocations) and error_perc (percentage of relocations removed by Bjoneraas2010 filter)

    collar_ID = strsplit(collar_file, split = "_")[[1]][1]

    beg_date = as.POSIXct(get_individual_info(collar_ID, individual_info_file, "date_pose"), tz="GMT", format="%d/%m/%Y %H:%M:%S")
    end_date = as.POSIXct(get_individual_info(collar_ID, individual_info_file, "date_retrait"), tz="GMT", format="%d/%m/%Y %H:%M:%S")
    day_prop = as.numeric(gsub(",", ".", get_individual_info(collar_ID, individual_info_file, "proportion_jour_allume"))) # proportion of day with collar switched on

    n_loc_theory = as.numeric(difftime(end_date, beg_date, units = "secs")) * day_prop / sampling_period
    print(paste("Working on", collar_ID, "from", beg_date, "to", end_date))

    indicators = data.frame(name = collar_ID)
    traject$ID <- collar_ID

    # Filter on dates and compute corresponding indicators
    traject <- date_filter(traject, beg_date, end_date)
        # If there is less than one day missing at the end of the time-series, we consider that the collar worked until the end (1)
    indicators$worked_until_end = ifelse(as.numeric(difftime(end_date, traject$date[nrow(traject)], units = "secs")) <= 24*3600, 1, 0)
    indicators$nloc = nrow(traject)

    # Filter on position and speed (Bjoneraas2010) and compute corresponding indicators
    traject <- position_filter(traject, medcrit=bjoneraas.medcrit, meancrit=bjoneraas.meancrit, spikesp=bjoneraas.spikesp, spikecos=bjoneraas.spikecos)
    indicators$R1error = sum(traject$R1error, na.rm=TRUE)
    indicators$R2error = sum(traject$R2error, na.rm=TRUE)
    indicators$localisation_rate = indicators$nloc/n_loc_theory
    indicators$error_perc = (indicators$R1error + indicators$R2error)/nrow(traject)

    traject <- traject[!traject$R1error & !traject$R2error,]
    traject <- traject[!is.na(traject$x), ]
    traject$alpage <- alpage_name
    traject$R1error <- NULL
    traject$R2error <- NULL
    traject$species <- get_individual_info(collar_ID, individual_info_file, "Espece")
    traject$race <- get_individual_info(collar_ID, individual_info_file, "Race")

    colnames(traject)[4] <- "time"

    save_append_replace_IDs(traject, file = output_rds_file)
    return(indicators)
}

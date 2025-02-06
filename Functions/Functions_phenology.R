load_and_perc_by_vegetation_period <- function(pheno, charge) {
        charge_by_season <- data.frame(pheno[, c("x", "y")])
        charge_by_season$growing = 0
        charge_by_season$plateau = 0
        charge_by_season$senesc = 0
        for (d in unique(charge$day)) {
            charge_jour = charge[charge$day==d, ]
            growing = (d >= pheno$t0 & d <= pheno$t1)
            plateau = (d >= pheno$t1 & d <= pheno$t2)
            charge_by_season$growing = charge_by_season$growing + charge_jour$Charge * growing
            charge_by_season$plateau = charge_by_season$plateau + charge_jour$Charge * plateau
            charge_by_season$senesc = charge_by_season$senesc + charge_jour$Charge * (!growing & !plateau)
        }
        charge_by_season$total = charge_by_season$growing + charge_by_season$plateau + charge_by_season$senesc

        charge_by_season = charge_by_season[charge_by_season$total > 0, ]

        charge_by_season$perc_growing = 100 * charge_by_season$growing / charge_by_season$total
        charge_by_season$perc_plateau = 100 * charge_by_season$plateau / charge_by_season$total
        charge_by_season$perc_senesc <- 100 - charge_by_season$perc_growing - charge_by_season$perc_plateau

        return(charge_by_season)
}

get_charge_by_period <- function(charge, breakdays) {
    # Computes the flock load by periods
    # INPUTS
    #    charge : a data.frame with x, y, day and Charge columns
    #    breakdays : the julian days defining the begining and end of periods (the first period being anything before the first breakday,
    #                 the last on anithing after the last breakday)
    # OUPUTS
    #    a data.frame with x, y, day and period columns, periods going from 1 to length(breakdays)+1

    charge %>%
        mutate(period = findInterval(day, breakdays)+1) %>%
        group_by(x, y, period) %>%
        summarize(charge = sum(Charge)) %>%
        as.data.frame() %>%
    return()
}

get_chargement_by_period_and_vegetation_units <- function(charge_by_period, vegetation_units, period_labs) {
    # Computes the number of sheep.days per vegetation units and per period
    # INPUTS
    #    charge_by_period : a data.frame with x, y, day and period columns, periods ranging from 1 to length(period_labs)
    #    vegetation_units : a terra SpatVector containing the vegetation polygons, with a vegetation_type and attribute
    #    period_labs : the name of the different periods, ordered by time
    # OUPUTS
    #    a data.frame with x, y, period, vegetation_type, charge and area (of polygon) columns 
    pixel_surface = get_pixel_surface(charge_by_period)
    test = charge_by_period %>%
        terra::rast() %>%
        terra::extract(y= vegetation_units, fun = sum, na.rm = T, ID = F) %>%
        mutate(vegetation_type=vegetation_units$vegetation_type) %>%
        mutate(area=terra::expanse(vegetation_units, unit="ha")) %>%
        pivot_longer(cols=3:ncol(charge_by_period)-2,
                        names_to='period',
                        values_to='charge') %>%
        mutate(period = period_labs[as.numeric(period)]) %>%
        mutate(period = factor(period, levels = period_labs[length(period_labs):1])) %>%
        mutate(charge = charge*pixel_surface/10000) %>%
        as.data.frame() %>%
    return()
}



###  FIGURES ###
#**************#

plot_presence_perc_by_period_and_habitat <- function(data, period_title, period_labs, habitat_labs, title="") {
    data = data %>%
                count(vegetation_type, period, wt = charge, name = "charge")
    print(ggplot(data, aes(y = vegetation_type, x = charge, fill = period)) +
        geom_col(position = "fill") +
        scale_x_continuous(labels = scales::percent_format()) +
        scale_fill_viridis_d(option = "viridis", direction = -1, period_title,
                            breaks=levels(data$period)[nlevels(data$period):1], labels=period_labs) +
        xlab("Flock presence") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        ggtitle(title))
}

plot_charge_by_period_and_habitat <- function(data, period_title, period_labs, habitat_labs, title="") {
        print(ggplot(data, aes(y = vegetation_type, x = charge/area, fill = period)) +
        geom_boxplot(outlier.shape = NA, varwidth = F) +
        scale_fill_viridis_d(option = "viridis", direction = -1, period_title,
                            breaks=levels(data$period)[nlevels(data$period):1]) +
        xlab("Flock load (sheep.days/ha)") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        ggtitle(title) +
        coord_cartesian(xlim = c(0, quantile(data$charge/data$area, 0.98, na.rm=T)))) # contrairement à scale_x_continuous(xlim =...), coord_cartesian ne supprime pas les données hors graphe avant de calculer les statistiques
}


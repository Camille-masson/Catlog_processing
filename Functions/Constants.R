### PLOTS ###
#-----------#

# momentuHMM’s color palette (modified order 3 1 2 4 5...)
pal <- c("#009E73", "#E69F00", "#56B4E9", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7")
hmm_scale = scale_colour_manual(values = c("Repos" = pal[1], "Paturage" = pal[2], "Deplacement" = pal[3]), drop=F, # drop=F to keep all of the factor values even if one is missing
                                name = "Comportement", guide = guide_legend(order = 1, override.aes = list(size = 1.2)))
hmm_scale_fill = scale_fill_manual(values = c("Repos" = pal[1], "Paturage" = pal[2], "Deplacement" = pal[3]), drop=F, # drop=F to keep all of the factor values even if one is missing
                                name = "Comportement", guide = guide_legend(order = 1, override.aes = list(size = 1.2)))


### GIS ###
#---------#
library(sp)
CRS_L93 <- "epsg:2154"
CRS_WSG84 = "+proj=longlat +datum=WGS84"
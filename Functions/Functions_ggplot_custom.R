library(ggplot2)
theme_set(theme_light()) # theme for ggplot2
library(ggnewscale)
library(magick)

### ADJUSTING LEGEND ###
#**********************#

colourbar_bottom <- function(length = "long") {
  if(length == "short") {
    return( theme(legend.position="bottom",
                  legend.key.height = unit(0.3, "cm"),
                  legend.box.spacing = unit(0, "cm"),
                  legend.title = element_text(vjust=1),
                  legend.text = element_text(angle=45, size=5)) )
  }
  else {
    return( theme(legend.position="bottom",
                  legend.key.height = unit(0.3, "cm"),
                  legend.key.width = unit(dev.size()[1] / 7, "cm"),
                  legend.box.spacing = unit(0, "cm"),
                  legend.title = element_text(vjust=1),
                  legend.text = element_text(size=6)) )
  }
}

colourbar_right <- function() {
    theme(legend.position="right",
          legend.key.height = unit(dev.size()[2] / 7, "cm"),
          legend.key.width = unit(0.3, "cm"),
          legend.box.spacing = unit(0, "cm"),
          legend.text = element_text(size=6) )
}

### SAVING PLOTS ###
#******************#

save_ggplot <- function(filename, plot = last_plot(), width = 7, height =  NA) {
    # Save a ggplot to an image while avoiding white spaces around the pannels when fixed aspect ratio are used in pannels (such as in maps)
    # If one only of width or height is NA, the other dimention will be adapted to avoid white space
    # If both width and height are NA, they will be set to ggsave’s default and there will be no adaptation to avoid white spaces
    # If both width and height are given, there will be no adaptation to avoid white spaces
    # https://stackoverflow.com/questions/59089928/how-to-save-a-ggplot2-graphic-with-the-proper-aspect-ratio/59440180#59440180

    if (is.na(height) & !is.na(width)) {
        height = 2*width
    }
    if (!is.na(height) & is.na(width)) {
        width = 2*height
    }
    ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
    m_png <- image_border(image_trim(image_read(filename)), "white", "30x30")
    image_write(m_png, filename)
}

ggGetAr <- function(p, default.ar=-1){
    # Return the aspect ratio of a NON FACETTED ggplot p
    # https://stackoverflow.com/questions/16422847/save-plot-with-a-given-aspect-ratio

    gb <- ggplot_build(p)
    # first check if theme sets an aspect ratio
    ar <- gb$plot$coordinates$ratio

    # second possibility: aspect ratio is set by the coordinates, which results in 
    # the use of 'null' units for the gtable layout. let's find out
    g <- ggplot_gtable(gb)
    nullw = grid::unitType(g$widths) == "null"
    nullh = grid::unitType(g$heights) == "null"
    # ugly hack to extract the aspect ratio from these weird units
    if(any(nullw))
        ar <- unlist(g$widths[nullw]) / unlist(g$heights[nullh])

    if(is.null(ar)) # if the aspect ratio wasn't specified by the plot
        ar <- default.ar

    ar[1]
}
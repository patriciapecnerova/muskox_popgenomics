if (dir.exists("rEEMSplots")){
    install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
        stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

setwd("muskox/eems/")

library(rEEMSplots)
require(grDevices)
require(graphics)
require(rgeos)
require(Rcpp)
require(RcppEigen)
require(raster)
library("rgdal")
library("rworldmap")
library("rworldxtra")
require(reemsplots2) 
require(ggplot2)

tiles2contours_standardize <- function(tiles, rates, seeds, marks, distm) {
    .Call('_rEEMSplots_tiles2contours_standardize', PACKAGE = 'rEEMSplots', tiles, rates, seeds, marks, distm)
}

tiles2contours <- function(tiles, rates, seeds, marks, distm) {
    .Call('_rEEMSplots_tiles2contours', PACKAGE = 'rEEMSplots', tiles, rates, seeds, marks, distm)
}

map <- rworldmap::getMap(resolution = "high")
map <- broom::tidy(map)

## muskox
eems_results <- file.path("chain1/")
name_figures <- file.path("chain1/eems_muskoxmap_chain1")

coord <- read.table(paste0("muskox.coord"))

eemsPlot <- make_eems_plots(mcmcpath = eems_results, longlat = FALSE, 
                               add_outline = FALSE, add_demes = TRUE, 
                               col_outline = gray, dpi = 300, 
                               add_grid = FALSE, prob_level = 0.9, m_colscale = c(-2.5, 2.5))

eemsPlot2 <- eemsPlot$mrates01  
#    geom_point(data = coord, x =coord$V1, y= coord$V2)
eemsPlot3 <- eemsPlot2 + geom_path(data = map, aes(x = long, y = lat, group = group),
                                         color = "#888888", size = 0.1) +
    coord_map(projection = "orthographic")

plotpath <- file.path(path.expand("chain1/"), "EEMS_muskoxmap_chain1")
ggsave(paste0(plotpath, "-mrates01.png"), eemsPlot3, dpi = 600,width = 6, height = 4)

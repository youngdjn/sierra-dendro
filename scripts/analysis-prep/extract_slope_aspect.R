library(tidyverse)
library(raster)
library(sf)

## Load plots

d = read.csv("data/compiled-for-analysis/plots.csv")
t = read.csv("data/compiled-for-analysis/trees.csv")


d_sf = st_as_sf(d,coords=c("x.plot","y.plot"), crs=3310)
t_sf = st_as_sf(t,coords=c("x","y"),crs=3310)

elev = raster("/home/derek/gis/DEM/CAmerged15_albers.tif")

slope = terrain(elev,opt="slope",unit="degrees")
aspect = terrain(elev,opt="aspect", unit="degrees")

plots_terrain = d_sf %>%
  dplyr::select(plot.id)

plots_terrain$slope = extract(slope,d_sf,method="bilinear")
plots_terrain$aspect = extract(aspect,d_sf,method="bilinear")

trees_terrain = t_sf %>%
  dplyr::select(tree.id)

trees_terrain$slope = extract(slope,t_sf,method="bilinear")
trees_terrain$aspect = extract(aspect,t_sf,method="bilinear")


st_geometry(plots_terrain) = NULL
st_geometry(trees_terrain) = NULL

write.csv(plots_terrain,"data/compiled-for-analysis/plots_terrain.csv", row.names=FALSE)
write.csv(trees_terrain,"data/compiled-for-analysis/trees_terrain.csv", row.names=FALSE)


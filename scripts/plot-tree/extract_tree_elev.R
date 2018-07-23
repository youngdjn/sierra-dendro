## this script adds elevation (if it doesn't already exist) to trees_loc.csv (in data\plot-and-tree\processed)

library(sf)
library(raster)

trees = read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE)

if(is.null(trees$elev)) {
  trees_sp = st_as_sf(trees, coords = c("x","y"), crs = 3310)
  dem = raster("data/non-synced/DEM/CAmerged12_albers2.tif")
  trees$elev = extract(dem,trees_sp,method="bilinear")
}

write.csv(trees,"data/plot-and-tree/processed/trees_loc.csv",row.names=FALSE)

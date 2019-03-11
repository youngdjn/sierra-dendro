setwd("C:/Users/DYoung/Dropbox/Research projects/PSME climate adaptation/Data-analysis")

library(raster)
library(ncdf4)


######################################################
##### Declare global variables and directories #######
######################################################

topowx.dir <- "~/UC Davis/GIS/CA abiotic layers/TopoWx/"
tmin.dir <- "tmin_annual/"
tmax.dir <- "tmax_annual/"

mo.num <- 1:12
mo.chr <- sprintf("%02d",mo.num)

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
ppt.dir <- "ppt_annual/"
td.dir <- "2010 normal means/tdmean/"

vic.dir <- "~/UC Davis/GIS/CA abiotic layers/VIC/"
rad.dir <- "~/UC Davis/GIS/CA abiotic layers/solar rad/rsun nldas adjusted/"
rad.hr.dir <- "~/UC Davis/GIS/CA abiotic layers/solar rad/rsun hours/"
wind.dir <-"~/UC Davis/GIS/CA abiotic layers/wind/"
dem.dir <- "~/UC Davis/GIS/CA abiotic layers/DEM/"
snow.dir <- "~/UC Davis/GIS/CA abiotic layers/snow/"

dobr.pub.dir <- "~/UC Davis/GIS/CA abiotic layers/DobrowskiPublished"





##############################################################
###### Open and stack all monthly-yearly climate layers ######
##############################################################

### TopoWx tmin and tmax monthly-yearly ###

tmax.yr.mo <- list()
tmin.yr.mo <- list()
for(year in years) {
  
  tmax.file <- paste(topowx.dir,tmax.dir,"tmax_",year,".nc",sep="")
  tmax.year <- brick(tmax.file)
  
  tmin.file <- paste(topowx.dir,tmin.dir,"tmin_",year,".nc",sep="")
  tmin.year <- brick(tmin.file)
  
  for(month in 1:12) {
    
    tmax.layer.index <- paste("tmax",year,mo.chr[month],sep=".")
    tmin.layer.index <- paste("tmin",year,mo.chr[month],sep=".")
    
    tmax.yr.mo[tmax.layer.index] <- tmax.year[[month]]
    names(tmax.yr.mo[tmax.layer.index]) <- tmax.layer.index
    tmin.yr.mo[tmin.layer.index] <- tmin.year[[month]]
    names(tmin.yr.mo[tmin.layer.index]) <- tmax.layer.index
  }
}


### VIC water balance ###

# open layers

vic.aet <- brick(paste(vic.dir,"ET.historical.191501-200612.monthly_tot.nc",sep=""))
yearname <- paste("vic.aet",1915:2006,sep=".")
yearname.rep <- as.vector(sapply(yearname,rep,times=12))
yearmoname <- paste(yearname.rep,mo.chr,sep=".")
names(vic.aet) <- yearmoname

vic.pet <- brick(paste(vic.dir,"PET3.historical.191501-200612.monthly_tot.nc",sep=""))
yearname <- paste("vic.pet",1915:2006,sep=".")
yearname.rep <- as.vector(sapply(yearname,rep,times=12))
yearmoname <- paste(yearname.rep,mo.chr,sep=".")
names(vic.pet) <- yearmoname

### Spring snow duration and Dobr pub wb (this is already annualized) ###

spring.snow <- brick(paste(snow.dir,"CA_snowpack_origres.grd",sep=""))

dobr.pub.aet <- list()
dobr.pub.def <- list()

dobr.years <- years[years %in% 1950:2009]

for(year in dobr.years) {
  
  aet.file <- paste(dobr.pub.dir,"/aet_sum_",year,".tif",sep="")
  def.file <- paste(dobr.pub.dir,"/def_sum_",year,".tif",sep="")
  
  aet.name <- paste("dobr.pub.aet.",year,sep="")
  def.name <- paste("dobr.pub.def.",year,sep="")
  
  dobr.pub.aet[[aet.name]] <- raster(aet.file) 
  dobr.pub.def[[def.name]] <- raster(def.file)
  
}

dobr.pub.aet.brick <- brick(dobr.pub.aet)
dobr.pub.def.brick <- brick(dobr.pub.def)



# make a stack of the VIC layers just from the years and months desired for analysis
years.vic <- 1915:2006
years.keep <- years[years %in% years.vic]
keepyearname <- paste("vic.aet",years.keep,sep=".")
keepyearname.rep <- as.vector(sapply(keepyearname,rep,times=12))
keepyearmoname <- paste(keepyearname.rep,mo.chr,sep=".") 
nkeepyearmo <- length(keepyearmoname)

vic.aet.list <- list()
for(i in 1:nkeepyearmo) {
  vic.aet.list[[keepyearmoname[i]]] <- vic.aet[[keepyearmoname[i]]]
}

keepyearname <- paste("vic.pet",years.keep,sep=".")
keepyearname.rep <- as.vector(sapply(keepyearname,rep,times=12))
keepyearmoname <- paste(keepyearname.rep,mo.chr,sep=".") 

nkeepyearmo <- length(keepyearmoname)

vic.pet.list <- list()
for(i in 1:nkeepyearmo) {
  vic.pet.list[[keepyearmoname[i]]] <- vic.pet[[keepyearmoname[i]]]
}


### PRISM precip monthly-yearly ###

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
ppt.dir <- "ppt_annual/"

ppt.yr.mo <- list()

for(year in years) {
  
    if(year > 1980) {
      ppt.yr.files <- paste(prism.dir,ppt.dir,"PRISM_ppt_stable_4kmM3_",as.character(year),mo.chr,"_bil.bil",sep="")
    } else {
      ppt.yr.files <- paste(prism.dir,ppt.dir,"PRISM_ppt_stable_4kmM2_",as.character(year),mo.chr,"_bil.bil",sep="")
    }
    
    ppt.yr <- stack(ppt.yr.files)
    
    for(month in 1:12) {
      
      layer.index <- paste("ppt",year,mo.chr[month],sep=".")
      ppt.yr.mo[layer.index] <- ppt.yr[[month]]
    }
}



### PRISM dewpoint monthly normals (do not have monthly-yearly) ###
td.normal.mo <- list()
td.files <- paste(prism.dir,td.dir,"tdmean_1981_2010_normal_",mo.chr,".tif",sep="")
td.normal <- stack(td.files)
for(month in 1:12) {
  layer.index <- paste("td.normal",mo.chr[month],sep=".")
  td.normal.mo[layer.index] <- td.normal[[month]]
}


### Radiation
rad.file <- paste(rad.dir,"glob_rad_monthly_dobr1_int.tif",sep="")
rad.normal.mo <- list()
rad.normal <- stack(rad.file)

for(month in 1:12) {
  layer.index <- paste("rad.normal",mo.chr[month],sep=".")
  rad.normal.mo[layer.index] <- rad.normal[[month]]
}


### Rad hours
rad.hr.file <- paste(rad.hr.dir,"rad_hr_monthly_30m.tif",sep="")
rad.hr.normal.mo <- list()
rad.hr.normal <- stack(rad.hr.file)

for(month in 1:12) {
  layer.index <- paste("rad.hr.normal",mo.chr[month],sep=".")
  rad.hr.normal.mo[layer.index] <- rad.hr.normal[[month]]
}



  
### Wind
wind.file <- paste(wind.dir,"NLDAS_wind_snclip_monthly.tif",sep="")
wind.normal.mo <- list()
wind.normal <- stack(wind.file)

for(month in 1:12) {
  layer.index <- paste("wind.normal",mo.chr[month],sep=".")
  wind.normal.mo[layer.index] <- wind.normal[[month]]
}

### elevation
dem.file <- paste(dem.dir,"CAmerged9_SNclip_nd.tif",sep="")
dem <- raster(dem.file)


### stack everything of the same resolution and extent into one stack
temp.800m.rast <- unlist(list(tmin.yr.mo,tmax.yr.mo))
ppt.4km.rast <- unlist(list(ppt.yr.mo))
td.4km.rast <- unlist(list(td.normal.mo))
rad.rast <- unlist(list(rad.normal.mo))
rad.hr.rast <- unlist(list(rad.hr.normal.mo))
wind.rast <- unlist(list(wind.normal.mo))
vic.aet.rast <- unlist(list(vic.aet.list))
vic.pet.rast <- unlist(list(vic.pet.list))
snow.rast <- spring.snow
dobr.pub.aet.rast <- dobr.pub.aet.brick
dobr.pub.def.rast <- dobr.pub.def.brick

##### Extract climate values for points

#### Load tree data with locations ####
trees <- read.csv("~/../Dropbox/Research projects/PSME climate adaptation/Data-analysis/Processed data/trees_loc5.csv",header=TRUE)
albers <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
trees <- SpatialPointsDataFrame(coords=cbind(trees$x,trees$y),data=data.frame(tree.id=trees[,"tree.id"]),proj4string=albers)


#### Interpolate weather/climate vals from grids ####

# define utility function #
extract.single <- function(layer, trees) {
  cat("\rInterpolating layer: ",names(layer))
  vals <- extract(layer,trees,method="bilinear")
  return(vals)
}

# perform bilinear interpolation
ppt.4km <- sapply(X = ppt.4km.rast,FUN = extract.single, trees = trees)
td.4km <- sapply(X = td.4km.rast,FUN = extract.single, trees = trees)
temp.800m <- sapply(X = temp.800m.rast,FUN = extract.single, trees = trees)
rad <- sapply(X = rad.rast,FUN = extract.single, trees = trees)
rad.hr <- sapply(X = rad.hr.rast,FUN = extract.single, trees = trees)
wind <- sapply(X = wind.rast,FUN = extract.single, trees = trees)
vic.aet <- sapply(X = vic.aet.rast,FUN = extract.single, trees = trees)
vic.pet <- sapply(X = vic.pet.rast,FUN = extract.single, trees = trees)
elev <- extract(dem,trees, method="bilinear")
snow <- extract(snow.rast,trees, method="bilinear")
dobr.pub.aet <- extract(dobr.pub.aet.rast,trees,method="bilinear")
dobr.pub.def <- extract(dobr.pub.def.rast,trees,method="bilinear")


# #########THIS PART IS FOR DOWNSCALING USING GIDS, IF DESIRED (WILL NEED TO SPEED UP THE TOPOWX DOWNSCALING)
#  #!##### need to extract point elevation to elev attribute
# 
#
# #!update topowx dem
# dem.prism.ppt <- raster("~/UC Davis/GIS/CA abiotic layers/PRISM/dem_clipped_prism4km.tif")
# dem.prism.td <- raster("~/UC Davis/GIS/CA abiotic layers/PRISM/dem_clipped_prism4km.tif")
# dem.topowx <- raster("~/UC Davis/GIS/CA abiotic layers/PRISM/dem_clipped_prism800m.tif")
# dem.prism.ppt <- projectRaster(dem.prism.ppt,ppt.4km.rast[[1]])
# dem.prism.td <- projectRaster(dem.prism.td,td.4km.rast[[1]])
# dem.topowx <- projectRaster(dem.topowx,temp.800m.rast[[1]])
# 
# downscale.single <- function(layer, plots, dem, res) {
#   cat("\rDownscaling layer: ",names(layer))
#   vals <- gids.spatial(layer,dem,plots,res)
#   return(vals)
# }
# 
# 
# ppt.4km <- sapply(X = ppt.4km.rast,FUN = downscale.single, plots = plots, dem=dem.prism.ppt, res=4000)
# td.800m <- sapply(X = td.800m.rast,FUN = downscale.single, plots = plots, dem=dem.prism.td, res=4000)
# temp.800m <- sapply(X = temp.800m.rast,FUN = downscale.single, plots = plots, dem=dem.topowx, res=800)
# 


## compute tmean
tmax.cols <- grep("tmax",colnames(temp.800m))
tmin.cols <- grep("tmin",colnames(temp.800m))
temp.yearmo <- substr(colnames(temp.800m[,tmax.cols]),6,12)
tmax <- temp.800m[,tmax.cols]
tmin <- temp.800m[,tmin.cols]

tmean <- (tmax + tmin) / 2
colnames(tmean) <- paste("tmean",temp.yearmo,sep=".")

## VIC deficit
vic.def <- vic.pet - vic.aet
colnames(vic.def) <- gsub("pet","def",colnames(vic.def))

## extract latitude
nad83geo <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
trees.geo <- spTransform(trees,nad83geo)
lat <- coordinates(trees.geo)[,2]



clim <- data.frame(temp.800m,tmean,ppt.4km,td.4km,rad,rad.hr,wind,elev,lat,vic.aet,vic.def,snow,dobr.pub.aet,dobr.pub.def)
trees.clim <- data.frame(trees,clim)

#! possibly want to restrict to complete.cases in case any of the climate variables is NA for any of the trees
# (to make fair comparison between methods)


write.csv(trees.clim,"Dendro/Climate correlations/tree_climate_extracted.csv",row.names=FALSE)




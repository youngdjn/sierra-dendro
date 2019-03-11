## could make this script export monthly water balance values by changing the parallelized function (removing the summation by year) and then during the re-conversion to raster produce one raster layer per year and month

source("Dendro/Climate correlations/R/utils.R")
source("Dendro/Climate correlations/R/thorn_wil_monthly.R")
library("snowfall")
library(reshape)

run.tw.wb <- function(clim.df) {

  
  temp.vars <- grep("tmean",names(clim.df))
  precip.vars <- grep("ppt",names(clim.df))
  
  if(length(temp.vars) != length(precip.vars)) stop("ThornWil WB not provided equal number of months fro temp and precip.")
  
  nyears <- length(temp.vars)/12
  if(round(nyears) != nyears) stop("Error in ThornWill WB: Must supply a multiple of 12 months of weather data.")
  
  
  run.thornwil <- function(i) { #i from 1:nrow(clim.df)
  
    clim.tree <- clim.df[i,]
    
    tree.precip <- t(clim.tree[,precip.vars])
    tree.temp <- t(clim.tree[,temp.vars])
    
    tree.df <- data.frame("precip"=tree.precip,"temp"=tree.temp)
    names(tree.df) <- c("precip","temp")
    lat <- clim.tree$lat
    
    tree.df$year <- as.numeric(substr(rownames(tree.df),start=5,stop=8))
    tree.df$month <- as.numeric(substr(rownames(tree.df),start=10,stop=11))
  
    
    #replicate first year of weather 5 times to spin up model
    first.year <- tree.df[1:12,]
    tree.df <- rbind(first.year,first.year,first.year,first.year,first.year,tree.df)
    
    wb.params <- list(T.m = tree.df$temp, P.m = tree.df$precip, month=tree.df$month, L = lat, S.max = 150, soil.start=0)

    tree.df$pet <- PET.Thorn(wb.params)
    tree.df$aet <- AET.Wil(tree.df$pet,wb.params,soil.start=150)
    tree.df$def <- tree.df$pet-tree.df$aet
    
    
    
    #remove the first 5 years (spin-up)
    tree.df <- tree.df[61:nrow(tree.df),]
    
    tree.df$tree.id <- clim.tree$tree.id
    
    # sum by water year
    nyear <- (nrow(tree.df)/12) -1 # subtract 1 because we need water years
    yearstarts <- (0:(nyear-1)) *12 + 1
    yearends <- yearstarts + 11
    yearmonths <- mapply(function(x, y) x:y, yearstarts, yearends)
    yearmonths <- yearmonths + 9 # add 9 so it starts with october, not january
    
    pet <- mapply(function(x,y) sum(tree.df$pet[x:y]),yearstarts,yearends)
    aet <- mapply(function(x,y) sum(tree.df$aet[x:y]),yearstarts,yearends)
    def <- mapply(function(x,y) sum(tree.df$def[x:y]),yearstarts,yearends)
    
    year.min <- min(tree.df$year)
    year.ending <- (year.min + 1) : (year.min + nyear)
  
    ID <- as.character(clim.tree$tree.id)
    
    tree.wb <- data.frame(tree.id=ID,year.ending,pet,aet,def)
    
    return(tree.wb)
  }
  
  
  i <- 1:nrow(clim.df)
  
  sfInit(parallel=TRUE,cpus=4)
  sfExport("clim.df","precip.vars","temp.vars")
  sfSource("Dendro/Climate correlations/R/utils.R")
  sfSource("Dendro/Climate correlations/R/thorn_wil_monthly.R")
  results <- sfLapply(i,fun=run.thornwil)
  sfStop()
  
  
  ## should take 3-4 min on 4 cores
  
  
  
  
  #save back to rasters
  trees.wb <- do.call("rbind",results)
  trees.wb <- as.data.frame(trees.wb)
  
  tw.aet <- cast(trees.wb,tree.id~year.ending,value=c("aet"))
  tw.def <- cast(trees.wb,tree.id~year.ending,value=c("def"))
  
  #remove tree.id column
  tw.aet <- tw.aet[,-1]
  tw.def <- tw.def[,-1]
  
  
  names(tw.aet) <- paste("tw.aet.",names(tw.aet),sep="")
  names(tw.def) <- paste("tw.def.",names(tw.def),sep="")
  
  
  
  return(list(tw.aet=tw.aet, tw.def=tw.def))
  
}


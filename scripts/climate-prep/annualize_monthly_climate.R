#### Annualizes monthly climate by water year

annualize.clim <- function(trees.clim) { 
  
  
  first.year <- years[1]
  last.year <- years[length(years)]
  
  nyears <- last.year - first.year + 1
  if(nyears < 3) {stop("You must supply at least 3 years.")}
  n.water.years <- nyears - 1
  inter.years <- years[2:(nyears-1)]
  inter.years.rep <- rep(inter.years,each=12)
  
  mo.chr <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  first.year.mos <- c("10","11","12")
  last.year.mos <- c("01","02","03","04","05","06","07","08","09")
  inter.mos <- rep(mo.chr,nyears-2)
  
  first.year.mo.yr <- paste(first.year,first.year.mos,sep=".")
  last.year.mo.yr <- paste(last.year,last.year.mos,sep=".")
  inter.mos.mo.yr <- paste(inter.years.rep,inter.mos,sep=".")
  all.mo.yr <- c(first.year.mo.yr,inter.mos.mo.yr,last.year.mo.yr)
  
  all.mo.yr.tmin <- paste("tmin",all.mo.yr,sep=".")
  all.mo.yr.tmax <- paste("tmax",all.mo.yr,sep=".")
  all.mo.yr.ppt <- paste("ppt",all.mo.yr,sep=".")
  all.mo.yr.vic.aet <- paste("vic.aet",all.mo.yr,sep=".")
  all.mo.yr.vic.def <- paste("vic.def",all.mo.yr,sep=".")

  tmin.col.nums <- sapply(all.mo.yr.tmin,function(x) {grep(x,names(trees.clim))})
  tmax.col.nums <- sapply(all.mo.yr.tmax,function(x) {grep(x,names(trees.clim))})
  ppt.col.nums <- sapply(all.mo.yr.ppt,function(x) {grep(x,names(trees.clim))})
  
  vic.aet.col.nums <- sapply(all.mo.yr.vic.aet,function(x) {grep(x,names(trees.clim))})
  vic.def.col.nums <- sapply(all.mo.yr.vic.def,function(x) {grep(x,names(trees.clim))})
  

  tmin.col.nums.mat <- matrix(tmin.col.nums,ncol=12,byrow=TRUE)
  tmax.col.nums.mat <- matrix(tmax.col.nums,ncol=12,byrow=TRUE)
  ppt.col.nums.mat <- matrix(ppt.col.nums,ncol=12,byrow=TRUE)
  vic.aet.col.nums.mat <- matrix(vic.aet.col.nums,ncol=12,byrow=TRUE)
  vic.def.col.nums.mat <- matrix(vic.def.col.nums,ncol=12,byrow=TRUE)

  #! ### This may produce year delayed by one?
  water.year.summary.annual <- matrix(nrow=nrow(trees.clim),ncol=0)
  for (i in 1:n.water.years) {
    
    water.year.ending <- years[i]+1
    
    tmin.col.nums.year <- tmin.col.nums.mat[i,]
    tmin.cols.year <- trees.clim[,tmin.col.nums.year]
    
    tmax.col.nums.year <- tmax.col.nums.mat[i,]
    tmax.cols.year <- trees.clim[,tmax.col.nums.year]
    
    ppt.col.nums.year <- ppt.col.nums.mat[i,]
    ppt.cols.year <- trees.clim[,ppt.col.nums.year]
    
    tmin.avg.year <- apply(tmin.cols.year,1,mean)
    tmax.avg.year <- apply(tmax.cols.year,1,mean)
    ppt.tot.year <- apply(ppt.cols.year,1,sum)
    
    tmean.avg.year <- (tmin.avg.year+tmax.avg.year)/2
    
    
    vic.aet.col.nums.year <- vic.aet.col.nums.mat[i,]
    vic.def.col.nums.year <- vic.def.col.nums.mat[i,]
    
    #only comput vic water year if data exist in that year
    if(length(unlist(vic.aet.col.nums.year)) == 12) {
      vic.aet.cols.year <- trees.clim[,unlist(vic.aet.col.nums.year)]
      vic.def.cols.year <- trees.clim[,unlist(vic.def.col.nums.year)]
      vic.aet.tot.year <- apply(vic.aet.cols.year,1,sum)
      vic.def.tot.year <- apply(vic.def.cols.year,1,sum)
    } else {
      vic.aet.tot.year <- rep(NA,nrow(trees.clim))
      vic.def.tot.year <- rep(NA,nrow(trees.clim))
    }
    

    ### merge into DF
    colname.prefixes <- c("tmin","tmax","tmean","ppt","vic.aet","vic.def")
    colnames <- paste(colname.prefixes,water.year.ending,sep=".")
    
    clim.year <- cbind(tmin.avg.year,tmax.avg.year,tmean.avg.year,ppt.tot.year,vic.aet.tot.year,vic.def.tot.year)
    colnames(clim.year) <- colnames
    
    water.year.summary.annual <- cbind(water.year.summary.annual,clim.year)
    
  }
  water.year.summary.annual <- as.data.frame(water.year.summary.annual)
  water.year.summary.annual$tree.id <- trees.clim$tree.id
  
  return(water.year.summary.annual)
  
}
  
  
if(FALSE) {
  
  ### pull out temp and precip
  temp <- t(water.year.summary.annual[,grep("tmean",colnames(water.year.summary.annual))])
  colnames(temp) <- trees.clim$ID
  temp <- as.data.frame(temp)
  temp$year <- as.numeric(substr(rownames(temp),7,10))
  
  precip <- t(water.year.summary.annual[,grep("ppt",colnames(water.year.summary.annual))])
  colnames(precip) <- trees.clim$ID
  precip <- as.data.frame(precip)
  precip$year <- as.numeric(substr(rownames(precip),5,8))
  
  vic.aet <- t(water.year.summary.annual[,grep("vic.aet",colnames(water.year.summary.annual))])
  colnames(vic.aet) <- trees.clim$ID
  vic.aet <- as.data.frame(vic.aet)
  vic.aet$year <- as.numeric(substr(rownames(vic.aet),9,12))
  
  vic.def <- t(water.year.summary.annual[,grep("vic.def",colnames(water.year.summary.annual))])
  colnames(vic.def) <- trees.clim$ID
  vic.def <- as.data.frame(vic.def)
  vic.def$year <- as.numeric(substr(rownames(vic.def),9,12))
}
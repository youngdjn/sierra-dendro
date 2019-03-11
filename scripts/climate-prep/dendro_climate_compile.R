setwd("C:/Users/DYoung/Dropbox/Research projects/PSME climate adaptation/Data-analysis")

################################
########## Set global variables
################################
years <- 1956:2014  # set years of analysis



###################################
######### 1) Extract monthly climate variables (and normals for variables that have no annuals)
###################################

source("Dendro/Climate correlations/R/dendro_extract_climate.R") ###!!! also other problems, but problem opening vic layers? (invalid layer names, only when year >2006)

write.csv(trees.clim,"Dendro/Climate correlations/tree_climate_extracted.csv",row.names=FALSE)

#If it was run on the most recent data already, just open it
trees.clim <- read.csv("Dendro/Climate correlations/tree_climate_extracted.csv",header=TRUE)




##################################
########## 2) Annualize monthly raw climate
##################################

source("Dendro/Climate correlations/R/annualize_monthly_climate.R")
trees.clim.ann <- annualize.clim(trees.clim) #adds a tree.id column

#add snowpack (already annualized)
all.yr.snow <- paste("snow.late",2001:2015,sep=".")
all.yr.snow.colnums <- sapply(all.yr.snow,function(x) {grep(x,names(trees.clim))})
all.yr.snow <- trees.clim[,all.yr.snow.colnums]

#add dobr pub aet, def
yr.dobr.pub <- years[years %in% 1950:2009]
all.yr.dobr.pub.aet <- paste("dobr.pub.aet",yr.dobr.pub,sep=".")
all.yr.dobr.pub.def <- paste("dobr.pub.def",yr.dobr.pub,sep=".")
all.yr.dobr.pub.aet.colnums <- sapply(all.yr.dobr.pub.aet,function(x) {grep(x,names(trees.clim))})
all.yr.dobr.pub.def.colnums <- sapply(all.yr.dobr.pub.def,function(x) {grep(x,names(trees.clim))})
all.yr.dobr.pub.aet <- trees.clim[,all.yr.dobr.pub.aet.colnums]
all.yr.dobr.pub.def <- trees.clim[,all.yr.dobr.pub.def.colnums]

trees.clim.ann <- cbind(trees.clim.ann,all.yr.snow,all.yr.dobr.pub.aet,all.yr.dobr.pub.def)


##################################
######### 3) Compute (including annualize) Dobrowski AET and Deficit from monthly raw climate
##################################

source("Dendro/Climate correlations/R/dobrowski_wb.R")
dobr.wb <- run_dob_wb_timeseries(trees.clim,years) #computes for all water years within the calendar year range provided (which must be subset of climate data provided)
dobr.aet <- dobr.wb$AET.Dobr
dobr.def <- dobr.wb$Deficit.Dobr



##################################
######### 4) Compute (including annualize) ThornWil AET and Deficit from monthly raw climate
##################################

source("Dendro/climate correlations/R/run_thorn_wil_monthly.R")
tw.wb <- run.tw.wb(trees.clim) #computes for all water years possible given the years of the data supplied
tw.aet <- tw.wb$tw.aet
tw.def <- tw.wb$tw.def



##################################
######### 6) Combine water balance values with annualized non-wb climate vals
##################################


trees.clim.ann.wb <- data.frame(trees.clim.ann,dobr.aet,dobr.def,tw.aet,tw.def)

write.csv(trees.clim.ann.wb,"Dendro/Climate correlations/trees_clim_ann_wb.csv",row.names=FALSE)

trees.clim.ann.wb <- read.csv("Dendro/Climate correlations/trees_clim_ann_wb.csv",header=TRUE)



##################################
######### 7) Create data frame with normal, sd, and lagged (raw, z, and anom)
##################################

#### First step, put all climate data into long-form DF
d <- melt(trees.clim.ann.wb,id.vars="tree.id")
d$variable <- as.character(d$variable)

## Separate the year and variable names

d$year <- substr(d$variable,start=nchar(d$variable)-3,stop=nchar(d$variable))
d$variable2 <- substr(d$variable,start=1,stop=nchar(d$variable)-5)

#see all variable names
unique(d$variable2)

#get rid of tmax and tmin
d <- d[(d$variable2 != "tmin") & (d$variable2 != "tmax"),]

#get rid of the original varibale (mashup) column
d <- d[,names(d) != "variable"]

#see all variable names
unique(d$variable2)

# convert back to wide format with one column for each weather variable
d.w <- cast(d,year+tree.id~variable2)

# compute normals for defined year range (note that not all datasets fully span the same date range)
norm.years <- 1975:2015
normal <- aggregate(d.w[d.w$year %in% norm.years,],by=list(d.w[d.w$year %in% norm.years,]$tree),FUN=mean,na.rm=TRUE)
normal <- normal[,-c(2,3)]
names(normal) <- c(names(normal)[1],paste(names(normal[,-1]),"normal",sep="."))
d.w.norm <- merge(d.w,normal,by.x="tree.id",by.y="Group.1",all.x=TRUE)


# compute sds for defined year range (note that not all datasets fully span the same date range)
sds <- aggregate(d.w[d.w$year %in% norm.years,],by=list(d.w[d.w$year %in% norm.years,]$tree),FUN=sd,na.rm=TRUE)
sds <- sds[,-c(2,3)]
names(sds) <- c(names(sds)[1],paste(names(sds[,-1]),"sd",sep="."))
d.w.base <- merge(d.w.norm,sds,by.x="tree.id",by.y="Group.1",all.x=TRUE)



# for each variable, compute anomalies and z-scores and proportion of normal
weather.vars <- c("dobr.aet","dobr.def","ppt","snow.late","tmean","tw.aet","tw.def","vic.aet","vic.def","dobr.pub.aet","dobr.pub.def")

for(var in weather.vars) {
 
  var.normal <- paste0(var,".normal")
  var.sd <- paste0(var,".sd")
  var.z <- paste0(var,".z")
  var.anom <- paste0(var,".a")
  var.prop <- paste0(var,".prop")
  var.logprop <- paste0(var,".logprop")
  
  d.w.base[,var.z] <- (d.w.base[,var] - d.w.base[,var.normal]) / d.w.base[,var.sd]
  d.w.base[,var.anom] <- (d.w.base[,var] - d.w.base[,var.normal])
  d.w.base[,var.prop] <- (d.w.base[,var] / d.w.base[,var.normal])
  d.w.base[,var.logprop] <- log(d.w.base[,var] / d.w.base[,var.normal])

  
}

# compute sds for proportion departure from normal (note that not all datasets fully span the same date range)
prop.dep.norm.cols <- grep("prop",names(d.w.base))
sds <- aggregate(d.w.base[d.w.base$year %in% norm.years,prop.dep.norm.cols],by=list(d.w.base[d.w.base$year %in% norm.years,]$tree),FUN=sd,na.rm=TRUE)
names(sds) <- c(names(sds)[1],paste(names(sds[,-1]),"sd",sep="."))
d.w.base <- merge(d.w.base,sds,by.x="tree.id",by.y="Group.1",all.x=TRUE)


# for each variable, compute z-score for proportion of normal #! be careful because this is based on log deviation from normal (average), not log deviation from the mean of all the deviations--what if deviation above normal are larger than deviations below normal?
weather.vars <- c("dobr.aet","dobr.def","ppt","snow.late","tmean","tw.aet","tw.def","vic.aet","vic.def","dobr.pub.def","dobr.pub.aet")

for(var in weather.vars) {
  
  var.logprop.sd <- paste0(var,".logprop.sd")
  var.logprop <- paste0(var,".logprop")
  var.logprop.z <- paste0(var,".logprop.z")
  
  var.prop.sd <- paste0(var,".prop.sd")
  var.prop <- paste0(var,".prop")
  var.prop.z <- paste0(var,".prop.z")
  
  d.w.base[,var.logprop.z] <- d.w.base[,var.logprop] / d.w.base[,var.logprop.sd]
  d.w.base[,var.prop.z] <- d.w.base[,var.prop] / d.w.base[,var.prop.sd]
  
}


  

# for each variable, for each tree id, grab lagged climate
weather.vars.z <- paste0(weather.vars,".z")
weather.vars.a <- paste0(weather.vars,".a")
weather.vars.prop <- paste0(weather.vars,".prop")
weather.vars.prop.z <- paste0(weather.vars,".prop.z")
weather.vars.logprop <- paste0(weather.vars,".logprop")
weather.vars.logprop.z <- paste0(weather.vars,".logprop.z")
weather.vars.ext <- c(weather.vars,weather.vars.z,weather.vars.a,weather.vars.prop,weather.vars.prop.z,weather.vars.logprop,weather.vars.logprop.z)
tree.ids <- unique(d.w.base$tree.id)

d.w.base <- d.w.base[order(d.w.base$tree.id,d.w.base$year),]

shifts <- 1:3

for(tree.id in tree.ids) {
  
  #find what rows apply to this tree
  tree.rows <- which(d.w.base$tree == tree.id)

  for(var in weather.vars.ext) {

    #values of variable unlagged
    vals.orig <- d.w.base[tree.rows,var]
    
    for(shift in shifts) {
    
      #values of variable lagged
      vals.lagged <- c(rep(NA,shift),vals.orig[1:(length(vals.orig)-shift)])
      
      #lag name
      var.lag <- paste0(var,shift)
      
      d.w.base[tree.rows,var.lag] <- vals.lagged
    }
  }
}

## add radiation
# get cols
rad.cols <- grep("rad.normal",names(trees.clim))
rad.tot <- rowSums(trees.clim[,rad.cols])
rad.03 <- trees.clim[,"rad.normal.03"]
rad.06 <- trees.clim[,"rad.normal.06"]
rad.hr.03 <- trees.clim[,"rad.hr.normal.03"]
rad.hr.06 <- trees.clim[,"rad.hr.normal.06"]
id <- trees.clim[,"tree.id"]
rad.df <- data.frame(rad.tot=rad.tot,rad.03=rad.03,rad.06=rad.06,rad.hr.03=rad.hr.03,rad.hr.06=rad.hr.06,tree.id=id)

d.w.base <- merge(d.w.base,rad.df,by="tree.id",all.x=TRUE)

write.csv(d.w.base,"Dendro/Climate correlations/tree_clim_full.csv",row.names=FALSE)

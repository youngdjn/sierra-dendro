#setwd("~/UC Davis/Research projects/Deficit project/For running at RSL")


### This is the script published by Dobrowski, with additions by Derek Young at the end



# This script contains 4 functions used to model ET0 and water balance:
# 1. 'snowmod' estimates snowfall and snowpack and net moisture input as a function of temperature, precip, and existing
#   snowpack.  It also outputs a vector of albedo values, generally 0.2 if there is no snow, or 0.8 if there is snow.
# 2. 'monthlyETO' for calculating monthly reference evapotranspiration
# 3. 'dailyET0' for calculating daily reference evapotranspiration
# 4. 'aetmod' estimates actual et, deficit, soil moisture and runoff as a function of moisture input, existing
#   soil moisture, and soil water capacity. 
#
# Author: Alan Swanson 2012
###############################################################################

snowmod <- function(tmean,ppt,radiation=NULL,snowpack_prev=NULL,albedo=0.23,albedo_snow=0.8){
  # This function computes monthly estimated snowfall and snowmelt. Output includes end-of-month snowpack,
  # water "input" (snowmelt plus rain), and albedo.
  # Arguments:
  #  tmean - vector of mean monthly temperatures
  #  radiation - vector of shortwave solar radiation in MJ/m^2/day.
  #  snowpack_prev - vector of snowpack at the beginning of the month.  If NULL this is 
  #    taken to be zero.
  #  albedo - a single value for albedo in the absence of snow cover.  
  #  albedo_snow - single value of albedo given snow cover
  #
  # Value:  dataframe with three columns for end-of-month snowpack, H2O input (rain plus snowmelt),
  #  and albedo.
  
  
  N <- length(tmean)
  if(is.null(snowpack_prev)) snowpack_prev <- rep(0,N)
  snowpack <- rep(NA,N)
  input <- rep(NA,N)
  
  # this is for radiation in MJ/m^2/day
  mf <- function(t,t0,t1) pmin(pmax((t-t0)/(t1-t0),0),1)
  linrmelt <- function(temp,radiation,b0,b1,b2) pmax((b0+temp*b1+radiation*b2),0)
  parvec <- c(-4.604,6.329,-398.4,81.75,25.05)
  mfsnow <- mf(tmean,parvec[1],parvec[2])
  mfmelt <- linrmelt(tmean,radiation,parvec[3],parvec[4],parvec[5])
  
  # calculate values
  snow <- (1-mfsnow)*ppt
  rain <- mfsnow*ppt	
  melt <- pmin(mfmelt,snow+snowpack_prev) 
  snowpack <- snowpack_prev+snow-melt 
  input <-rain+melt
  
  # make vector of albedo values
  albedo <- rep(albedo,N)
  albedo[snowpack>0 | (snowpack_prev>0)] <- albedo_snow
  
  return(data.frame(snowpack=snowpack,input=input,albedo=albedo))
}


monthlyET0 <- function(radiation,tmax,tmin,wind,lat,elev,dpt,tmean_prev,albedo=0.23,month){
  # This function runs Reference ET estimates for monthly timesteps using methods based on
  # the Penman-Montieth equation as presented in Allen et al (1998).  It incorporates a 
  # modification which adjusts stomatal conductance downwards at temperatures below 5 C.
  #
  # Arguments:
  # radiation: vector of monthly average shortwave radiation in MJ/m^2/day
  # tmax, tmin: vectors of monthly average maximum and minimum temperatures in C, 
  # wind: vector of monthly average wind speed in m/s at 10m above ground, 
  # tmean_prev: vector of mean temp for the previous month, 
  # lat: vector of latitude in degrees 
  # elev: vector of elevation in meters, 
  # dpt: vector of dewpoint temperature in C.
  # tmean_prev: vector of mean temp of previous month in C
  # albedo: vector or scalar of albedo values, 
  # month: scalar 1-12.
  
  #
  # Value: 
  # Returns a vector of ET0 values.
  
  t0<-unclass(Sys.time())	
  daysinmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
  d2=c(31,59,90,120,151,181,212,243,273,304,334,365)
  d1=c(1,32,60,91,121,152,182,213,244,274,305,335)
  DoY <- (d1[month]+d2[month])/2 # use middle day of month to represent monthly average. 
  n_days <- daysinmonth[month]
  
  # calculate soil heat flux (total for the month) using change in temperature from previous month
  tmean <- (tmax+tmin)/2 
  G <- 0.14*(tmean-tmean_prev) # fixed from previous version
  
  # convert to wind height at 2m
  hw=10 # height of wind measurements 
  wind <- wind*(4.87/log(67*hw-5.42))  # convert to wind height at 2m
  
  # stomatal conductance adjustment for low temperatures
  sr=100 # stomatal resistance sec/m
  ks_min=.01 # minimum value for temps below T1
  Tl=-10       # minimum temp (sc goes to ks_min below this temp)
  T0=5		# optimal temp
  Th=100     # maximum temp (sc goes to zero above this)
  thresh=5   # temperature threshold below which to apply Jarvis equation (ks=1 above this temp)
  b4 <- (Th-T0)/(Th-Tl)
  b3 <- 1/((T0-Tl)*(Th-T0)^b4)
  ks <- pmax(pmin(b3*(tmean-Tl)*(Th-tmean)^b4,1),ks_min)
  ks[is.na(ks)] <- ks_min
  ks[tmean>=thresh] <- 1
  
  # convert to stomatal resistance.
  sr  <- sr/ks
  
  # ra is aerodynamic resistance, rs is bulk surface resistance
  ra  <- 208/wind #(log((2-2/3*0.12)/(0.123*0.12))*log((2-2/3*0.12)/(0.1*0.123*0.12)))/(0.41^2*wind) # equal to 208/wind for hh=hw=2.
  #ra <- 208/wind
  rs <- sr/(0.5*24*0.12) # value of 70 when sr=100
  
  # Saturation vapor pressure , 
  es <- 0.6108*exp(tmin*17.27/(tmin+237.3))/2+0.6108*exp(tmax*17.27/(tmax+237.3))/2     
  ea <- 0.6108*exp((dpt)*17.27/((dpt)+237.3))
  vpd <- es - ea
  vpd[vpd<0] <- 0    # added because this can be negative if dewpoint temperature is greater than mean temp (implying vapor pressure greater than saturation).
  
  # delta - Slope of the saturation vapor pressure vs. air temperature curve at the average hourly air temperature 
  delta  <- (4098 * es)/(tmean + 237.3)^2  
  
  P <- 101.3*((293-0.0065*elev)/293)^5.26  # Barometric pressure in kPa
  lambda <- 2.501-2.361e-3*tmean # latent heat of vaporization    
  cp  <- 1.013*10^-3 # specific heat of air
  gamma <- cp*P/(0.622*lambda) # Psychrometer constant (kPa C-1)
  pa <- P/(1.01*(tmean+273)*0.287) # mean air density at constant pressure
  
  # Calculate potential max solar radiation or clear sky radiation	
  GSC=0.082      # MJ m -2 min-1 (solar constant)
  phi <- pi*lat/180 
  dr <- 1+0.033*cos(2*pi/365*DoY)      
  delt <- 0.409*sin(2*pi/365*DoY-1.39)     
  omegas <- acos(-tan(phi)*tan(delt)) 
  Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas))    # Daily extraterrestrial radiation
  Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)
  
  
  # radfraction is a measure of relative shortwave radiation, or of
  # possible radiation (cloudy vs. clear-sky)
  radfraction <- radiation/Rso
  radfraction[radfraction>1] <- 1
  
  # longwave  and net radiation
  longw <- 4.903e-9*n_days*((tmax+273.15)^4+(tmin+273.15)^4)/2*(.34-.14*sqrt(ea))*(1.35*radfraction-.35)     
  netrad <- radiation*n_days*(1-albedo)-longw     
  
  # ET0
  et0 <- .408*((delta*(netrad-G))+(pa*cp*vpd/ra*3600*24*n_days))/(delta+gamma*(1+rs/ra)) # the *2 is for comparison
  return(et0)
} 




aetmod <- function(et0,input,awc,soil_prev=NULL){
  # This function computes AET given ET0, H2O input, soil water capacity, and beginning-of-month soil moisture
  # Arguments:
  # et0: vector of monthly reference evapotranspiration in mm
  # input: vector of monthly water input to soil in mm
  # awc: vector of soil water capacity in mm
  # soil_prev: vector of soil water content for the previous month (mm).  If left NULL this is assigned to be zero.
  #
  # Value:
  # returns a data frame with columns for AET, deficit, end-of-month soil moisture, and runoff.
  
  N <- length(et0)
  runoff <- def <-  aet <- soil <- rep(NA,N) # 
  if(is.null(soil_prev)) soil_prev <- rep(0,N)
  
  deltasoil <- input-et0 # positive=excess H2O, negative=H2O deficit
  
  # Case when there is a moisture surplus:
  Case <- deltasoil>=0
  if(sum(Case)>0){
    aet[Case] <- et0[Case]
    def[Case] <- 0
    soil[Case] <- pmin(soil_prev[Case]+deltasoil[Case],awc[Case])	# increment soil moisture, but not above water holding capacity
    runoff[Case] <- pmax(soil_prev[Case]+deltasoil[Case]-awc[Case],0) # when awc is exceeded, send the rest to runoff
  }
  
  # Case where there is a moisture deficit:  soil moisture is reduced
  Case <- deltasoil<0
  if(sum(Case)>0){
    soildrawdown <- soil_prev[Case]*(1-exp(-(et0-input)[Case]/awc[Case]))	# this is the net change in soil moisture (neg)
    aet[Case] <- pmin(input[Case] + soildrawdown,et0[Case])
    def[Case] <- et0[Case] - aet[Case]
    soil[Case] <- soil_prev[Case]-soildrawdown
    runoff[Case] <- 0
  }
  
  return(data.frame(aet=aet,def=def,soil=soil,runoff=runoff))
  
}


monthlengths <- function()
{
  return(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
}





############## From here down are additions by Derek Young ##########





##########################################
### Run the Dobrowski methods


### This function combines the separate functions published by Dobrowski and returns AET, Def, and PET
### This function is called by the function "run_dobr_wb". This function does not have to be called directly.
dobr_wb <- function(snowpack_prev_init,tmean_prev_init,soil_prev_init,monthlist,T.m,P.m,R.m,Tmin.m,Tmax.m,wind.m,Tdew.m,L,E,S.max) {
  
  snowpack_prev_cur <- snowpack_prev_init
  tmean_prev_cur <- tmean_prev_init
  soil_prev_cur <- soil_prev_init
  
  input <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  pet <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  aet <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  def <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  snowpack <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  soilwater <- matrix(NA,nrow=nrow(T.m),ncol=ncol(T.m))
  
  for(i in 1:ncol(T.m)) {
    
    month <- monthlist[i]
    
    input.df <- snowmod(tmean=T.m[,i],ppt=P.m[,i],radiation=R.m[,i],snowpack_prev=snowpack_prev_cur,albedo=0.23,albedo_snow=0.8)
    input[,i] <- input.df$input
    pet[,i] <- monthlyET0(radiation=R.m[,i],tmax=Tmax.m[,i],tmin=Tmin.m[,i],wind=wind.m[,i],lat=L,elev=E,dpt=Tdew.m[,i],tmean_prev=tmean_prev_cur,albedo=input.df$albedo,month=month)
  
    pet[pet < 0] <- 0
    
    aet.df <- aetmod(et0=pet[,i],input=input.df$input,awc=S.max,soil_prev=soil_prev_cur)
    
    aet[,i] <- aet.df$aet
    def[,i] <- aet.df$def
    snowpack[,i] <- input.df$snow
    soilwater[,i] <- aet.df$soil
    
    snowpack_prev_cur <- input.df$snow
    tmean_prev_cur <- T.m[,i]
    soil_prev_cur <- aet.df$soil
    
    
    
  }
  
  wb <- list(aet=aet,def=def,pet=pet,snowpack=snowpack,soilwater=soilwater)
  
  return(wb)
  
}





## This function takes the 12 months from Jan to Dec of each monthly clim var.
## The data are intended to represent climate normals
## The script replicates the 12-months time series three times in order to provide sufficient spin-up
run_dob_wb <- function(T.m,P.m,R.m,Tmin.m,Tmax.m,wind.m,Tdew.m,L,E,S.max) {
  
  if(ncol(T.m) != 12) { error("You must supply 12 months of climate data.")}
  
  T.m <- cbind(T.m,T.m,T.m)
  P.m <- cbind(P.m,P.m,P.m)
  R.m <- cbind(R.m,R.m,R.m)
  Tmin.m <- cbind(Tmin.m,Tmin.m,Tmin.m)
  Tmax.m <- cbind(Tmax.m,Tmax.m,Tmax.m)
  Tdew.m <- cbind(Tdew.m,Tdew.m,Tdew.m)
  wind.m <- cbind(wind.m,wind.m,wind.m)
  #wind.m <- cbind(wind.m,wind.m,wind.m)
  
  snowpack_prev_init <- rep(0,nrow(T.m))
  tmean_prev_init <- T.m[,1]
  soil_prev_init <- S.max
  monthlist <- rep(c(1,2,3,4,5,6,7,8,9,10,11,12),3)
  
  
  
  wb <- dobr_wb(snowpack_prev_init,tmean_prev_init,soil_prev_init,monthlist,T.m,P.m,R.m,Tmin.m,Tmax.m,wind.m,Tdew.m,L,E,S.max)
    
  ### Take the last year (12 months) of AET and Def (after model has stabilized) and sum to get annual
  
  aet.last <- wb[["aet"]][,25:36]
  def.last <- wb[["def"]][,25:36]
  pet.last <- wb[["pet"]][,25:36]
  
  aet.ann <- apply(aet.last,1,sum)
  def.ann <- apply(def.last,1,sum)
  pet.ann <- apply(pet.last,1,sum)
  
  wb.df <- data.frame(PET.Dobr=pet.ann,AET.Dobr=aet.ann,Deficit.Dobr=def.ann)
  return(wb.df)
}

## This function takes any number of consecutive months of each monthly clim var (temp and precip).
## Radiation and dewpoint temp are climate normals
## You should provide at least 12 monts of data prior to period of interest for sufficient spinup
run_dob_wb_timeseries <- function(clim.df,years) {
  
  Tmin.m <- clim.df[,grep("tmin",names(clim.df))]
  Tmax.m <- clim.df[,grep("tmax",names(clim.df))]
  T.m <- clim.df[,grep("tmean",names(clim.df))]
  P.m <- clim.df[,grep("ppt",names(clim.df))]
  
  nyears <- ncol(Tmin.m)/12
  if(round(nyears) != nyears) {stop("Dobr data not supplied in multiple of 12 months")}
    
  if((ncol(Tmin.m) != ncol(Tmax.m)) | (ncol(Tmin.m) != ncol(T.m)) | (ncol(Tmin.m) != ncol(P.m))) stop("Dobr WB not supplied with equal number of months of Tmin, Tmax, Tmean, or Ppt.")

  Tdew.normal <- clim.df[,grep("td.normal",names(clim.df))]
  Tdew.m <- Tdew.normal[,rep(1:ncol(Tdew.normal),nyears)]
  
  rad.normal <- clim.df[,grep("rad.normal",names(clim.df))]
  rad.normal <- rad.normal *.0036 # convert rad to correct units for Dobr (from watt-hour/sq m day to MJ/sq m day)
  R.m <- rad.normal[,rep(1:ncol(rad.normal),nyears)]
  
  wind.normal <- clim.df[,grep("wind.normal",names(clim.df))]
  wind.m <- rad.normal[,rep(1:ncol(wind.normal),nyears)]
  
  L <- clim.df$lat
  E <- clim.df$elev
  S.max <- rep(150,length(L))
  

  snowpack_prev_init <- rep(0,nrow(T.m))
  tmean_prev_init <- T.m[,1]
  soil_prev_init <- S.max
  monthlist <- rep(c(1,2,3,4,5,6,7,8,9,10,11,12),nyears)
  

  wb <- dobr_wb(snowpack_prev_init,tmean_prev_init,soil_prev_init,monthlist,T.m,P.m,R.m,Tmin.m,Tmax.m,wind.m,Tdew.m,L,E,S.max)
  
  ### Summarize annual wb
  
  aet <- wb[["aet"]]
  def <- wb[["def"]]
  pet <- wb[["pet"]]
  
  aet.ann <- matrix(NA,nrow=nrow(aet),ncol=nyears-1)
  def.ann <- matrix(NA,nrow=nrow(aet),ncol=nyears-1)
  pet.ann <- matrix(NA,nrow=nrow(aet),ncol=nyears-1)
  
  rownames(aet.ann) <- rownames(def.ann) <- rownames(pet.ann) <- clim.df$tree.id
  colnames(aet.ann) <- colnames(def.ann) <- colnames(pet.ann) <- years[-1]
  
  colnames(aet.ann) <- paste("dobr.aet",colnames(aet.ann),sep=".")
  colnames(def.ann) <- paste("dobr.def",colnames(def.ann),sep=".")
  
  for(i in 1:(nyears-1)) {
    
    month.indices <- (12*(i-1)+1):(12*(i-1)+12)
    month.indices <- month.indices + 9 # shift this so water year starts in Oct (Oct-Sept water year)
    
    aet.year <- aet[,month.indices]
    aet.ann[,i] <- apply(aet.year,1,sum)
    
    def.year <- def[,month.indices]
    def.ann[,i] <- apply(def.year,1,sum)
    
    pet.year <- pet[,month.indices]
    pet.ann[,i] <- apply(pet.year,1,sum)
    
  }
    
  wb.df <- list(PET.Dobr=pet.ann,AET.Dobr=aet.ann,Deficit.Dobr=def.ann)
  return(wb.df)
}

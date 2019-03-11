#
# utils:
#     Some utility functions for calculating water balance parameters.
#

daylength <- function(latitude,month)
{
	# Day length from Forthsith (Gavin paper)
	# phi: constant based on julian date.  These are for the 15th of the month
	phi <- c(	-0.374329268,
			-0.231959485,
			-0.042826351,
			 0.165449779,
			 0.326011171,
			 0.406337883,
			 0.378109740,
			 0.249923394,
			 0.058515387,
			-0.143269334,
			-0.318990215,
			-0.405457073)
	
	phi.m <- phi[month]
	
	
	return(24 - (24 / pi) * acos((sin(latitude) * sin(phi.m)) / (cos(latitude) * cos(phi.m))))
}

monthlengths <- function(months)
{
  
  ml <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	return(ml[months])
}

monthinits <- function()
{
	return(c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
}

as.radians <- function(x)
{
	return(x * pi / 180)
}

as.degrees <- function(x)
{
	return(x * 180 / pi)
}

# temp in C
SatVapPresSlope.old <- function(temp)
{
	return(((2508.3 / (temp + 237.3)^2)) * exp(17.3 * temp / (temp + 237.3)))
}

SatVapPresSlope <- function(temp) {
  sl <- (4098*0.6108*exp(17.27*temp/(temp+237.3)))/((temp+237.3)^2)
  return(sl)
}

#Bernoulli formula
ElevToPress <- function(elev) {
  
  101325*exp(-9.80665*0.0289644*elev/(8.31447*288.15))/1000
  
}


WindTo2m <- function(wind,elev) {
  U <- wind*4.87/log(67.8*elev-5.42)
  return(U)
}


PET.Thorn <- function(WB.params)
{
  
  T.m     <- WB.params$T.m
  month <- WB.params$month
  L <- as.radians(WB.params$L)
  
  nmonths <- length(month)
  nyears <- nmonths/12
  if(round(nyears) != nyears) stop("You did not supply a multiple of 12 months.")

  T.m.frz <- T.m < 0
  T.m.hot <- T.m >= 26.5
    
  T.m.i <- ifelse(T.m < 0,0,T.m) #because the I calculations below cannot accept a T.m less than 0

  I <- sum((T.m.i / 5)^1.514)/nmonths
  a <- 6.75e-7 * I^3 - 7.71e-5 * I^2 + 1.79e-2 * I + 0.49
  
  PET.m          <- 16 * (10 * T.m / I)^a                            #    0-26.5 degrees
  PET.m[T.m.frz] <- 0                                                #  < 0 degrees
  PET.m[T.m.hot] <- (-415.85 + 32.24 * T.m - 0.43 * T.m^2)[T.m.hot]  # >= 26.5 degrees
  
  return(PET.m * (monthlengths(month) / 30) * (daylength(L,month) / 12))
}


AET.Wil <- function(PET.m, WB.params, soil.start)
{
  T.m   <- WB.params$T.m
  P.m   <- WB.params$P.m
  S.max <- WB.params$S.max
  month <- WB.params$month
  

  
  # Set initial guesses for W.s and W
  W.s.0 <- 0
  W.0   <- soil.start
  W.s   <- 1
  W     <- 0

    # In this method, everything is calculated on a quasi-daily basis
    # So temp, precip, PET, etc. are all expanded to reflect this
    nmonths <- length(T.m)
    nyears <- nmonths/12
    ndays <- nyears*365
    monthlengths.all <- monthlengths(month)
  
    T <- rep(T.m, monthlengths.all)
    P.s <- P.r <- unlist(mapply(function(x, y) rep(x, y) / y, P.m, monthlengths.all))
    P.s[T >= -1] <- 0
    P.r[T <  -1] <- 0
    E.0 <- unlist(mapply(function(x, y) rep(x, y) / y, PET.m, monthlengths.all))
    
    # Amount that melts each day
    M <- 2.63 + (2.55 * T) + (0.0912 * T * P.r)
    M[M < 0] <- 0
    
    # Amount of water stored as snow each day
    W.s <- c(W.s.0, rep(NA, ndays-1))
    for(i in 2:ndays)
    {
      snow.cover <- W.s[i - 1] + P.s[i]
      if(M[i] > snow.cover)
        M[i] <- snow.cover
      
      W.s[i] <- W.s[i - 1] + P.s[i] - M[i]
    }
    
    # Evaporative demand: (+) is recharging, (-) is a demand
    D <- M + P.r - E.0
    
    # Evapotranspiration proportion: evapotranspiration slows down as the soil water becomes depleted
    beta.d <- rep(1, ndays)
    if(D[1] < 0)
      beta.d[1] <- 1 - exp(-6.68)
    
    # Amount of water in the soil and surplus each day
    W <- c(W.0, rep(NA, ndays-1))
    S <- rep(0, ndays)
    for(i in 2:ndays)
    {
      if(D[i] < 0)
        beta.d[i] <- 1 - exp(-6.68 * W[i - 1] / S.max)
      W[i] <- W[i - 1] + beta.d[i]*D[i]
      if(W[i] > S.max)
      {
        S[i] <- W[i] - S.max
        W[i] <- S.max
      }
    }
    
#     # Amount of water in the soil on the last day of each month
#     W.30s <- W[cumsum(monthlengths())]
#     # Difference between months
#     delta.W <- c(W.30s[1] - W.30s[12], W.30s[2:12] - W.30s[1:11])
#     
#     # Dec 31 is the new initial condition for Jan 1
#     W.s.0 <- W.s[365]
#     W.0   <- W[365]
#     ctr   <- ctr + 1
#     if(ctr > 10)
#       stop("Error in Wilmott AET")
#   }
    

    # Amount of water in the soil on the last day of each month
    W.30s <- W[cumsum(monthlengths.all)]
    # Difference between months
    delta.W <- c(0, W.30s[2:nmonths] - W.30s[1:(nmonths-1)])


    monthends <- cumsum(monthlengths.all)
    monthstarts <- c(1,1+monthends)[1:nmonths]
    monthdays <- mapply(function(x, y) x:y, monthstarts, monthends)
    
    # Reduce daily variables to monthly variables
    P.r <- sapply(monthdays, function(x) sum(P.r[x]))
    M   <- sapply(monthdays, function(x) sum(M[x]))
    S   <- sapply(monthdays, function(x) sum(S[x]))
    
    
    AET.m <- P.r + M - delta.W - S
    AET.m[AET.m > PET.m] <- PET.m[AET.m > PET.m]
  

  return(AET.m)

}

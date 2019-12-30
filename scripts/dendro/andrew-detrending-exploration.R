# Test out detrending methods on the raw ring widths. Loads a data set saved from the script andrew-prelim-eda.R

library(dplR)

dat.rwl <- csv2rwl("./working-data/tree_rings_out.csv")

n.missings <- apply(dat.rwl, 2, f<-function(x) return(sum(is.na(x))))
rwl.report(dat.rwl[,-(which(n.missings>40))])

plot(dat.rwl[,301:350], plot.type="spag")

# Try negative exponential detrending on some of the data sets 
dat.rwi <- detrend.series(dat.rwl[,720], method = "ModNegExp", verbose=TRUE)
names(dat.rwl)[720]

# Try stiff spline detrending on some of the data sets 
dat.rwi <- detrend.series(dat.rwl[,720], method = "Spline", nyrs=30, f=0.5, verbose=TRUE)
names(dat.rwl)[720]
# Try stiff spline detrending on some of the data sets 
dat.rwi <- detrend.series(dat.rwl[,738], method = "Spline", nyrs= 30, f=0.5, verbose=TRUE)
names(dat.rwl)[720]
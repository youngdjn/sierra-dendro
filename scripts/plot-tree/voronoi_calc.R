setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro")

library(sp)
library(rgeos)

deg.to.rad <- function(x) {
  return((x*3.14159)/180)
}


d <- read.csv("data/plot-and-tree/field-data/adult_trees_merged 2015 10 13.csv",header=TRUE,stringsAsFactors=FALSE)
head(d)

d$foc <- ifelse((d$foc == "") | is.na(d$foc),FALSE,TRUE)
d$cmp <- ifelse((d$cmp == "") | is.na(d$cmp),FALSE,TRUE)
d$backbrg <- ifelse((d$backbrg == "") | is.na(d$backbrg),FALSE,TRUE)
d$dist.to.ref <- ifelse(d$dist.to.ref == 0,0.1,d$dist.to.ref) #if distance listed as 0, it really just means it was very close

#turn backbearings to regular bearings
d$brg.to.ref <- as.numeric(d$brg.to.ref)
d$brg.comp <- ifelse(d$backbrg,((d$brg.to.ref + 180) %% 360),d$brg.to.ref)
d$brg.perp <- (d$brg.comp + 90) %% 360


# add the focal tree ID to each competitor tree
# go line-by-line

current.focal <- NA
for(i in 1:nrow(d)) {
  
  if(d[i,]$foc) {
    current.focal <- d[i,]$tree.id
  } else {
    d[i,]$tree.id <- current.focal
  }

}

## for each focal tree, compute its voronoi
trees <- unique(d$tree.id)

voronoi.area <- data.frame()
for(tree in trees) {
  
  cat("Running tree",tree,"\n")
  
  d.tree.cmp <- d[(d$tree.id == tree) & (d$cmp),]
  
  #if tree has no competitors listed, skip
  if(nrow(d.tree.cmp)==0) next()
  
  #if any of the distance columns are NA, then there was no clear competitor, so skip. alternatively could assume the competitor is at X distance
  if(sum(is.na(d.tree.cmp$dist.to.ref))>0) next()
  
  #if any of the bearing columns are NA, then there was an error in recording, so skip. alternatively could assume the competitor is at X distance
  if(sum(is.na(d.tree.cmp$brg.to.ref))>0) next()
  
  cmp.loc <- data.frame(x=0,y=0,type="Focal") #start it with the focal point
  
  lines <- list()
  for(i in 1:nrow(d.tree.cmp)) {
    brg.rad <- deg.to.rad(d.tree.cmp[i,]$brg.comp)
    brg.perp.rad <- deg.to.rad(d.tree.cmp[i,]$brg.perp)
    dist <- d.tree.cmp[i,]$dist.to.ref
    
    # get the competitor coords
    x.tree <- cos(brg.rad) * dist
    y.tree <- sin(brg.rad) * dist
    
    #get the midpoint of the line between focal and competitor
    x <- cos(brg.rad) * dist/2
    y <- sin(brg.rad) * dist/2
    
    
    
    #get endpoints of a line extending out 50m perpendicularly to that line
    perp.x1 <- x + cos(brg.perp.rad) * 50
    perp.y1 <- y + sin(brg.perp.rad) * 50
    
    perp.x2 <- x - cos(brg.perp.rad) * 50
    perp.y2 <- y - sin(brg.perp.rad) * 50
    
    line.coords<- cbind(c(perp.x1,perp.x2),c(perp.y1,perp.y2))
    
    line <- Line(line.coords)
    lines[[i]] <- line

    loc <- data.frame(x=x.tree,y=y.tree,type="Competitor")
    cmp.loc <- rbind(cmp.loc,loc)
  }
  
  lines2 <- Lines(lines,"a")
  
  lines.sp <- SpatialLines(list(lines2))
  cmp.sp <- SpatialPointsDataFrame(cmp.loc[,c("x","y")],data=cmp.loc)
  
  #plot them
#   plot(lines.sp)
#   points(cmp.sp)
  
  
  
  # make a circle polygon around the central tree with radius 20m, and chop it using the lines
  focal.point <- SpatialPoints(cbind(0,0))
  focal.poly <- gBuffer(focal.point,width=40)
  lpi <- gIntersection(focal.poly, lines.sp)
  blpi <- gBuffer(lpi, width = 0.000001)
  dpi <- gDifference(focal.poly, blpi)
  dpi <- disaggregate(dpi)
  
  #find the polygon under the focal tree
  focal.poly.index <- focal.point %over% dpi
  focal.poly.1 <- dpi[focal.poly.index]
  
  
  plot(lines.sp)
  plot(focal.poly.1,col="slateblue",add=TRUE)
  points(cmp.sp,col=cmp.loc$type,pch=19)
  text(20,20,labels=tree)
  
  
  #check if any of the polygons extends all the way to the outside of the starting circle (error in tree coords)
  focal.poly.larger <- gBuffer(focal.poly,width=0.5)
  outer.ring <- gDifference(focal.poly.larger,focal.poly)
  outer.ring.larger <- gBuffer(outer.ring,width=.01)
  
  #see if the focal voronoi polygon touches the ring (i.e. does not have an outer boundary)
  a <- over(focal.poly.1,outer.ring.larger)
  if(!is.na(a[1])) {
    cat("Tree",tree,"polygon extends infinitely; skipping\n")
    next()
  }
  

  area <- focal.poly.1@polygons[[1]]@Polygons[[1]]@area 
  a <- data.frame(tree.id=tree,voronoi.area=area)
  voronoi.area <- rbind(voronoi.area,a)

}

write.csv(voronoi.area,"data/plot-and-tree/processed/voronoi_area.csv",row.names=FALSE)

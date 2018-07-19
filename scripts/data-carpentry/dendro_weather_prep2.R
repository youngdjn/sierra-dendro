setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro")

library(dplyr)

###################################
##### Open and merge datasets #####
###################################

# get list of trees with plots, and other ancillary data
# this file comes from:
trees <- read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE,stringsAsFactors=FALSE)

#drop trees that were resampled under a different ID
trees$tree.id <- toupper(trees$tree.id)
trees$former.id <- toupper(trees$former.id)

former.ids <- trees$former.id
trees <- trees[!(trees$tree.id %in% former.ids),]


#do any trees have duplicate samples?
ndup <- sum(duplicated(trees$tree.id))
if(ndup > 0) {warning("Some tree ids duplicated in trees_loc spreadsheet.")}


#### open ring width data
source("scripts/data-carpentry/dendro-functions/read_pos_extended_bai.R")

chron <- open.chron(samples=trees$tree.id,samples.secondary=trees$former.id,unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE)

length(unique(trees$tree.id)) # were all the requested series provided?
ncol(chron$chron) # no; why not? because the chronology drops colums for the cores it could not find

# clean it based on externally-recorded records of correlation
chron.clean <- clean.chron(chron) #fix error in colSums  #this function does:

# remove A, B, Z, or T from core names to get tree.ids
names(chron.clean$chron) <- gsub("a?A?b?B?t?T?z?Z?","",names(chron.clean$chron))
names(chron$chron) <- gsub("a?A?b?B?t?T?z?Z?","",names(chron$chron))


chron.rwi <- spline.na.rm(chron.clean$chron) #! why errors here?
chron.rwi$year <- rownames(chron.rwi)

radii <- chron.clean$radius
radii <- data.frame(tree.id=names(radii),radius=radii)


### calculate basal area using internal radius
bai.ba.results <- ba.bai.calc(chron.clean$chron,chron.clean$radius)
ba <- bai.ba.results$ba
bai <- bai.ba.results$bai
bai.ba <- bai.ba.results$bai.ba
ba.prev <- bai.ba.results$ba.prev
names(ba) <- names(bai) <- names(bai.ba) <- names(ba.prev) <- names(chron.clean$chron)
ba$year <- rownames(ba)
bai$year <- rownames(bai)
bai.ba$year <- rownames(bai.ba)
ba.prev$year <- rownames(ba.prev)


### calculate basal area using external radius

# get radius from tree survey data
#radius.external <- as.numeric(trees$dbh)/2*10

# estimate the "internal" radius from the external radius based on the observed empirical relationship
radius.external <- 26.08439 + as.numeric(trees$dbh)*3.48909

names(radius.external) <- trees$tree.id

### exclude trees that did not have ring measurements to the outside (measured diam/radius not representative of outermost ring)
outer.rings <- chron.clean$chron[c("2012","2013","2014"),]
trees.meas.tobark <- colnames(outer.rings[(colSums(outer.rings,na.rm=TRUE) > 0)])
chron.bacalc <- chron.clean$chron[,trees.meas.tobark]

bai.ba.external.results <- ba.bai.calc(chron.bacalc,radius.external)
ba.ext <- bai.ba.external.results$ba
bai.ext <- bai.ba.external.results$bai
bai.ba.ext <- bai.ba.external.results$bai.ba
ba.ext.prev <- bai.ba.external.results$ba.prev
names(ba.ext) <- names(bai.ext) <- names(bai.ba.ext) <- names(ba.ext.prev) <- names(chron.bacalc)
ba.ext$year <- rownames(ba.ext)
bai.ext$year <- rownames(bai.ext)
bai.ba.ext$year <- rownames(bai.ba.ext)
ba.ext.prev$year <- rownames(ba.ext.prev)


#### open climate data
# trees.clim <- read.csv("data/non-synced/climate-extracted/tree_clim_full.csv",header=TRUE)
# 
# ## for now, thin to a few columns
# trees.clim = trees.clim %>%
#   select(tree.id,year,ppt,tmean,ppt.z,tmean.z,ppt,
#          ppt1,tmean1,ppt.z1,tmean.z1,ppt1,
#          ppt2,tmean2,ppt.z2,tmean.z2,ppt2,
#          ppt3,tmean3,ppt.z3,tmean.z3,ppt3,
#          rad.tot,rad.03,rad.06) %>%
#   mutate_at(vars(-tree.id,-year),funs(signif))
# 
# ## write simplified file
# write.csv(trees.clim,"data/plot-and-tree/processed/trees_climate.csv",row.names=FALSE)

trees.clim <- read.csv("data/plot-and-tree/processed/trees_climate.csv",header=TRUE,stringsAsFactors=FALSE)


#### merge chron, ba, radius, and climate data and tree survey data

# melt chron
library(reshape)
rwi.melt <- melt(chron.rwi,id.vars="year",variable_name="tree.id")
names(rwi.melt)[3] <- "rwi"

#internal radius
ba.melt <- melt(ba,id.vars="year",variable_name="tree.id")
bai.melt <- melt(bai,id.vars="year",variable_name="tree.id")
bai.ba.melt <- melt(bai.ba,id.vars="year",variable_name="tree.id")
ba.prev.melt <- melt(ba.prev,id.vars="year",variable_name="tree.id")
names(ba.melt)[3] <- "ba"
names(bai.melt)[3] <- "bai"
names(bai.ba.melt)[3] <- "bai.ba"
names(ba.prev.melt)[3] <- "ba.prev"

#external radius
ba.ext.melt <- melt(ba.ext,id.vars="year",variable_name="tree.id")
bai.ext.melt <- melt(bai.ext,id.vars="year",variable_name="tree.id")
bai.ba.ext.melt <- melt(bai.ba.ext,id.vars="year",variable_name="tree.id")
ba.ext.prev.melt <- melt(ba.ext.prev,id.vars="year",variable_name="tree.id")
names(ba.ext.melt)[3] <- "ba.ext"
names(bai.ext.melt)[3] <- "bai.ext"
names(bai.ba.ext.melt)[3] <- "bai.ba.ext"
names(ba.ext.prev.melt)[3] <- "ba.ext.prev"

# merge climate and rwi data
trees.clim$tree.id <- toupper(trees.clim$tree.id)
rwi.melt$tree.id <- toupper(rwi.melt$tree.id)
clim.rwi <- merge(trees.clim,rwi.melt,by=c("tree.id","year")) # only get climate for trees and years that have RWI data
clim.rwi <- merge(clim.rwi,ba.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,bai.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,bai.ba.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,ba.prev.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,ba.ext.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,bai.ext.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,bai.ba.ext.melt,by=c("tree.id","year"),all.x=TRUE)
clim.rwi <- merge(clim.rwi,ba.ext.prev.melt,by=c("tree.id","year"),all.x=TRUE)

# make a composite ba, bai, bai.ba column that uses internal when available and external when not

clim.rwi$ba.comb <- ifelse(is.na(clim.rwi$ba),clim.rwi$ba.ext,clim.rwi$ba)
clim.rwi$bai.comb <- ifelse(is.na(clim.rwi$bai),clim.rwi$bai.ext,clim.rwi$bai)
clim.rwi$bai.ba.comb <- ifelse(is.na(clim.rwi$bai.ba),clim.rwi$bai.ba.ext,clim.rwi$bai.ba)
clim.rwi$ba.prev.comb <- ifelse(is.na(clim.rwi$ba.prev),clim.rwi$ba.ext.prev,clim.rwi$ba.prev)
clim.rwi <- clim.rwi[clim.rwi$bai.ba.comb != Inf,]

  

# for export purposes (e.g. comparing with DBH,) internal radius only counts for trees that were measured to the outermost ring
radii <- radii[radii$tree.id %in% trees.meas.tobark,]

# merge in the internal radius
clim.rwi <- merge(clim.rwi,radii,by="tree.id",all.x=TRUE)

# merge in tree survey data
trees$tree.id <- toupper(trees$tree.id)
trees.clim.rwi <- merge(clim.rwi,trees,by="tree.id")

# merge in competition data
voronoi <- read.csv("data/plot-and-tree/processed/voronoi_area.csv",header=TRUE,stringsAsFactors=FALSE)
voronoi$tree.id <- toupper(voronoi$tree.id)
trees.clim.rwi <- merge(trees.clim.rwi,voronoi,by="tree.id",all.x=TRUE)


trees.clim.rwi[(trees.clim.rwi$dbh == 0) | (trees.clim.rwi$dbh == ""),]$dbh <- NA # <------ This data frame ready for stats

trees.clim.rwi = trees.clim.rwi %>%
  mutate_at(vars(ppt:radius),funs(signif))

write.csv(trees.clim.rwi,"data/compiled-for-analysis/trees_clim_rwi.csv",row.names=FALSE)


#### output year-level data
clim.rwi.out = clim.rwi %>%
  filter(year > 1959) %>%
  mutate_at(vars(-tree.id,-year),funs(signif)) # truncate to 6 digits
write.csv(clim.rwi.out,"data/compiled-for-analysis/chronology_climate.csv",row.names=FALSE)

#### output tree-level data
trees.out.pre = left_join(trees,voronoi,by="tree.id")
trees.out = left_join(trees.out.pre,radii,by="tree.id")
write.csv(trees.out,"data/compiled-for-analysis/trees.csv",row.names=FALSE)


trees.clim.rwi <- read.csv("data/compiled-for-analysis/trees_clim_rwi.csv",header=TRUE,stringsAsFactors=FALSE)


# calculate plot-level climate normals (averaging tree-level normals)
tree.dup <- duplicated(trees.clim.rwi$tree.id) # only keep one record per tree (because the normals are replicated across all years of each tree)
normal.cols <- grep("normal",names(trees.clim.rwi),value=TRUE)
#tree.norm <- trees.clim.rwi[!tree.dup,c("plot.id","tree.id",normal.cols,"rad.tot","rad.03","rad.06","x","y","dbh","radius","rad.hr.03","rad.hr.06","species","voronoi.area")]
tree.norm <- trees.clim.rwi[!tree.dup,c("plot.id","tree.id",normal.cols,"rad.tot","rad.03","rad.06","x","y","dbh","radius","species","voronoi.area")]

tree.norm$radius <- as.numeric(tree.norm$radius)
tree.norm$dbh <- as.numeric(tree.norm$dbh)

tree.norm$cluster <- NA
tree.norm[which(tree.norm$y > 150000),"cluster"] <- "Plumas"
tree.norm[which((tree.norm$y < 150000) & (tree.norm$y > 0)),"cluster"] <- "Tahoe"
tree.norm[which((tree.norm$y < 0) & (tree.norm$y > -70000)),"cluster"] <- "Yose"
tree.norm[which((tree.norm$y < -70000)),"cluster"] <- "Sierra"



plot.norm <- aggregate(tree.norm,by=list(tree.norm$plot.id),FUN=mean)
plot.norm$plot.id <- plot.norm$Group.1
plot.norm <- plot.norm[,-1]

# assign plots to plumas, tahoe, yose, sierra
plot.norm$cluster <- NA
plot.norm[plot.norm$y > 150000,"cluster"] <- "Plumas"
plot.norm[(plot.norm$y < 150000) & (plot.norm$y > 0),"cluster"] <- "Tahoe"
plot.norm[(plot.norm$y < 0) & (plot.norm$y > -70000),"cluster"] <- "Yose"
plot.norm[(plot.norm$y < -70000),"cluster"] <- "Sierra"




#### Determine 































# see how exterior radius relates to internal radius
# 
# # only include trees that were measured all the way to outside
# tree.norm.radius <- tree.norm[!(tree.norm$tree.id %in% c("2528","1231","1250","2534","2531")),]
# tree.norm.radius <- tree.norm.radius[tree.norm.radius$tree.id %in% trees.meas.tobark,]
# 
# 
# # add cluster
# tree.norm.radius <- merge(tree.norm.radius,plot.norm[,c("plot.id","cluster")],by="plot.id")
# 
# 
# plot(radius~dbh,data=tree.norm.radius)
# #identify(tree.norm.radius$dbh,tree.norm.radius$radius,labels=tree.norm.radius$tree.id)
# cor(tree.norm.radius$radius,tree.norm.radius$dbh,use="complete.obs")^2
# m1 <- lm(radius~dbh,data=tree.norm.radius)
# m2 <- lm(radius~dbh+ppt.normal,data=tree.norm.radius)
# m3 <- lm(radius~dbh+cluster,data=tree.norm.radius)
# 
# AIC(m1,m2,m3) # BIC prefers model with only dbh; may want to add cluster and ppt.norm to predictive model?
# 
# summary(m1)
# # internal radius = 26.08439 + 3.48909 * dbh




# make a tree-level normal df
tree.dup <- duplicated(trees.clim.rwi$tree.id) # only keep one record per tree (because the normals are replicated across all years of each tree)
normal.cols <- grep("normal",names(trees.clim.rwi),value=TRUE)
tree.norm <- trees.clim.rwi[!tree.dup,c("tree.id","plot.id","species","dbh",normal.cols,"rad.tot","rad.03","rad.06","x","y","bai.ba","bai.ba.comb","bai.ba.ext","voronoi.area")]

# assign plots to plumas, tahoe, yose, sierra
tree.norm$cluster <- NA
tree.norm[which(tree.norm$y > 150000),"cluster"] <- "Plumas"
tree.norm[which((tree.norm$y < 150000) & (tree.norm$y > 0)),"cluster"] <- "Tahoe"
tree.norm[which((tree.norm$y < 0) & (tree.norm$y > -70000)),"cluster"] <- "Yose"
tree.norm[which((tree.norm$y < -70000)),"cluster"] <- "Sierra"

tree.norm$dbh <- as.numeric(tree.norm$dbh)




if(FALSE) {




  
  
  
  
  
  ################################
  ######### Stats. This has been moved to a separate script ################
  ################################
  library(lme4)
  source("Dendro/Climate correlations/R/rsquared_glmm.R")
  
  #restrict to 30-year window
  #trees.clim.rwi <- subset(trees.clim.rwi,(year > 1976) & (year < 2007))
  #trees.clim.rwi <- subset(trees.clim.rwi,(year > 1960) & (year < 2007))
  #trees.clim.rwi <- subset(trees.clim.rwi,(year > 2000)) # alternate for snowpack data
  
  
  # thin to VIC subset
  #trees.clim.rwi <- subset(trees.clim.rwi,(year > 1965) & (year < 2007))
  
  #thin to subset of years
  #trees.clim.rwi <- subset(trees.clim.rwi,(year > 2002) & (year < 2015))
  #trees.clim.rwi <- subset(trees.clim.rwi,year < 2007)
  
  
  
  rsq.list <- list()
  slope.list <- list()
  stderr.list <- list()
  pval.list <- list()
  
  ### only necessary for fitting with BAI: remove trees with no BAI and with BAI less than 20000
  trees.clim.rwi <- trees.clim.rwi[!is.na(trees.clim.rwi$ba),]
  trees.clim.rwi <- trees.clim.rwi[trees.clim.rwi$ba > 20000,]
  
  #trees.clim.rwi <- trees.clim.rwi[!is.na(trees.clim.rwi$ba.ext),]
  #trees.clim.rwi <- trees.clim.rwi[trees.clim.rwi$ba > 25000,]
  
  
  
  for(tree in unique(trees.clim.rwi$tree.id)) {
    
    tree.clim.rwi <- subset(trees.clim.rwi,tree.id==tree)
    if(sum(!(is.na(tree.clim.rwi$rwi))) < 25) next() # skip if tree has less than 25 rings in the focal period
    
    #m <- lm(rwi ~ dobr.def.z + dobr.def.z1 + dobr.def.z2, data=tree.clim.rwi)
    m <- lm(rwi ~ ppt.z, data=tree.clim.rwi)
    #m <- lm(rwi ~ snow.late.z + snow.late.z1 + snow.late.z2, data=tree.clim.rwi)
    #m <- lm(rwi ~ vic.def.z + vic.def.z1 + vic.def.z2, data=tree.clim.rwi)
    
    rsq <- summary(m)$r.squared
    slope <- coef(m)[2]
    stderr <- summary(m)$coefficients[2,2]
    pval <- summary(m)$coefficients[2,4]
    
    rsq.list[[tree]] <- rsq
    slope.list[[tree]] <- slope
    stderr.list[[tree]] <- stderr
    pval.list[[tree]] <- pval
    
  }
  
  mod.fits <- data.frame(rsq=unlist(rsq.list),slope=unlist(slope.list),stderr=unlist(stderr.list),pval=(unlist(pval.list)))
  mod.fits$tree.id <- rownames(mod.fits)
  
  # get species column
  trees.fits.norm <- merge(tree.norm,mod.fits,by="tree.id")
  plot(rsq~ppt.normal,data=trees.fits.norm)
  plot(slope~ppt.normal,data=trees.fits.norm)
  
  trees.fits.norm$ppt.norm.std <- (trees.fits.norm$ppt.normal - mean(trees.fits.norm$ppt.normal))/(trees.fits.norm$ppt.normal)
  #trees.fits.norm.sp <- subset(trees.fits.norm,species=="PSME")
  trees.fits.norm.sp <- trees.fits.norm
  m <- lmer(slope~ppt.norm.std+(1+ppt.norm.std|plot.id),data=trees.fits.norm.sp)
  r.squared(m)
  
  
  # divide into PSME and OTHER
  #trees.fits.norm.sp$species <- ifelse(trees.fits.norm.sp$species=="PSME","PSME","OTHER")
  trees.fits.norm.sp <- trees.fits.norm.sp[trees.fits.norm.sp$species=="PSME",]
  
  #compute plot-level averages for these vars
  
  plots.fits.norm.sp <- aggregate(trees.fits.norm.sp,by=list(trees.fits.norm.sp$plot.id,trees.fits.norm.sp$species),FUN=mean)
  plots.fits.norm.sp$plot.id <- plots.fits.norm.sp$Group.1
  plots.fits.norm.sp$species <- plots.fits.norm.sp$Group.2
  
  plot(rsq~ppt.normal,data=plots.fits.norm.sp)
  plot(slope~ppt.normal,data=plots.fits.norm.sp)
  
  m <- lm(slope~rad.03,data=plots.fits.norm.sp)
  summary(m)
  
  
  ### add plot-level variables (most importantly, coordinates)
  plots.fits <- merge(plots.fits.norm.sp,plot.norm[,c("plot.id","x","y")],by="plot.id",all.x=TRUE)
  
  ### make spatial data frame
  library(sp)
  albers <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  plots.fits.sp <- SpatialPointsDataFrame(coords=data.frame(x=plots.fits$x.x,y=plots.fits$y.x),data=plots.fits,proj4string=albers)
  
  write.csv(plots.fits.sp,"Dendro/Climate correlations/plots_fits_sp.csv",row.names=FALSE)
  
  #### Designate as plumas, tahoe, yose, or sierra
  plots.fits.sp$cluster <- NA
  plots.fits.sp[plots.fits.sp$y > 150000,"cluster"] <- "Plumas"
  plots.fits.sp[(plots.fits.sp$y < 150000) & (plots.fits.sp$y > 0),"cluster"] <- "Tahoe"
  plots.fits.sp[(plots.fits.sp$y < 0) & (plots.fits.sp$y > -70000),"cluster"] <- "Yose"
  plots.fits.sp[(plots.fits.sp$y < -70000),"cluster"] <- "Sierra"
  
  plots.fits <- as.data.frame(plots.fits.sp)
  
  library(ggplot2)
  
  
  #thin to only PSME and PIPO for now
  #plots.fits <- subset(plots.fits,species %in% c("PSME","PIPO"))
  
  library(gridExtra)
  
  plot.levels <- unique(plots.fits$plot.id[order(plots.fits$cluster,-plots.fits$ppt.normal)])
  plots.fits$plot.id <- factor(plots.fits$plot.id,levels=plot.levels)
  
  ggplot(plots.fits,aes(x=plot.id,y=slope,colour=ppt.normal)) +
    geom_point(aes(size=rsq),colour="black") +
    geom_errorbar(aes(ymin=slope-stderr,ymax=slope+stderr),colour="black")+
    facet_grid(. ~ cluster,scales="free_x",space="free_x") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    geom_rug(size=4,sides="b")+
    scale_colour_gradient(low="#993300",high="#99FFFF")
  
  ggplot(plots.fits,aes(x=ppt.normal,y=slope,colour=cluster)) +
    geom_point(size=3) +
    geom_smooth(method="lm",se=TRUE)
  
  
  m <- glm(slope~ppt.normal +cluster,data=plots.fits,family=Gamma)
  
  
  
  #### Designate as plumas, tahoe, yose, or sierra
  trees.fits.norm.sp$cluster <- NA
  trees.fits.norm.sp[trees.fits.norm.sp$y > 150000,"cluster"] <- "1Plumas"
  trees.fits.norm.sp[(trees.fits.norm.sp$y < 150000) & (trees.fits.norm.sp$y > 0),"cluster"] <- "2Tahoe"
  trees.fits.norm.sp[(trees.fits.norm.sp$y < 0) & (trees.fits.norm.sp$y > -70000),"cluster"] <- "3Yose"
  trees.fits.norm.sp[(trees.fits.norm.sp$y < -70000),"cluster"] <- "4Sierra"
  
  
  ggplot(trees.fits.norm.sp,aes(x=ppt.normal,y=slope,colour=cluster)) +
    geom_point(size=3) +
    geom_smooth(method="lm",se=TRUE)
  
  
  m <- lmer(slope~ppt.normal+(1|cluster/plot.id),data=trees.fits.norm.sp)
  
  
  
  
  
  
  
  
  
  
  
  
  
  summary(m)
  
  
    geom_errorbar(aes(ymin=slope-stderr,ymax=slope+stderr),colour="black") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    scale_colour_gradient(low="#993300",high="#99FFFF")
  
  
  
  
  
  ggplot(plots.fits,aes(x=plot.id,y=1,fill=ppt.normal)) +
    geom_tile(stat="identity") +
    facet_grid(. ~ cluster,scales="free_x",space="free_x")
  
  
  
  
  
  #### Test synchrony ####
  trees.rwi.sp <- trees.clim.rwi[trees.clim.rwi$species == "PSME",] #Thin to PSME only
  
  plot.ids <- unique(trees.rwi.sp$plot.id)
  
  rs.plt <- NA
  
  for(plot.id in plot.ids) {
  
    trees.rwi.sp.plt <- trees.rwi.sp[trees.rwi.sp$plot.id==plot.id,]
    tree.ids <- unique(trees.rwi.sp.plt$tree.id)   
    
    rs.trees <- matrix(NA,nrow=length(tree.ids),ncol=length(tree.ids),dimnames=list(tree.ids,tree.ids))
    
    for(tree1 in tree.ids) { #for each tree
      for(tree2 in tree.ids) { #for each other tree
        
        tree1.df <- trees.rwi.sp.plt[trees.rwi.sp.plt$tree.id==tree1,]
        tree2.df <- trees.rwi.sp.plt[trees.rwi.sp.plt$tree.id==tree2,]
        
        a <- try(r <- cor(tree1.df$rwi,tree2.df$rwi,use="complete.obs"),silent=TRUE)
        if(class(a) == "try-error") {
          #cat("No complete pairs for",tree1,tree2,"\n")
          next()
        } else {
          rs.trees[tree1,tree2] <- r
        }
      }
    }
    
    #average the corr coefs
    r.plt.avg <- mean(rs.trees[lower.tri(rs.trees)],na.rm=TRUE)
    
    #store in rs.plt
    rs.plt[plot.id] <- r.plt.avg
  }
  
  rs.plt <- data.frame(rs.plt)
  rs.plt$synchrony <- rs.plt$rs.plt
  rs.plt$plot.id <- rownames(rs.plt)
  rs.plt <- rs.plt[,c("plot.id","synchrony")]
  
  
  plots.fits.sync <- merge(plots.fits,rs.plt,by="plot.id",all.x=TRUE)
  
  
  plot.levels <- unique(plots.fits.sync$plot.id[order(plots.fits.sync$cluster,-plots.fits.sync$ppt.normal)])
  
  plots.fits.sync$plot.id <- factor(plots.fits.sync$plot.id,levels=plot.levels)
  
  
  ggplot(plots.fits.sync,aes(x=plot.id,y=synchrony,colour=ppt.normal)) +
    geom_point() +
    facet_grid(. ~ cluster,scales="free_x",space="free_x") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
  
  m <- lm(synchrony~rad.03,data=plots.fits.sync)
  summary(m)
  
  
  
  
  
  
  
  ###################################################################
  ####### Get each plot's AIC drop for all candidate models #########
  ###################################################################
  
  # we have trees.clim.rwi to work with
  plots <- unique(trees.clim.rwi$plot.id)
  
  m.aics <- data.frame()
  
  
  for(plot in plots) {
    
    d.plot <- subset(trees.clim.rwi, (plot.id == plot) & (species == "PSME") )
  
    m <- list()
    
    m[["n"]] <- try(lmer(rwi ~ 1 + (1|tree.id),data=d.plot))
    m[["p1"]] <- try(lmer(rwi ~ ppt.z + (ppt.z|tree.id),data=d.plot))
    m[["p2"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + (ppt.z + ppt.z1|tree.id),data=d.plot))
    m[["p3"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + ppt.z2 + (ppt.z + ppt.z1 + ppt.z2|tree.id),data=d.plot))
    m[["t1"]] <- try(lmer(rwi ~ tmean.z + (tmean.z|tree.id),data=d.plot))
    m[["t2"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + (tmean.z + tmean.z1|tree.id),data=d.plot))
    m[["t3"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + tmean.z2 + (tmean.z + tmean.z1 + tmean.z2|tree.id),data=d.plot))
    m[["tp1"]] <- try(lmer(rwi ~ tmean.z + ppt.z + (tmean.z + ppt.z|tree.id),data=d.plot))
    m[["tp2"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + (tmean.z + ppt.z + tmean.z1 + ppt.z1 |tree.id),data=d.plot))
    m[["tp3"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + tmean.z2 + ppt.z2 + (tmean.z + ppt.z + tmean.z1 + ppt.z1 + tmean.z2 + ppt.z2 |tree.id),data=d.plot))
    m[["s1"]] <- try(lmer(rwi ~ snow.late.z + (snow.late.z|tree.id),data=d.plot))
    m[["s2"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + (snow.late.z + snow.late.z1|tree.id),data=d.plot))
    m[["s3"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + snow.late.z2 + (snow.late.z + snow.late.z1 + snow.late.z2|tree.id),data=d.plot))
    m[["dc.d1"]] <- try(lmer(rwi ~ dobr.def.z + (dobr.def.z|tree.id),data=d.plot))
    m[["dc.d2"]] <- try(lmer(rwi ~ dobr.def.z + dobr.def.z1 + (dobr.def.z+ dobr.def.z1|tree.id),data=d.plot))
    m[["dc.a1"]] <- try(lmer(rwi ~ dobr.aet.z + (dobr.aet.z|tree.id),data=d.plot))
    m[["dc.a2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.aet.z1 + (dobr.aet.z+ dobr.aet.z1|tree.id),data=d.plot))
    m[["dc.da1"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + (dobr.aet.z+ dobr.def.z|tree.id),data=d.plot))
    m[["dc.da2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + dobr.aet.z1 + dobr.def.z1 + (dobr.aet.z + dobr.def.z +dobr.aet.z1 + dobr.def.z1|tree.id),data=d.plot))
    m[["dp.d1"]] <- try(lmer(rwi ~ dobr.pub.def.z + (dobr.pub.def.z|tree.id),data=d.plot))
    m[["dp.d2"]] <- try(lmer(rwi ~ dobr.pub.def.z + dobr.pub.def.z1 + (dobr.pub.def.z+ dobr.pub.def.z1|tree.id),data=d.plot))
    m[["dp.a1"]] <- try(lmer(rwi ~ dobr.pub.aet.z + (dobr.pub.aet.z|tree.id),data=d.plot))
    m[["dp.a2"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.aet.z1 + (dobr.pub.aet.z+ dobr.pub.aet.z1|tree.id),data=d.plot))
    m[["dp.da1"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.def.z + (dobr.pub.aet.z+ dobr.pub.def.z|tree.id),data=d.plot))
    m[["dp.da2"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.def.z + dobr.pub.aet.z1 + dobr.pub.def.z1 + (dobr.pub.aet.z + dobr.pub.def.z +dobr.pub.aet.z1 + dobr.pub.def.z1|tree.id),data=d.plot))
    m[["v.d1"]] <- try(lmer(rwi ~ vic.def.z + (vic.def.z|tree.id),data=d.plot))
    m[["v.d2"]] <- try(lmer(rwi ~ vic.def.z + vic.def.z1 + (vic.def.z+ vic.def.z1|tree.id),data=d.plot))
    m[["v.a1"]] <- try(lmer(rwi ~ vic.aet.z + (vic.aet.z|tree.id),data=d.plot))
    m[["v.a2"]] <- try(lmer(rwi ~ vic.aet.z + vic.aet.z1 + (vic.aet.z+ vic.aet.z1|tree.id),data=d.plot))
    m[["v.da1"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + (vic.aet.z+ vic.def.z|tree.id),data=d.plot))
    m[["v.da2"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + vic.aet.z1 + vic.def.z1 + (vic.aet.z + vic.def.z +vic.aet.z1 + vic.def.z1|tree.id),data=d.plot))
    
    AIC.if <- function(x) {
      
      if(class(x) == "try-error") {
        return(NA)
      }
      
      
      if(length(x@optinfo$conv$lme4$messages > 0)) { #then there was no convergence
        return(NA)
      } else {
        return(AIC(x))
      }
    }
  
    m.aics.plot <- unlist(lapply( m, FUN=AIC.if )  )
    m.aics.plot <- m.aics.plot - m.aics.plot["n"] # get difference from null
  
    for(i in 1:length(m.aics.plot)){
     df <- data.frame(plot.id=plot,model=names(m.aics.plot[i]),AIC=m.aics.plot[i])
     m.aics <- rbind(m.aics,df)
    }
  
  }
  
  
  m.aics.norm <- merge(m.aics,plot.norm,by="plot.id")
  
  m.aics.norm$cluster <- NA
  m.aics.norm[m.aics.norm$y > 150000,"cluster"] <- "1Plumas"
  m.aics.norm[(m.aics.norm$y < 150000) & (m.aics.norm$y > 0),"cluster"] <- "2Tahoe"
  m.aics.norm[(m.aics.norm$y < 0) & (m.aics.norm$y > -70000),"cluster"] <- "3Yose"
  m.aics.norm[(m.aics.norm$y < -70000),"cluster"] <- "4Sierra"
  
  plot.levels <- unique(m.aics.norm$plot.id[order(m.aics.norm$cluster,-m.aics.norm$ppt.normal)])
  
  m.aics.norm$plot.id <- factor(m.aics.norm$plot.id,levels=plot.levels)
  
  ggplot(m.aics.norm,aes(x=plot.id,y=model,fill=AIC)) +
    geom_tile() +
    facet_grid(. ~ cluster,scales="free_x",space="free_x") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
  
  
  
  ###########################################################################################
  ####### Get each plot's AIC drop for all candidate models (random intercept only) #########
  ###########################################################################################
  
  # we have trees.clim.rwi to work with
  plots <- unique(trees.clim.rwi$plot.id)
  
  m.aics <- data.frame()
  
  
  for(plot in plots) {
    
    d.plot <- subset(trees.clim.rwi, (plot.id == plot) )
    
    m <- list()
    
    m[["n"]] <- try(lmer(rwi ~ 1 + (1|tree.id),data=d.plot))
    m[["p1"]] <- try(lmer(rwi ~ ppt.z + (1|tree.id),data=d.plot))
    m[["p2"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + (1|tree.id),data=d.plot))
    m[["p3"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + ppt.z2 + (1|tree.id),data=d.plot))
    m[["t1"]] <- try(lmer(rwi ~ tmean.z + (1|tree.id),data=d.plot))
    m[["t2"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + (1|tree.id),data=d.plot))
    m[["t3"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + tmean.z2 + (1|tree.id),data=d.plot))
    m[["tp1"]] <- try(lmer(rwi ~ tmean.z + ppt.z + (1|tree.id),data=d.plot))
    m[["tp2"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + (1|tree.id),data=d.plot))
    m[["tp3"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + tmean.z2 + ppt.z2 + (1|tree.id),data=d.plot))
    m[["s1"]] <- try(lmer(rwi ~ snow.late.z + (1|tree.id),data=d.plot))
    m[["s2"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + (1|tree.id),data=d.plot))
    m[["s3"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + snow.late.z2 + (1|tree.id),data=d.plot))
    m[["dc.d1"]] <- try(lmer(rwi ~ dobr.def.z + (1|tree.id),data=d.plot))
    m[["dc.d2"]] <- try(lmer(rwi ~ dobr.def.z + dobr.def.z1 + (1|tree.id),data=d.plot))
    m[["dc.a1"]] <- try(lmer(rwi ~ dobr.aet.z + (1|tree.id),data=d.plot))
    m[["dc.a2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.aet.z1 + (1|tree.id),data=d.plot))
    m[["dc.da1"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + (1|tree.id),data=d.plot))
    m[["dc.da2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + dobr.aet.z1 + dobr.def.z1 + (1|tree.id),data=d.plot))
    m[["dp.d1"]] <- try(lmer(rwi ~ dobr.pub.def.z + (1|tree.id),data=d.plot))
    m[["dp.d2"]] <- try(lmer(rwi ~ dobr.pub.def.z + dobr.pub.def.z1 + (1|tree.id),data=d.plot))
    m[["dp.a1"]] <- try(lmer(rwi ~ dobr.pub.aet.z + (1|tree.id),data=d.plot))
    m[["dp.a2"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.aet.z1 + (1|tree.id),data=d.plot))
    m[["dp.da1"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.def.z + (1|tree.id),data=d.plot))
    m[["dp.da2"]] <- try(lmer(rwi ~ dobr.pub.aet.z + dobr.pub.def.z + dobr.pub.aet.z1 + dobr.pub.def.z1 + (1|tree.id),data=d.plot))
    m[["v.d1"]] <- try(lmer(rwi ~ vic.def.z + (1|tree.id),data=d.plot))
    m[["v.d2"]] <- try(lmer(rwi ~ vic.def.z + vic.def.z1 + (1|tree.id),data=d.plot))
    m[["v.a1"]] <- try(lmer(rwi ~ vic.aet.z + (1|tree.id),data=d.plot))
    m[["v.a2"]] <- try(lmer(rwi ~ vic.aet.z + vic.aet.z1 + (1|tree.id),data=d.plot))
    m[["v.da1"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + (1|tree.id),data=d.plot))
    m[["v.da2"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + vic.aet.z1 + vic.def.z1 + (1|tree.id),data=d.plot))
    
    AIC.if <- function(x) {
      
      if(class(x) == "try-error") {
        return(NA)
      }
      
      
      if(length(x@optinfo$conv$lme4$messages > 0)) { #then there was no convergence
        return(NA)
      } else {
        return(AIC(x))
      }
    }
    
    m.aics.plot <- unlist(lapply( m, FUN=AIC.if )  )
    m.aics.plot <- m.aics.plot - m.aics.plot["n"] # get difference from null
    m.aics.plot[m.aics.plot > 0] <- 0
    
    for(i in 1:length(m.aics.plot)){
      df <- data.frame(plot.id=plot,model=names(m.aics.plot[i]),AIC=m.aics.plot[i])
      m.aics <- rbind(m.aics,df)
    }
    
  }
  
  
  m.aics.norm <- merge(m.aics,plot.norm,by="plot.id")
  
  m.aics.norm$cluster <- NA
  m.aics.norm[m.aics.norm$y > 150000,"cluster"] <- "1Plumas"
  m.aics.norm[(m.aics.norm$y < 150000) & (m.aics.norm$y > 0),"cluster"] <- "2Tahoe"
  m.aics.norm[(m.aics.norm$y < 0) & (m.aics.norm$y > -70000),"cluster"] <- "3Yose"
  m.aics.norm[(m.aics.norm$y < -70000),"cluster"] <- "4Sierra"
  
  plot.levels <- unique(m.aics.norm$plot.id[order(m.aics.norm$cluster,-m.aics.norm$ppt.normal)])
  
  m.aics.norm$plot.id <- factor(m.aics.norm$plot.id,levels=plot.levels)
  
  ggplot(m.aics.norm,aes(x=plot.id,y=model,fill=AIC)) +
    geom_tile() +
    facet_grid(. ~ cluster,scales="free_x",space="free_x") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
  
  
  
  
  
  ###########################################################################################
  ####### Get AIC drop for all candidate models across all plots (random slopes by plot) #########
  ###########################################################################################
  
  #(does not converge)
  
  # we have trees.clim.rwi to work with
  plots <- unique(trees.clim.rwi$plot.id)
  
  m.aics <- data.frame()
  
  
  
  d.plot <- trees.clim.rwi
  
  m <- list()
  
  
  m[["n"]] <- try(lmer(rwi ~ 1 + (1|plot.id),data=d.plot))
  m[["p1"]] <- try(lmer(rwi ~ ppt.z + (ppt.z|plot.id),data=d.plot))
  m[["p2"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + (ppt.z + ppt.z1|plot.id),data=d.plot))
  m[["p3"]] <- try(lmer(rwi ~ ppt.z + ppt.z1 + ppt.z2 + (ppt.z + ppt.z1 + ppt.z2|plot.id),data=d.plot))
  m[["t1"]] <- try(lmer(rwi ~ tmean.z + (tmean.z|plot.id),data=d.plot))
  m[["t2"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + (tmean.z + tmean.z1|plot.id),data=d.plot))
  m[["t3"]] <- try(lmer(rwi ~ tmean.z + tmean.z1 + tmean.z2 + (tmean.z + tmean.z1 + tmean.z2|plot.id),data=d.plot))
  m[["tp1"]] <- try(lmer(rwi ~ tmean.z + ppt.z + (tmean.z + ppt.z|plot.id),data=d.plot))
  m[["tp2"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + (tmean.z + ppt.z + tmean.z1 + ppt.z1 |plot.id),data=d.plot))
  m[["tp3"]] <- try(lmer(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1 + tmean.z2 + ppt.z2 + (tmean.z + ppt.z + tmean.z1 + ppt.z1 + tmean.z2 + ppt.z2 |plot.id),data=d.plot))
  m[["s1"]] <- try(lmer(rwi ~ snow.late.z + (snow.late.z|plot.id),data=d.plot))
  m[["s2"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + (snow.late.z + snow.late.z1|plot.id),data=d.plot))
  m[["s3"]] <- try(lmer(rwi ~ snow.late.z + snow.late.z1 + snow.late.z2 + (snow.late.z + snow.late.z1 + snow.late.z2|plot.id),data=d.plot))
  m[["dc.d1"]] <- try(lmer(rwi ~ dobr.def.z + (dobr.def.z|plot.id),data=d.plot))
  m[["dc.d2"]] <- try(lmer(rwi ~ dobr.def.z + dobr.def.z1 + (dobr.def.z+ dobr.def.z1|plot.id),data=d.plot))
  m[["dc.a1"]] <- try(lmer(rwi ~ dobr.aet.z + (dobr.aet.z|plot.id),data=d.plot))
  m[["dc.a2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.aet.z1 + (dobr.aet.z+ dobr.aet.z1|plot.id),data=d.plot))
  m[["dc.da1"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + (dobr.aet.z+ dobr.def.z|plot.id),data=d.plot))
  m[["dc.da2"]] <- try(lmer(rwi ~ dobr.aet.z + dobr.def.z + dobr.aet.z1 + dobr.def.z1 + (dobr.aet.z + dobr.def.z +dobr.aet.z1 + dobr.def.z1|plot.id),data=d.plot))
  m[["v.d1"]] <- try(lmer(rwi ~ vic.def.z + (vic.def.z|plot.id),data=d.plot))
  m[["v.d2"]] <- try(lmer(rwi ~ vic.def.z + vic.def.z1 + (vic.def.z+ vic.def.z1|plot.id),data=d.plot))
  m[["v.a1"]] <- try(lmer(rwi ~ vic.aet.z + (vic.aet.z|plot.id),data=d.plot))
  m[["v.a2"]] <- try(lmer(rwi ~ vic.aet.z + vic.aet.z1 + (vic.aet.z+ vic.aet.z1|plot.id),data=d.plot))
  m[["v.da1"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + (vic.aet.z+ vic.def.z|plot.id),data=d.plot))
  m[["v.da2"]] <- try(lmer(rwi ~ vic.aet.z + vic.def.z + vic.aet.z1 + vic.def.z1 + (vic.aet.z + vic.def.z +vic.aet.z1 + vic.def.z1|plot.id),data=d.plot))
  
  
  
  
  
  AIC.if <- function(x) {
    
    if(class(x) == "try-error") {
      return(NA)
    }
    
    
    if(length(x@optinfo$conv$lme4$messages > 0)) { #then there was no convergence
      return(NA)
    } else {
      return(AIC(x))
    }
  }
  
  m.aics.plot <- unlist(lapply( m, FUN=AIC.if )  )
  m.aics.plot <- m.aics.plot - m.aics.plot["n"] # get difference from null
  m.aics.plot[m.aics.plot > 0] <- 0
  
  
  
  
  
  ####################################################################
  #### Calc each tree's sens to an increase of 1 sd in predictors ####
  ####################################################################
  library(stats)
  
  trees <- unique(trees.clim.rwi$tree.id)
  
  
  ### precip
  sens <- data.frame()
  for(tree in trees) {
    
    tree.clim.rwi <- trees.clim.rwi[trees.clim.rwi$tree.id == tree,]
    
    if(sum(!is.na(tree.clim.rwi$rwi))< 25) next() #skip trees with < 25 years of data
    
    
    m <- lm(rwi~ dobr.def.z + dobr.def.z1 + dobr.def.z2,data=tree.clim.rwi)
    newdata <- data.frame(dobr.def.z=1,dobr.def.z1=0,dobr.def.z2=0)
    
    #m <- lm(rwi~ ppt.z,data=tree.clim.rwi)
    #newdata <- data.frame(ppt.z=1)
    
    pred <- predict(m,newdata,se.fit=FALSE,interval="confidence",level=0.95)
    
    sens.tree <- data.frame(tree.id=tree,pred)
    sens <- rbind(sens,sens.tree)
  
  }
  
  
  
  # get species column
  trees.sens.norm <- merge(tree.norm[,c("tree.id","plot.id","species")],sens,by="tree.id")
  
  
  
  
  
  
  # divide into PSME and OTHER
  #trees.fits.norm.sp$species <- ifelse(trees.fits.norm.sp$species=="PSME","PSME","OTHER")
  #trees.fits.norm.sp <- trees.fits.norm.sp[trees.fits.norm.sp$species=="PSME",]
  
  #compute plot-level averages for these vars
  
  plots.sens.norm <- aggregate(trees.sens.norm,by=list(trees.sens.norm$plot.id,trees.sens.norm$species),FUN=mean)
  plots.sens.norm$plot.id <- plots.sens.norm$Group.1
  plots.sens.norm$species <- plots.sens.norm$Group.2
  
  #get coordinates
  plots.sens.norm <- merge(plots.sens.norm,plot.norm,by="plot.id")
  
  
  #### Designate as plumas, tahoe, yose, or sierra
  plots.sens.norm$cluster <- NA
  plots.sens.norm[plots.sens.norm$y > 150000,"cluster"] <- "1Plumas"
  plots.sens.norm[(plots.sens.norm$y < 150000) & (plots.sens.norm$y > 0),"cluster"] <- "2Tahoe"
  plots.sens.norm[(plots.sens.norm$y < 0) & (plots.sens.norm$y > -70000),"cluster"] <- "3Yose"
  plots.sens.norm[(plots.sens.norm$y < -70000),"cluster"] <- "4Sierra"
  
  plots.fits <- as.data.frame(plots.sens.norm)
  
  library(ggplot2)
  
  #thin to only PSME and PIPO for now
  #plots.fits <- subset(plots.fits,species %in% c("PSME","PIPO"))
  plots.sens.norm <- subset(plots.sens.norm,species=="PSME")
  
  library(gridExtra)
  
  plot.levels <- unique(plots.sens.norm$plot.id[order(plots.sens.norm$cluster,plots.sens.norm$dobr.def)])
  
  plots.sens.norm$plot.id <- factor(plots.sens.norm$plot.id,levels=plot.levels)
  
  ggplot(plots.sens.norm,aes(x=plot.id,y=fit,colour=dobr.def.normal)) +
    geom_point(aes(size=rsq)) +
    geom_errorbar(aes(ymin=lwr,ymax=upr))+
    facet_grid(. ~ cluster,scales="free_x",space="free_x") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
    geom_rug(sides="b",size=4)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  m1 <- lm(rwi ~ ppt.z + ppt.z1 + ppt.z2, data=trees.clim.rwi)
  m2 <- lm(rwi ~ tmean.z + tmean.z1 + tmean.z2, data=trees.clim.rwi)
  m3 <- lm(rwi ~ dobr.aet.z + dobr.aet.z1 + dobr.aet.z2 + dobr.def.z + dobr.def.z1 + dobr.def.z2, data=trees.clim.rwi)
  m4 <- lm(rwi ~ dobr.def.z + dobr.def.z1 + dobr.def.z2, data=trees.clim.rwi)
  m5 <- lm(rwi ~ vic.aet.z + vic.aet.z1 + vic.aet.z2, data=trees.clim.rwi)
  m6 <- lm(rwi ~ vic.def.z + vic.def.z1 + vic.def.z2, data=trees.clim.rwi)
  m7 <- lm(rwi ~ tw.aet.z + tw.aet.z1 + tw.aet.z2, data=trees.clim.rwi)
  m8 <- lm(rwi ~ tw.def.z + tw.def.z1 + tw.def.z2, data=trees.clim.rwi)
  m9 <- lm(rwi ~ snow.late.z + snow.late.z1 + snow.late.z2, data=trees.clim.rwi)
  
  AIC(m1,m2,m3,m4,m7,m8)
  
  summary(m1)
  summary(m2)
  summary(m3)
  summary(m4)
  summary(m5)
  summary(m6)
  summary(m9)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##################################
  ### Initial explorations (old) ###
  ##################################
  
  m <- lm(rwi ~ tmean.z + ppt.z + tmean.z1 + ppt.z1, data=trees.clim.rwi)
  m <- lm(rwi ~ vic.aet.z + vic.def.z + vic.aet.z1 + vic.def.z1, data=trees.clim.rwi) # temp, ppt, vic, dob all perform similarly; tw is bad.
  m <- lm(rwi ~ ppt.z + ppt.z1, data=trees.clim.rwi) # temp, ppt, vic, dob all perform similarly; tw is bad.
  m.plot <- lmer(rwi0 ~ ppt.z + (ppt.z|plot.id), data=trees.clim.rwi, start=c(0.2,0.2,0.2))
  m.plot <- lmer(rwi0 ~ snow.late.z + (snow.late.z|plot.id), data=trees.clim.rwi, start=c(0.2,0.2,0.2))
  m.plot.tree <- lmer(rwi ~ ppt.z + (ppt.z|plot.id/tree.id), data=trees.clim.rwi, start=c(0.2,0.2,0.2,0,0,0))
  m.plot.tree <- lmer(rwi0 ~ snow.late.z + (snow.late.z|plot.id/tree.id), data=trees.clim.rwi, start=c(0.2,0.2,0.2,0,0,0))
  
  m.tree <- lmer(rwi0 ~ ppt.z + (ppt.z|tree.id), data=trees.clim.rwi, start=c(0.2,0.2,0.2))
  
  
  #r.squared(m)
  
  m.coef <- coef(m.plot.tree)$plot.id
  names(m.coef) <- paste0("b.",names(m.coef))
  m.coef$plot.id <-rownames(m.coef)
  
  #merge the plot-level random effect of the coefficient and the plot-level normal climate
  plot.coef.clim <- merge(m.coef,plot.norm,by="plot.id")
  
  plot(b.ppt.z~ppt.normal,data=plot.coef.clim)
  
  m.l2 <- lm(b.ppt.z~ppt.normal,data=plot.coef.clim)
  summary(m.l2)
  
  
  
  
  #### test whether absolute precip amount, or absolute precip anom, or proportional precip anom, or precip z-score
  
  
  
  
  plot(rwi ~ (ppt.prop), data=trees.clim.rwi)
  
  
  hist(trees.clim.rwi$ppt.a) #long right tail?
  hist(trees.clim.rwi$ppt.z) # long right tail? maybe not
  hist(trees.clim.rwi$ppt.prop) # long right tail? but might be the best of them all
  hist(trees.clim.rwi$ppt.logprop) # long left tail?
  hist(trees.clim.rwi$ppt.logprop.z) # long left tail?
  
  
  
  #! explore distribution of prop departure and log prop departure for individual cores: symmetrical?
  
  d.tree <- trees.clim.rwi[trees.clim.rwi$tree.id=="4007",]
  plot(rwi ~ ppt.logprop.z,data=d.tree)
  
  
  
  
  ## see if one species correlates better
  sp.clim.rwi <- trees.clim.rwi[trees.clim.rwi$species=="PSME",]
  #sp.clim.rwi <- trees.clim.rwi
  
  m.abs <- lm(rwi ~ ppt + ppt1 + ppt2, data=sp.clim.rwi)
  m.a.abs <- lm(rwi ~ ppt.a + ppt.a1 + ppt.a2, data=sp.clim.rwi)
  m.prop <- lm(rwi ~ ppt.prop + ppt.prop2, data=sp.clim.rwi)
  m.logprop <- lm(rwi ~ ppt.logprop + ppt.logprop2, data=sp.clim.rwi)
  m.z <- lm(log(rwi) ~ ppt.z + ppt.z1 + ppt.z2, data=sp.clim.rwi)
  m.logprop.z <- lm(log(rwi) ~ ppt.logprop.z + ppt.logprop.z1 + ppt.logprop.z2, data=sp.clim.rwi)
  
  summary(m.abs)
  summary(m.a.abs)
  summary(m.prop)
  summary(m.logprop)
  summary(m.z)          ## winner!
  summary(m.logprop.z)  ## also winner!
  
  # repeat for vic def
  m.abs <- lm(rwi ~ vic.def, data=sp.clim.rwi)
  m.a.abs <- lm(rwi ~ vic.def.a + vic.def.a1 + vic.def.a2, data=sp.clim.rwi)
  m.prop <- lm(rwi ~ vic.def.prop + vic.def.prop2 + vic.def.prop3, data=sp.clim.rwi)
  m.logprop <- lm(rwi ~ vic.def.logprop + vic.def.logprop2 + vic.def.logprop3, data=sp.clim.rwi)
  m.z <- lm(rwi ~ vic.def.z + vic.def.z1 + vic.def.z2 + vic.aet.z + vic.aet.z1 + vic.aet.z2, data=sp.clim.rwi)
  m.logprop.z <- lm(rwi ~ vic.def.logprop.z + vic.def.logprop.z1 + vic.def.logprop.z2, data=sp.clim.rwi)
  
  summary(m.a.abs)
  summary(m.prop)
  summary(m.logprop)
  summary(m.z)         
  summary(m.logprop.z)  
  
  # repeat for dobr def
  m.abs <- lm(rwi ~ dobr.def, data=sp.clim.rwi)
  m.a.abs <- lm(rwi ~ dobr.def.a + dobr.def.a1 + dobr.def.a2, data=sp.clim.rwi)
  m.prop <- lm(rwi ~ dobr.def.prop + dobr.def.prop2 + dobr.def.prop3, data=sp.clim.rwi)
  m.logprop <- lm(rwi ~ dobr.def.logprop + dobr.def.logprop2 + dobr.def.logprop3, data=sp.clim.rwi)
  m.z <- lm(rwi ~ dobr.def.z + dobr.def.z1 + dobr.def.z2 + dobr.aet.z + dobr.aet.z1 + dobr.aet.z2, data=sp.clim.rwi)
  m.logprop.z <- lm(rwi ~ dobr.def.logprop.z + dobr.def.logprop.z1 + dobr.def.logprop.z2, data=sp.clim.rwi)
  
  summary(m.a.abs)
  summary(m.prop)
  summary(m.logprop)
  summary(m.z)         
  summary(m.logprop.z)  
  
  
  ##### Summary of results so far: ppt is by far the best variable (with dobr or vic coming close for PSME), with the exception that for PILA, AET is better
  ## For representing PPT departure from normal, z-score or logprop.z seems to be the best
}

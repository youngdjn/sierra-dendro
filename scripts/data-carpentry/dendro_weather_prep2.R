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

trees$dbh <- as.numeric(as.character(trees$dbh))


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

# remove A, B, Z, or T from core names to get tree.ids (this should not be necessary because the names should already be the tree IDs, not the core ids)
names(chron.clean$chron) <- gsub("a?A?b?B?t?T?z?Z?","",names(chron.clean$chron))
names(chron$chron) <- gsub("a?A?b?B?t?T?z?Z?","",names(chron$chron))

chron.rwi <- spline.na.rm(chron.clean$chron)
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


# estimate the "internal" radius from the external radius based on the observed empirical relationship (see below for linear model)
radius.external <- 31.2315 + as.numeric(trees$dbh)*3.3216

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

inf.rows = which(clim.rwi$bai.ba.comb == Inf)
clim.rwi = clim.rwi[-inf.rows,]


# for export purposes (e.g. comparing with DBH,) internal radius only counts for trees that were measured to the outermost ring
radii <- radii[radii$tree.id %in% trees.meas.tobark,]
# 
# # merge in the internal radius
# clim.rwi <- merge(clim.rwi,radii,by="tree.id",all.x=TRUE)

# merge in tree survey data
trees$tree.id <- toupper(trees$tree.id)
trees.clim.rwi <- merge(clim.rwi,trees,by="tree.id")
## add internal radius
trees.clim.rwi <- left_join(trees.clim.rwi,radii,by="tree.id")

# merge in competition data
voronoi <- read.csv("data/plot-and-tree/processed/voronoi_area.csv",header=TRUE,stringsAsFactors=FALSE)
voronoi$tree.id <- toupper(voronoi$tree.id)
trees.clim.rwi <- merge(trees.clim.rwi,voronoi,by="tree.id",all.x=TRUE)

trees.clim.rwi = trees.clim.rwi %>%
  mutate_at(vars(ppt:ba.prev.comb),funs(signif))

#### output year-level data
clim.rwi.out = clim.rwi %>%
  filter(year > 1959) %>%
  mutate_at(vars(-tree.id,-year),funs(signif)) %>% # truncate to 6 digits
  select(-(rad.tot:rad.06))
write.csv(clim.rwi.out,"data/compiled-for-analysis/years.csv",row.names=FALSE)

#### output tree-level data
trees.out.pre.pre = trees %>%  ## filter to only trees for which we have ring data
  filter(tree.id %in% clim.rwi$tree.id)
trees.out.pre = left_join(trees.out.pre.pre,voronoi,by="tree.id")
trees.out = left_join(trees.out.pre,radii,by="tree.id")

## add the external-radius
radius.ext <- data.frame(tree.id = names(radius.external),radius.external = radius.external)
trees.out <- left_join(trees.out,radius.ext,by="tree.id")

write.csv(trees.out,"data/compiled-for-analysis/trees.csv",row.names=FALSE)


### determine how internal radius is related to DBH
m <- lm(radius~dbh,data=trees.out[trees.out$dbh < 100,])
summary(m)
plot(radius~dbh,data=trees.out[trees.out$dbh < 100,])



#### Calculate plot-level variables, including average plot location and normal climate averaged across 1981 to 2010 and all trees in a plot ####

### Normal climate and location (plot-level)
plot.normal.climate <- trees.clim.rwi %>%
  filter(year %in% 1981:2010) %>%
  group_by(plot.id) %>%
  dplyr::summarize(ppt.norm = mean(ppt),
         tmean.norm = mean(tmean),
         x.plot = mean(x),
         y.plot = mean(y),
         rad.tot = mean(rad.tot),
         rad.03 = mean(rad.03),
         rad.06 = mean(rad.06),
         voronoi.area.mean = mean(voronoi.area,na.rm=TRUE))


### Number of trees per plot by species
plot.ntrees <- trees.out %>%
  group_by(plot.id) %>%
  dplyr::summarize(n.trees = n(),
            n.psme = sum(species == "PSME"),
            n.pipo = sum(species == "PIPO"),
            n.abco = sum(species == "ABCO"),
            n.pila = sum(species == "PILA"))

### Merge plot-level data
plot.data = left_join(plot.normal.climate,plot.ntrees,by="plot.id")


### Determine which cluster
plot.data$cluster <- NA
plot.data[which(plot.data$y.plot > 150000),"cluster"] <- "Plumas"
plot.data[which((plot.data$y.plot < 150000) & (plot.data$y.plot > 0)),"cluster"] <- "Tahoe"
plot.data[which((plot.data$y.plot < 0) & (plot.data$y.plot > -70000)),"cluster"] <- "Yose"
plot.data[which((plot.data$y.plot < -70000)),"cluster"] <- "Sierra"

### Export data
write.csv(plot.data,"data/compiled-for-analysis/plots.csv",row.names=FALSE)



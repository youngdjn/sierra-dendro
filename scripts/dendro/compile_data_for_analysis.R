setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro") # Derek on Derek's computer
setwd("~/Research projects/Sierra dendro/sierra-dendro") # Derek on Latim-GIS-S




# load the summarize.cluster() function and several other convenience functions
source("scripts/dendro/dendro-functions/summarize_cluster.R")


library(dplyr)
library(purrr)
library(stringr)
library(reshape)

leads <- function(var, n=10){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(lead(!!var, !!.x)) ) %>% 
    set_names(sprintf("lead_%s_%02d", rlang::quo_text(var), indices))
  
}


tree.ring.folder = "data/dendro/coorecorder-measurements" # on repository





#### Read in tree data ####

# get list of trees with plots, and other ancillary data
trees <- read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE,stringsAsFactors=FALSE)

##drop trees that were resampled under a different ID
trees$tree.id <- toupper(trees$tree.id)
trees$former.id <- toupper(trees$former.id)
former.ids <- trees$former.id
trees <- trees[!(trees$tree.id %in% former.ids),]

trees$dbh <- as.numeric(as.character(trees$dbh))

#do any trees have duplicate samples?
ndup <- sum(duplicated(trees$tree.id))
if(ndup > 0) {warning("Some tree ids duplicated in trees_loc spreadsheet.")}



#### Read in competition data ####
# merge in competition data
voronoi <- read.csv("data/plot-and-tree/processed/voronoi_area.csv",header=TRUE,stringsAsFactors=FALSE)
voronoi$tree.id <- toupper(voronoi$tree.id)




#### Read in and clean ring width data ####


chron_rwi = NULL
chron_raw = NULL
truncs = NULL
radii = NULL
ages = NULL

focal_clusters = c("SL") # c("NL","NH", "SL","SH") #,"NL","NH")

for(focal_cluster in focal_clusters) {
    
  chron_out = open_cluster_chron_truncate_window(focal_cluster=focal_cluster,
                                                 trunc_unjustified=TRUE,
                                                 window_width=12,
                                                 cor_threshold = 0.1,
                                                 max_bad_years = 15,
                                                 min_good_years = 10)
  
  if(is.null(chron_rwi)) {
    chron_rwi = chron_out$chron_rwi
    chron_raw = chron_out$chron_raw
  } else {
    chron_rwi = full_join(chron_rwi, chron_out$chron_rwi,by="year")
    chron_raw = full_join(chron_raw, chron_out$chron_raw,by="year")
  }
  truncs = c(truncs,chron_out$truncs)
  radii = rbind(radii,chron_out$radii)
  ages = c(ages,chron_out$ages)

}

radii2 = radii$radius
names(radii2) = radii$tree.id
radii.old = radii
radii = radii2


#### Compute basal area / BAI ####

### Using internal radius, where available
bai.ba.results <- ba.bai.calc(chron_raw %>% select(-year),radii)
ba <- bai.ba.results$ba
bai <- bai.ba.results$bai
bai.ba <- bai.ba.results$bai.ba
ba.prev <- bai.ba.results$ba.prev
names(ba) <- names(bai) <- names(bai.ba) <- names(ba.prev) <- names(chron_raw %>% select(-year)) 
ba$year <- rownames(ba)
bai$year <- rownames(bai)
bai.ba$year <- rownames(bai.ba)
ba.prev$year <- rownames(ba.prev)

### Using external radius
# estimate the "internal" radius from the external radius based on the observed empirical relationship (see below for linear model)
radius.external <- 31.2315 + as.numeric(trees$dbh)*3.3216
names(radius.external) <- trees$tree.id
## exclude trees that did not have ring measurements to the outside (measured diam/radius not representative of outermost ring)
outer.rings <- chron_raw[chron_raw$year %in% c(2012,2013,2014),] %>% select(-year)
trees.meas.tobark <- colnames(outer.rings[(colSums(outer.rings,na.rm=TRUE) > 0)])
chron.bacalc <- chron_raw[,trees.meas.tobark]
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


#### Add climate data and merge into chron, ba, radius data and tree survey data ####

trees.clim <- read.csv("data/plot-and-tree/processed/trees_climate.csv",header=TRUE,stringsAsFactors=FALSE)

# melt chron
library(reshape)
rwi.melt <- melt(chron_rwi,id.vars="year",variable_name="tree.id")
names(rwi.melt)[3] <- "rwi"

raw.melt = melt(chron_raw,id.vars="year",variable_name="tree.id")
names(raw.melt)[3] <- "raw_width"

#internal radius BA data
ba.melt <- melt(ba,id.vars="year",variable_name="tree.id")
bai.melt <- melt(bai,id.vars="year",variable_name="tree.id")
bai.ba.melt <- melt(bai.ba,id.vars="year",variable_name="tree.id")
ba.prev.melt <- melt(ba.prev,id.vars="year",variable_name="tree.id")
names(ba.melt)[3] <- "ba"
names(bai.melt)[3] <- "bai"
names(bai.ba.melt)[3] <- "bai.ba"
names(ba.prev.melt)[3] <- "ba.prev"

#external radius BA data
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
clim.rwi <- merge(clim.rwi,raw.melt,by=c("tree.id","year"))
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
if(length(inf.rows)>0) {
  clim.rwi = clim.rwi[-inf.rows,]
}


# for export purposes (e.g. comparing with DBH,) internal radius only counts for trees that were measured to the outermost ring
# it's ok to use it for internal BAI-BA calc, but there's no other use
radii = data.frame(tree.id = names(radii),radius=radii)
radii <- radii[radii$tree.id %in% trees.meas.tobark,]

# merge in tree survey data
trees$tree.id <- toupper(trees$tree.id)
trees.clim.rwi <- merge(clim.rwi,trees,by="tree.id")
## add internal radius
trees.clim.rwi <- left_join(trees.clim.rwi,radii,by="tree.id")

# merge in competition data
voronoi$tree.id <- toupper(voronoi$tree.id)
trees.clim.rwi <- merge(trees.clim.rwi,voronoi,by="tree.id",all.x=TRUE)

trees.clim.rwi = trees.clim.rwi %>%
  mutate_at(vars(ppt:ba.prev.comb),funs(signif))


#### output year-level data
clim.rwi.out = clim.rwi %>%
  filter(year > 1800) %>%
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


#### Compute plot-level data ####
#Calculate plot-level variables, including average plot location and normal climate averaged across 1981 to 2010 and all trees in a plot 

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


### Determine which latitudinal cluster
plot.data$cluster <- NA
plot.data[which(plot.data$y.plot > 150000),"cluster"] <- "Plumas"
plot.data[which((plot.data$y.plot < 150000) & (plot.data$y.plot > 0)),"cluster"] <- "Tahoe"
plot.data[which((plot.data$y.plot < 0) & (plot.data$y.plot > -70000)),"cluster"] <- "Yose"
plot.data[which((plot.data$y.plot < -70000)),"cluster"] <- "Sierra"

### Export data
write.csv(plot.data,"data/compiled-for-analysis/plots.csv",row.names=FALSE)





#### Relate internal radius to DBH ####
m <- lm(radius~dbh,data=trees.out[trees.out$dbh < 100,])
summary(m)
plot(radius~dbh,data=trees.out[trees.out$dbh < 100,])



#### Summary and diagnostic statistics ####

## For each core, get its start year, end year, nyears, trunc reason, plot, and species

# which cores have no data?
cores.nodata = clim.rwi.out %>%
  group_by(tree.id) %>%
  summarize(non_na = sum(!is.na(raw_width))) %>%
  filter(non_na == 0)

tree.ids.nodata = cores.nodata$tree.id

cores.summary = clim.rwi.out %>%
  filter(!is.na(raw_width)) %>%
  group_by(tree.id) %>%
  summarize(last_year = max(year),
            first_year = min(year),
            nyears = last_year-first_year) 

if(length(tree.ids.nodata) > 0) {
  trees.nodata.summary = data.frame(tree.id=tree.ids.nodata,last_year=NA,first_year=NA,nyears=0)
  cores.summary = bind_rows(cores.summary,trees.nodata.summary)
}


## pull in cores that were truncated to no rings due to poor align throughout, mult potential alignments.

truncs.df = data.frame(tree.id = names(truncs),trunc=truncs)

cores_summary = full_join(cores.summary,truncs.df) %>%
  mutate(nyears = ifelse(is.na(nyears),-1,nyears)) %>%
  mutate(nyears = nyears+1)

trees_for_summary = trees %>%
  select(tree.id,species, dbh, plot.id, cluster)

cores_summary = left_join(cores_summary,trees_for_summary,by="tree.id") %>%
  select(cluster,plot.id,species,tree.id,everything()) %>%
  arrange(cluster,plot.id,species,tree.id)

## now make a plot summary

plot_summary = cores_summary %>%
  group_by(cluster,plot.id,species) %>%
  summarize(over_15y = sum(nyears>15),
            over_30y = sum(nyears>30),
            over_50y = sum(nyears>50))




setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro/") # Derek

## this script merges tree-level with plot-level data to compute a summary data table with tree, species, plot, and cluster, 

library(tidyverse)

trees = read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE)

lat.cutoff <- 50000
trees.s <- trees[trees$y < lat.cutoff,]
trees.n <- trees[trees$y >= lat.cutoff,]
n.mid.elev <- mean(c(min(trees.n$elev),max(trees.n$elev)))
s.mid.elev <- mean(c(min(trees.s$elev),max(trees.s$elev)))

trees$cluster <- NA
trees[(trees$y < lat.cutoff) & (trees$elev < s.mid.elev),]$cluster <- "SL"
trees[(trees$y < lat.cutoff) & (trees$elev >= s.mid.elev),]$cluster <- "SH"
trees[(trees$y >= lat.cutoff) & (trees$elev < n.mid.elev),]$cluster <- "NL"
trees[(trees$y >= lat.cutoff) & (trees$elev >= n.mid.elev),]$cluster <- "NH"

write.csv(trees,"data/plot-and-tree/processed/trees_loc.csv",row.names=FALSE)

trees = trees %>%
  dplyr::select(tree.id,former.id,species,dbh,x,y,elev,plot.id,cluster) %>%
  mutate_at(.vars =vars(x,y,elev),.funs=funs(round))

write.csv(trees,"data/plot-and-tree/processed/trees_loc_simple.csv",row.names=FALSE)

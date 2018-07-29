source("scripts/dendro/dendro-functions/read_pos_extended_bai.R")
#source("scripts/dendro/dendro-functions/read_google_sheet.R")

list.plots <- function() {
  
  trees <- read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE,stringsAsFactors=FALSE)
  #cat(unique(trees$plot.id))
  return(unique(trees$plot.id))
  
}



summarize.cluster <- function(cluster,type="cluster",name=cluster,clean.ref.chron=FALSE) {

  name = paste(name,collapse="") # in case a vector of core names was provided for the cluster and there was no additional name provided
  
  trees <- read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE,stringsAsFactors=FALSE)
  trees$tree.id <- toupper(trees$tree.id)
  trees$former.id <- toupper(trees$former.id)
  
  ###define clusters
  #open DEM
  
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
  
  
  if(type != "tree" && cluster=="ALL") {
    # open all trees
    trees.cluster <- trees
  # } else if (cluster == "OTHER") {
  #   # open all cores that have a POS file but are not in the trees datasheet
  #   files <- list.files(path="CooRecorder measurements",pattern=".pos",full.names=FALSE)
  #   files <- toupper(files)
  #   # take the filename before the t or . or A or B or Z
  #   files.split <- sapply(files,function(x) strsplit(x,"[ABZT.]+")[[1]][1])
  #   unmatched.indexes <- which(!((files.split %in% trees$tree.id) | (files.split %in% toupper(trees$former.id))))
  #   unmatched.pos.ids <- (files.split[unmatched.indexes])
  #   unmatched.pos.ids <- unmatched.pos.ids[!duplicated(unmatched.pos.ids)]
  #   trees.cluster <- trees[FALSE,]
  #   trees.cluster <- data.frame(tree.id = as.character(unmatched.pos.ids),former.id="",species=NA,ref.chron=NA)
  #   trees.cluster$tree.id <- as.character(trees.cluster$tree.id)
  } else {
    # open a chronology of all trees in the selected cluster or plot
    
    # in case the script is being run to crossdate specific trees or cores or clusters, check to see what trees and cores match
    trees.cluster <- trees[(trees$cluster %in% cluster) | (trees$plot.id %in% cluster) | (trees$tree.id %in% cluster),]
  }
  
  
  ## look up the core numbers corresponding to those trees

  tree.id.patterns = paste(trees.cluster$tree.id,"A?B?Z?t?$",sep="")
  former.id.patterns = paste(trees.cluster$former.id,"A?B?Z?t?$",sep="")
  
  patterns =cbind(tree.id.patterns,former.id.patterns)
  
  patterns[(patterns == "NAA?B?Z?t?$") | (patterns == "A?B?Z?t?$")] <- NA
  
  patterns = toupper(patterns)
  patterns = as.vector(patterns)
  patterns = patterns[!is.na(patterns)]
  pattern = paste(patterns,collapse="|")
  
  files = list.files(path=tree.ring.folder,pattern=".pos$",full.names=FALSE,ignore.case=TRUE)
  
  ## take the extension and "t" off of all of them
  core.names <- toupper(sapply(files,function(x) strsplit(x,"[tT]*.pos$")[[1]][1]))
  
  ## find which core names (from the pos files) match trees
  cores.match.trees = grepl(pattern,core.names)
  cores.cluster <- core.names[cores.match.trees]
  #cores.cluster <- core.names
    
  open.chron.results <- open.chron(cores.cluster,cores.cluster,unc.stop=TRUE,core.name=TRUE)
  
  ##!! NOTE that this script does not clean the ring-width series based on the crossdating progress records (e.g., areas of poor correlation to a reference chronology)
  
  if(clean.ref.chron == TRUE) {
    open.chron.results = clean.chron(open.chron.results)
    chron <- open.chron.results[["chron"]]
    trunc <- open.chron.results[["trunc"]]
    nrings <- open.chron.results[["nyears"]]
  } else {
    chron <- open.chron.results[["chron"]]
    trunc <- open.chron.results[["trunc"]]
    nrings <- open.chron.results[["nrings"]]
  }
  
  chron.detrended <- spline.na.rm(chron)
  
  
  #### make a list stating whether dated
  # for each core, get the years in which there are rings
  chron.years <- ifelse(!is.na(chron.detrended),rownames(chron.detrended),NA)
  chron.years <- chron.years < 2013
  chron.dated <- apply(chron.years,2,any,na.rm=TRUE)

  dated.cores.list <- (as.data.frame(chron.dated))
  names(dated.cores.list) <- "dated"
  dated.cores.list$sample <- rownames(dated.cores.list)
  dated.cores.list$tree.id <- toupper(sapply(dated.cores.list$sample,function(x) strsplit(x,"[ABZTabzt]+")[[1]][1]))
  
  #if there is a new tree ID for a former tree, use that instead
  dated.cores.list.newnames <- suppressWarnings(merge(dated.cores.list,trees.cluster[,c("former.id","tree.id")],by.x="tree.id",by.y="former.id",all.x=TRUE))
  names(dated.cores.list.newnames) <- c("tree.id","dated","sample","new.id")
  dated.cores.list.newnames$tree.id <- ifelse(is.na(dated.cores.list.newnames$new.id),dated.cores.list.newnames$tree.id,dated.cores.list.newnames$new.id)
  
  
  # look up cluster for each tree
  dated.cores.list.2 <- merge(dated.cores.list.newnames,trees.cluster,by.x="tree.id",by.y="tree.id",all.x=TRUE)
  
  #write.csv(dated.cores.list.2,"Ancillary data/dated_cores.csv",row.names=FALSE)
  #dated.cores.list.2 <- read.csv("Ancillary data/dated_cores.csv",header=TRUE)
  
  undated.list <- as.character(dated.cores.list.2[dated.cores.list.2$dated==FALSE,]$sample)
  dated.list <- as.character(dated.cores.list.2[dated.cores.list.2$dated==TRUE,]$sample)
  
  chron.dated <- chron.detrended[,dated.list]
  chron.dated <- chron.dated[as.numeric(rownames(chron.dated)) < 2015,] # remove the 3000's chrons (which means they don't have a date)
  
  corr <- corr.rwl.seg(chron.dated,seg.length=20,bin.floor=14,prewhiten=TRUE,pcrit=0.1,label.cex=.7,main=paste("All dated cores in cluster",name))
  
  # use only cores that had overall rho > 0.4
  corr.overall <- as.data.frame(corr$overall)
  ref.cores <- rownames(corr.overall[corr.overall$rho > 0.4,])
  
  if(length(ref.cores) < 4) stop("Too few cores in reference chronology. That is, too few cores that correlate well with the mean across all cores in the cluster.")
  
  chron.dated.ref <- chron.dated[,ref.cores]
  corr.ref <- corr.rwl.seg(chron.dated.ref,seg.length=20,bin.floor=14,pcrit=0.10,prewhiten=TRUE,label.cex=.7,main=paste("Cores in reference chronology for",name))
  
  # add a col to the dated cores list for whether it is included in master chron (i.e., rho > 0.5)
  dated.cores.list.2$ref.chron <- ifelse(dated.cores.list.2$sample %in% ref.cores,TRUE,FALSE)
  
# add a col to the dated cores list giving the length of the chron to truncation point
  trunc.nrings <- merge(trunc,nrings,by="row.names",all=TRUE)
  names(trunc.nrings) <- c("sample","truncated","nrings")
  dated.cores.list.3 <- merge(dated.cores.list.2,trunc.nrings,by="sample",all=TRUE)
  dated.cores.list.4 <- dated.cores.list.3[,c("tree.id","sample","species","ref.chron","dated","truncated","nrings")]
  dated.cores.list.4$sample.clean <- gsub("T","",dated.cores.list.4$sample)

  # load in comments from crossdating progress records
  
  crossdating.records <- read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv")

  
  if(length(crossdating.records) == 0) {
    
    cat("!!! Could not connect to Google spreadsheet.\n!!! Previously-processed cores cannot be displayed in cluster summary.\n")
    crossdating.records <- data.frame("Core number"=NA,Person=NA,Timestamp=NA)
    names(crossdating.records) <- c("Core number","Person","Timestamp")
  
  } else {
    crossdating.records <- crossdating.records[rowSums(is.na(crossdating.records)) != ncol(crossdating.records),] # remove empty rows

    # take only the most recent entry for each core
    crossdating.records <- crossdating.records[nrow(crossdating.records):1,] # reverse order
    crossdating.records.dups <- duplicated(crossdating.records$Core.number)
    crossdating.records <- crossdating.records[!crossdating.records.dups,]

  }

  crossdating.records$Core.number = toupper(crossdating.records$Core.number)

  cluster.summary <- merge(dated.cores.list.4,crossdating.records[,c("Core.number","Person","Timestamp","Core.status","Good.alignment.through.year","Other.notes","Derek.double.check.needed.")],by.x="sample.clean",by.y="Core.number",all.x=TRUE)
  names(cluster.summary) <- c("sample.clean","tree.id","sample","species","ref.chron","dated","truncated","nrings","completed.by","completed.date","previous crossdating results","Alignment good through","previous crossdating notes","previous crossdating notes2")



  #if more than one entry per core, use most recent (lowest down)
  cluster.summary <- cluster.summary[!duplicated(cluster.summary$sample,fromLast=TRUE),]



  # sort cluster summary
  if("ab" %in% cluster.summary$truncated)  cluster.summary$truncated <- relevel(cluster.summary$truncated,"ab")
  if("unc" %in% cluster.summary$truncated)  cluster.summary$truncated <- relevel(cluster.summary$truncated,"unc")
  cluster.summary <- cluster.summary[order(-cluster.summary$ref.chron,-cluster.summary$dated,cluster.summary$truncated,-cluster.summary$nrings),]
    
    


  # write cluster summary
  a <- suppressWarnings(try(write.csv(cluster.summary,paste("data/dendro/crossdating-summaries/",name,".csv",sep="",na=" "),row.names=FALSE),silent=TRUE))
  if(class(a) == "try-error") {
    stop("You must close the cluster_summary spreadsheet before running this function.\n")
  }
  
  ### make mean value reference chronology
  chron.dated.ref.mean <- chron(chron.dated.ref,prefix=substr(cluster,1,3),biweight=TRUE,prewhiten=FALSE)

  if(type == "plot" | type == "tree") {
    chron.dated.ref.mean <- chron.dated.ref.mean[chron.dated.ref.mean$samp.depth > 3,] #must have at least 4 cores to compare against
  } else   {
    chron.dated.ref.mean <- chron.dated.ref.mean[chron.dated.ref.mean$samp.depth > 8,] #must have at least 8 cores to compare against
  }
  
  
  
  
  if(nrow(chron.dated.ref.mean) == 0) stop("Too few cores in reference chronology.")

  # export mean value ref chron to CDendro
  filename <- paste("data/dendro/reference-chronologies/ref_chron_",name,".rwl",sep="")

  write.df <- data.frame(chron=chron.dated.ref.mean[,1])
  rownames(write.df) <- rownames(chron.dated.ref.mean)
  write.rwl(write.df,filename,format="tucson")
  
  
  ## see how far back the ref chron goes
  oldest.ref.chron <- min(as.numeric(rownames(chron.dated.ref.mean)))
  nyears.ref.chron <- nrow(chron.dated.ref.mean)
  cat("--------------------------------------------\n")
  cat("Reference chronology for ",name," goes to ",oldest.ref.chron," (",nyears.ref.chron," years)\n",sep="")
  cat("--------------------------------------------\n")
  
  
  
  #### Inspecting individual cores ####
  # core to inspect
  #focal.core <- "2006"
  
  #### Correlogram plotting ####
  # !! is it ok that the focal core may be a member of the reference series? yes because the ref series must have at least 5 members
  #corr.series.seg(chron.dated.ref.mean,series=chron[,focal.core],series.yrs=as.numeric(rownames(chron)),seg.length=20,bin.floor=0,prewhiten=FALSE)
  
  #### CCF plotting ####
  #ccf.series.rwl(chron.dated.ref.mean,series=chron[,focal.core],series.yrs=as.numeric(rownames(chron)),seg.length=20,bin.floor=0,prewhiten=FALSE)

  #CDendro-type plot
  # series.rwl.plot(chron.dated,series=names(chron.dated)[5],seg.length=10,bin.floor=0)
  
  
}





####NEXT STEPS

# More accurate way to tell if core is dated


# none of the old YOSE cores are making it into the best-correlated. Are they dated correctly?
# note that can only crossdate out to end of ref chron: could date a core but part of the chron could be bad if there were few cores to crossdate against
# will resolve unc's in ref chron to some historic date, then once those are resolved,  resolve the next set of uncs in reference chron. Then could either date undated chrons and then resolve their uncs, or attempt to fix poorly-correlating cores (dated cores that are not in the master chron)
# make sure the old YOSE cores were dated to 2012
#chron.dated.ref[rownames(chron.dated.ref) %in% c(2012,2013),]

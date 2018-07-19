library(plyr)
library(dplR)


pith.problems <- data.frame(sample=NA,problem=NA)
pos.problems <- data.frame(sample=NA,problem=NA)
unc.problems <- data.frame(sample=NA,problem=NA)
dated.problems <- data.frame(sample=NA,problem=NA)




split.line <- function(line) {
  
  split <- strsplit(line,"  | #|,")[[1]]
  
  xcoord1 <- split[1]
  ycoord1 <- split[2]
  
  #if the third element is a number with a decimal, then this is a dual point situation
  if(grepl("[0-9].[0-9]",split[3])) {
    xcoord2 <- split[3]
    ycoord2 <- split[4]
    comment <- split[5]
  } else {
    xcoord2 <- NA
    ycoord2 <- NA
    comment <- split[3]
  }
  
  xcoord1 <- as.numeric(xcoord1)
  ycoord1 <- as.numeric(ycoord1)
  xcoord2 <- as.numeric(xcoord2)
  ycoord2 <- as.numeric(ycoord2)
  comment <- as.character(comment)
  
  ret <- data.frame(xcoord1=xcoord1,ycoord1=ycoord1,xcoord2=xcoord2,ycoord2=ycoord2,comment=comment,stringsAsFactors=FALSE)
  return(ret)
  
}



read.pos <- function(sample,sample.secondary, unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE) {

  
  #cat("Reading pos for cores",sample,sample.secondary,"\n")
  
  ## for testing
  #sample <- "2146"
  
  
  ### !!!!  need to allow for multiple cores per tree when looking up filename
  #             how to find best one? or open both? store inferior one in other folder?
  
  ### !!! make sure working properly with point offsets
  
  ### !!!! need to add core ID for first-ring comment check

  pattern <- paste(sample,"A?B?Z?t?.pos$",sep="")
  filename <- list.files(path="data/dendro/coorecorder-measurements",pattern=pattern,full.names=TRUE,ignore.case=TRUE)
  
  if(length(filename)>1) {
    problem.string <- paste("More than one pos file matching:",sample,". Using first one")
    problem.row <- data.frame(sample=sample,problem=problem.string)
    pos.problems <<- rbind(pos.problems,problem.row)
    filename <- filename[1]
  }
  
  if(length(filename)<1) {
    
    if(!is.na(sample.secondary) & (sample.secondary != "")) {
      
      pattern.secondary <- paste(sample.secondary,"A?B?Z?t?.pos$",sep="")
      filename.secondary <- list.files(path="data/dendro/coorecorder-measurements",pattern=pattern.secondary,full.names=TRUE,ignore.case=TRUE)
      
      if(length(filename.secondary)<1) {
        problem.string <- paste("No pos files matching primary:",sample,"or secondary:",sample.secondary)
        cat(problem.string,"\n")
        problem.row <- data.frame(sample=sample,problem=problem.string)
        pos.problems <<- rbind(pos.problems,problem.row)
        return(NULL)
      } else {  # no pos file for the main sample ID, but there is one for the "former.id" column, so open it
        filename <- filename.secondary[1]
      }

    } else {
      
      problem.string <- paste("No pos files matching primary:",sample,"and there is no secondary.")
      cat(problem.string,"\n")
      problem.row <- data.frame(sample=sample,problem=problem.string)
      pos.problems <<- rbind(pos.problems,problem.row)
      return(NULL)
    }

  }
  
  
  file.core <- filename[1]
  pos.core <- strsplit(file.core,"/",fixed=TRUE)[[1]][4]
  name.core <- strsplit(pos.core,".",fixed=TRUE)[[1]][1]
  name.tree <- gsub("t","",name.core,fixed=TRUE)
  name.tree <- gsub("T","",name.tree,fixed=TRUE)
  
  pos <- scan(filename[1],character(0),sep="\n",quiet=TRUE)
  
  #extract distance to pith (if exists)
  pith.dist.line <- grep("^#C.+DistanceToPith=",pos)
  
  #from comments, extract the pith comment
  pith.comment <- NA
  pith.comment.line <- grep("^#C.+pith",pos)
  if((length(pith.comment.line) == 0) ) {
    
    if(length(pith.dist.line) == 0) {
      problem.string <- paste("No pith comment and no pith measurement for core",name.core)
      problem.row <- data.frame(sample=name.core,problem=problem.string)
      pith.problems <<- rbind(pith.problems,problem.row)
    }
    
    # no pith comment on its own is OK as long as there is a measurement (means pith present but forgot to measure)
    

  } else {
    pith.comment.pre <- pos[pith.comment.line]
    pith.comment <- substr(pith.comment.pre,4,100000)
  
    pith.comments.options <- c("pith offset","hit pith","no pith","pith hit")
    if(!(pith.comment %in% pith.comments.options)) {
      problem.string <- paste("Unrecognized pith comment for core",name.core,"; comment is:",pith.comment)
      problem.row <- data.frame(sample=name.core,problem=problem.string)
      pith.problems <<- rbind(pith.problems,problem.row)
    }
  }
  
  #if measured distance to pith, extract distance
  if(length(pith.dist.line) == 0) {
    if(pith.comment %in% c("pith offset","hit pith","pith hit")) {
      problem.string <- paste("Comment specifies pith present, but DTP not measured for core:",name.core)
      problem.row <- data.frame(sample=name.core,problem=problem.string)
      pith.problems <<- rbind(pith.problems,problem.row)
    }
  }
  
  pith.dist <- NA
  radius <- NA
  if(length(pith.dist.line) == 1) {
    pith.dist.line.text <- substr(pos[pith.dist.line],start=4,stop=10000)
    
    pith.dist.parts <- strsplit(pith.dist.line.text,"; ",fixed=TRUE)[[1]]
    
    pith.dist.index <- grep("^DistanceToPith",pith.dist.parts)
    #radius.index <- grep("^Radius",pith.dist.parts)
    
    pith.dist <- as.numeric(strsplit(pith.dist.parts[pith.dist.index],"=")[[1]][2])
    #radius <- as.numeric(strsplit(pith.dist.parts[radius.index],"=")[[1]][2])
  }
  

  # for now, remove the earlywood-latewood boundary points
  annual.rings <- grep("^[0-9]",pos)
  coords <- pos[annual.rings]
  
  
  #extract youngest ring date
  dated.text <- grep("#C DATED [0-9]{4}",pos,value=TRUE)
  dated.year <- as.numeric(substr(dated.text,10,13))
  if(length(dated.year)==0){
    cat("No date assigned to core:",name.core,"\n",sep=" ")
    return(NULL)
  }
  
  if(dated.year > 2500){
    problem.string <- paste("Core still dated to 3000:",name.core,sep=" ")
    problem.row <- data.frame(sample=name.core,problem=problem.string)
    dated.problems <<- rbind(dated.problems,problem.row)
  }
  
  #get coords and convert to ring widths
  
  first.coord <- 1 # first line with coordinates
  last.coord <- length(coords) #assume last line of file contains last coord
  
  d <- rbind.fill(lapply(coords,split.line))
  
  #if(!is.na(d[1,]$comment)) {
  # cat("Comment in first ring of core: ",sample,"\n")
  # return(NULL)
  #}
  
  
  # wherever the second coord is NA (i.e. there was not a second set of coords for that ring), replace with the values of the first coord
  d$xcoord2 <- ifelse(is.na(d$xcoord2),d$xcoord1,d$xcoord2)
  d$ycoord2 <- ifelse(is.na(d$ycoord2),d$ycoord1,d$ycoord2)
  
  # convert coords to widths
  xdist.ends <- c(d$xcoord1,0) - c(0,d$xcoord2) # offset and subtract
  ydist.ends <- c(d$ycoord1,0) - c(0,d$ycoord2)
  diagdist.ends <- sqrt(xdist.ends^2 + ydist.ends^2)    
  widths <- diagdist.ends[2:(length(diagdist.ends)-1)] #youngest ring is first

  nrings <- length(widths)

  years <- (dated.year:(dated.year-(nrings-1)))
  
  core <- data.frame(year=years,width=widths,comment=d$comment[1:nrings],stringsAsFactors=FALSE)
  
  
  ### calculate radius as sum of all widths plus DTP
  if(!is.na(pith.dist)) {
    radius <- sum(widths,na.rm=TRUE) + pith.dist
  } else {
    radius <- NA
  }
  
  
  # remove widths from rings with cracks, if specified
  if(cr.del) {
    #core[5,]$comment <- "cr" #for testing
    crack.rows <- grep("^cr",core$comment)
    core[crack.rows,]$width <- rep(NA,length(crack.rows))
    
  }
  
  #! could add code here for skipping ab, if desired (add the specified number of ab rings with NA width immediately following ab comment)
  
    
  #estimate tree age if pith distance measured
  #! user must understand that it is probably wrong if core had uncs or abs
  
  age <- NA
  if((!is.na(pith.dist))) {
    
    rings.with.width <- which(!is.na(core$width))
    nrings.with.width <- length(rings.with.width)
    last.3.rings.with.width <- rings.with.width[(nrings.with.width-2):nrings.with.width]
    last.rings <- core$width[last.3.rings.with.width]
    
    last.width.avg <- mean(last.rings)
    nrings.pith <- floor(pith.dist / last.width.avg)
    age <- nrings + nrings.pith
  }
  
  
  # remove uncs and abs, if specified
  truncated <- "no" #for recording whether the core was truncated due to comments
  first.unc <- NA
  first.ab <- NA
  
  unc.lines <- grep("^unc",core$comment)
  if(length(unc.lines)>0) { # if there was at least one match
    
    first.unc <- unc.lines[1]
    
    if(unc.stop) {
      core <- core[1:(first.unc-1),]
      truncated <- "unc"
      age <- NA
      problem.string <- paste("Truncated core ",name.core," at ring ",first.unc-1," due to unc.")
      problem.row <- data.frame(sample=name.core,problem=problem.string)
      unc.problems <<- rbind(unc.problems,problem.row)
      
      if(nrow(d) < 4) {
        cat("Truncated at unc in first 4 rings of core: ",name.core," -- skipping core.\n")
        return(NULL)
      }
    }
  }
  
  
  ab.rows <- grep("^ab",core$comment)
  if(length(ab.rows) > 0) {
    first.ab <- min(ab.rows)
    
    if(ab.stop) {
      core <- core[1:(first.ab-1),]
      truncated <- "ab"
      age <- NA
      cat("Truncated core ",name.core," at ring ",first.ab-1," due to ab.\n")
      
      if(nrow(d) < 4) {
        cat("Truncated at ab in first 4 rings of core: ",name.core," -- skipping core.\n")
        return(NULL)
      }
    }
  }
  
  name.tree <- toupper(name.tree)
  
  
  series <- data.frame(year=core$year,sample=core$width)
  names(series) <- c("year",sample)
  
  nrings <- nrow(series)
  
  names(truncated) <- sample
  names(nrings) <- sample
  names(age) <- sample
  names(radius) <- sample
  names(name.tree) <- sample
  
  return(list(series=series,truncated=truncated,nrings=nrings,age=age,radius=radius,core=name.tree))
}

merge.all <- function(x,y) {
  merge(x,y,all=TRUE,by="year")
}

open.chron <- function(samples,samples.secondary,unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE) {
  
  seriess.trunc <- list(rep(NA,length(samples)))
  
  for(i in 1:length(samples)) {
    sample <- samples[i]
    sample.secondary <- samples.secondary[i]
    results <- read.pos(sample,sample.secondary,unc.stop=unc.stop,cr.del=cr.del,ab.stop=ab.stop)
    seriess.trunc[[i]] <- results
  }
  
  seriess <- lapply(seriess.trunc,function(x) x[["series"]])
  trunc <- lapply(seriess.trunc,function(x) x[["truncated"]])
  nrings <- lapply(seriess.trunc,function(x) x[["nrings"]])
  age <- lapply(seriess.trunc,function(x) x[["age"]])
  radius <- lapply(seriess.trunc,function(x) x[["radius"]])
  core <- lapply(seriess.trunc,function(x) x[["core"]])
  
  seriess <- seriess[!sapply(seriess,is.null)]
  chron <- Reduce(merge.all,seriess)
  yearss <- chron$year
  chron.noyr <- chron[,colnames(chron) != "year"]
  chron <- as.data.frame(chron.noyr) # this change might have created a bug?
  rownames(chron) <- yearss  
  
  trunc.all <- unlist(trunc)
  trunc.all <- trunc.all[!is.null(trunc.all)]
  
  nrings.all <- unlist(nrings)
  nrings.all <- nrings.all[!is.null(nrings.all)]
  
  age.all <- unlist(age)
  age.all <- age.all[!is.null(age.all)]
  
  radius.all <- unlist(radius)
  radius.all <- radius.all[!is.null(radius.all)]
  
  core.all <- unlist(core)
  core.all <- core.all[!is.null(core.all)]
  
  return(list(chron=chron,trunc=trunc.all,nrings=nrings.all,age=age.all,radius=radius.all,core=core.all))
}


clean.chron <- function(chron) {
 
  widths <- chron$chron
  truncs <- chron$trunc
  ages <- chron$age
  radii <- chron$radius
  cores <- chron$core
  
  sample.ids <- names(widths)
  
  #### for each sample, apply data from crossdating progress records
  cd.records <- read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv")
  
  ####!!!! sort by record date
  
  
  #get most recent record for each core
  cd.mostrecent <- cd.records[!duplicated(cd.records$Core.number,fromLast=TRUE),]
  
  #!assume that if there was no response, there was no problem
  cd.mostrecent$Was.alignment.limited.by.reference.chronology.length.[cd.mostrecent$Was.alignment.limited.by.reference.chronology.length. == ""] <- "No"
  
  cd.mostrecent$Core.number <- toupper(cd.mostrecent$Core.number)
  
  for(i in 1:length(cores)) {
    
    core <- cores[i]
    sample <- names(widths)[i]
    
    core <- toupper(core)
    sample <- toupper(sample)
    
    cd.mostrecent.samp <- cd.mostrecent[cd.mostrecent$Core.number == core,]
    
    if(nrow(cd.mostrecent.samp) == 0) {
      cat("No crosdating records found for core",core,"; dropping it.\n")
      widths[,sample] <- NA
      ages[sample] <- NA
      radii[sample] <- NA
      
      next()
    }
    
    
    widths.samp <- widths[,sample] #may not need this if modifying whole DF directly
    years <- as.numeric(rownames(widths))
    
    # get most recent year of core
    widths.value <- which(!is.na(widths.samp)) #which rows of data frame have a value for width
    years.value <- years[widths.value]
    dated.samp <- max(years.value)
    
    if(dated.samp > 2020) {
      cat("Series not dated for core",core,"; dropping it.")
      widths[,sample] <- NA
      ages[sample] <- NA
      radii[sample] <- NA
      
      next()
    }
    
    
    if(cd.mostrecent.samp$Was.alignment.limited.by.reference.chronology.length == "Yes") {
      if(is.na(cd.mostrecent.samp$Good.alignment.through.year)) {
        cat("Alignment limited by ref chron but last year of alignment not specified for core: ",core,"; using most recent 40 years.\n")
          
        # clip to most recent 40 years (set widths to NA after 40 years)
        yr.40 <- dated.samp - 39
        widths[(years<yr.40),sample]
      
      }
      
      ####!!! otherwise, clip to good year
      
    }
    
    if(cd.mostrecent.samp$Core.status == "Alignment good through year X") {
      if(is.na(cd.mostrecent.samp$Good.alignment.through.year)) {
        cat("Alignment good through year X, but year X not specified for core: ",core,"; using most recent 40 years.\n")

        # clip to most recent 40 years (set widths to NA after 40 years)
        yr.40 <- dated.samp - 39
        widths[(years<yr.40),sample] <- NA
      
      } else {
        
        #drop ring widths from poor-alignment section
        last.good.alignment <- cd.mostrecent.samp$Good.alignment.through.year
        widths[(years<last.good.alignment),sample] <- NA
        
      }
    }
    
    if(cd.mostrecent.samp$Core.status %in% c("Poor alignment  throughout","Poor alignment throughout")) {
      cat("Core",core,"has poor alignment throughout; dropping ring width data, and age estimate may be inaccurate.\n")
      widths[,sample] <- NA
      
      next()
    }
    
    
    status.options <- c("Poor alignment throughout","Poor alignment  throughout","Good alignment throughout","Good alignment through year X")
    if(!(cd.mostrecent.samp$Core.status %in% status.options)) {
      cat("Core alignment status unrecognized for core",core,"; proceeding anyway.\n")      
    }

    
    # look up and truncate at most recent unclear modification date if one exists. This is the most recent date at which, in order to get good alignment, the analyst had to mark a ring that did not make sense in the image
    
    ##!!@@ need to do this using crossdating progress records too (?)
    
    unc.records <- read.csv("data/dendro/crossdating-records/Measurement checking_cleaned.csv")
    unc.records$ID <- toupper(unc.records$ID)
    unc.year <- unc.records[which(unc.records$ID == core),]$Most.recent.year.of.unclear.modification
    if(length(unc.year) > 1) unc.year <- max(unc.year) # if there is more then one answer due to some fluke, use the most recent year
    if(length(unc.year) != 0) {
      if(!is.na(unc.year)) {
        cat("Truncating at most recent year of modification unsupported by image (",unc.year,") for core:",sample,"\n")
        widths[(years<=unc.year),sample] <- NA
      }
    }
    
    #make a flag for truncated? probably not necessary
    
    
  }
  
  
  # any cores with widths that are all NAs: remove them from data returned.
  # also remove cores with less than 3 ring widths
  nyears <- colSums(!is.na(widths))
  widths.ret <- widths[,nyears>3]
  nyears <- nyears[nyears > 3]
  ret <- list(chron=widths.ret,nyears=nyears,trunc=truncs,age=ages,radius=radii)
  
  return(ret)
  
}
  
interior.na <- function(vect) {
  
  nas <- which(is.na(vect))
  non.nas <- which(!is.na(vect))
  interior.nas <- rep(FALSE,length(vect))  
  
  if(length(non.nas) > 0) {
    first.num <- min(non.nas)
    last.num <- max(non.nas)
    interior.nas.index <- nas[(nas > first.num) & (nas < last.num)]
    interior.nas[interior.nas.index] <- TRUE
  }
  
  names(interior.nas) <- names(vect)
  
  return(interior.nas)
  
}  

#replace interior na values (indicated by the value -987) with the average of the nearest two surrounding non-missing values
interior.na.replace <- function(vect) {
  
  newvect <- vect
  interior.na.indeces <- which(vect == -987)
  
  for(interior.na.index in interior.na.indeces) {
    
    #find the nearest neighbor elements that are real numbers
    
    prev.index <- interior.na.index - 1
    while(vect[prev.index] == -987) {
      prev.index <- prev.index-1
    }
    
    next.index <- interior.na.index + 1
    while(vect[next.index] == -987) {
      next.index <- next.index +1
    }
    
    # average them and set the focal element to that value
    
    prev.width <- vect[prev.index]
    next.width <- vect[next.index]
    dummy.width <- mean(c(prev.width,next.width))
    newvect[interior.na.index] <- dummy.width
    
    
  }
  
  return(newvect)
  
}



spline.na.rm <- function(chron) {
  
  #identify all interior nas
  interior.nas <- apply(chron, 2, interior.na)
  chron[interior.nas] <- -987
  
  #fill them with the average of the two adjacent nas
  chron.nas.replaced <- apply(chron,2,interior.na.replace)
  chron.nas.replaced <- as.data.frame(chron.nas.replaced)
  
  #calculate rwi
  chron.rwi <- detrend(chron.nas.replaced,make.plot=TRUE,method="Spline",nyrs=30)

  #change the rings that were previously NA back to NA
  chron.rwi[interior.nas] <- NA
  
  #return the updated detrended spline
  return(chron.rwi)
  
}



ba.bai.calc <- function(widths,radius) {
  
  bas <- widths
  bai <- widths
  bai.ba <- widths
  ba.prev <- widths
  bas[] <- NA
  bai[] <- NA
  bai.ba[] <- NA
  ba.prev[] <- NA
  
  sample.ids <- names(widths)
  
  for(sample in sample.ids) {
    
    widths.samp <- widths[,sample]
    radius.samp <- radius[sample]
    names(widths.samp) <- rownames(widths)
    
    if(sum(!is.na(widths.samp))==0) { # if all the ring widths are NA
      
      next()
      
    } else {
    
      if(!is.na(radius.samp)) {
        
        if(radius.samp < sum(widths.samp,na.rm=TRUE)) { #if due to some fluke of measurement the radius is less than the sum of the widths, set radius equal to sum of widths
          radius.samp <- sum(widths.samp,na.rm=TRUE)
        }
        
        interior.nas <- interior.na(widths.samp)
        widths.samp[interior.nas] <- -987
        widths.samp <- interior.na.replace(widths.samp)
        widths.samp <- rev(widths.samp)
        
        exterior.nas <- is.na(widths.samp) # because all the interior ones were already replaced
        widths.samp[exterior.nas] <- 0
        radius.subtract <- c(0,cumsum(widths.samp))
        radii <- radius.samp - radius.subtract
        radii <- radii[-length(radii)]
        names(radii) <- names(widths.samp) #this is the radius at the end of that year
        bas.samp <- (radii^2)*3.141593
        
        bai.samp <- c(0,bas.samp) - c(bas.samp,0) 
        bai.samp <- bai.samp[-1]
        
        #set last value of bai to NA (it is meaningless)
        bai.samp[length(bai.samp)] <- NA
        
        # get bai/ba
        bai.ba.samp <- c(0,bai.samp) / c(bas.samp,1)
        bai.ba.samp <- bai.ba.samp[-1]
        
        ba.prev.samp <- c(bas.samp[-1],NA)
        
        bas.samp[exterior.nas] <- NA
        bai.samp[exterior.nas] <- NA
        bai.ba.samp[exterior.nas] <- NA
        ba.prev.samp[exterior.nas] <- NA
        
        bas.samp[interior.nas] <- NA
        bai.samp[interior.nas] <- NA
        bai.ba.samp[interior.nas] <- NA
        ba.prev.samp[interior.nas] <- NA
        
        bas[,sample] <- rev(bas.samp)
        bai[,sample] <- rev(bai.samp)
        bai.ba[,sample] <- rev(bai.ba.samp)
        ba.prev[,sample] <- rev(ba.prev.samp)
      
      } else { #we did not have radius data for this tree
        
        next()
        
      }
    }

  }
  
  return(list(ba=bas,bai=bai,bai.ba=bai.ba,ba.prev=ba.prev))
  
}




# read.chron.trees <- function(file) {
#   d <- read.csv(file,header=TRUE)
#   sites <- unique(d$Site)
#   sitetrees <- NULL
#   for(site in sites) {
#     sitetrees[[site]] <- as.character(d[d$Site==site,]$Tree)
#   }
#   return(sitetrees)
# }
# 
# 
# open.all.chrons <- function(file) {
#   chrons.trees <- read.chron.trees(file)
#   chrons.widths <- lapply(chrons.trees,open.chron)
#   return(chrons.widths)
# }
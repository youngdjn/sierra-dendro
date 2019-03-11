
## Normalize a series of numbers by getting each number's proportional difference from the average of the two surrounding it
norm = function(x) {
  
  x_lag = lag(x)
  x_lead = lead(x)
  neighbor_mean = ((x_lag+x_lead)/2)
  
  x_prop = (x-neighbor_mean)/neighbor_mean
  
  return(x_prop)
  
}


## Create n new columns with leads of different lengths (1:n) from the specified variable
leads <- function(var, n=10){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(lead(!!var, !!.x)) ) %>% 
    set_names(sprintf("lead_%s_%02d", rlang::quo_text(var), indices))
  
}




get_tree_name = function(core_names) {
  
  tree_names = str_replace(core_names,pattern="[aAbBzZ]",replacement="")
  
  new_names = NULL
  for(i in 1:length(tree_names)) {
    tree_name = tree_names[i]
    new_names[i] = ifelse(tree_name %in% trees$former.id,trees[trees$former.id==tree_name,]$tree.id,tree_name)
  }
  
  return(new_names)
}



pith.problems <- data.frame(sample=NA,problem=NA)
pos.problems <- data.frame(sample=NA,problem=NA)
unc.problems <- data.frame(sample=NA,problem=NA)
dated.problems <- data.frame(sample=NA,problem=NA)



# function to take a line in a .pos file and return the x and y coordinates
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


# funciton to take a core number ("sample") and read the .pos file and output the ring-width series and some related metrics like tree age.
read.pos <- function(sample,sample.secondary, unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE,core.name=FALSE) {

  
  #cat("Reading pos for cores",sample,sample.secondary,"\n")
  
  ## for testing
  #sample <- "2146"
  
  
  ### !!!!  need to allow for multiple cores per tree when looking up filename
  #             how to find best one? or open both? store inferior one in other folder?
  
  ### !!! make sure working properly with point offsets
  
  ### !!!! need to add core ID for first-ring comment check

  if(!core.name) { #if core.name is true, then look for an exact match
    pattern <- paste(sample,"A?B?Z?t?.pos$",sep="")
  } else {
    pattern <- paste(sample,"t?.pos$",sep="")
  }
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
  name.core <- gsub("t","",name.core,fixed=TRUE)
  name.core <- gsub("T","",name.core,fixed=TRUE)
  name.core <- toupper(name.core)
  
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
  
  #see if there is a comment saying the core is rotated
  rotated.text <- grep("#C CORE ROTATED",pos,value=TRUE,ignore.case=TRUE)
  if(length(rotated.text) > 0){
    cat("Core rotated:",name.core,", dropping radius and age\n",sep=" ")
    rotated=TRUE
  } else {
    rotated=FALSE
  }
  

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
  
  
  
  
  series <- data.frame(year=core$year,sample=core$width)
  names(series) <- c("year",sample)
  
  nrings <- nrow(series)
  
  names(truncated) <- sample
  names(nrings) <- sample
  names(age) <- sample
  names(radius) <- sample
  
  if(rotated==TRUE) {
    age = NA
    radius = NA
  }
  
  return(list(series=series,truncated=truncated,nrings=nrings,age=age,radius=radius,core=name.core,tree.id=sample))
}

merge.all <- function(x,y) {
  merge(x,y,all=TRUE,by="year")
}

# funciton to open multiple .pos files and compile them into a chronology
open.chron <- function(samples,samples.secondary=NA,unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE, core.name=FALSE) {
  
  seriess.trunc <- list(rep(NA,length(samples)))
  
  for(i in 1:length(samples)) {
    sample <- samples[i]
    sample.secondary <- samples.secondary[i]
    results <- read.pos(sample,sample.secondary,unc.stop=unc.stop,cr.del=cr.del,ab.stop=ab.stop,core.name=core.name)
    seriess.trunc[[i]] <- results
  }
  
  seriess <- lapply(seriess.trunc,function(x) x[["series"]])
  trunc <- lapply(seriess.trunc,function(x) x[["truncated"]])
  nrings <- lapply(seriess.trunc,function(x) x[["nrings"]])
  age <- lapply(seriess.trunc,function(x) x[["age"]])
  radius <- lapply(seriess.trunc,function(x) x[["radius"]])
  core <- lapply(seriess.trunc,function(x) x[["core"]])
  tree.id <- lapply(seriess.trunc,function(x) x[["tree.id"]])
  
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
  
  tree.id.all <- unlist(tree.id)
  tree.id.all <- tree.id.all[!is.null(tree.id.all)]
  
  return(list(chron=chron,trunc=trunc.all,nrings=nrings.all,age=age.all,radius=radius.all,core=core.all,tree.id=tree.id.all))
}


clean.chron <- function(chron,trunc.unjustified.by.image = TRUE,clean_using_alignment = FALSE, drop_poor_align_throughout = TRUE, drop_mult_potential_ways = TRUE) {
  
  widths <- chron$chron
  truncs <- chron$trunc
  ages <- chron$age
  radii <- chron$radius
  cores <- chron$core
  trees.id <- chron$tree.id
  
  sample.ids <- names(widths)
  
  #### for each sample, apply data from crossdating progress records
  #### NOTE that this must already be sorted by modification date
  cd.records <- read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv")
  cd.records$Core.number <- toupper(cd.records$Core.number)
  
  ## aggregate it to determine for each core the most recent date of modifications that did not make sense
  mod.justification <- cd.records %>%
    dplyr::group_by(Core.number) %>%
    dplyr::summarize(most.recent.unjustified.mod = max(`If.modification.did.not.make.sense.in.image...`,na.rm=TRUE)) %>%
    filter(most.recent.unjustified.mod > 0)
  

  #get most recent record for each core
  cd.mostrecent <- cd.records[!duplicated(cd.records$Core.number,fromLast=TRUE),]
  
  #!assume that if there was no response, there was no problem
  cd.mostrecent$Was.alignment.limited.by.reference.chronology.length.[cd.mostrecent$Was.alignment.limited.by.reference.chronology.length. == ""] <- "No"
  
  cd.mostrecent$Core.number <- toupper(cd.mostrecent$Core.number)
  
  for(i in 1:length(cores)) {
    
    core <- cores[i]
    sample <- names(widths)[i]
    trunc = truncs[i]
    
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
    
    
    if(clean_using_alignment == TRUE) {
    
      
      if(cd.mostrecent.samp$Was.alignment.limited.by.reference.chronology.length == "Yes") {
            if(is.na(cd.mostrecent.samp$Good.alignment.through.year)) {
              cat("Alignment limited by ref chron but last year of alignment not specified for core: ",core,"; using most recent 40 years.\n")
                
              # clip to most recent 40 years (set widths to NA after 40 years)
              yr.40 <- dated.samp - 39
              widths[(years<yr.40),sample]
            
            }
            
      
            
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
    }
    
    if(drop_poor_align_throughout == TRUE) {
      if(cd.mostrecent.samp$Core.status %in% c("Poor alignment  throughout","Poor alignment throughout")) {
        cat("Core",core,"has poor alignment throughout; dropping ring width data, and age estimate may be inaccurate.\n")
        widths[,sample] <- NA
        truncs[sample] = "poor_align_throughout"
        next()
      }
    }
    
    if(drop_mult_potential_ways == TRUE) {
      
      if(cd.mostrecent.samp$Were.there.multiple.potential.ways.to.get.good.alignment. == "Multiple potential ways") {
        cat("Core",core,"has multiple potential alignments; dropping ring width data, and age estimate may be inaccurate.\n")
        widths[,sample] <- NA
        truncs[sample] = "multiple_alignments"
      
        next()
      }
      
    }
    

    if(trunc.unjustified.by.image == TRUE) {
    
      # look up and truncate at most recent unclear modification date if one exists. This is the most recent date at which, in order to get good alignment, the analyst had to mark a ring that did not make sense in the image
      
      unjustified.year <- mod.justification[mod.justification$Core.number == core,]$most.recent.unjustified.mod
      
      unc.year <- unjustified.year
      
      
      if(length(unc.year) > 1) unc.year <- max(unc.year,na.rm=TRUE) # if there is more then one answer due to some fluke, use the most recent year
      if(length(unc.year) != 0) {
        if(!is.na(unc.year)) {
          
          cat("Truncating at most recent year of modification unsupported by image (",unc.year,") for core:",core,"\n")
          widths[(years<=unc.year),sample] <- NA
          truncs[sample] = "unjustified_mod"
        }
      }
      

    }
    
  }
  
  
  # # any cores with widths that are all NAs: remove them from data returned.
  # # also remove cores with less than 3 ring widths
  nyears <- colSums(!is.na(widths))
  widths.ret <- widths[,nyears>3]
  nyears <- nyears[nyears > 3]
   truncs = truncs[nyears > 3]
   ages = ages[nyears > 3]
   radii = radii[nyears > 3]
  
  
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


# function to take a series of ring withds and a tree radius and calculate basal area and basal area increment for each year
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


list.plots <- function() {
  
  trees <- read.csv("data/plot-and-tree/processed/trees_loc.csv",header=TRUE,stringsAsFactors=FALSE)
  #cat(unique(trees$plot.id))
  return(unique(trees$plot.id))
  
}



summarize.cluster <- function(cluster,type="cluster",name=cluster,clean.ref.chron=TRUE,drop_mult_potential_ways=TRUE, drop_poor_align_throughout=TRUE, trunc_unjustified =TRUE,...) {
  
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
  
  files = list.files(path="data/dendro/coorecorder-measurements",pattern=".pos$",full.names=FALSE,ignore.case=TRUE)
  
  ## take the extension and "t" off of all of them
  core.names <- toupper(sapply(files,function(x) strsplit(x,"[tT]*.pos$")[[1]][1]))
  
  ## find which core names (from the pos files) match trees
  cores.match.trees = grepl(pattern,core.names)
  cores.cluster <- core.names[cores.match.trees]
  #cores.cluster <- core.names
  
  open.chron.results <- open.chron(cores.cluster,cores.cluster,unc.stop=TRUE,core.name=TRUE)
  
  if(clean.ref.chron == TRUE) {
    open.chron.results = clean.chron(open.chron.results, drop_mult_potential_ways=drop_mult_potential_ways, drop_poor_align_throughout = drop_poor_align_throughout, trunc.unjustified.by.image = trunc_unjustified)
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
  
  corr <- corr.rwl.seg(chron.dated,seg.length=20,bin.floor=14,prewhiten=FALSE,pcrit=0.1,label.cex=.7,main=paste("All dated cores in cluster",name))
  
  # use only cores that had overall rho > 0.4
  corr.overall <- as.data.frame(corr$overall)
  ref.cores <- rownames(corr.overall[corr.overall$rho > 0.5,])
  
  if(length(ref.cores) < 4) stop("Too few cores in reference chronology. That is, too few cores that correlate well with the mean across all cores in the cluster.")
  
  chron.dated.ref <- chron.dated[,ref.cores]
  corr.ref <- corr.rwl.seg(chron.dated.ref,seg.length=20,bin.floor=14,pcrit=0.10,prewhiten=FALSE,label.cex=.7,main=paste("Cores in reference chronology for",name))
  
  # add a col to the dated cores list for whether it is included in master chron (i.e., rho > 0.5)
  dated.cores.list.2$ref.chron <- ifelse(dated.cores.list.2$sample %in% ref.cores,TRUE,FALSE)
  
  # add a col to the dated cores list giving the length of the chron to truncation point
  trunc.nrings <- merge(trunc,nrings,by="row.names",all=TRUE)
  names(trunc.nrings) <- c("sample","truncated","nrings")
  dated.cores.list.3 <- merge(dated.cores.list.2,trunc.nrings,by="sample",all=TRUE)
  dated.cores.list.4 <- dated.cores.list.3[,c("tree.id","sample","species","ref.chron","dated","truncated","nrings")]
  dated.cores.list.4$sample.clean <- gsub("T","",dated.cores.list.4$sample)
  
  # # load in comments from crossdating progress records
  # 
  # crossdating.records <- read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv")
  # 
  # 
  # if(length(crossdating.records) == 0) {
  #   
  #   cat("!!! Could not connect to Google spreadsheet.\n!!! Previously-processed cores cannot be displayed in cluster summary.\n")
  #   crossdating.records <- data.frame("Core number"=NA,Person=NA,Timestamp=NA)
  #   names(crossdating.records) <- c("Core number","Person","Timestamp")
  # 
  # } else {
  #   crossdating.records <- crossdating.records[rowSums(is.na(crossdating.records)) != ncol(crossdating.records),] # remove empty rows
  # 
  #   # take only the most recent entry for each core
  #   crossdating.records <- crossdating.records[nrow(crossdating.records):1,] # reverse order
  #   crossdating.records.dups <- duplicated(crossdating.records$Core.number)
  #   crossdating.records <- crossdating.records[!crossdating.records.dups,]
  # 
  # }
  # 
  # crossdating.records$Core.number = toupper(crossdating.records$Core.number)
  # 
  # cluster.summary <- merge(dated.cores.list.4,crossdating.records[,c("Core.number","Person","Timestamp","Core.status","Good.alignment.through.year","Other.notes","Derek.double.check.needed.")],by.x="sample.clean",by.y="Core.number",all.x=TRUE)
  # names(cluster.summary) <- c("sample.clean","tree.id","sample","species","ref.chron","dated","truncated","nrings","completed.by","completed.date","previous crossdating results","Alignment good through","previous crossdating notes","previous crossdating notes2")
  
  cluster.summary <- dated.cores.list.4
  names(cluster.summary) <- c("tree.id","sample","species","ref.chron","dated","truncated","nrings","sample.clean")
  
  
  
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
  
  
  cluster_summary = list()
  cluster_summary$cores_in_cluster = dated.cores.list.4$sample
  cluster_summary$ref_chron = chron.dated.ref.mean
  
  return(cluster_summary)
  
  #### Inspecting individual cores ####
  # core to inspect
  #focal.core <- "2006"
  
  #### Correlogram plotting ####
  # !! is it ok that the focal core may be a member of the reference series? yes because the ref series must have at least 5 members
  #corr.series.seg(chron.dated.ref.mean,series=chron[,focal.core],series.yrs=as.numeric(rownames(chron)),seg.length=20,bin.floor=0,prewhiten=FALSE)
  
  #### CCF plotting ####
  #ccf.series.rwl(chron.dated.ref.mean,series=chron[,focal.core],series.yrs=as.numeric(rownames(chron)),seg.length=20,bin.floor=0,prewhiten=FALSE)
  
  #CDendro-type plot
  # series.rwl.plot(chron.dated,series="2006",seg.length=10,bin.floor=0)
  
  
  
}




open_cluster_chron_truncate_window = function(focal_cluster,trunc_unjustified,window_width, cor_threshold, max_bad_years, min_good_years) {
  
  
  ## Summarize the cluster.
  ## In doing so, clean it, specifically: truncate where unjustified by image, where poor alignment throughout, and where multiple potential ways
  
  
  cluster_summary = summarize.cluster(cluster = focal_cluster,type = "cluster",clean.ref.chron=TRUE, clean_using_alignment=FALSE, drop_mult_potential_ways=TRUE, drop_poor_align_throughout=TRUE, trunc_unjustified =TRUE) # clean.ref.chron (which is true by default) truncates the core when corssdating records indicate: a ring removal or addition did not make sense in the image; crossdating was limited by reference chronology length 
  
  
  ## Open all cores, clean, compare against the ref, find where correlation goes bad
  chron = open.chron(cluster_summary$cores_in_cluster,cluster_summary$cores_in_cluster,unc.stop=TRUE,cr.del=TRUE,ab.stop=TRUE, core.name=FALSE)
  
  ## convert the core names to tree ids
  names(chron$chron) = get_tree_name(names(chron$chron))
  names(chron$trunc) = get_tree_name(names(chron$trunc))
  names(chron$nrings) = get_tree_name(names(chron$nrings))
  names(chron$age) = get_tree_name(names(chron$age))
  names(chron$radius) = get_tree_name(names(chron$radius))
  names(chron$core) = get_tree_name(names(chron$core))
  
  chron2 = clean.chron(chron, drop_mult_potential_ways=TRUE, drop_poor_align_throughout = TRUE, trunc.unjustified.by.image = TRUE, clean_using_alignment=FALSE)
  chron3 = spline.na.rm(chron2$chron)
  chron_rwi = chron3
  chron_rwi$year = row.names(chron_rwi)
  chron_raw = chron2$chron
  chron_raw$year = row.names(chron_raw)
  truncs = chron2$trunc
  ages = chron2$age
  radii = chron2$radius
  radii <- data.frame(tree.id=names(radii),radius=radii)
  
  ref = cluster_summary$ref_chron
  ref$year = row.names(ref)
  names(ref) = c("width","samp.depth","year")
  
  for(i in seq_along(names(chron_rwi))) {
    
    core = names(chron_rwi)[i]
    if(core == "year") next()
    
    series = chron_rwi[,c(core,"year")]
    names(series) = c("width","year")
    
    # plot(norm(ref$width)~ref$year,type="l")
    # lines(norm(series$width)~series$year,col="red")
    # 
    # 
    # ###tem;
    # ref2 = ref
    # names(ref2) = c("ref","samp","year")
    # 
    # series2 = series
    # names(series) = c("core","year")
    # 
    # series.ref = left_join(series2,ref2)
    # 
    # row.names(series.ref) = series.ref$year
    # 
    # series.ref = series.ref %>%
    #   select(width,ref)
    #   
    # detrended = spline.na.rm(series.ref)
    # detrended$year= row.names(detrended)
    # 
    # plot(detrended$width~detrended$year,type="l",col="red")
    # lines(detrended$ref~detrended$year,col="black")
    
    ## run a moving window of correlations
    first_year = min(series[!is.na(series$width),]$year) %>% as.numeric
    last_year = max(series[!is.na(series$width),]$year) %>% as.numeric
    first_center = first_year + window_width/2
    last_center = last_year - window_width/2
    window_centers = last_center:first_center
    
    ##make sure core has at least window_width number of years
    if((last_year-first_year) < (window_width+2)) {
      next()
    }
    
    
    cor_windows = data.frame()
    for(center in window_centers) {
      
      window_low = center - window_width/2
      window_high = center + window_width/2
      window = window_low:window_high
      
      series_window = series[series$year %in% window,]$width
      ref_window = ref[ref$year %in% window,]$width
      
      if(length(ref_window) < window_width+1) next() #reached the end of the reference chronology
      
      window_cor = cor(norm(series_window),norm(ref_window),use="complete.obs",method="pearson")
      
      cor_window = data.frame(year=center,cor=window_cor)
      
      cor_windows = rbind(cor_windows,cor_window)  
      
    }
    
    #plot(cor~year,data=cor_windows,type="l",main=core,ylim=c(-1,1))
    
    # encode the cor as either good or bad
    cor_windows = cor_windows %>%
      mutate(quality = ifelse(cor > cor_threshold,"good","bad"))
    
    quality_runs = rle(cor_windows$quality)
    run_lengths = rep(quality_runs$lengths,times=quality_runs$lengths)
    cor_windows$qual_num = ifelse(cor_windows$quality=="good",run_lengths,NA)
    
    cor_windows = cor_windows %>%
      mutate(!!!leads(qual_num,max_bad_years))
    
    lead_cols = cor_windows %>% select(starts_with("lead_qual_num")) %>% names()
    
    cor_windows = cor_windows %>%
      mutate(max_past_qual = pmax(!!!rlang::syms(lead_cols),na.rm=TRUE)) %>%
      mutate(max_past_qual = ifelse(is.na(max_past_qual),0,max_past_qual)) %>%
      dplyr::select(-starts_with("lead_qual_num"))
    
    
    ### find the first year where the quality is bad and the max_past_qual is less then the min_good_length
    bad_rows = cor_windows %>%
      filter(quality == "bad" &
               max_past_qual < min_good_years)    
    most_recent_bad_year = max(bad_rows$year)
    truncate_year = most_recent_bad_year + window_width/2
    
    #remove the ring width data from the bad year and older
    if(truncate_year != -Inf) {
      chron_rwi[chron_rwi$year<truncate_year,core] = NA
      chron_raw[chron_raw$year<truncate_year,core] = NA
      truncs[core] = paste0("alignment_drops, ",truncs[core])
    }
    
    
  }
  
  return(list(chron_rwi=chron_rwi,chron_raw=chron_raw,truncs=truncs,radii=radii,ages=ages))
}





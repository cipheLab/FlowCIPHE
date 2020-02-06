#@ `export
read.FCS.CIPHE <- function(fcs){
  out <- tryCatch({
    read.FCS(fcs,emptyValue=FALSE)
  },
  error = function(cond){
    return(NULL)
  })
}

#@ `export
catch.create.FCS.from.CSV <- function(csv){
  if(dim(csv)[2]<2){return(NULL)}
  out <- tryCatch({
    createFCSfromCSV(csv)
  },
  error = function(cond){
    return(NULL)
  })
}

#@ `export
found.spill.CIPHE <- function(fcs){
  comp <- "NULL"
  dim <- "NULL"
  desc <- fcs@description
  id1 <- unlist(lapply(desc, function(i){return(dim(i)[1])}))
  if(!is.null(id1)){
    comp <- names(id1)
    s <- fcs@description[[names(id1)]]
    if(dim(s)[1]==dim(s)[2]){
      dim <- dim(s)[1]
    }
  }
  return(c(comp,dim))
}

#@ `export

compensate.CIPHE <- function(flow.frame, spill=NULL){
  if(is.null(spill)){
    spill <- foundSpillCIPHE(flow.frame)[1]
    if(is.null(spill)){
      warning("No compensation apply and/or found")
      return(flow.frame)
    }
  }
}


#@ `export
delete.column.FCS.CIPHE <- function(fcs, marker, spill=NULL){
  id <- which(colnames(fcs)==marker)
  data <- exprs(fcs)
  data <- data[,-id]
  exprs(fcs) <- data
  desc <- fcs@description
  new.desc <- list()
  new.name <- list()
  for(i in c(1:length(desc))){
    if(length(grep(paste0("P",id),names(desc[i])))>0){
      print(names(desc)[i])
    } else {
      new.name <- c(new.name,names(desc)[i])
      new.desc <- c(new.desc,desc[i])
    }
  }
  names(new.desc) <- new.name
  fcs@description <- new.desc
  if(!is.null(spill)){
    new.spill <- fcs@description[[spill]]
    id.spill <- which(colnames(new.spill)==marker)
    new.spill <- new.spill[-id.spill,-id.spill]
    fcs@description[[spill]] <- new.spill
  }
  return(fcs)
}

#@ `export
enrich.FCS.CIPHE <- function(original, new.column){
  new_p <- parameters(original)[1,]
  ## Now, let's change it's name from $P1 to $P26 (or whatever the next new number is)
  new_p_number <- as.integer(dim(original)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))
  ## Now, let's combine the original parameter with the new parameter
  library('BiocGenerics') ## for the combine function
  allPars <- combine(parameters(original), new_p)
  ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
  new_p_name <- colnames(new.column)
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_name
  new_exprs <- cbind(original@exprs, new.column)
  new_kw <- original@description
  new_kw["$PAR"] <- as.character(new_p_number)
  new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
  new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
  new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
  new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P1R"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- new_kw["flowCore_$P1Rmin"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- new_kw["flowCore_$P1Rmax"]
  ## Now, let's just combine it into a new flowFrame
  new_fcs <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)
  return(new_fcs)
}

#@ `export
logicle.CIPHE <- function(flow.frame, value = NULL, markers = NULL){

  if(is.null(markers)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- markers
  }

  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)

  if(is.null(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])){
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      )
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }

  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5

  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }

  return(flow.frame)
}

#@ `export
invers.logicle.CIPHE <- function(flow.frame, value = NULL, markers = NULL){
  if(is.null(markers)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- markers
  }

  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)

  if(is.null(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      )
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }

  w.values <- (4.5-log10(262144/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5

  flow.frame.inv <- flow.frame

  for(t in 1:length(markers.transform)){
    invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
    flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
  }

  return(flow.frame.inv)
}

#@ `export
arcsinh.CIPHE <- function(flow.frame, marker=NULL, arg){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  if(is.null(marker) || length(marker)<1){
    marker <- colnames(flow.frame)
  }
  mat <- mat[,marker]
  colnames(mat) <- marker
  mat <- asinh(mat/arg)
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

#@ `export
invers.arcsinh.CIPHE <- function(flow.frame, marker_untrans, arg){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  marker_untrans_index <- which(colnames(flow.frame)%in%marker_untrans)
  mat <- mat[,-marker_untrans_index]
  marker <- colnames(mat)
  mat <- sinh(mat)*arg
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

#@ `export
decompensate.CIPHE <- function(x, spillover) {
  if(!is.null(spillover)){
    cols <- colnames(spillover)
    sel <- cols %in% colnames(x)
    if(!all(sel)) {
      stop(keyword(x)[["FILENAME"]], "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
           paste(cols[!sel], collapse=", "), call.=FALSE)
    }
    e <- exprs(x)
    e[, cols] <- e[, cols] %*% spillover
    exprs(x) = e
    return(x)
  } else {
    return(x)
  }
}

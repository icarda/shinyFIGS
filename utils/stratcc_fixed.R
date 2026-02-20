### From archived ccChooser package -- fixed allocc function
D2 <-
  function(size, newSize, grpIDs, grpSize, groups, data){
    sumGrpDist <- unlist( lapply( 1:length(grpIDs), function(i){
      tmp<- mean( daisy(data[groups==grpIDs[i],]) * grpSize[i] )
      if(is.nan(tmp))
        tmp <- 0
      return(tmp)
    } ))
    
    sumDist <- sum(sumGrpDist)
    newGrpSizes<- unlist (lapply( 1: length(grpIDs), function(i){
      newGrpSize <- round( newSize * (sumGrpDist[i] / sumDist) )
      if(newGrpSize == 0)
        newGrpSize <- 1
      return(newGrpSize)
    }) )
    
    return(newGrpSizes)
  }

D3 <-
  function(size, newSize, grpIDs, grpSize, groups, data){
    sumGrpDist <- unlist( lapply( 1: length(grpIDs), function(i){
      tmp <- mean( daisy( data[groups==grpIDs[i] , ] ) )* log(grpSize[i]) * grpSize[i]
      if(is.nan(tmp))
        tmp <- 0
      return(tmp)
    }))
    
    sumDist <- sum(sumGrpDist)
    
    newGrpSizes <- unlist( lapply( 1: length(grpIDs), function(i){
      newGrpSize <- round(newSize * (sumGrpDist[i] / sumDist))
      if(newGrpSize == 0)
        newGrpSize <- 1
      return(newGrpSize)
    }) )
    
    return(newGrpSizes)
  }

LOG <-
  function(groups, size, newSize, grpIDs, grpSize){
    sumGrpLog<-unlist( lapply( 1:length(grpIDs), function(i, grpSize){
      return(log(grpSize[i]) * grpSize[i])
    }, grpSize) )
    sumLog <- sum(sumGrpLog)
    
    newGrpSizes<-unlist( lapply( 1: length(grpIDs), function(i, newSize, sumGrpLog, sumLog){
      newGrpSize <- round(newSize * ( sumGrpLog[i] / sumLog ))
      if(newGrpSize == 0)
        newGrpSize <- 1
      return(newGrpSize)
    }, newSize, sumGrpLog, sumLog) )
    
    return(newGrpSizes)
  }

PRO <-
  function(groups, size, newSize, grpIDs){
    newGrpSizes<-unlist(lapply(1:length(grpIDs), function(i, groups){
      newGrpSize <- round(newSize*(sum(groups==i)/size))
      if(newGrpSize == 0)
        newGrpSize <- 1
      return(newGrpSize)
    }, groups))
    
    return(newGrpSizes)
  }

allocc <- function (x, groups, fraction = 0.1, method = "Pro") 
{
  if (!any(method %in% c("Pro", "Log", "D2", "D3"))) 
    stop("Bad name for allocation method!")
  size <- nrow(x)
  if (size != length(groups)) 
    stop("Number of groups must be equal number of rows in x")
  newSize <- size * fraction
  grpIDs <- sort(unique(groups))
  grpSize <- unlist(lapply(grpIDs, function(i, groups) {
    sum(groups == i)
  }, groups))
  if (method == "Pro") {
    cat("Proportional allocation method\n")
    newGrpSizes <- PRO(groups, size, newSize, grpIDs)
  }
  else if (method == "Log") {
    cat("Logarytmic allocation method\n")
    newGrpSizes <- LOG(groups, size, newSize, grpIDs, grpSize)
  }
  else if (method == "D2") {
    cat("D2 allocation method\n")
    newGrpSizes <- D2(size, newSize, grpIDs, grpSize, groups, 
                      x)
  }
  else if (method == "D3") {
    cat("D3 allocation method\n")
    newGrpSizes <- D3(size, newSize, grpIDs, grpSize, groups, 
                      x)
  }
  results <- cbind(grpIDs, newGrpSizes, deparse.level = 0)
  cat(paste("Number of accessions in core colection:", sum(as.numeric(results[, 
                                                                              2])), "\n", sep = " "))
  results <- results[order(results[, 1]), , drop = FALSE]
  colnames(results) <- c("GroupID", "NewSize")
  invisible(results)
}

stratcc <- function (x, groups, alloc = "Pro", fraction = 0.1, clustering = FALSE, 
                     cluster_method = "ward") 
{
  allocated <- allocc(x, groups, fraction, method = alloc)
  grpIDs <- unique(groups)
  if (clustering) {
    CC <- lapply(grpIDs, function(grpID) {
      grpData <- subset(x, groups == grpID)
      if (nrow(grpData) == 1) {
        return(grpData)
      }
      else {
        tree <- stats::hclust(cluster::daisy(grpData, stand = TRUE), 
                              method = cluster_method)
        
        clGrpID <- stats::cutree(tree, k = allocated[grpID, 2])
        tmp <- lapply(1:allocated[allocated[, 1] == grpID, 
                                  2], function(alloID) {
                                    alloData <- grpData[clGrpID == alloID, ]
                                    resID <- sample(1:nrow(alloData), 1)
                                    return(alloData[resID, ])
                                  })
      }
      tmp <- as.data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
      return(tmp)
    })
  }
  else {
    CC <- lapply(grpIDs, function(grpID) {
      grpData <- subset(x, groups == grpID)
      allocID <- sample(1:nrow(grpData), allocated[allocated[, 
                                                             1] == grpID, 2])
      newData <- grpData[allocID, ]
      return(newData)
    })
  }
  CC <- as.data.frame(do.call(rbind, CC), stringsAsFactors = FALSE)
  return(CC)
}
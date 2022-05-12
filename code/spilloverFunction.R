

calcBetaMuVar <- function(alpha, beta) {
  mu <- alpha / (alpha+beta)
  var <- alpha * beta / ((alpha+beta)^2 * (alpha+beta+1))
  return(list(mu=mu, var=var))
}


estimateSpilloverZNB <- function(spilloverPerArea, population, startYear=NULL, endYear=NULL, minYear=NULL, 
                                 zeroYr=T,use_adm1=T,useRange=T, numReplicates=1e3){
  require(glmmTMB)
  
  # check input parametes
  stopifnot(c("spillovers", "YEAR","GID_1") %in% colnames(spilloverPerArea))
  
  rownames(population) <- population$GID_1
  pop2015 <- population[spilloverPerArea$GID_1,"pop2015"]
  spillovers <- spilloverPerArea

  # fit negative binomial parameters describing spillover
  spT=transform(spilloverPerArea,yr=scale(YEAR))
  if(zeroYr==T){
    #Include year term for zero-inflation process
    
    if(use_adm1==T){
      nb1fit=glmmTMB(spillovers~1+(1|yr)+(1|GID_0:GID_1)+(1|GID_0)+offset(log(pop2015)),
                     data=spT,
                     ziformula = ~1|yr,
                     family=nbinom2)        
    }else{
      nb0fit=glmmTMB(spillovers~1+(1|yr)+(1|GID_0)+offset(log(pop2015)),
                     data=spT,
                     ziformula = ~1|yr,
                     family=nbinom2)        
    }
  }else{
    nb0fit=glmmTMB(spillovers~1+(1|yr)+(1|GID_0:GID_1)+(1|GID_0)+offset(log(pop2015)),
                   data=spT,
                   ziformula = ~1,
                   family=nbinom2)  
  }
  
  print(AIC(nb0fit))

  spT$est.mu=predict(nb0fit,newdata=spT,type="conditional")
  spT$est.theta=predict(nb0fit,newdata=spT,type="disp") 
  spT$est.pr0=predict(nb0fit,newdata=spT,type="zprob")
  
  return(spT)
}


# Generate spillovers per area
estimateSpillover <- function(spilloverPerArea, population, alpha, startYear=NULL, endYear=NULL, minYear=NULL, useRange=T, numReplicates=1e3){
  require(mgcv)
  
  # check input parametes
  stopifnot(c("spillovers", "YEAR","GID_1") %in% colnames(spilloverPerArea))
  stopifnot(c("pop2015") %in% colnames(population))
  stopifnot(length(alpha) == 2)
  rownames(population) <- population$GID_1

  # get spillover variation symptomatic rate
  scale <- sum(alpha)
  
  # get number of year
  if(is.null(startYear))
    startYear <- unlist(lapply(split(spilloverPerArea$YEAR, spilloverPerArea$GID_1), min))
  if(is.null(endYear))
    endYear <- unlist(lapply(split(spilloverPerArea$YEAR, spilloverPerArea$GID_1), max))
  population$startYear <- startYear[population$GID_1]
  population$endYear <- endYear[population$GID_1]
  if(useRange){
    population$spillover.nYears <- population$endYear - population$startYear
    years <- min(startYear):max(endYear)
  }else{
    numYears <- unlist(lapply(split(spilloverPerArea$YEAR, spilloverPerArea$GID_1), length))
    population$spillover.nYears <- numYears[population$GID_1]
    years <- sort(unique(spilloverPerArea$YEAR))
  }
  # fit negative binomial parameters describing spillover
  params.spillover = lapply(population$GID_1, function(adm1){
    message(adm1)
    obsSpillovers <- spilloverPerArea[spilloverPerArea$GID_1==adm1,c("YEAR","spillovers")]
    rownames(obsSpillovers) <- obsSpillovers$YEAR
    estSpillover <- do.call(rbind, lapply(years, function(yy){
      sp <- obsSpillovers[as.character(yy),"spillovers"]
      sp <- ifelse(is.na(sp),0,sp)
      max.spillovers = max(round(c(scale, scale*sp),0))
      spillovers = seq(from=max(0, sp), 
                       to=min(max.spillovers, ceiling(population[adm1,"pop2015"])),
                       by=1)
      # estimate probability to have a spillover of x cases
      Pr.spillovers = ddirmnom(x=cbind(sp, spillovers - sp),
                               size=spillovers,
                               alpha=alpha)
      spillovers <- sample(spillovers, size=numReplicates, prob=Pr.spillovers, replace=T)
      data.frame(spillover=spillovers, years=yy)
    }))
    if(length(years)<10){
      gamCases <- gam(spillover ~ s(years, k=length(years)), data=estSpillover, method = "REML", family=nb())
    }else{
      gamCases <- gam(spillover ~ s(years), data=estSpillover, method = "REML", family=nb())
    }
    predictCases <- as.vector(predict(gamCases, newdata=data.frame(years=years), type="response"))
    predictMean <- mean(predictCases)
    theta <- gamCases$family$getTheta(TRUE)
    return(list(predictedMeans=predictCases, theta=theta))
  })

  means <- do.call(rbind, lapply(params.spillover, function(pp){ pp$predictedMeans }))
  colnames(means) <- paste0("spillover.mean_", years)
  population <- cbind(population, means)
  population$spillover.theta = unlist(lapply(params.spillover, function(pp){ pp$theta }))
  return(population)
}



# Generate spillovers per area
generateSpillover <- function(spilloverEstPerArea, alpha, spilloversEmp, nonZeroRatio=NULL, numReplicates=1e3){

  if(is.null(nonZeroRatio)){
    stopifnot(length(grep("spillover.nonZero", colnames(spilloverEstPerArea)))>0)
  }
  
  # check input parametes
  stopifnot(c("spillover.theta") %in% colnames(spilloverEstPerArea))
  stopifnot(length(alpha) == 2)
  stopifnot(length(grep("spillover.mean", colnames(spilloverEstPerArea)))>0)
            
  meanIdx <- grep("spillover.mean", colnames(spilloverEstPerArea))
  names(meanIdx) <- sub("spillover.mean_", "", colnames(spilloverEstPerArea)[meanIdx])

  # spillover process
  spilloverLst <- lapply(1:nrow(spilloverEstPerArea), function(i){
    yearReplicates <- do.call(cbind, lapply(meanIdx, function(yy){
      spill_over <- rnbinom(n = numReplicates,
                            mu=spilloverEstPerArea[i,yy],
                            size=spilloverEstPerArea$spillover.theta[i])
      spill_over <- rdirmnom(numReplicates, 
                             spill_over, 
                             alpha)[,1]
    }))
    colnames(yearReplicates) <- names(meanIdx)
    return(yearReplicates)
  })

  # zero-inflation process
  if(nonZeroRatio<1 || is.null(nonZeroRatio)){
    for(rr in 1:numReplicates){
      # mix year randomly to mix up zero years
      #spilloversEmp <- sample(spilloversEmp)
      oId <- order(spilloversEmp, decreasing=F)
      for(aa in seq_along(spilloverLst)){
        if(is.null(nonZeroRatio)){
          nonZeroRatio <- spilloverEstPerArea$spillover.nonZero[aa]
        }
        zeroSpillover <- sort(rbinom(length(meanIdx), size=1, prob=nonZeroRatio))
        spilloverLst[[aa]][rr,oId[zeroSpillover==0]] <- 0
      }
    }
  }
  return(spilloverLst)
}

sseEcdf <- function(x, empric, simulated){
  diff <- ecdf(empric)(x) - ecdf(simulated)(x)
  sum(diff^2)
}  


estimateZeroInflationRatio <- function(spilloverEstPerArea, spilloversEmp, alpha, numReplicates=10){
  
  # check input parameters
  stopifnot(c("spillover.theta") %in% colnames(spilloverEstPerArea))
  stopifnot(length(alpha) == 2)
  stopifnot(length(grep("spillover.mean", colnames(spilloverEstPerArea)))>0)
  
  meanIdx <- grep("spillover.mean", colnames(spilloverEstPerArea))
  names(meanIdx) <- sub("spillover.mean_", "", colnames(spilloverEstPerArea)[meanIdx])

  spilloverLst <- lapply(1:nrow(spilloverEstPerArea), function(i){
    # spillover process
    yearReplicates <- do.call(cbind, lapply(meanIdx, function(yy){
      spill_over <- rnbinom(n = numReplicates,
                            mu=spilloverEstPerArea[i,yy],
                            size=spilloverEstPerArea$spillover.theta[i])
      spill_over <- rdirmnom(numReplicates, 
                             spill_over, 
                             alpha)[,1]
    }))
    colnames(yearReplicates) <- names(meanIdx)
    return(yearReplicates)
  })

  # mix year randomly to mix up zero years
  oId <- order(spilloversEmp, decreasing=F)
  sseRes <- matrix(NA, nrow=0, ncol=2)
  zi <- c(seq(0,1,10^-1))
  for(it in 1:10){
    message(it)
    fit <- unlist(lapply(zi, function(zz){
      yearReplicates <- lapply(1:numReplicates, function(rr){
        singleYear <- do.call(rbind, lapply(seq_along(spilloverLst), function(aa){
          zeroSpillover <- sort(rbinom(length(meanIdx), size=1, prob=zz))
          spilloverLst[[aa]][rr,oId[zeroSpillover==0]] <- 0
          return(spilloverLst[[aa]][rr,])
        }))
        colSums(singleYear)
      })
      spilloversSim <- yearReplicates[[1]]
      sum(unlist(lapply(yearReplicates, function(yr){
        sseEcdf(seq(from=0, to=1, by=10^-1), spilloversEmp, yr)
      })))
    }))
    sseRes <- rbind(sseRes, cbind(sse=fit, zi=zi))
    
    if(it<=2){
      ziIdx <- order(fit, decreasing = F)
      ziRange <- range(zi[ziIdx][1:3], zi[ziIdx][1]+0.05, zi[ziIdx][1]-0.05)
      zi <- seq(max(0,ziRange[1]), min(1, ziRange[2]), 10^-2)
    }
  }
  return(sseRes)
}


### OUTBREAK SIMULATION FUNCTIONS
simulateOutbreakNetwork = function(
  total.cases=5000,
  index.cases,
  R,
  R.dispersion,
  shape1.infec,
  shape2.infec,
  shape.incub,
  scale.incub,
  shape.death,
  scale.death)
{
  # simulate index cases, which is imported
  # limit spillover to the total.cases numbers
  if(index.cases>total.cases)
    index.cases <- total.cases+1
  edges = matrix(NA, nrow=index.cases, ncol=2)
  colnames(edges) <- c("Parent","Offspring")
  edges[,"Parent"] <- 0
  edges[,"Offspring"] <- 1:index.cases
  total.cases_count = index.cases
  time.infec = sort(rbeta(index.cases, shape1=shape1.infec, shape2=shape2.infec)*365)
  time.incub = time.infec + rgamma(index.cases, shape=shape.incub, scale=scale.incub)
  if(!is.na(shape.death) & !is.na(scale.death)){
    time.death = time.incub + rgamma(index.cases, shape=shape.death, scale=scale.death)
  }else{
    time.death = rep(NA, index.cases)
  }
  # specify that this first cases needs to be visited and offspring simulated
  tovisit = 1:index.cases
  
  # continue to simulate cases until we have enough or gets extinct
  while(total.cases_count < total.cases & length(tovisit) > 0){
    if(length(tovisit) > 0){
      # if we have local cases who need to have their offspring simulated
      # simulate offspring and their status
      parent = tovisit[1]
      if(is.null(R.dispersion) | is.na(R.dispersion)){
        number.offspring <- rpois(n=1, lambda=R)
      }else{
        # An alternative parametrization (often used in ecology) is by the mean mu (see above), 
        # and size, the dispersion parameter, where prob = size/(size+mu). 
        # The variance is mu + mu^2/size in this parametrization.
        number.offspring <- rnbinom(n=1, mu=R, size=R.dispersion)
      }
      
      if(number.offspring > 0){
        # if there are any offspring
        # add them to the network edge list and list of nodes to visit
        newToVisit = max(edges)+(1:number.offspring)
        tovisit = c(tovisit, newToVisit)
        edges = rbind(
          edges,
          cbind(
            rep(parent, number.offspring),
            max(edges)+(1:number.offspring)))
        time.infec = c(time.infec, time.incub[parent] + (time.death[parent] - time.incub[parent]) * runif(number.offspring))
        time.incub = c(time.incub, time.infec[newToVisit] + rgamma(number.offspring, shape=shape.incub, scale=scale.incub))
        if(!is.na(shape.death) & !is.na(scale.death)){
          time.death = c(time.death, time.incub[newToVisit] + rgamma(number.offspring, shape=shape.death, scale=scale.death))
        }else{
          time.death = rep(NA, number.offspring)
        }
        total.cases_count = total.cases_count + number.offspring
      }
      
      # remove the node we just visited from the list of nodes to visit
      tovisit = tovisit[-1]
    } 
  }
  
  # mark not visited cases, resp. cases of which no offsprings where simulated
  incomplete = rep(0, dim(edges)[1])
  incomplete[tovisit] <- 1
  return(list(edges=edges, times=cbind(time.infec,time.incub,time.death), 
              incomplete=incomplete))
}



simulateVaccination = function(ntw, 
                               thresholdCases=10, 
                               thresholdTimeFrame=28, 
                               vacFrac_sp=NA, 
                               vacFrac_h2h=NA, 
                               vacDelay=NA, 
                               vacDuration=0, 
                               vacEfficacy=NA, 
                               protDelay=10, 
                               vac2Delay=NULL, 
                               vac2Efficacy=NA){
  require(igraph)
  
  # at vac times 
  ntw$times <- cbind(ntw$times, time.vac1=NA, time.vac2=NA, time.prot=NA)
  ntw$toPrun <- rep(F, dim(ntw$times)[1])
  
  # vac only if number of cases exceeds threshold and cases occur within a month
  ntw$VacCampainStarted <- checkThresholdCriteria(ntw$times[,"time.incub"], thresholdCases, thresholdTimeFrame)


  if(!is.na(ntw$VacCampainStarted)){
    startVac <- vacDelay + ceiling(ntw$VacCampainStarted)
    # selected at time of vac start not yet infected cases
    vacToday = which(floor(ntw$times[,"time.incub"])>=startVac)

    if(length(vacToday) > 0){
      # vacFrac is probability of vaccination
      vacYes = ifelse(ntw$edges[,1][vacToday]==0, 
                      rbinom(length(vacToday), 1, vacFrac_sp), 
                      rbinom(length(vacToday), 1, vacFrac_h2h))
      # determine cases which will be vac
      vacToday = vacToday[vacYes==1]
      
      # set vaccination1 time
      ntw$times[vacToday, "time.vac1"] = runif(length(vacToday), min=startVac, max=startVac+vacDuration)
      # determine protection cases
      vacProtect = rbinom(length(vacToday), 1, vacEfficacy)
      
      # set vaccination2 time
      if(!is.na(vac2Delay)){
        ntw$times[vacToday, "time.vac2"] <- ntw$times[vacToday, "time.vac1"]+vac2Delay
        effInc <- (vac2Efficacy-vacEfficacy)/(1-vacEfficacy)
        vac2Protect = rbinom(length(vacToday[vacProtect==0]), 1, effInc)
        vacProtect[vacProtect==0] <- 2*vac2Protect
      }
      
      if(sum(vacProtect) > 0){
        # set protection time
        protToday = vacToday[vacProtect==1]
        ntw$times[protToday, "time.prot"] = ntw$times[protToday, "time.vac1"] + protDelay
        protToday = vacToday[vacProtect==2]
        ntw$times[protToday, "time.prot"] = ntw$times[protToday, "time.vac2"] + protDelay
        
        # prune protected cases from the outbreak
        protToday = vacToday[vacProtect>0]
        protCases <- protToday[ntw$times[protToday,"time.prot"] < ntw$times[protToday,"time.infec"]]
        if(length(protCases) > 0){
          g <- graph_from_edgelist(apply(ntw$edges,2, as.character), directed=TRUE)
          toDeleteCases <- lapply(as.character(protCases), function(dd){
            subcomponent(g, v=dd, mode="out")
          })
          toDeleteCases <- as.integer(unique(names(unlist(toDeleteCases))))
        }else{ 
          toDeleteCases = NULL
        }
        ntw$toPrun = ntw$edges[,2] %in% toDeleteCases
      }
    }
  }
  return(ntw)
}


checkThresholdCriteria <- function(times, thresholdCases=10, thresholdTimeFrame=28){
  times <- sort(times)
  startCampain <- unlist(lapply(seq_along(times), function(i){
    slideWindow <- times>=times[i] & times<=(times[i]+thresholdTimeFrame)
    numCases <- sum(slideWindow)
    if(numCases>=thresholdCases)
      return(times[which(slideWindow)[thresholdCases]])
    else
      return(NA)
  }))
  if(any(!is.na(startCampain)))
    return(min(startCampain, na.rm=T))
  else
    return(NA)
}

evalCampainSimulation <- function(ntwLstAA){

    lapply(ntwLstAA, function(ntwLst){
      do.call(rbind, lapply(seq_along(ntwLst), function(i){
        # SKIP IF THERE IS NO SPILLOVER
        if(is.null(ntwLst[[i]]$VacCampainStarted)){
          return(c(outbreakSize=0, 
                   numSpillover=0,
                   # numH2h=NA,
                   numAverted=NA,
                   numAvertedHCW=NA,
                   numVacOnce=NA,
                   numVacTwice=NA,
                   numProtected=NA,
                   numDeath=NA,
                   numDeathAverted=NA,
                   numDeathAvertedHCW=NA,
                   vacCampainStarted=NA))
        }else{
          if(is.null(ntwLst[[i]])) stop(sprintf("Network should not be NULL, i=%i and j=%s",i,j))
          # VACCINATION RESPONSE
          return(c(outbreakSize=dim(ntwLst[[i]]$edges)[1], 
                   numSpillover=sum(ntwLst[[i]]$edges[,1]==0),
                   # numH2h=sum(ntwLst[[i]]$edges[,1]!=0),
                   numAverted=sum(ntwLst[[i]]$toPrun),
                   numAvertedHCW=sum(ntwLst[[i]]$toPrun & ntwLst[[i]]$edges[,1]!=0),
                   numVacOnce=sum(!is.na(ntwLst[[i]]$times[,"time.vac1"])),
                   numVacTwice=sum(!is.na(ntwLst[[i]]$times[,"time.vac2"])),
                   numProtected=sum(!is.na(ntwLst[[i]]$times[,"time.prot"])),
                   numDeath=sum(ntwLst[[i]]$death),
                   numDeathAverted=sum(ntwLst[[i]]$death & ntwLst[[i]]$toPrun),
                   numDeathAvertedHCW=sum(ntwLst[[i]]$death & ntwLst[[i]]$toPrun & ntwLst[[i]]$edges[,1]!=0),
                   vacCampainStarted=ntwLst[[i]]$VacCampainStarted))
        }
      }))
    })

}

statsSummary <- function(lst){
  stats <- do.call(cbind, lapply(lst, function(ll){
    qt <- quantile(ll, probs=c(0.025, 0.25, 0.75, 0.975), na.rm=T)
    rbind(Min=min(ll, na.rm=T),
          '2.5%'=qt[1],
          '25%'=qt[2],
          Median=median(ll, na.rm=T),
          Mean=mean(ll, na.rm=T),
          '75%'=qt[3],
          '97.5%'=qt[4],
          Max=max(ll, na.rm=T))
  }))
  colnames(stats) <- colnames(lst)
  return(stats)
}

getSecondRing <- function(ids, ntw.vac){
  sIds <- ntw.vac$edges[ntw.vac$edges[,1] %in% ids,"Offspring"]
  if(length(sIds)>0){
    idx <- ntw.vac$time[sIds,"time.incub"]>=ntw.vac$VacCampainStarted
    return(sum(idx) + getSecondRing(sIds, ntw.vac))
  }else{
    return(0)
  }
}


# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load source
source("code/outbreakFunction.R")

library(igraph)

##############
# Read Scenarios
scenarios <- read.csv("data/mers_scenarios.csv", header=T, stringsAsFactors = F)

# Read catchment area information
fnLst <- list.files("results", pattern="mers_networksSimulation")
postfixLst <- sub("mers_networksSimulation_","", fnLst)
postfixLst <- sub(".RDS","", postfixLst)

# Read health care worker information
hcwTab <- read.csv("data/HealthCareWorkers.csv", stringsAsFactors = F)
hcwPer10000 <- rowSums(hcwTab[,c("MD_Per10.000population","Nurse_Per10000population")], na.rm=T)
names(hcwPer10000) <- hcwTab$ISO3

lapply(postfixLst, function(postfix){
  # READ IN SPILLOVER INFORMATION
  clusters <- read.csv(sprintf("results/mers_spillovers_%s.csv", postfix), stringsAsFactors = F)
  useCol <- grep("^SO[0-9]", colnames(clusters), invert=F)
  spill_over <- clusters[,useCol]
  if(postfix == "adm1"){ 
    population <- clusters$pop2015 
  }else{
    population <- clusters$population/clusters$numCatchAreas
  }
  hcw <- hcwPer10000[clusters$GID_0]*population/10000
  
  # Read outbreak information
  ntwLstAll <- readRDS(file=sprintf("results/mers_networksSimulation_%s.RDS", postfix))
  
  # RUN Vaccination
  avertedMat <- lapply(1:dim(scenarios)[1], function(ss){
  #lapply(1:14, function(ss){
    set.seed(1) # set random value generator to same starting value
    message(scenarios[ss,]$Scenarios)
    
    scen <- scenarios[ss,]
    scen[scen == "NULL"] <- NA
    
    ntwLst <- ntwLstAll[[as.character(scen$R0)]]
    
    # PARAMETER SET UP
    thresholdCases = scen$thresholdCases # number of cases needed to activate vaccination campaign
    thresholdTimeFrame = scen$thresholdTimeFrame
    vacEfficacy = scen$vacEfficacy #(0.35) efficacy of vaccine
    vac2Efficacy = scen$vac2Efficacy #(0.7) efficacy of vaccine
    vacFrac_sp = scen$vacFrac_sp # vaccination coverage
    vacFrac_h2h = scen$vacFrac_h2h # vaccination coverage
    vacDelay = scen$vacDelay # number of days needed to start vaccination campaign after thresholdCases is reached
    vac2Delay = scen$vac2Delay # number of days needed to start vaccination campaign after thresholdCases is reached
    vacDuration = 0 # number of days needed for vaccination campaign
    protDelay = scen$protDelay #(7) number of days needed to get protection after vaccination
    
    # Simulate vaccination campaign
    # By replicates
    res <- lapply(1:dim(spill_over)[2], function(j){
      # By cluster
      do.call(rbind, lapply(1:dim(spill_over)[1], function(i){
        #message(paste(ss, j, i, sep=" "))
        ##SKIP IF THERE IS NO SPILLOVER
        #if(spill_over[i,j] == 0){
        if(is.null(ntwLst[[i]][[j]])){
          return(c(outbreakSize=0, 
                   numSpillover=0,
                   numH2H=NA,
                   numCasesAfterStart=NA,
                   numRingIndex=NA,
                   numRingIndexH2H=NA,
                   numAverted=NA,
                   numAvertedH2H=NA,
                   numVacOnce=NA,
                   numVacTwice=NA,
                   numProtected=NA,
                   numDeath=NA,
                   numDeathAverted=NA,
                   numDeathAvertedH2H=NA,
                   vacCampaignStarted=NA,
                   incomplete=F))
        }else{
          ##VACCINATION RESPONSE
          #set.seed(1) # set random value generator to same starting value
          ntw.vac <- simulateVaccination(ntwLst[[i]][[j]], 
                                         thresholdCases=thresholdCases, thresholdTimeFrame=thresholdTimeFrame,
                                         vacFrac_sp=vacFrac_sp, vacFrac_h2h=vacFrac_h2h,
                                         vacDelay=vacDelay, vac2Delay=vac2Delay, vacDuration=vacDuration,
                                         vacEfficacy=vacEfficacy, vac2Efficacy=vac2Efficacy, protDelay=protDelay)
          if(!is.na(ntw.vac$vacCampaignStarted)){
            numCasesAfterStart <- sum(ntw.vac$time[,"time.incub"]>=ntw.vac$vacCampaignStarted)
            numRingIndex <- sum((ntw.vac$edges[,1]==0) & (ntw.vac$time[,"time.incub"]>=ntw.vac$vacCampaignStarted))
            # get offspring of index cases with infectious period prior to start time
            idx <- (ntw.vac$edges[,1]==0) & (ntw.vac$time[,"time.incub"]<ntw.vac$vacCampaignStarted)
            numRingIndexH2H <- getSecondRing(ntw.vac$edges[idx,"Offspring"], ntw.vac)
          }else{
            numCasesAfterStart=0
            numRingIndex=0
            numRingIndexH2H=0
          }
          return(c(outbreakSize=dim(ntw.vac$edges)[1], 
                   numSpillover=sum(ntw.vac$edges[,1]==0),
                   numH2H=sum(ntw.vac$edges[,1]!=0),
                   numCasesAfterStart=numCasesAfterStart,
                   numRingIndex=numRingIndex,
                   numRingIndexH2H=numRingIndexH2H,
                   numAverted=sum(ntw.vac$toPrun),
                   numAvertedH2H=sum(ntw.vac$toPrun & ntw.vac$edges[,1]!=0),
                   numVacOnce=sum(!is.na(ntw.vac$times[,"time.vac1"])),
                   numVacTwice=sum(!is.na(ntw.vac$times[,"time.vac2"])),
                   numProtected=sum(!is.na(ntw.vac$times[,"time.prot"])),
                   numDeath=sum(ntw.vac$death),
                   numDeathAverted=sum(ntw.vac$death & ntw.vac$toPrun),
                   numDeathAvertedH2H=sum(ntw.vac$death & ntw.vac$toPrun & ntw.vac$edges[,1]!=0),
                   vacCampaignStarted=ntw.vac$vacCampaignStarted,
                   incomplete=any(ntw.vac$incomplete==1)))
        }
      }))
    })
    saveRDS(res, file=sprintf("results/mers_valuesNtwVac_Scenario_%02i_%s.RDS", scen$Scenarios, postfix))
    # get averted per year
    tab <- do.call(rbind, lapply(seq_along(res), function(j){
      osize <- res[[j]][,"outbreakSize"]
      av <- ifelse(!is.na(res[[j]][,"vacCampaignStarted"]),
                   res[[j]][,"numAverted"],
                   NA)
      avH2H <- ifelse(!is.na(res[[j]][,"vacCampaignStarted"]),
                   res[[j]][,"numAvertedH2H"],
                   NA)
      numCampaigns <- !is.na(res[[j]][,"vacCampaignStarted"])
      numIndex <- res[[j]][!is.na(res[[j]][,"vacCampaignStarted"]),"outbreakSize"]
      numRingCasesAfterStart <- res[[j]][,"numCasesAfterStart"]
      numRingIndex <- res[[j]][,"numRingIndex"]
      numRingIndexH2H <- res[[j]][,"numRingIndexH2H"]
      spsize <- res[[j]][,"numSpillover"]
      numPopulation <- population[!is.na(res[[j]][,"vacCampaignStarted"])]
      numHCW <- hcw[!is.na(res[[j]][,"vacCampaignStarted"])]
      incomplete <- sum(res[[j]][,"incomplete"])
      return(data.frame(TotalAverted=sum(av, na.rm=T), 
                        TotalAvertedH2H=sum(avH2H, na.rm=T), 
                        TotalCases=sum(osize, na.rm=T), 
                        TotalSpillover=sum(spsize, na.rm=T),
                        Campaigns=sum(numCampaigns, na.rm=T),
                        Regimens_Index=sum(numIndex, na.rm=T),
                        Regimens_RingCasesAfterStart=sum(numRingCasesAfterStart, na.rm=T),
                        Regimens_RingIndex=sum(numRingIndex, na.rm=T),
                        Regimens_RingH2H=sum(numRingIndexH2H, na.rm=T),
                        Regimens_RingTotal=sum(numRingIndex, na.rm=T)+sum(numRingIndexH2H, na.rm=T),
                        Regimens_Population=sum(numPopulation, na.rm=T),
                        Regimens_HCW=sum(numHCW, na.rm=T),
                        incomplete=incomplete))
    }))
    return(tab)
  })
  names(avertedMat) <- paste("Scenario", scenarios$Scenarios)
  saveRDS(avertedMat, file=sprintf("results/mers_averted_Scenario_all_%s.RDS", postfix))
  # avertedMat <- readRDS(file=sprintf("results/mers_averted_Scenario_all_%s.RDS", postfix))
  
  ##########################
  # stats of number averted
  sumAv <- do.call(rbind, lapply(seq_along(avertedMat), function(ii){
    tab <- statsSummary(avertedMat[[ii]])
    tab <- cbind(names(avertedMat)[ii], rownames(tab), tab)
  }))
  rownames(sumAv) <- NULL
  write.csv(sumAv, sprintf("results/mers_statistics_%s.csv", postfix), row.names = F, quote=T)

  ##########
  numScen <- 3
  ntwLst <- readRDS(file=sprintf("results/mers_valuesNtwVac_Scenario_%02i_%s.RDS", numScen, postfix))
  numCampaigns <- colSums(do.call(rbind, lapply(ntwLst, function(ntw){
    #startCampaigns <- ntw[,"vacCampaignStarted"]
    numCampaigns <- !is.na(ntw[,"vacCampaignStarted"])
  })))
  
  clusters <- read.csv(sprintf("results/mers_spillovers_%s.csv", postfix), stringsAsFactors = F)
  useCol <- grep("^SO[0-9]", colnames(clusters), invert=F)
  clusters <- clusters[,-useCol]
  write.csv(cbind(clusters, numCampaigns), file=sprintf("results/mers_campaign_%s.csv", postfix), row.names=F, quote=T)
  
})



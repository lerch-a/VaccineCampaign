
# setwd("~/Repos/VaccineCampaign/")
rm(list=ls())

# load source
source("code/outbreakFunction.R")

##############
# Read Scenarios
scenarios <- read.csv("data/mers_scenarios.csv", header=T, stringsAsFactors = F)
R <- unique(scenarios$R0)
thresholdCases = scenarios$thresholdCases # number of cases needed to activate vaccination campain
thresholdTimeFrame = scenarios$thresholdTimeFrame
vacEfficacy = scenarios$vacEfficacy #(0.35) efficacy of vaccine
vac2Efficacy = scenarios$vac2Efficacy #(0.7) efficacy of vaccine
vacFrac_sp = scenarios$vacFrac_sp # vaccination coverage
vacFrac_h2h = scenarios$vacFrac_h2h # vaccination coverage
vacDelay = scenarios$vacDelay # number of days needed to start vaccination campain after thresholdCases is reached
vac2Delay = scenarios$vac2Delay # number of days needed to start vaccination campain after thresholdCases is reached
vacDuration = 0 # number of days needed for vaccination campain
protDelay = scenarios$protDelay #(7) number of days needed to get protection after vaccination


# Read parameters
param <- read.csv("data/mers_parameter_estimates.csv", header=F, stringsAsFactors = F, row.names = 1)
param <- as.data.frame(t(param))


# Read catchment area information
fnLst <- list.files("results", pattern="mers_spillovers_a")
postfixLst <- sub("mers_spillovers_","", fnLst)
postfixLst <- sub(".csv","", postfixLst)
postfixLst <- postfixLst[postfixLst != "year"]

lapply(postfixLst, function(postfix){
  
  # Read spillover information
  clusters <- read.csv(sprintf("results/mers_spillovers_%s.csv", postfix), stringsAsFactors = F)
  useCol <- grep("^SO[0-9]", colnames(clusters), invert=F)
  clusters <- clusters[,useCol]
  
  
  ######################
  # RUN Outbreaks and check response trigger
  set.seed(1) # set random value generator to same starting value
  
  # PARAMETER SET UP
  total.cases=15000 # maximum outbreak size
  
  # i=1
  # j=1
  ntwLst <- lapply(R, function(rr){
    lapply(1:dim(clusters)[1], function(i){
      #message(paste(rr, i, sep="-"))
      spill_over <- clusters[i,]
      
      lapply(seq_along(spill_over), function(j){
        ##SKIP IF THERE IS NO SPILLOVER
        if(spill_over[j] == 0){
          return(NULL)
        }else{
          # simulate network
          ntw <- simulateOutbreakNetwork(total.cases = total.cases, 
                                         index.cases = as.numeric(spill_over[j]), 
                                         R = rr,
                                         R.dispersion = param$R0_dispersion[1],
                                         shape1.infec = param$timing_shape1[1],
                                         shape2.infec = param$timing_shape2[1],
                                         shape.incub = param$incubation_shape[1],
                                         scale.incub = 1/param$incubation_rate[1],
                                         shape.death = param$duration_shape[1],
                                         scale.death = 1/param$duration_rate[1])
  
          return(ntw)
        }
      })
    })
  })
  names(ntwLst) <- R
  saveRDS(ntwLst, file=sprintf("results/mers_networksSimulation_%s.RDS", postfix))

})


# setwd("~/Repos/VaccineCampaign/")
rm(list=ls())

library(extraDistr)
library(tidyverse)
library(glmmTMB)
source("code/helperFunctions.R")
source("code/spilloverFunction.R")

set.seed(1234)

# selected countries
countries <- c("OMN","QAT","SAU","JOR","ARE","YEM","KWT","BHR")

# get population
population <- read.csv("data/pop2015_africa_asia_adm1.csv", stringsAsFactors = F)
population <- subset(population, subset=population$GID_0 %in% countries, select=c("GID_0","NAME_0","GID_1","NAME_1","pop2015","Long","Lat"))

# load case data
cases <- read.csv("data/mers_case_reports.csv", stringsAsFactors = F)
cases$GID_1[is.na(cases$GID_1)] <- ""
cases$Cases <- rowSums(cases[,c("Active_Unknown","Recovered","Deceased")], na.rm=T)
cases$Secondary[is.na(cases$Secondary)] <- 0
cases$Cases <- cases$Cases - cases$Secondary
cases <- cases[grep("secondary", cases$Case.classification, invert = T),]
cases <- cases[,c("Country","GID_0","Region","GID_1","Year","Cases", "Case.classification")]

years <- seq(range(cases$Year)[1],range(cases$Year)[2])
years <- seq(2010,range(cases$Year)[2])
n.years = diff(range(cases$Year))+1

# aggregate spillovers per year
sp_per_year <- data.frame(aggregate(cases$Cases, by=list(YEAR=cases$Year, GID_1=cases$GID_1, GID_0=cases$GID_0), FUN=sum, na.rm=T))
colnames(sp_per_year)[colnames(sp_per_year)=="x"] <- "spillovers"

sp_per_year <- do.call(rbind, lapply(years, function(yy){
  do.call(rbind, lapply(c(population$GID_1,""), function(gg){
    entry <- sp_per_year[sp_per_year$YEAR==yy & sp_per_year$GID_1==gg,]
    if(nrow(entry)>0)
      return(entry)
    else
      return(data.frame(YEAR=yy,GID_1=gg,GID_0=substr(gg,1,3),spillovers=0))
  }))
}))
sp_per_year$YEAR <- as.integer(sp_per_year$YEAR)

sp_per_year <- data.frame(aggregate(sp_per_year$spillovers, by=list(YEAR=sp_per_year$YEAR, GID_1=sp_per_year$GID_1, GID_0=sp_per_year$GID_0), FUN=sum, na.rm=T))
colnames(sp_per_year)[colnames(sp_per_year)=="x"] <- "spillovers"

sp_per_year <- sp_per_year[sp_per_year$GID_0 != "",]
write.csv(sp_per_year, file="results/mers_spillovers_year.csv", row.names=F, quote=T)

# remove case data without adm1
sp_per_year <- sp_per_year[sp_per_year$GID_1 != "",]


##########
# estimate spillover process parameter
spilloversEst <- estimateSpilloverZNB(sp_per_year, population=population, zeroYr=F, startYear=2010, endYear=2020, useRange=T, numReplicates=1e3)
write.csv(spilloversEst, file="results/mers_spilloversEst_ZNB.csv", row.names=F, quote=T)

# estimte non-zero spillover process parameter
spilloversEmp <- read.csv(file="results/mers_spillovers_year.csv", stringsAsFactors=F)
spilloversEmp <- aggregate(spilloversEmp$spillovers, by=list(YEAR=spilloversEmp$YEAR), FUN=sum, na.rm=TRUE)
rownames(spilloversEmp) <- spilloversEmp$YEAR

# generate spillovers
spill_l=list()
for(ii in 1:1000){
  spill_z=rbinom(n=nrow(spilloversEst),size=1,prob=1-spilloversEst$est.pr0)
  spill_nb=rnbinom(n=nrow(spilloversEst),size=spilloversEst$est.theta,mu=spilloversEst$est.mu)
  spill_nb=spill_nb*spill_z
  spill_l[[ii]]=data.frame(spillovers=spill_nb)
}
spill_est=do.call("cbind",spill_l)
names(spill_est)=paste0("spillovers_",1:1000)
spill_est$year=spilloversEst$YEAR
spill_est$GID_0=spilloversEst$GID_0
spill_est$GID_1=spilloversEst$GID_1

spill_est_l = spill_est %>% pivot_longer(cols=starts_with("spillovers_"),names_to="sim",values_to="spillovers")
spill_yr = spill_est_l %>% group_by(year,sim) %>% summarize(spillovers=sum(spillovers))

pdf("plots/distSpilloverPerYear_ecdf_mers_glmmZNB.pdf", height=4)
plot(ecdf(spilloversEmp[,2]), cex=0.5, main="MERS", ylab="ECDF", xlab="Cases", las=1, col=1, xlim=c(0,1500))
#lines(ecdf(sample(spilloversSimYear, size = length(spilloversEmp))), col=2, cex=0.5)
lines(ecdf(spill_yr$spillovers), col=2, cex=0.5)
dev.off()


# select current years
spill_est_recent=spill_est_l %>% filter(year>2015)
## Get a random 1000 to reduce total sims to 1000
spill_est_recent = spill_est_recent %>% group_by(GID_1) %>% 
  sample_n(1000) %>% 
  mutate(sim_num=1:1000) %>%
  ungroup() %>% select(-c(year,sim))
## Reorg table
spillover_sims=spill_est_recent %>% pivot_wider(names_from=sim_num,values_from = spillovers,names_prefix="SO")

write.csv(spillover_sims, file="results/mers_spillovers_adm1.csv", row.names=F, quote=T)

pdf("plots/distSpilloverLastYears_ecdf_mers_glmmZNB.pdf", height=4)
plot(ecdf(spilloversEmp[as.character((2020-4):2020),2]), cex=0.5, main="MERS", ylab="ECDF", xlab="Cases", las=1, col=1, xlim=c(0,600))
lines(ecdf(spill_yr$spillovers[spill_yr$year>2015]), col=2, cex=0.5)
dev.off()


#######################
# Distribute spillovers to hospital
### Admin1
clusters <- read.csv("data/Hostpital_1k_admin1.csv")
clusters <- clusters[clusters$GID_0 %in% countries,]

sp <- matrix(NA,nrow=nrow(clusters), ncol=1000)
for(i in 1:nrow(spillover_sims)) {
  idxC <- which(clusters$GID_1==as.character(spillover_sims$GID_1[i]))
  for(j in 1:ncol(sp)) {
    ##Allocate spillovers by adm1_pct for each cluster
    sp_size=pull(spillover_sims[i,j+2])
    sp[idxC,j] <- rmultinom(1, size=sp_size, prob=1/clusters$numCatchAreas[idxC])
  }
}
colnames(sp) <- colnames(spillover_sims)[-(1:2)]
write.csv(cbind(clusters, sp), file="results/mers_spillovers_adm1_Hosp.csv", row.names=F, quote=T)

### Admin2
clusters <- read.csv("data/Hostpital_10k_admin2.csv")
clusters <- clusters[clusters$GID_0 %in% countries,]

# Allocate spillovers to clusters
sp <- matrix(NA,nrow=nrow(clusters), ncol=1000)
for(i in 1:nrow(spillover_sims)) {
  idxC <- which(clusters$GID_1==as.character(spillover_sims$GID_1[i]))
  for(j in 1:ncol(sp)) {
    ##Allocate spillovers by adm1_pct for each cluster
    sp_size=pull(spillover_sims[i,j+2])
    sp[idxC,j] <- rmultinom(1, size=sp_size, prob=1/clusters$numCatchAreas[idxC])
  }
}
colnames(sp) <- colnames(spillover_sims)[-(1:2)]
write.csv(cbind(clusters, sp), file="results/mers_spillovers_adm2.csv", row.names=F, quote=T)




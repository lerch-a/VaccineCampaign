
# functions for transforming probabilities
logit <- function(x){ log(x/(1-x)) }
logit.inv <- function(x){ exp(x)/(1+exp(x)) }

convert2doy <- function(year=NULL, month=NULL, week=NULL, day=NULL, checkRange=F){
  
}

week2day <- function(week=NULL){
  doy <- (week*7)-4
  return(doy)
}

month2day <- function(month=NULL){
  doy <- as.numeric(difftime(as.Date(paste('2016', month, 15, sep='-')),
                            as.Date('2015-12-31')))
  return(doy)
}
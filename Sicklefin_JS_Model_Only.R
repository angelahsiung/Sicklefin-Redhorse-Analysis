rm(list=ls())

# Set working directory
basedirectory <- "C:/Users/solit/Documents/GitHub/Sicklefin-Redhorse-Analysis"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "RMark","lubridate", "tidyverse", "rjags", "coda")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

###################################################
## Operational years (occasion) for each gear type
## Seine: 2014-2016 (1-3)
## Fyke: 2016-2018 (3-5)
## PIT: 2017-2018 (4-5)
###################################################

# Read in capture history
ms.cap.hist<-read.csv("Sicklefin Capture History.csv")
ms.cap.hist<-ms.cap.hist[-278,]
Sex<-ms.cap.hist$Sex

js.CH<-ms.cap.hist
js.CH$Sex<-NULL

# Convert capture history to data for JS model
js.CH[is.na(js.CH)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(js.CH)[1]), js.CH)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL
head(CH.du)
nz<-1000 # augmented individuals

ms.js.CH.aug<-rbind(CH.du, matrix(1, ncol=dim(CH.du)[2], nrow=nz)) # data agumentation

# Jags data
# Known z values
z.known<-CH.du
n1<-n2<-rep(NA, nrow(CH.du))
for(i in 1:nrow(CH.du)){
  n1<-min(which(CH.du[i,]>1))
  n2<-max(which(CH.du[i,]>1))
  z.known[i,n1:n2]<-2
  for(j in 1:ncol(CH.du)){
    z.known[i,j]<-ifelse(z.known[i,j]==1, NA, z.known[i,j])
  }
}
z.known<-rbind(z.known, matrix(NA, ncol=dim(CH.du)[2], nrow=nz))


# Calculating new and recaps
# get.first<-function(x) min(which(x>1))
# f<-apply(ms.cap.hist, 1, get.first)
# Recap<-NewCap<-rep(NA, 4) # recaps and new caps for 2015-2018
#
# for(n in 1:4){
#   NewCap[n]<-length(f[f==(n+1)]) # number of newly captured individuals 2015-2018
#   Recap[n]<-nrow(ms.cap.hist[!is.na(ms.cap.hist[,n+1])&ms.cap.hist[,n+1]>1&!is.na(ms.cap.hist[,n]),]) # number of recaptured individuals
# }

# Effor data
fykeEffort <- c(0, 0, 0, 2.5, 7, 4) # number of sampling days
##pitEffort <- c(0, 0, 0, 0, 1, 5.25) # ratio of number of sampling days from all antenna in 2017 and 2018, 66 and 347 respectively (Tried to use number of sampling days but model would not run)
pitEffort <- c(0, 0, 0, 0, 66, 347)

# Bundle data
dat3<-list(y = ms.js.CH.aug, z=z.known, n.occ = dim(ms.js.CH.aug)[2], Sex=c(Sex, rep(1, nz/2), rep(2, nz/2)), M=dim(ms.js.CH.aug)[1], fykeEffort=fykeEffort, pitEffort=pitEffort)

# Initial values

js.multistate.init <- function(ch, nz){
  ch[ch==1] <- NA
  state <- ch
  for (i in 1:nrow(ch)){
    n1 <- min(which(ch[i,]>1))
    n2 <- max(which(ch[i,]>1))
    state[i,n1:n2] <- 2
  }
  state[state==0] <- NA
  get.first <- function(x) min(which(!is.na(x)))
  get.last <- function(x) max(which(!is.na(x)))
  f <- apply(state, 1, get.first)
  l <- apply(state, 1, get.last)
  for (i in 1:nrow(ch)){
    state[i,1:f[i]] <- 1
    if(l[i]!=ncol(ch)) state[i, (l[i]+1):ncol(ch)] <- 2
    state[i, f[i]] <- 2
  }
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  return(state)
}
z.init<-js.multistate.init(CH.du, nz)
z.init[z.known==2]<-NA


inits <- function() {list(z=cbind(rep(NA, nrow(z.init)),z.init[,-1]), phi= runif(2, 0, 1), gamma=runif(5, 0, 1), p0pit=runif(1, 0, 0.01))}

# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "phi", "gamma", "Nsuper", "N", "B", "psi")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 3



sr.ms.js.jm1<-jags.model(data=dat3, inits = inits,
                         file = "ms_js_phiSex_pGearEffort.jags",
                         n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc1 <- coda.samples(sr.ms.js.jm1, params, n.iter=ni)

out.jc1<-summary(sr.ms.js.jc1) # model output

jc1.stats<-out.jc1$statistics # parameter estimates
jc1.quants<-out.jc1$quantiles # 95% confidence intervals

plot(sr.ms.js.jc1, ask=TRUE)

sr.ms.js.jc2 <- coda.samples(sr.ms.js.jm1, params, n.iter=1000)


plot(sr.ms.js.jc2, ask=TRUE)





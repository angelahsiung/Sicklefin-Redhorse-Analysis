rm(list=ls())

# Set working directory
basedirectory <- "C:/Users/solit/Documents/GitHub/Sicklefin-Redhorse-Analysis"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "lubridate", "tidyverse", "rjags","jagsUI", "coda","parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

load("SRH_js_phi0_gamTime_pGearEffort_100000.gzip")

###################################################
## Operational years (occasion) for each gear type
## Seine: 2014-2016 (2-4)
## Fyke: 2016-2018 (4-6)
## PIT: 2017-2018 (5-6)
###################################################


#### Creating matrix of individual and tag type ("1" for HDX tags, and "0" for FDX-B tags)

# Read in capture data
raw.dat<-read.csv("Sicklefin_Capture_Data.csv")
raw.dat<-raw.dat[raw.dat$Sex!="Unknown",] # remove rows where Sex is "Unknown"
raw.dat$CollectDate<-strptime(raw.dat$CollectDate, format="%d-%b-%y") # convert date format
raw.dat$CollectYear.1<-ifelse(is.na(raw.dat$CollectDate), raw.dat$CollectYear, year(raw.dat$CollectDate))

# Convert tag type
raw.dat$Tag_Type<-ifelse(raw.dat$PIT_Type=="FDX-B",0,1) # if fish was tagged wtih FDX-B, cannot be detected by PIT antenna. If tagged with HDX tag, can be detected by antenna

# Create table with individual and tag detectability over time ("0" for FDX-B tags, "1" for HDX tags)
tag.dat<-data.frame(with(raw.dat, table(raw.dat$Individual_ID, raw.dat$Tag_Type, raw.dat$CollectYear.1)))
tag.dat<-tag.dat[tag.dat$Freq>0,]
tag.dat<-spread(tag.dat, Var3, Var2)

coalesce_by_column <- function(df) {
  return(coalesce(df[1], df[2]))
}

tag.dat<-tag.dat%>%
  group_by(Var1) %>%
  summarise_all(coalesce_by_column) 

tag.dat<-data.frame(tag.dat[,3:7])

get.first<-function(x){min(which(!is.na(x)))}
f<-apply(tag.dat, 1, get.first)

for(i in 1:nrow(tag.dat)){
  for(j in (f[i]+1):ncol(tag.dat)){
    tag.dat[i,j]<-ifelse(tag.dat[i,j-1]==0, 0, 1)
    }
}        # create matrix of individual detectability by PIT antenna based on the tag they received


# prepare tag data for model
tag.dat<-tag.dat[,1:5]
colnames(tag.dat)<-NULL
tag.dat[is.na(tag.dat)]<-0
tag.dat[6,5]<-1 #new HDX tag added to this individual in 2018
tag.dat<-tag.dat[-278,] #remove fish that was supposedly captured by seine in Valley River
tag.dat <- apply(tag.dat, 2, as.numeric)

##### Preparing data for model

# capture history
ms.cap.hist<-read.csv("Sicklefin Capture History.csv")
ms.cap.hist<-ms.cap.hist[-278,] #remove fish that was supposedly captured by seine in Valley River

# prepare data for model
# Sex<-ifelse(ms.cap.hist$Sex==1, 0, 1) #0 is female, 1 is male

js.CH<-ms.cap.hist
js.CH$Sex<-NULL

# Convert CJS capture history to data for JS model
js.CH[is.na(js.CH)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(js.CH)[1]), js.CH)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL
head(CH.du)
nz<-2000 # augmented individuals

# augmenting individuals to original capture history
ms.js.CH.aug<-rbind(CH.du, matrix(1, ncol=dim(CH.du)[2], nrow=nz)) 

#augment PIT tag data
tag.dat.aug<-as.matrix(cbind(rep(0,dim(tag.dat)[1]), tag.dat)) 
tag.dat.aug<-rbind(tag.dat.aug, matrix(0, ncol=dim(tag.dat.aug)[2], nrow=nz))

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


# Effor data
#effort<-read.csv("SRH_Effort_v2.csv", header=T)
#fykeEffort <- tapply(effort$UB_FYKE, effort$Year, sum)
fykeEffort<- c(0, 0, 0, 3, 5, 4) # number of sampling days
#pitEffort <- aggregate(effort[,4:7], by=list(effort$Year), sum)
pitEffort <- c(0, 0, 0, 0, 66, 347)

# Bundle data
dat3<-list(y = ms.js.CH.aug, n.occ = dim(ms.js.CH.aug)[2], M=dim(ms.js.CH.aug)[1], fykeEffort=fykeEffort, pitEffort=pitEffort, tag.dat.aug=as.matrix(tag.dat.aug)) #Sex=c(Sex, rep(NA, nz))

# Initial values
js.multistate.init <- function(ch, nz){
  ch[ch==1] <- NA
  state <- ch
  for (i in 1:nrow(ch)){
    n1 <- min(which(ch[i,]!=1))
    n2 <- max(which(ch[i,]!=1))
    state[i,n1:n2] <- 2
  }
  #state[state==0] <- NA
  get.first <- function(x) min(which(!is.na(x)))
  get.last <- function(x) max(which(!is.na(x)))   
  f <- apply(state, 1, get.first)
  l <- apply(state, 1, get.last)
  for (i in 1:nrow(ch)){
    state[i,1:f[i]] <- 1
    if(l[i]!=ncol(ch)) state[i, (l[i]+1):ncol(ch)] <- 3
    state[i, f[i]] <- 2
  }   
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  return(state)
}

z.init<-js.multistate.init(CH.du, nz)



inits <- function() {list(z=cbind(rep(NA, nrow(z.init)),z.init[,-1]), mean.phi=runif(1, 0, 1), p0seine=runif(1, 0, 1), p0fyke=runif(1, 0, 1), p0pit=runif(1, 0, 0.01), gamma=runif(5, 0, 1))} #ER=runif(1, 0, 1))} #Sex=c(rep(NA,length(Sex)), rbinom(nz ,1, 0.5)),

# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "mean.phi", "B", "gamma", "Nsuper", "N", "psi") #, "alpha", "beta") 

# MCMC settings
ni <- 500
nt <- 1
nb <- 100
nc <- 1

## Run models

#  nz=2000, nb=100, nt=1, ni=500, nc=1
ptm <- proc.time()
sr.ms.js.jm1 <- jags.model(data=dat3, inits = inits,
                            file = "srh_js_phi0_gamTime_pGearEffort_multistate.jags",
                            n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc1 <- coda.samples(sr.ms.js.jm11, params, n.iter=ni)

proc.time() - ptm

#ER=500, nz=1800, nb=10000, nt=1, ni=200000, nc=3
ptm <- proc.time()

sr.ms.js.jm14 <- autojags(dat3, inits, params,
                          model.file = "srh_js_phi0_gam0_pGearEffort.jags",
                          n.chains = nc, n.adapt = 2000, iter.increment=20000,
                          n.burnin=nb,n.thin=nt,
                          parallel=TRUE,n.cores=3,Rhat.limit=1.1, max.iter=ni, verbose=TRUE)

proc.time() - ptm



summary(sr.ms.js.jc13)
plot(sr.ms.js.jc13, ask=TRUE)
out.jc13<-summary(sr.ms.js.jc13)
stats.jc13<-out.jc13$statistics
quants.jc13<-out.jc13$quantiles

save(sr.ms.js.jc13, file="SRH_js_phi0_gamTime_pGearEffort_100000.gzip")

# plotting results

plot(stats.jc13[2:5,1], type="l", ylim=c(0, 1200))
lines(quants.jc13[2:5, 1], lty=2)
lines(quants.jc13[2:5, 5], lty=2)

plot(stats.jc13[6:10], type="l", ylim=c(1000, 2200))
lines(quants.jc13[6:10, 1], lty=2)
lines(quants.jc13[6:10, 5], lty=2)



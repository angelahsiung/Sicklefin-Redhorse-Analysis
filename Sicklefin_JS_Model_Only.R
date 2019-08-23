rm(list=ls())

# Set working directory
#basedirectory <- "C:/Users/solit/Documents/GitHub/Sicklefin-Redhorse-Analysis"
basedirectory <- "C:/Users/brang/Documents/GitHub/Sicklefin-Redhorse-Analysis"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "lubridate", "tidyverse", "rjags","jagsUI", "coda","parallel", "doParallel", "foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

#load("SRH_js_phi0_gamTime-pGearEffort_100000.gzip")

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


#Pulling out Brasstown data
brasstown<-raw.dat[raw.dat$Waterbody=="Brasstown Creek",]

# Create table with individual and tag detectability over time ("0" for FDX-B tags, "1" for HDX tags)
tag.dat<-data.frame(with(brasstown, table(brasstown$Individual_ID, brasstown$Tag_Type, brasstown$CollectYear.1)))
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
tag.dat[6,5:ncol(tag.dat)]<-1 #new HDX tag added to this individual in 2018
tag.dat <- apply(tag.dat, 2, as.numeric)

##### Preparing data for model
# organize data into pivot table
ms.dat<-dcast(brasstown, brasstown$Individual_ID+brasstown$CollectYear.1~brasstown$PriMethod)
ms.dat<-ms.dat[rowSums(ms.dat[,3:5])>0,] # keep rows where individual was captured at least once

for(i in 1:nrow(ms.dat)){
  for(j in 3:ncol(ms.dat)){
    if(ms.dat[i,j]>1){
      ms.dat[i,j]<-1
    }
  }
} # transform capture records into 1s and 0s

# Turn capture record into gear used to capture individual (e.g. if caught by PIT array, denoted as "P")
ms.dat$`PIT Array`<-ifelse(ms.dat$`PIT Array`==1, "P", 0)
ms.dat$`Seining`<-ifelse(ms.dat$`Seining`==1, "S", 0)
ms.dat$`Trapping (Weir or Fyke)`<-ifelse(ms.dat$`Trapping (Weir or Fyke)`==1, "T", 0)

# Assign simpler column names
colnames(ms.dat)<-c("ID", "Year", "PIT", "Sein", "Trap")

# Using paste function to denote "state" of captured individual based on gear(s) with which it was captured
for(i in 1:nrow(ms.dat)){
  for(j in 1:ncol(ms.dat)){
    ms.dat[i,j]<-ifelse(ms.dat[i,j]==0, '', ms.dat[i,j])
  }
  ms.dat$state[i]<-paste(ms.dat[i,3:5], collapse='')
}

# Turn table into wide shape
ms.dat<-ms.dat%>%
  spread(Year, state)

# extract capture history from table
ms.cap.hist<-ms.dat[,c(1,5:9)]

# combine multiple rows of the same individual
coalesce_by_column <- function(df) {
  return(coalesce(df[1], df[2]))
}

ms.cap.hist<-ms.cap.hist%>%
  group_by(ID) %>%
  summarise_all(coalesce_by_column)

for (i in 1:nrow(ms.cap.hist)){
  for (j in f[i]:ncol(ms.cap.hist)){
    ms.cap.hist[i,j]<-ifelse(is.na(ms.cap.hist[i,j]), 1, ifelse(ms.cap.hist[i,j]=="S", 2, ifelse(ms.cap.hist[i,j]=="T", 3, ifelse(ms.cap.hist[i,j]=="P", 4,ifelse(ms.cap.hist[i,j]=="PT", 5, ifelse(ms.cap.hist[i,j]=="ST", 6, ms.cap.hist[i,j]))))))
  }
}

brasstown.CH<-ms.cap.hist[,-1]
brasstown.CH<-apply(brasstown.CH, 2, as.numeric)
brasstown.CH<-as.matrix(brasstown.CH)

# Convert CJS capture history to CH for JS model
brasstown.CH[is.na(brasstown.CH)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(brasstown.CH)[1]), brasstown.CH)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL
head(CH.du)
nz<-1500# augmented individuals

# augmenting individuals to original capture history
brasstown.CH.aug<-rbind(CH.du, matrix(1, ncol=dim(CH.du)[2], nrow=nz)) 

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
effort<-read.csv("SRH_Effort_v2.csv", header=T)
fykeEffort <- tapply(effort$UB_FYKE, effort$Year, sum)
fykeEffort<- c(0, 0, 0, 3, 5, 4) # number of sampling days
pitEffort <- aggregate(effort[,4:5], by=list(effort$Year), sum)
pitEffort <- c(0, 0, 0, 0, 66, 186)

# Bundle data
dat3<-list(y = brasstown.CH.aug, n.occ = dim(brasstown.CH.aug)[2], M=dim(brasstown.CH.aug)[1], fykeEffort=fykeEffort, pitEffort=pitEffort, tag.dat.aug=as.matrix(tag.dat.aug)) #Sex=c(Sex, rep(NA, nz))

# Initial values

js.multistate.init <- function(ch, nz){
  ch[ch==1] <- NA
  state <- ch
  for (i in 1:nrow(ch)){
    n1 <- min(which(ch[i,]>1))
    state[i,n1:ncol(ch)] <- 1
    state[i, 1:(n1-1)]<- 0
  }
  # state[state==0] <- NA
  # get.first <- function(x) min(which(!is.na(x)))
  # get.last <- function(x) max(which(!is.na(x)))
  # f <- apply(state, 1, get.first)
  # l <- apply(state, 1, get.last)
  # for (i in 1:nrow(ch)){
  #   state[i,1:f[i]] <- 1
  #   if(l[i]!=ncol(ch)) state[i, (l[i]+1):ncol(ch)] <- 2
  #   state[i, f[i]] <- 2
  # }
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  return(state)
}
z.init<-js.multistate.init(CH.du, nz)



inits <- function() {list(z=cbind(rep(NA, nrow(z.init)),z.init[,-1]), mean.phi=runif(1, 0, 1)
                          ,p0seine=runif(1, 0, 1), p0fyke=runif(1, 0, 1), p0pit=runif(1, 0, 0.01) 
                          #,gamma=runif(5, 0, 1) 
                          ,ER=runif(1, 0, 1))}

# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "mean.phi", "ER", "Nsuper", "N", "psi") 
#codaOnly<-c("Nsuper", "N", "B", "psi")

# MCMC settings
ni <- 500
nt <- 1
nb <- 100
nc <- 3

## Run models

# ER=500, nz=1500, nb=10000, nt=1, ni=100000, nc=3
ptm <- proc.time()
sr.ms.js.jm13 <- jags.model(data=dat3, inits = inits,
                            file = "srh_js_phi0_gam0_pGearEffort.jags",
                            n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc13 <- coda.samples(sr.ms.js.jm13, params, n.iter=ni)

proc.time() - ptm

#ER=500, nz=1800, nb=10000, nt=1, ni=200000, nc=3
ptm <- proc.time()

sr.ms.js.jm14 <- autojags(dat3, inits, params,
                          model.file = "srh_js_phi0_gam0_pGearEffort.jags",
                          n.chains = nc, n.adapt = 100, iter.increment=100,
                          n.burnin=nb,n.thin=nt,
                          parallel=TRUE,n.cores=3,Rhat.limit=1.1, max.iter=ni, verbose=TRUE)

proc.time() - ptm



summary(sr.ms.js.jc13)
plot(sr.ms.js.jc13, ask=TRUE)
out.jc13<-summary(sr.ms.js.jc13)
stats.jc13<-out.jc13$statistics
quants.jc13<-out.jc13$quantiles

summary(sr.ms.js.jc14)
plot(sr.ms.js.jm14, ask=TRUE)
out.jc14<-summary(sr.ms.js.jc14)
stats.jc14<-out.jc14$statistics
quants.jc14<-out.jc14$quantiles

#save(sr.ms.js.jc13, file="SRH_js_phi0_gamTime_pGearEffort_100000.gzip")

# plotting results







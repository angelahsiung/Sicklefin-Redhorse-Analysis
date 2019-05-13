rm(list=ls())

# Set working directory
basedirectory <- "C:/Users/solit/Documents/GitHub/Sicklefin-Redhorse-Analysis"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "RMark","lubridate", "tidyverse", "rjags", "coda")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

load("SRH_JS_phiSex_pGearEffort_100000.gzip")

###################################################
## Operational years (occasion) for each gear type
## Seine: 2014-2016 (1-3)
## Fyke: 2016-2018 (3-5)
## PIT: 2017-2018 (4-5)
###################################################


#### Creating matrix of individual and tag type

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
tag.dat<-as.matrix(tag.dat[,1:5])
colnames(tag.dat)<-NULL
tag.dat[is.na(tag.dat)]<-0
tag.dat[6,5]<-1 #new HDX tag added to this individual in 2018
tag.dat<-tag.dat[-278,] #remove fish that was supposedly captured by seine in Valley River

##### Preparing data for model

# capture history
ms.cap.hist<-read.csv("Sicklefin Capture History.csv")
ms.cap.hist<-ms.cap.hist[-278,] #remove fish that was supposedly captured by seine in Valley River
Sex<-ms.cap.hist$Sex

js.CH<-ms.cap.hist
js.CH$Sex<-NULL

# Convert capture history to data for JS model
js.CH[is.na(js.CH)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(js.CH)[1]), js.CH)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL
head(CH.du)
nz<-3500 # augmented individuals

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
effort<-read.csv("SRH_Effort_v2.csv", header=T)
#fykeEffort <- tapply(effort$UB_FYKE, effort$Year, sum)
fykeEffort<- c(0, 0, 0, 3, 5, 4) # number of sampling days
#pitEffort <- aggregate(effort[,4:7], by=list(effort$Year), sum)
pitEffort <- c(0, 0, 0, 0, 66, 347)

# Bundle data
dat3<-list(y = ms.js.CH.aug, z=z.known, n.occ = dim(ms.js.CH.aug)[2], Sex=c(Sex, rep(1, nz/2), rep(2, nz/2)), M=dim(ms.js.CH.aug)[1], fykeEffort=fykeEffort, pitEffort=pitEffort, tag.dat.aug=tag.dat.aug)

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
params <- c("pFyke", "pSeine", "pPit", "phi", "Nsuper", "N", "B", "psi") #"gamma")

# MCMC settings
ni <- 50000
nt <- 1
nb <- 100
nc <- 3


# Run models
## 3000 augmented individuals, 100000 iterations
sr.ms.js.jm1 <- jags.model(data=dat3, inits = inits,
                         file = "ms_js_phiSex_pGearEffort.jags",
                         n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc1 <- coda.samples(sr.ms.js.jm1, params, n.iter=ni)

## 3500 augmented individuals, 50000 iterations
sr.ms.js.jm2 <- jags.model(data=dat3, inits = inits,
                           file = "ms_js_phiSex_pGearEffort.jags",
                           n.chains = nc, n.adapt = nb, quiet = F)
sr.ms.js.jc2 <- coda.samples(sr.ms.js.jm2, params, n.iter=ni)


out.jc1<-summary(sr.ms.js.jc1) # model output
out.jc2<-summary(sr.ms.js.jc2)

jc1.stats<-out.jc1$statistics # parameter estimates
jc1.quants<-out.jc1$quantiles # 95% confidence intervals
jc2.stats<-out.jc2$statistics # parameter estimates
jc2.quants<-out.jc2$quantiles
plot(sr.ms.js.jc1, ask=TRUE)

#save(sr.ms.js.jc1, file="SRH_JS_phiSex_pGearEffort.gzip")
#save(sr.ms.js.jc1, file="SRH_JS_phiSex_pGearEffort_100000.gzip")


# plotting results
pop.size<-data.frame(cbind(Year=c(2014:2018),pop.est=jc2.stats[6:10, 1], lower=jc2.quants[6:10, 1], upper=jc2.quants[6:10, 5]))
pop.size.plot<-ggplot(pop.size, aes(Year, y=pop.est, ymin=lower, ymax=upper))
pop.size.plot+geom_pointrange(size=1) + labs(y="Population Size")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))+ylim(0,3100)

recruits<-data.frame(cbind(Year=c(2015:2018),pop.est=jc1.stats[2:5, 1], lower=jc1.quants[2:5, 1], upper=jc1.quants[2:5, 5]))
recruits.plot<-ggplot(recruits, aes(Year, y=pop.est, ymin=lower, ymax=upper))
recruits.plot+geom_pointrange(size=1) + labs(y="Number of Recruits")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))+ylim(0,1000)


rm(list=ls())

# Set working directory
basedirectory <- "C:/Users/solit/Documents/GitHub/Sicklefin-Redhorse-Analysis"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "lubridate", "tidyverse", "rjags","jagsUI", "coda","parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)



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
Sex<-ifelse(ms.cap.hist$Sex==1, 0, 1) #0 is female, 1 is male

js.CH<-ms.cap.hist
js.CH$Sex<-NULL

# Convert CJS capture history to data for JS model
js.CH[is.na(js.CH)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(js.CH)[1]), js.CH)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL
head(CH.du)
nz<-2200 # augmented individuals

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



inits <- function() {list(z=cbind(rep(NA, nrow(z.init)),z.init[,-1]), mean.phi=runif(1, 0, 1), p0seine=runif(1, 0, 1), p0fyke=runif(1, 0, 1), p0pit=runif(1, 0, 0.01), gamma=runif(5, 0, 1))} #ER=runif(1, 0, 1))} #Sex=c(rep(NA,length(Sex)), rbinom(nz ,1, 0.5)),

# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "mean.phi", "B", "gamma", "Nsuper", "N", "psi") #, "alpha", "beta") 
#codaOnly<-c("Nsuper", "N", "B", "psi")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 3

## Run models

#ER=200, nz=1000, nb=100, nt=1, ni=300
# sr.ms.js.jm10 <- jags.model(data=dat3, inits = inits,
#                             file = "ms_js_phiSex_gam0_pGearEffort.jags",
#                             n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.ms.js.jc10 <- coda.samples(sr.ms.js.jm10, params, n.iter=ni)

# ER=200, nz=1000, nb=100, nt=1, ni=300
# 
# sr.ms.js.jm11 <- jags.model(data=dat3, inits = inits,
#                             file = "srh_js_phi0_gam0_pGearEffort.jags",
#                             n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.ms.js.jc11 <- coda.samples(sr.ms.js.jm11, params, n.iter=ni)


# ER=500, nz=1100, nb=100, nt=1, ni=1000

# sr.ms.js.jm12 <- jags.model(data=dat3, inits = inits,
#                             file = "srh_js_phi0_gam0_pGearEffort.jags",
#                             n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.ms.js.jc12 <- coda.samples(sr.ms.js.jm12, params, n.iter=ni)

# ER=500, nz=1800, nb=100, nt=1, ni=500, nc=1
ptm <- proc.time()
sr.ms.js.jm13 <- jags.model(data=dat3, inits = inits,
                            file = "srh_js_phi0_gamTime_pGearEffort.jags",
                            n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc13 <- coda.samples(sr.ms.js.jm13, params, n.iter=ni)

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


# load("SRH_js_phi0_gam0_pGearEffort_7000.gzip")
# plot(sr.ms.js.jc13, ask=TRUE)
# 
# out.1<-sr.ms.js.jc13[[1]]
# out.2<-sr.ms.js.jc13[[1]]
# out.3<-sr.ms.js.jc13[[1]]
# out.4<-sr.ms.js.jc13[[1]]
# out.5<-sr.ms.js.jc13[[1]]
# out.6<-sr.ms.js.jc13[[1]]
# out.all<-rbind(out.1,out.2,out.3,out.4,out.5, out.6)
# 
# hist(out.all[,24], breaks=seq(0,1,0.01))
# plot(out.all[,6],type="l")
# mean<-quants2.5<-quants97.5<-rep(NA,ncol(out.all))
# 
# for(i in 1:ncol(out.all)){
#   mean[i]<-mean(out.all[,i])
#   quants2.5[i]<-quantile(out.all[,i],probs=0.025)
#   quants97.5[i]<-quantile(out.all[,i],probs=0.975)
# }
# 
# estimates<-data.frame(cbind(Param=colnames(out.all), Mean=mean, Lower=quants2.5, Upper=quants97.5))
# estimates$Mean<-as.numeric(as.character(estimates$Mean))
# estimates$Lower<-as.numeric(as.character(estimates$Lower))
# estimates$Upper<-as.numeric(as.character(estimates$Upper))
# 
# pop.size<-cbind(Year=c(2015:2018),estimates[3:6,])
# pop.size.plot<-ggplot(pop.size, aes(Year, y=Mean, ymin=Lower, ymax=Upper))
# pop.size.plot+ coord_cartesian(ylim = c(0, 1500))+ geom_pointrange(size=1) + labs(y="Population Size")+theme(plot.margin=margin(.75,.75,.75,.75,"cm"), axis.text=element_text(size=20), axis.title=element_text(size=25, face="bold"), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)))
# 
# est.rec<-cbind(Recruits="Recruits", estimates[1,])
# est.rec.plot<-ggplot(est.rec, aes(Recruits, y=Mean, ymin=Lower, ymax=Upper))
# est.rec.plot+coord_cartesian(ylim = c(0, 450)) + geom_pointrange(size=1.25) + labs(y="Expected Annual Recruits")+theme(plot.margin=margin(.75,.75,.75,.75,"cm"), axis.text=element_text(size=20), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=25, face="bold"), axis.title.x=element_blank(), axis.text.x=element_blank())+scale_y_continuous(expand = c(0, 0))
# 
# 
# est.surv<-cbind(Survival="Survival", estimates[13,])
# 
# est.surv.plot<-ggplot(est.surv, aes(Survival, y=Mean, ymin=Lower, ymax=Upper))
# est.surv.plot+geom_pointrange(size=1.15) + coord_cartesian(ylim = c(0, 1)) + labs(y="Estimated Annual Survival")+theme(plot.margin=margin(.75,.75,.75,.75,"cm"), axis.text=element_text(size=20), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0), size=25, face="bold"), axis.title.x=element_blank(), axis.text.x=element_blank())+scale_y_continuous(expand = c(0, 0))
# 
# 
# det.prob<-estimates[c(16:19, 22:25,28:31),]
# 
# for(i in 1:nrow(det.prob)){
#   for(j in 2:ncol(det.prob)){
#     det.prob[i,j]<-ifelse(det.prob[i,j]==0, NA, det.prob[i,j])
#   }
# }
# det.prob<-det.prob[,2:4]
# 
# det.prob<-cbind(det.prob, Gear=c(rep("Fyke",4), rep("Antenna", 4), rep("Seine", 4)), Year=c(rep(2015:2018, 3)))
# 
# det.prob.plot<-ggplot(det.prob, aes(y=Mean, x=Year, group=Gear, ymin=Lower, ymax=Upper))
# 
# det.prob.plot+ geom_pointrange(aes(colour=Gear), position=position_dodge(width=c(0.2, 0)), size = 1, alpha = 0.7)+
#   labs(y="Detection Probability")+
#   theme(plot.margin=margin(.75,.75,.75,.75,"cm"), axis.text=element_text(size=20), axis.title=element_text(size=25, face="bold"), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)), legend.title=element_text(size=16),legend.text=element_text(size=15))+ 
#   expand_limits(y=c(0,1))+
#   scale_y_continuous(expand = c(0, 0))
# 
# 


## 3000 augmented individuals, 100000 iterations
# sr.ms.js.jm1 <- jags.model(data=dat3, inits = inits,
#                          file = "ms_js_phiSex_pGearEffort.jags",
#                          n.chains = nc, n.adapt = 100, quiet = F)
# 
# sr.ms.js.jc1 <- coda.samples(sr.ms.js.jm1, params, n.iter=300)

## 3500 augmented individuals, 50000 iterations
# sr.ms.js.jm2 <- jags.model(data=dat3, inits = inits,
#                            file = "ms_js_phiSex_pGearEffort.jags",
#                            n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.ms.js.jc2 <- coda.samples(sr.ms.js.jm2, params, n.iter=ni)

## 4000 augmented individuals, 50000 iterations
# sr.ms.js.jm3<-jags.model(data=dat3, inits = inits,
#                          file = "ms_js_phiSex_pGearEffort.jags",
#                          n.chains = nc, n.adapt = nb, quiet = F)
# sr.ms.js.jc3 <- coda.samples(sr.ms.js.jm3, params, n.iter=ni)








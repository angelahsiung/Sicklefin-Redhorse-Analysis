rm(list=ls())

############################################################
#### Potential additions/modifications
#### 1. use robust design to more accurately estimate
#### detection probability (account for effort based on gear)
#### 2. expand multi-state model to add spatial component
#### (how fishes move between different sections of the river)
#### 3. Add environmental data (e.g. stream temps)
#### 4. Incorporate antenna functional periods
############################################################

###################################
## Operational years (occasion)
## Seine: 2014-2016 (1-3)
## Fyke: 2016-2018 (3-5)
## PIT: 2017-2018 (4-5)
###################################

# Set working directory
basedirectory <- "C:/Users/solit/Desktop/Dissertation/Sicklefin"
setwd(basedirectory)

# Prepare packages
list.of.packages <- c("stringr", "reshape2", "RMark","lubridate", "tidyverse", "rjags")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

# Read in data
dat<-read.csv("Sicklefin_Capture_Data.csv")

dat<-dat[dat$Sex!="Unknown",] # remove rows where Sex is "Unknown"
dat$CollectDate<-strptime(dat$CollectDate, format="%d-%b-%y") # convert date format
dat$CollectYear.1<-ifelse(is.na(dat$CollectDate), dat$CollectYear, year(dat$CollectDate)) # create Year column

seine.dat<-dat[dat$PriMethod=="Seining",]
with(seine.dat, table(CollectYear.1, CollectDate))


##############################
##### summary statistics #####
##############################
length(unique(dat$Individual_ID)) # number of individual fish in the data set
nrow(dat) # total captures
nrow(dat[dat$Sex=="Male",])
nrow(dat[dat$Sex=="Female",])
sex.table<-data.frame(with(dat, table(Individual_ID, Sex)))
sex.table<-sex.table[sex.table$Freq>0,]

nrow(sex.table[sex.table$Sex=="Male",]) # number of males captured
nrow(sex.table[sex.table$Sex=="Female",]) # number of females captured
with(dat, table(dat$PriMethod)) # number of captures by each method
with(dat, table(dat$CollectYear.1)) # number of captures by year
with(dat, table(CollectYear.1, PriMethod)) # number of captures by year and method
ind.gear.year<-with(dat, table(dat$Individual_ID, PriMethod, CollectYear.1))

## number of captures per year
caps.by.year<-with(dat, table(Individual_ID, CollectYear.1))
caps.by.year<-data.frame(ifelse(caps.by.year>1, 1, caps.by.year))
caps.by.year<-colSums(caps.by.year)


## captures by fyke
fyke.2016<-ind.gear.year[,"Trapping (Weir or Fyke)", "2016"]
fyke.2017<-ind.gear.year[,"Trapping (Weir or Fyke)", "2017"]
fyke.2018<-ind.gear.year[,"Trapping (Weir or Fyke)", "2018"]
fyke.caps.2016<-length(fyke.2016[fyke.2016>0]) # number of individuals captured by fyke in 2016
fyke.caps.2017<-length(fyke.2017[fyke.2017>0]) # number of individuals captured by fyke in 2017
fyke.caps.2018<-length(fyke.2018[fyke.2018>0]) # number of individuals captured by fyke in 2018

## captures by seine
seine.2014<-ind.gear.year[,"Seining", "2014"]
seine.2015<-ind.gear.year[,"Seining", "2015"]
seine.2016<-ind.gear.year[,"Seining", "2016"]
seine.caps.2014<-length(seine.2014[seine.2014>0]) # number of individuals captured by seine in 2014
seine.caps.2015<-length(seine.2015[seine.2015>0]) # number of individuals captured by seine in 2015
seine.caps.2016<-length(seine.2016[seine.2016>0]) # number of individuals captured by seine in 2016

## captures by antenna
PIT.2017<-ind.gear.year[,"PIT Array", "2017"]
PIT.2018<-ind.gear.year[,"PIT Array", "2018"]
PIT.caps.2017<-length(PIT.2017[PIT.2017>0]) # number of individuals captured by PIT in 2017
PIT.caps.2018<-length(PIT.2018[PIT.2018>0]) # number of individuals captured by PIT in 2018

###########################
##### Data formatting #####
###########################

# organize data into pivot table
ms.dat<-dcast(dat, dat$Individual_ID+dat$CollectYear.1~dat$PriMethod)
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

#ms.cap.hist<-ms.cap.hist[,2:6]
ms.cap.hist<-read.csv("Sicklefin Capture History.csv")[-1]
# Creating capture history for multi-state analysis
get.first<-function(x) min(which(!is.na(x)))
f<-apply(ms.cap.hist, 1, get.first)

for (i in 1:nrow(ms.cap.hist)){
  for (j in f[i]:ncol(ms.cap.hist)){
    ms.cap.hist[i,j]<-ifelse(is.na(ms.cap.hist[i,j]), 1, ifelse(ms.cap.hist[i,j]=="S", 2, ifelse(ms.cap.hist[i,j]=="T", 3, ifelse(ms.cap.hist[i,j]=="P", 4,ifelse(ms.cap.hist[i,j]=="PT", 5, ifelse(ms.cap.hist[i,j]=="ST", 6, ms.cap.hist[i,j]))))))
  }
}

ms.cap.hist<-as.matrix(ms.cap.hist) # Capture history
# write.csv(ms.cap.hist, "Sicklefin Capture History.csv")

# Extracting sex covariate
ms.dat.sex<-dcast(dat, dat$Individual_ID~dat$Sex)
for (i in 1:nrow(ms.dat.sex)){
  ms.dat.sex$Female[i]<-ifelse(ms.dat.sex$Female[i]!=0, "F", NA)
  ms.dat.sex$Male[i]<-ifelse(ms.dat.sex$Male[i]!=0, "M", NA)
  #ms.dat.sex$Unknown[i]<-ifelse(ms.dat.sex$Unknown[i]!=0, "U", NA)
}
ms.dat.sex$Sex <- coalesce(ms.dat.sex$Female, ms.dat.sex$Male)
ms.dat.sex$Sex <- ifelse(ms.dat.sex$Sex=="F", 1, 2) 

### Calculating new and recaps 
Recap<-NewCap<-rep(NA, 4) # recaps and new caps for 2015-2018

for(n in 1:4){
  NewCap[n]<-length(f[f==(n+1)]) # number of newly captured individuals 2015-2018
  Recap[n]<-nrow(ms.cap.hist[!is.na(ms.cap.hist[,n+1])&ms.cap.hist[,n+1]>1&!is.na(ms.cap.hist[,n]),]) # number of recaptured individuals
} 


#################################
##### Multi-state CJS Model #####
#################################
nind<-dim(ms.cap.hist)[1]
n.occ<-dim(ms.cap.hist)[2]

sink("sicklefin_redhorse_ms_cjs.jags")
cat("
    model{
    #phi ~ dunif(0,1) # Annual survival prob
    alpha0 ~ dnorm(0,0.37)
    alpha1 ~ dnorm(0,0.37)

    for(t in 1:2){
      pFyke[t]<-0 # det prob zero during years when fyke was not used
    }                                        
    for(t in 3:n.occ){
      p.d[(t-2)] ~ dunif(0,1) # Pr(catch in fyke net)      #### FIRST 3 PLACES IN p.d RESERVED FOR RANDOM pFyke VALUES
      pFyke[t] <- p.d[(t-2)]
    }
    
    for(t in 1:3){
      p.d[(t+3)] ~ dunif(0,1)                              #### SLOTS 4-6 IN p.d RESERVED FOR RANDOM pSeine VALUES
      pSeine[t] <- p.d[(t+3)]
    }
    for(t in 4:n.occ){
      pSeine[t] <-0 # det prob zero during years when seine was not used
    }
    
    for(t in 1:3){
      pPit[t]<-0 # det prob zero during years when PIT antenna was not set up
    }
    for(t in 4:n.occ){
      p.d[(t+3)] ~ dunif(0,1)                              #### SLOTS 7-8 IN p.d RESERVED FOR RANDOM pPit VALUES
      pPit[t] <- p.d[(t+3)]
    }

    phi.male <- 1/(1+exp(-(alpha0+alpha1)))   # back-transformed male survival rates
    phi.female <- 1/(1+exp(-(alpha0)))   # back-transformed female survival rates
    
    for (i in 1:nind){
    logit(phi[i]) <- alpha0 + alpha1*Sex[i]
    }

 
    for(t in 1:n.occ){
    cap.probs[1,t,1] <- (1-pSeine[t])*(1-pFyke[t])*(1-pPit[t]) # Pr(not detected|alive)
    cap.probs[1,t,2] <- pSeine[t]*(1-pFyke[t])*(1-pPit[t]) # Pr(detected in seine, but not by other gear)
    cap.probs[1,t,3] <- pFyke[t]*(1-pSeine[t])*(1-pPit[t])
    cap.probs[1,t,4] <- pPit[t]*(1-pFyke[t])*(1-pSeine[t])
    cap.probs[1,t,5] <- pFyke[t]*pPit[t]*(1-pSeine[t])
    cap.probs[1,t,6] <- pSeine[t]*pFyke[t]*(1-pPit[t])
    cap.probs[1,t,7] <- pSeine[t]*pPit[t]*(1-pFyke[t])
    cap.probs[1,t,8] <- pSeine[t]*pFyke[t]*pPit[t]
    cap.probs[2,t,1] <- 1 # If dead, can't be detected
    cap.probs[2,t,2] <- 0
    cap.probs[2,t,3] <- 0
    cap.probs[2,t,4] <- 0
    cap.probs[2,t,5] <- 0
    cap.probs[2,t,6] <- 0
    cap.probs[2,t,7] <- 0
    cap.probs[2,t,8] <- 0
    }

    for (i in 1:nind){
    z[i, f[i]] <- 1 # Known to be alive in first capture year
    for (t in (f[i]+1):n.occ){
      z[i,t] ~ dbern(z[i,t-1]*phi[i]) # Actual alive/dead state (partially observed)
      state[i,t] <- ifelse(z[i,t]<1, 2, 1)
      y[i,t] ~ dcat(cap.probs[state[i,t],t,1:8])
    }
    }
    
    for(t in 2:n.occ){
      PrNewCap[t] <- 1 - (1-pFyke[t])*(1-pSeine[t])
      PrRecap[t] <- 1 - (1-pFyke[t])*(1-pSeine[t])*(1-pPit[t]) ## wouldn't this also include new caps since this means prob of being captured by any gear whether it's new or recap?
      N[t] <- NewCaps[t-1]/PrNewCap[t] + Recaps[t-1]/PrRecap[t]
    }
    }
    
    ",fill = TRUE)
sink()


known.state.ms<-function(ms, notseen){
  # notsee: label for 'not seen"
  state<-ms
  state[state==notseen]<-NA
  for (i in 1:nind){
    m<-min(which(!is.na(state[i,])))
    state[i,m]<-NA
  }
  return(state)
}

# initial values
zi <- matrix(NA, nrow(ms.cap.hist), ncol(ms.cap.hist))
for(i in 1:nrow(ms.cap.hist)) {
  if(f[i]>4)
    next
  zi[i,(f[i]+1):5] <- 1
}


inits <- function() list(z=zi, 
                         p.d = runif(8, 0, 1),
                         alpha0 = runif(1, -10, 10), alpha1 = runif(1, -10, 10))


# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "alpha0", "alpha1", "phi.male", "phi.female", "N",  "PrNewCap", "PrRecap") #",pdot")

# Data
in.data <- list(y = ms.cap.hist, f = f, nind = dim(ms.cap.hist)[1], n.occ = dim(ms.cap.hist)[2], z = known.state.ms(ms.cap.hist, 1), Sex=ms.dat.sex$Sex, NewCaps=NewCap, Recaps=Recap)

# MCMC settings
ni <- 5000
nt <- 1
nb <- 100
nc <- 3


in.data$y <- data.matrix(in.data$y)
zdat <- data.matrix(in.data$z)
in.data$z <- NULL


sr.ms.jm1<-jags.model(data=in.data, inits = inits, file = "sicklefin_redhorse_ms_cjs.jags",
                      n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.jc1 <- coda.samples(sr.ms.jm1, params, n.iter=ni)


plot(sr.ms.jc1, ask=TRUE)


save(sr.ms.jc1, file="SRH_MultiState_Results.gzip") # Save MCMC results

# Model results and plots
out.jc1<-summary(sr.ms.jc1) # model output

jc1.stats<-out.jc1$statistics # parameter estimates
jc1.quants<-out.jc1$quantiles # 95% confidence intervals

### detection probability by gear
p.table<-data.frame(rbind(jc1.stats[3:7,1], jc1.stats[8:12, 1], jc1.stats[13:17, 1]))
p.table<-rbind(p.table, colSums(p.table))
colnames(p.table)<-c(2014:2018)
rownames(p.table)<-c("Fyke", "PIT Array", "Seine", "Any")
write.csv(p.table, "SRH_Det_Probs_by_Gear.csv")

### Plotting male and female survival rates over time
phi.sex<-data.frame(cbind(Sex=c("Female", "Male"),phi.est=round(jc1.stats[26:27, 1],2), lower=round(jc1.quants[26:27, 1],2), upper=round(jc1.quants[26:27, 5],2)))

phi.sex.plot<-ggplot(data=phi.sex, aes(Sex, y=phi.est))
phi.sex.plot+geom_pointrange(aes(ymin=lower, ymax=upper))+labs(y="Survival Probability")+theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold"))+expand_limits(y=c(0,1))

### Plotting detection probability by gear
det.prob<-cbind(jc1.stats[7:21, 1], jc1.quants[7:21, 1], jc1.quants[7:21, 5])

par(mar=c(5,5,5,5))
plot(c(NA, jc1.stats[paste0("pFyke", "[", c(3:5),"]"),1]), ylim=c(0,1), ylab="Detection probability", xlab=NA, xaxt='n', pch=16, cex.lab=2, cex.axis=1.5)
segments(2:4, jc1.quants[paste0("pFyke", "[", c(3:5),"]"),1], 2:4, jc1.quants[paste0("pFyke", "[", c(3:5),"]"),5], lwd=2)
axis(1, at=1:4, labels=c("2015", "2016", "2017", "2018"), cex.axis=1.5)
points(c(rep(NA, 2),jc1.stats[paste0("pPit", "[", c(4:5),"]"),1]), col="red", pch=16)
segments(3:4, jc1.quants[paste0("pPit", "[", c(4:5),"]"),1], 3:4, jc1.quants[paste0("pPit", "[", c(4:5),"]"),5], col="red", lwd=2)
points(x=1:4+0.02, c(jc1.stats[paste0("pSeine", "[", c(2:3),"]"),1], NA, NA), col="blue", pch=16)
segments(1:2+0.02, jc1.quants[paste0("pSeine", "[", c(2:3),"]"),1], 1:2+0.02, jc1.quants[paste0("pSeine", "[", c(2:3),"]"),5], col="blue", lwd=2)
legend(1, 1, c("Fyke", "Antenna", "Seine"), col=c("black", "red", "blue"), pch=16, cex=1.5)

### detection probability by any gear
overall.det.prob<-data.frame(cbind(Year=2015:2018, det.prob=jc1.stats[paste0("pdot", "[", 2:5,"]"),1], lower=jc1.quants[paste0("pdot", "[", 2:5,"]"),1], upper=jc1.quants[paste0("pdot", "[", 2:5,"]"),5]))
overall.p.plot<-ggplot(overall.det.prob, aes(Year, y=det.prob, ymin=lower, ymax=upper))
overall.p.plot+geom_pointrange() + labs(y="Detection Probability")+theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold"))

### Plotting population size
pop.size<-data.frame(cbind(Year=c(2015:2018),pop.est=jc1.stats[1:4, 1], lower=jc1.quants[1:4, 1], upper=jc1.quants[1:4, 5]))
pop.size.plot<-ggplot(pop.size, aes(Year, y=pop.est, ymin=lower, ymax=upper))
pop.size.plot+geom_pointrange() + labs(y="Population Size")+theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold"))


####################################
##### MS CJS model with effort #####
####################################

## Effort
effort<-read.csv("SRH_Effort_v1.csv", header=T)
fyke.eff<-tapply(effort$UB_FYKE, effort$Year, sum)
pit.eff<-aggregate(effort[,4:7], by=list(effort$Year), sum)
rowSums(pit.eff[,2:5])
## Data
dat2 <- in.data
#dat2$y <- data.matrix(caphist.in)

dat2$fykeEffort <- c(0, 0, 2.5, 7, 4)
dat2$pitEffort <- c(0, 0, 0, 1, 5.25)

str(dat2)

# Initial values
inits2 <- function() list(z = zi, phi = runif(2, 0.5, 1))

# Parameters monitored
params2 <- c("pFyke", "pSeine", "pPit", "N", "phi", "PrNewCap", "PrRecap")


sr.ms.jm2 <- jags.model(data=dat2, inits = inits2, file = "mscjs_phiSex_pGearEffort.jag",
                        n.chains = 3, n.adapt = 1000, quiet = F)

sr.ms.jc2 <- coda.samples(sr.ms.jm2, params2, n.iter=10000)

plot(sr.ms.jc2, ask=TRUE)

(ss2 <- summary(sr.ms.jc2))

plot(1:4, ss2$quant[1:4,3], ylim=c(0, 2000))
segments(1:4, ss2$quant[1:4,1], 1:4, ss2$quant[1:4,5])




# Model results and plots
ss2.stats<-ss2$statistics # parameter estimates
ss2.quants<-ss2$quantiles # 95% confidence intervals

### detection probability by gear
p.table<-data.frame(rbind(ss2.stats[13:17,1], ss2.stats[18:22, 1], ss2.stats[23:27, 1]))
colnames(p.table)<-c(2014:2018)
rownames(p.table)<-c("Fyke", "PIT Array", "Seine")
# write.csv(p.table, "SRH_Det_Probs_by_Gear.csv")

### Plotting male and female survival rates over time
phi.sex<-data.frame(cbind(Sex=c("Female", "Male"),phi.est=round(ss2.stats[28:29, 1],2), lower=round(ss2.quants[28:29, 1],2), upper=round(ss2.quants[28:29, 5],2)))

phi.sex.plot<-ggplot(phi.sex, aes(x=Sex, y=phi.est, ymin=lower, ymax=upper))
phi.sex.plot+geom_pointrange(size=1)+labs(y="Survival Probability")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))+expand_limits(y=c(0,1))

### Plotting detection probability by gear
det.prob<-data.frame(cbind(ss2.stats[c(14:17, 19:22,24:27), 1], ss2.quants[c(14:17, 19:22,24:27), 1], ss2.quants[c(14:17, 19:22,24:27), 5]))
for(i in 1:nrow(det.prob)){
  for(j in 1:ncol(det.prob)){
    det.prob[i,j]<-ifelse(det.prob[i,j]==0, NA, det.prob[i,j])
  }
}
row.names(det.prob)<-NULL
det.prob<-cbind(det.prob, c(rep("Fyke",4), rep("Antenna", 4), rep("Seine", 4)), c(rep(2015:2018, 3)))
colnames(det.prob)<-c("Estimate", "Lower", "Upper", "Gear", "Year")


det.prob.plot<-ggplot(det.prob, aes(y=Estimate, x=Year, group=Gear, ymin=Lower, ymax=Upper))
det.prob.plot+geom_pointrange(aes(colour=Gear), position=position_dodge(width=c(0.2, 0)), size = 1, alpha = 0.7)+labs(y="Detection Probability")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))+expand_limits(y=c(0,1))+scale_y_continuous(expand = c(0, 0))+theme(plot.margin=margin(1,1,1,1,"cm"), legend.title=element_text(size=16), legend.text=element_text(size=15))


# par(mar=c(5,5,5,5))
# plot(c(NA, ss2.stats[paste0("pFyke", "[", c(3:5),"]"),1]), ylim=c(0,1), ylab="Detection probability", xlab=NA, xaxt='n', pch=16, cex.lab=2, cex.axis=1.5)
# segments(2:4, ss2.quants[paste0("pFyke", "[", c(3:5),"]"),1], 2:4, ss2.quants[paste0("pFyke", "[", c(3:5),"]"),5], lwd=2)
# axis(1, at=1:4, labels=c("2015", "2016", "2017", "2018"), cex.axis=1.5)
# points(c(rep(NA, 2),ss2.stats[paste0("pPit", "[", c(4:5),"]"),1]), col="red", pch=16)
# segments(3:4, ss2.quants[paste0("pPit", "[", c(4:5),"]"),1], 3:4, ss2.quants[paste0("pPit", "[", c(4:5),"]"),5], col="red", lwd=2)
# points(x=1:4+0.02, c(ss2.stats[paste0("pSeine", "[", c(2:3),"]"),1], NA, NA), col="blue", pch=16)
# segments(1:2+0.02, ss2.quants[paste0("pSeine", "[", c(2:3),"]"),1], 1:2+0.02, ss2.quants[paste0("pSeine", "[", c(2:3),"]"),5], col="blue", lwd=2)
# legend(1, 1, c("Fyke", "Antenna", "Seine"), col=c("black", "red", "blue"), pch=16, cex=1.5)

### New capture probability
new.cap.prob<-data.frame(cbind(Year=2015:2018, det.prob=ss2.stats[paste0("PrNewCap", "[", 2:5,"]"),1], lower=ss2.quants[paste0("PrNewCap", "[", 2:5,"]"),1], upper=ss2.quants[paste0("PrNewCap", "[", 2:5,"]"),5]))
overall.p.plot<-ggplot(new.cap.prob, aes(Year, y=det.prob, ymin=lower, ymax=upper))
overall.p.plot+geom_pointrange() + labs(y="Detection Probability")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))+expand_limits(y=c(0,1))
### Recapture probability
overall.det.prob<-data.frame(cbind(Year=2015:2018, det.prob=ss2.stats[paste0("PrRecap", "[", 2:5,"]"),1], lower=ss2.quants[paste0("PrRecap", "[", 2:5,"]"),1], upper=ss2.quants[paste0("PrRecap", "[", 2:5,"]"),5]))
overall.p.plot<-ggplot(overall.det.prob, aes(Year, y=det.prob, ymin=lower, ymax=upper))
overall.p.plot+geom_pointrange(size=1) + labs(y="Recapture Probability")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))

### Plotting population size
pop.size<-data.frame(cbind(Year=c(2015:2018),pop.est=ss2.stats[1:4, 1], lower=ss2.quants[1:4, 1], upper=ss2.quants[1:4, 5]))
pop.size.plot<-ggplot(pop.size, aes(Year, y=pop.est, ymin=lower, ymax=upper))
pop.size.plot+geom_pointrange(size=1) + labs(y="Population Size")+theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"))



#######################################
####### Multistate J-S Model ##########
#######################################

## Convert capture data for CJS model to data for JS model
ms.cap.hist[is.na(ms.cap.hist)]<-1 # if NA, assign 1 for "not seen"

CH.du<-as.matrix(cbind(rep(1, dim(ms.cap.hist)[1]), ms.cap.hist)) # add extra (dummy) sampling occ in the beginning (Kerry and Schaub 2012)
colnames(CH.du)<-NULL

nz<-1000 # augmented individuals

ms.js.CH.aug<-rbind(CH.du, matrix(1, ncol=dim(CH.du)[2], nrow=nz)) # data agumentation


dat3<-list(y = ms.js.CH.aug, n.occ = dim(ms.js.CH.aug)[2], Sex=c(ms.dat.sex$Sex, rep(1, nz/2), rep(2, nz/2)), M=dim(ms.js.CH.aug)[1]) #,NewCaps=NewCap, Recaps=Recap)
dat3$fykeEffort <- c(0, 0, 0, 2.5, 7, 4)
dat3$pitEffort <- c(0, 0, 0, 0, 1, 5.25)

zi<-ms.js.CH.aug
for(i in 1:nrow(zi)){
  for(j in 1:ncol(zi)){
    zi[i,j]<-ifelse(zi[i,j]>1, 2, zi[i,j])
  }
}
inits <- function() {list(z=cbind(rep(NA, dim(ms.js.CH.aug)[1]), zi[, -1]), phi= runif(2, 0.5, 1))}
                         
# Parameters monitored
params <- c("pFyke", "pSeine", "pPit", "phi", "gamma") #, "Nsuper", "N", "B", "psi")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 1000
nc <- 3



sr.ms.js.jm1<-jags.model(data=dat3, inits = inits, file = "ms_js_phiSex_pGearEffort.jags",
                         n.chains = nc, n.adapt = nb, quiet = F)

sr.ms.js.jc1 <- coda.samples(sr.ms.jm1, params, n.iter=ni)


plot(sr.ms.jc1, ask=TRUE)













#######################################################################################################



########################################
####### Jolly-Seber Robust Design ######
########################################
# robust.dat<-dcast(cp.dat, cp.dat$Individual_ID~CollectDate, sum)
# 
# robust.dat<-robust.dat[,2:ncol(robust.dat)]
# 
# robust.dat[robust.dat>0]<-1
# 
# n.occasions<-3
# 
# nz<-1000
# 
# aug<-matrix(0, nz, ncol(robust.dat))
# 
# CH.aug<-rbind(as.matrix(robust.dat), aug)
# 
# nss<-c(4,7,7)
# cnss<-c(0, 4, 11)
# 
# 
# sink("sicklefin_redhorse_js_robust.jags")
# cat("
#     model
#     {
#     ## Jolly-Seber - restricted dynamic occupancy formulation
#     
#     ## Priors & constraints
#     for (i in 1:M) {
#     for (t in 1:(n.occasions-1)) {
#     phi[i, t] <- mean.phi 
#     }
#     for (t in 1:n.occasions) {
#     p[i, t] <- mean.p 
#     }
#     }
#     mean.phi ~ dunif(0, 1)        ## mean.phi is not updated if phi.vary=1!
#     mean.p ~ dunif(0, 1)          ## mean.p is not updated if p.vary=1!
#     
#     for (t in 1:n.occasions) {
#     gamma[t] ~ dunif(0, 1)
#     }
#     
#     
#     
#     ## Likelihood
#     for (i in 1:M){
#     z[i, 1] ~ dbern(gamma[1])        # First occasion, state process
#     
#     # First occasion, observation process
#     for(j in 1:nss[1]){
#     mu1[i,j]<-z[i,1]*p[i,j]
#     y[i,j]<-dbern(mu1[i,j])
#     }
#     
#     # Subsequent occasions
#     for (t in 2:n.occasions) {
#     # State process
#     q[i, t-1] <- 1-z[i, t-1]									              # Availability for recruitment
#     mu2a[i, t] <- phi[i, t-1] * z[i, t-1]
#     mu2b[i, t] <- gamma[t] * prod(q[i,1:(t-1)])				      # Prob entering pop given available
#     mu2[i, t] <- mu2a[i, t] + mu2b[i, t]
#     z[i, t] ~ dbern(mu2[i, t])
#     # Observation process
#     for(j in 1:nss[t]){
#     mu3[i,(cnss[t]+j)]<-z[i,t]*p[i, (cnss[t]+j)]
#     y[i,(cnss[t]+j)]~dbern(mu3[i, (cnss[t]+j)])
#     }
#     }
#     }
#     
#     ## Calculate derived population parameters
#     #   Entry probability
#     for (t in 1:n.occasions) {  qgamma[t] <- 1-gamma[t]  }
#     cprob[1] <- gamma[1]
#     for (t in 2:n.occasions) {  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])   }
#     psi <- sum(cprob[])                                    # Inclusion probability
#     for (t in 1:n.occasions) {  b[t] <- cprob[t] / psi  }  # Entry probability
#     
#     #   Population size and recruitment
#     for (i in 1:M) {
#     recruit[i, 1] <- z[i, 1]
#     for (t in 2:n.occasions) {  recruit[i, t] <- (1-z[i, t-1]) * z[i, t]  }
#     }
#     for (t in 1:n.occasions) {
#     N[t] <- sum(z[1:M, t])						# Actual population size
#     B[t] <- sum(recruit[1:M, t])			# Number of entries
#     f[t] <- B[t]/N[t]                 # Per capita entry probability
#     }
#     
#     for (t in 2:n.occasions) { lamda[t] <- N[t]/N[t-1]}  # Population growth rate
#     
#     #   Superpopulation size
#     for (i in 1:M) {
#     Nind[i] <- sum(z[i, 1:n.occasions])
#     Nalive[i] <- 1-equals(Nind[i], 0)
#     }
#     Nsuper <- sum(Nalive[])
#     }
#     ",fill = TRUE)
# sink()
# 
# z.init <- CH.aug
# z.init[z.init==0] <- 1
# 
# in.data <- list(y = CH.aug, n.occasions = n.occasions, M = dim(CH.aug)[1], nss=nss, cnss=cnss)
# 
# # Initial values
# inits <- function() list(z = z.init, mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))
# 
# # Parameters monitored
# params <- c("mean.phi", "mean.p", "time.phi", "time.p", "N", "B", "Nsuper", "psi", "gamma", "b")
# 
# # MCMC settings
# ni <- 5000
# nt <- 1
# nb <- 2000
# nc <- 3 
# 
# 
# # Burn in the model
# sr.robust.jm1<-jags.model(data=in.data, inits = inits, file = "sicklefin_redhorse_js_robust.jags",
#                           n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.robust.jc1<-coda.samples(model=sr.robust.jm1, variable.names=params, n.iter=ni, thin=nt)

########################################
####### Closed-population model ########
########################################

# Formatting data
# cp.dat<-dat[dat$PriMethod=="Trapping (Weir or Fyke)",] # Extract Fyke captures
# with(cp.dat, table(CollectDate, CollectYear.1))
# # Capture history function
# closed.pop.CH<-function(dat, year){
#   CH<-dat[dat$CollectYear.1==year,]
#   CH<-dcast(CH, CH$Individual_ID~CH$CollectDate)
#   CH[is.na(CH)]<-0
#   CH[CH==year]<-1
#   CH<-as.matrix(CH[, 2:ncol(CH)])
#   colnames(CH)<-NULL
#   return(CH)
# }
# 
# # apply cap hist function to each year
# cp.CH.2016<-closed.pop.CH(cp.dat, 2016)
# cp.CH.2017<-closed.pop.CH(cp.dat, 2017)
# cp.CH.2018<-closed.pop.CH(cp.dat, 2018)
# 
# rowSums(cp.CH.2016)
# rowSums(cp.CH.2017)
# rowSums(cp.CH.2018)
# 
# # Data augmentation
# nz<-1000
# 
# yaug<-rbind(cp.cap.hist, matrix(0, nz, ncol(cp.cap.hist)))
# yaug.2016<-rbind(as.matrix(cp.CH.2016), matrix(0, nz, ncol(cp.CH.2016)))
# yaug.2017<-rbind(as.matrix(cp.CH.2017), matrix(0, nz, ncol(cp.CH.2017)))
# yaug.2018<-rbind(as.matrix(cp.CH.2018), matrix(0, nz, ncol(cp.CH.2018)))
# 
# # Closed-population model
# sink("sicklefin_redhorse_closed_pop.jags")
# cat("
#     model{
#     #priors
#     omega ~ dunif(0, 1)
#     p ~dunif(0, 1)
#     
#     #likelihood
#     for (i in 1:M){
#     z[i] ~ dbern(omega)
#     for(j in 1:T){
#     yaug[i,j] ~ dbern(p.eff[i,j])
#     p.eff[i,j] <- z[i]*p
#     }
#     }
#     N <- sum(z[])
#     }
#     ", fill=TRUE)
# sink()
# 
# 
# in.data<-list(yaug=yaug.2017, M=nrow(yaug.2017), T=ncol(yaug.2017))
# inits<-function()list(z=rep(1, nrow(yaug.2017)), p=runif(1, 0, 1))
# params<-c("N", "p", "omega")
# 
# ni<-100000
# nt<-1
# nb<-5000
# nc<-3
# 
# sr.cp.jm1.2017<-jags.model(data=in.data, inits = inits, file = "sicklefin_redhorse_closed_pop.jags",
#                       n.chains = nc, n.adapt = nb, quiet = F)
# 
# sr.cp.jc1.2017 <- coda.samples(sr.cp.jm1.2017, params, n.iter=ni)
# 
# save(sr.cp.jc1.2016, file="SRH_Closed_Pop_Results_2016.gzip")
# save(sr.cp.jc1.2017, file="SRH_Closed_Pop_Results_2017.gzip")
# save(sr.cp.jc1.2018, file="SRH_Closed_Pop_Results_2017.gzip")
# 
# 
# plot(sr.cp.jc1.2017, ask=TRUE)
# 
# # model results
# out.sr.cp.jc1.2016<-summary(sr.cp.jc1.2016)
# out.sr.cp.jc1.2016
# 
# out.sr.cp.jc1.2017<-summary(sr.cp.jc1.2017)
# out.sr.cp.jc1.2017
# 
# stats.2016<-out.sr.cp.jc1$statistics
# quants.2016<-out.sr.cp.jc1$quantiles
# 
# stats.2017<-out.sr.cp.jc1.2017$statistics
# quants.2017<-out.sr.cp.jc1.2017$quantiles
# 
# stats.2018<-out.sr.cp.jc1$statistics
# quants.2018<-out.sr.cp.jc1$quantiles
# 
# N.stats<-c(stats.2016[1,1], stats.2017[1,1], stats.2018[1,1])
# N.lower<-c(quants.2016[1,1], quants.2017[1,1], quants.2018[1,1])
# N.upper<-c(quants.2016[1,5], quants.2017[1,5], quants.2018[1,5])
# 
# cbind(N.stats, N.lower, N.upper)
# 
# plot(N.stats, ylim=c(0,1100), ylab="Population size", xlab="Year", xaxt='n')
# segments(1:3, N.lower, 1:3, N.upper)
# axis(1, at=1:3, labels=c("2016", "2017", "2018"))





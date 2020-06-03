## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = NA,
  fig.width=6,
  fig.height=5
)

## ----message=FALSE------------------------------------------------------------
library(Grym)

## -----------------------------------------------------------------------------
SeasonDate <- as.Date("2015-12-01")

## -----------------------------------------------------------------------------
SurveyDate <- as.Date("2016-04-20")
surveyI <- as.numeric(SurveyDate-SeasonDate)+c(0,1)
surveyN <- c(29.87,90.53,1177.61,rep(0,7))
surveyB <- 3900.8173

## -----------------------------------------------------------------------------
SpawnDate <- as.Date("2016-11-30")
spawnI <- as.numeric(SpawnDate-SeasonDate) + c(0,1)

## -----------------------------------------------------------------------------
icefishPr <- function(M,F,Catch=0,surveyN,surveyB,surveyI,spawnI,
                      VB.t0=0.06671238,VB.K=0.36842178,VB.Linf=489.73706791,
                      WLa=9.157E-10,WLb=3.316,
                      age.selectivity=approxfun(c(0,2.5,3),c(0,0,1),rule=2),
                      Fmax=2.5) {
    
  ## Ensure 0 included in test fishing mortalities
  F <- sort(union(0,F))
  
  ## Two year projections of 10 age classes with a daily time step
  n.yr <- 2
  n.inc <- 365
  Ages <- 1:10
  Days <- seq(0,1,length=n.inc+1)
  
  ## Matrices of ages, lengths and weights for each day and age class
  as <- outer(Days,Ages,FUN="+")  
  ls <- vonBertalanffyAL(as,t0=VB.t0,K=VB.K,Linf=VB.Linf)
  ws <- powerLW(ls,a=WLa,b=WLb)

  ## Constant intra-annual natural mortality
  ms <- matrix(1,n.inc+1,length(Ages))
  Ms <- ctrapz(ms,1/n.inc)
  MMs <- M*Ms
  
  ## Within year fishing mortality is determined by an age based selectivity 
  fs <- array(age.selectivity(as),dim(as))
  Fs <- ctrapz(fs,1/n.inc)
  
  ### Projection to end of year from survey data
  if(Catch>0) {
    ## Adjust with-year fishing mortality for post-survey Catch
    fs0 <- rep.int(c(0,1),c(max(surveyI),n.inc+1-max(surveyI)))
    fs0 <- fs0/trapz(fs0,1/n.inc)*fs  
    Fs0 <- ctrapz(fs0,1/n.inc)
    pr0 <- projectC(ws,MMs,Fs0,fs0,Catch,surveyN,surveyI,surveyB,surveyI,yield=1,Fmax=Fmax)
    if(pr0$F==Fmax) warning("Target catch could not be recovered")
  } else {
    pr0 <- project(ws,MMs,0,0,surveyN,surveyI,surveyB,surveyI,yield=0)
    pr0$F <- 0
  }
  SSB0 <- meanStock(pr0$B,1,spawnI)

   
  ## Annual cohort totals
  d <- data.frame(Year=c(rep(0:n.yr,length(F))),F=0,Nf=0,Bf=0,Y=0,SSN=0,SSB=0,Escapement=0)
  k <- 0
  
  ## Project forward for prescribed fishing mortalities.
  for(Fk in F) {
    ## Reset to survey year
    pr <- pr0
    d[k <- k+1,] <- data.frame(Year=0,F=Fk,Nf=sum(final(pr$N)),Bf=sum(final(pr$B)),Y=sum(pr$Y),
                               SSN=meanStock(pr$N,1,spawnI),SSB=SSB0,Escapement=1)
    for(yr in seq_len(n.yr)) {
      ## Project
      N0 <- advance(pr$N)
      pr <- project(ws,MMs,Fk*Fs,Fk*fs,N0,yield=1)
      SSB <- meanStock(pr$B,1,spawnI)
      d[k <- k+1,] <- data.frame(Year=yr,F=Fk,Nf=sum(final(pr$N)),Bf=sum(final(pr$B)),Y=sum(pr$Y),
                                 SSN=meanStock(pr$N,1,spawnI),SSB=SSB,Escapement=SSB/SSB0)
    }
  }
  d
}

## -----------------------------------------------------------------------------
d <- icefishPr(M=0.4,F=seq(0.12,0.18,0.001),Catch=0,surveyN,surveyB,surveyI,spawnI)
head(d)

## ----fig.show='hold'----------------------------------------------------------
d2 <- d[d$Year==2,]
d2$RelEscape <- d2$Escapement/d2$Escapement[1]
F75 <- approx(d2$RelEscape[-1],d2$F[-1],0.75)$y
opar <- par(mar = c(5,4,4,4)+0.1)
plot(F~RelEscape,data=d2[-1,],type="l",col="grey60",
     xlab="Relative Escapement",ylab="Fishing Mortality")
abline(v=0.75,h=F75,col="dodgerblue")
points(F~RelEscape,data=d2[-1,],pch=16,col="grey60")
axis(3,at=0.75,col="dodgerblue")
axis(4,at=F75,labels=round(F75,4),col="dodgerblue",las=2)
par(opar)

## -----------------------------------------------------------------------------
d <- icefishPr(M=0.4,F=F75,Catch=0,surveyN,surveyB,surveyI,spawnI)
d

## -----------------------------------------------------------------------------
d <- icefishPr(M=0.4,F=seq(0.12,0.18,0.001),Catch=200,surveyN,surveyB,surveyI,spawnI)
head(d)

## ----fig.show='hold'----------------------------------------------------------
d2 <- d[d$Year==2,]
d2$RelEscape <- d2$Escapement/d2$Escapement[1]
F75 <- approx(d2$RelEscape[-1],d2$F[-1],0.75)$y
opar <- par(mar = c(5,4,4,5)+0.1)
plot(F~RelEscape,data=d2[-1,],type="l",col="grey60",
     xlab="Relative Escapement",ylab="Fishing Mortality")
abline(v=0.75,h=F75,col="dodgerblue")
points(F~RelEscape,data=d2[-1,],pch=16,col="grey60")
axis(3,at=0.75,col="dodgerblue")
axis(4,at=F75,labels=round(F75,5),col="dodgerblue",las=2)
par(opar)

## -----------------------------------------------------------------------------
d <- icefishPr(M=0.4,F=F75,Catch=200,surveyN,surveyB,surveyI,spawnI)
d

## -----------------------------------------------------------------------------
icefishRE <- function(target,M,F,
                      Catch=0,surveyN,surveyB,surveyI,spawnI,
                      VB.t0=0.06671238,VB.K=0.36842178,VB.Linf=489.73706791,
                      WLa=9.157E-10,WLb=3.316,
                      age.selectivity=approxfun(c(0,2.5,3),c(0,0,1),rule=2),
                      Fmax=2.5,tol=1.0E-6) {

  ## Extract summary data from a projection
  annualSummary <- function(yr,F,pr) {
  }

  ## Ensure 0 included in test fishing mortalities
  F <- sort(union(0,F))
  
  ## Two year projections of 10 age classes with a daily time step
  n.yr <- 2
  n.inc <- 365
  Ages <- 1:10
  Days <- seq(0,1,length=n.inc+1)
  
  ## Matrices of ages, lengths and weights for each day and age class
  as <- outer(Days,Ages,FUN="+")  
  ls <- vonBertalanffyAL(as,t0=VB.t0,K=VB.K,Linf=VB.Linf)
  ws <- powerLW(ls,a=WLa,b=WLb)

  ## Constant intra-annual natural mortality
  ms <- matrix(1,n.inc+1,length(Ages))
  Ms <- ctrapz(ms,1/n.inc)
  MMs <- M*Ms
  
  ## Within year fishing mortality is determined by an age based selectivity 
  fs <- array(age.selectivity(as),dim(as))
  Fs <- ctrapz(fs,1/n.inc)
  
  ### Projection to end of year from survey data
  if(Catch>0) {
    ## Adjust with-year fishing mortality for post-survey Catch
    fs0 <- rep.int(c(0,1),c(max(surveyI),n.inc+1-max(surveyI)))
    fs0 <- fs0/trapz(fs0,1/n.inc)*fs  
    Fs0 <- ctrapz(fs0,1/n.inc)
    pr0 <- projectC(ws,MMs,Fs0,fs0,Catch,surveyN,surveyI,surveyB,surveyI,yield=1,Fmax=Fmax)
    if(pr0$F==Fmax) warning("Target catch could not be recovered")
  } else {
    pr0 <- project(ws,MMs,0,0,surveyN,surveyI,surveyB,surveyI,yield=0)
    pr0$F <- 0
  }
  ## Numbers at end of survey year - no recruitment
  N0Survey <- advance(pr0$N)
  SSB0 <- meanStock(pr0$B,1,spawnI)
  
  ## Project ahead and return final SSB 
  ProjectSSB <- function(F,target=0) {
    ## Project and compute SSB for final year
    N0 <- N0Survey
    for(yr in seq_len(n.yr)) {
      pr <- project(ws,MMs,F*Fs,F*fs,N0,yield=0)
      N0 <- advance(pr$N)
    }
    SSB <- meanStock(pr$B,1,spawnI)
    SSB-target
  }
  
  SSB1 <- ProjectSSB(0)
  r <- uniroot(ProjectSSB,F,target=target*SSB1)
  F <- c(0,r$root)
  
  ## Annual cohort totals
  d <- data.frame(Year=c(rep(0:n.yr,length(F))),F=0,Nf=0,Bf=0,Y=0,SSN=0,SSB=0,Escapement=0)
  k <- 0
  
  ## Project forward for prescribed fishing mortalities.
  for(Fk in F) {
    ## Reset to survey year
    pr <- pr0
    d[k <- k+1,] <- data.frame(Year=0,F=pr$F,Nf=sum(final(pr$N)),Bf=sum(final(pr$B)),Y=sum(pr$Y),
                               SSN=meanStock(pr$N,1,spawnI),SSB=SSB0,Escapement=1)
    for(yr in seq_len(n.yr)) {
      ## Project
      N0 <- advance(pr$N)
      pr <- project(ws,MMs,Fk*Fs,Fk*fs,N0,yield=1)
      SSB <- meanStock(pr$B,1,spawnI)
      d[k <- k+1,] <- data.frame(Year=yr,F=Fk,Nf=sum(final(pr$N)),Bf=sum(final(pr$B)),Y=sum(pr$Y),
                                 SSN=meanStock(pr$N,1,spawnI),SSB=SSB,Escapement=SSB/SSB0)
    }
  }
  d
}

## -----------------------------------------------------------------------------
d <- icefishRE(target=0.75,M=0.4,F=c(0,0.5),Catch=0,surveyN,surveyB,surveyI,spawnI)
d$RelEscapement <- d$SSB/d$SSB[1:3]
d

## -----------------------------------------------------------------------------
d <- icefishRE(target=0.75,M=0.4,F=c(0,0.5),Catch=200,surveyN,surveyB,surveyI,spawnI)
d$RelEscapement <- d$SSB/d$SSB[1:3]
d


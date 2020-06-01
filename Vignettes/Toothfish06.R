## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = NA,
  fig.width=6,
  fig.height=5
)

## -----------------------------------------------------------------------------
library(Grym)

## -----------------------------------------------------------------------------
survey.df <- read.csv(textConnection("
Survey,Year,Frac,Age,Density,SE,Area,ObsTotal,ExpTotal
1,1989,0.49,3,0.01,0.01,53383.16,70.32,74.57
1,1989,0.49,4,30.56,8.96,53383.16,70.32,74.57
1,1989,0.49,5,6.83,7.13,53383.16,70.32,74.57
1,1989,0.49,6,0.01,0.01,53383.16,70.32,74.57
1,1989,0.49,7,0.01,0.01,53383.16,70.32,74.57
2,1992,0.77,3,8.01,8.97,53383.16,67.54,85.22
2,1992,0.77,4,27.06,12.9,53383.16,67.54,85.22
2,1992,0.77,5,0.01,0.01,53383.16,67.54,85.22
2,1992,0.77,6,16.8,19.26,53383.16,67.54,85.22
2,1992,0.77,7,5.66,21.84,53383.16,67.54,85.22
3,1998,0.33,3,25.85,7.63,80660.77,373.59,371.54
3,1998,0.33,4,0.01,0.01,80660.77,373.59,371.54
3,1998,0.33,5,85.13,65.51,80660.77,373.59,371.54
3,1998,0.33,6,174.83,104.99,80660.77,373.59,371.54
3,1998,0.33,7,0.01,0.01,80660.77,373.59,371.54
3,1998,0.33,8,66.34,31.68,80660.77,373.59,371.54
4,2000,0.48,3,27.32,8.31,85693.96,198.46,200.63
4,2000,0.48,4,5.8,15.56,85693.96,198.46,200.63
4,2000,0.48,5,59.59,35.74,85693.96,198.46,200.63
4,2000,0.48,6,32.98,47.78,85693.96,198.46,200.63
4,2000,0.48,7,29.64,30.16,85693.96,198.46,200.63
5,2001,0.48,3,14.4,9.37,85693.96,207.12,206.07
5,2001,0.48,4,47.26,17.19,85693.96,207.12,206.07
5,2001,0.48,5,0.01,0.01,85693.96,207.12,206.07
5,2001,0.48,6,101.72,42.56,85693.96,207.12,206.07
5,2001,0.48,7,9.3,37.05,85693.96,207.12,206.07
6,2002,0.42,3,24.57,10.36,42063.96,142.77,140.1
6,2002,0.42,4,28.16,23.4,42063.96,142.77,140.1
6,2002,0.42,5,18.55,30.15,42063.96,142.77,140.1
6,2002,0.42,6,56.89,21.35,42063.96,142.77,140.1
7,2003,0.43,3,0.01,0.01,85123.46,234.65,231.64
7,2003,0.43,4,102.51,28.86,85123.46,234.65,231.64
7,2003,0.43,5,24.19,66,85123.46,234.65,231.64
7,2003,0.43,6,54.69,74.47,85123.46,234.65,231.64
8,2004,0.43,3,0.01,0.01,85693.96,240.42,241.79
8,2004,0.43,4,0.01,0.01,85693.96,240.42,241.79
8,2004,0.43,5,168.88,29.37,85693.96,240.42,241.79
8,2004,0.43,6,20.36,29.24,85693.96,240.42,241.79
9,2005,0.47,3,0.01,0.01,85693.96,173.09,175.94
9,2005,0.47,4,52.75,11.17,85693.96,173.09,175.94
9,2005,0.47,5,0.01,0.01,85693.96,173.09,175.94
9,2005,0.47,6,99.76,18.49,85693.96,173.09,175.94"),header=T)

## -----------------------------------------------------------------------------
SurveyAdjust <- function(survey.df,Ms,M,rec.age) {
  ## Scale density to abundance
  r <- survey.df$Area*survey.df$ObsTotal/survey.df$ExpTotal
  ab.mn <- r*survey.df$Density
  ab.se <- r*survey.df$SE

  ## Compute log survival adjustment 
  inc <- ceiling((nrow(Ms)-1)*survey.df$Frac)+1
  S <- surveySurvival(survey.df$Year,survey.df$Age,inc,inc,Ms,M,rcls=rec.age)
  
  ## Weight rescaled abundance by 1/cv^2 
  ab.wt <- (ab.mn/ab.se)^2
  ## Compute the year of "recruitment" to the target age class
  rc.yr <- survey.df$Year-survey.df$Age+rec.age
  yr <- seq.int(min(rc.yr),max(rc.yr))
  rec.yf <- factor(rc.yr,yr)
  ## Compute the weighted geometric means
  data.frame(
    Year=yr,
    Rec=exp(tapply(ab.wt*(log(ab.mn)-log(S)),rec.yf,sum)/tapply(ab.wt,rec.yf,sum)))
}

## -----------------------------------------------------------------------------
## Constant intra-annual natural mortality
nsteps <- 24
Ages <- 4:35
Days <- seq(0,1,length=nsteps+1)
h <- 1/nsteps
ms <- matrix(1,nsteps+1,length(Ages))
Ms <- ctrapz(ms,h)
M <- 0.13
recruit.df <- SurveyAdjust(survey.df,Ms,M,4)
recruit.df

## -----------------------------------------------------------------------------
catch.df <- read.csv(textConnection("
Year,Catch
1986,0
1987,0
1988,0
1989,0
1990,0
1991,0
1992,0
1993,0
1994,0
1995,3000000
1996,9044000
1997,7915000
1998,3974000
1999,4720000
2000,4984000
2001,6245000
2002,4356000
2003,3501000
2004,3048000
2005,2696000"),header=T)

## -----------------------------------------------------------------------------
historic.df <- merge(recruit.df,catch.df)

## -----------------------------------------------------------------------------
length.df <- read.csv(textConnection("
Age, Length
0, 197.56
1, 251.01
2, 307.54
3, 367.28
4, 430.40
5, 497.03
6, 547.46
7, 594.75
8, 641.07
9, 686.46
10, 730.91
11, 774.47
12, 817.13
13, 858.93
14, 899.88
15, 940.00
16, 979.29
17, 1017.79
18, 1055.51
19, 1092.46
20, 1128.65
21, 1164.11
22, 1198.85
23, 1232.88
24, 1266.22
25, 1298.88
26, 1330.87
27, 1362.22
28, 1392.92
29, 1423.00
30, 1452.47
31, 1481.34
32, 1509.62
33, 1537.33
34, 1564.47
35, 1591.06"),header=T)
length.age <- approxfun(length.df$Age,length.df$Length,rule=2)
plot(length.age,0,55)

## -----------------------------------------------------------------------------
mkSelectivity <- function(ages,ls) {
  approxArray <- function(x,y,arr,rule=2) array(approx(x,y,arr,rule=rule)$y,dim(arr))
  select5pt <- function(x) approxArray(x,c(0,0,1,1,0),ages)

  ## Age based selectivity
  ss <- list()
  ## 1986-1994 - age based selectivity
  ss[[1]] <- approxArray(x=c(0.0,4.1,4.9,5.8,7.0,8.4,9.8,13.7,14.9,16.1,17.3,18.4),
                         y=c(0.0,0.0,0.14,0.5,0.8,0.9,1.0,1.0,0.9,0.85,0.4,0.3),
                         ages)
  ## 1995 - length based selectivity
  ss[[2]] <- rampOgive(ls,670,250)
  ## 1996 - age based selectivity
  ss[[3]] <- select5pt(c(0.0,5.8,7.0,8.2,8.4))
  ## 1997 - age based selectivity
  ss[[4]] <- select5pt(c(0.0,4.9,5.8,11.1,13.7))
  ## 1998 - age based selectivity
  ss[[5]] <- select5pt(c(0.0,5.3,5.8,14.9,17.3))
  ## 1999-2004 - age based selectivity
  ss[[6]] <- select5pt(c(0.0,4.1,8.4,16.1,17.3))

  ## Year 1 = 1986
  list(index=setNames(c(rep(1,9),2:5,rep(6,6),1),1986:2005),ss=ss)
}

## -----------------------------------------------------------------------------
ToothfishProjection <- function(Catches,historic.df,length.df,n.years=35) {


  approxArray <- function(x,y,arr,rule=2) array(approx(x,y,arr,rule=rule)$y,dim(arr))
  
  ## Daily time steps with 8 age classes
  nsteps <- 24
  Ages <- 4:35
  Days <- seq(0,1,length=nsteps+1)
  h <- 1/nsteps

  ## Spawning and monitoring interval
  spawnI <- 13
  monitorI <- 1:25

  ## Ages, length at age and weight at age
  ages <- outer(Days,Ages,FUN="+")
  ls <- approxArray(length.df$Age,length.df$Length,ages)
  ws <- powerLW(ls,2.59E-9,3.2064)

  ## Build selectivity matrices
  sel <- mkSelectivity(ages,ls)
  current.sel <- -1

  ## Constant intra-annual natural mortality
  ms <- matrix(1,nsteps+1,length(Ages))
  Ms <- ctrapz(ms,h)
  Msf <- final(Ms)
  M <- 0.13
  MMs <- M*Ms


  ## Length based maturity 
  gs <- rampOgive(ls,930,300)
  
  ## Within year fishing pattern 
  fwy <- double(nsteps+1)
  fwy[] <- 1
  fwy <- fwy/mean(fwy)


  ## Calculate recruitment parameters from historic data
  ## By method of moments
  rlsd <- sqrt(log(1+(sd(historic.df$Rec)/mean(historic.df$Rec))^2))
  rlmn <- log(mean(historic.df$Rec))-rlsd^2/2
  ## By maximum likelihood
  #rlmn <- mean(log(historic.df$Rec))
  #rlsd <- sd(log(historic.df$Rec))
  
  ## This function performs the a projection for each prescibed gamma.
  function(run) {

    
    ## Median spawning biomass estimated from 1000 samples
    ##R <- matrix(prRecruits(1000*length(Msf),ps),1000,length(Msf))
    ##ssb0 <- spawningB0S(R,gs,ws,Ms,M,spawn=spawnI)

    ## Stochastic initial age structure in the absence of fishing
    N0 <- ageStructureS(rlnorm(length(Msf),rlmn,rlsd),Msf,M)

    ## Recruitment for projection period
    Rs <- rlnorm(n.years,rlmn,rlsd)
    
    ## Annual summary quantities
    n <- nrow(historic.df)+n.years*length(Catches)
    Year <- integer(n)
    Target <- R <- N <- B <- SSN <- SSB <- Catch <- F <- double(n)
    k <- 1

    ## Project over historic period
    for(yr in 1:nrow(historic.df)) {    

      ## Compute fishing mortality when selectivity changes 
      if(sel$index[yr]!=current.sel) {
        current.sel <- sel$index[yr]
        ss <- sel$ss[[current.sel]]
        fs <- fwy*ss
        Fs <- ctrapz(fs,h)
      }
      
      ## Project over year
      Year[k] <- historic.df$Year[yr]
      Target[k] <- historic.df$Catch[yr]
      R[k] <- historic.df$Rec[yr]
      if(yr > 1) N0 <- advance(pr$N,R[k])
      pr <- projectC(ws,MMs,Fs,fs,Target[k],Nref=N0,yield=1,Fmax=5,tol=1.0E-8)

      ## Collate annual summaries
      N[k] <- sum(final(pr$N))
      B[k] <- sum(final(pr$B))
      SSN[k] <- spawningStock(pr$N,ss,spawnI)
      SSB[k] <- spawningStock(pr$B,ss,spawnI)
      Catch[k] <- sum(pr$Y)
      F[k] <- pr$F
      k <- k+1
    }
    
    ## Record pre-projection state
    pr0 <- pr
    Year0 <- historic.df$Year[yr]

    ## Set projection selectivity
    ss <- sel$ss[[1]]
    fs <- fwy*ss
    Fs <- ctrapz(fs,h)


    ## Project for each catch
    for(catch in Catches) {
      ## Reset to pre-projection state
      pr <- pr0

      for(yr in 1:n.years) {
        ## Project over year
        N0 <- advance(pr$N,Rs[yr])
        pr <- projectC(ws,MMs,Fs,fs,catch,Nref=N0,yield=1,Fmax=5,tol=1.0E-8)
        #if(pr$F==5) return(NULL)

        ## Collate annual summaries
        Year[k] <- yr+Year0
        Target[k] <- catch
        R[k] <- Rs[yr]
        N[k] <- sum(final(pr$N))
        B[k] <- sum(final(pr$B))
        SSN[k] <- spawningStock(pr$N,ss,spawnI)
        SSB[k] <- spawningStock(pr$B,ss,spawnI)
        Catch[k] <- sum(pr$Y)
        F[k] <- pr$F
        k <- k+1
      }
    }
    data.frame(Run=run,M=M,Year=Year,Target=Target,
               R=R,N=N,B=B,SSN=SSN,SSB=SSB,Catch=Catch,F=F)
  }
}

## -----------------------------------------------------------------------------
sim <- ToothfishProjection(Catches=c(2.8e6,2.84e6,2.85e6,2.9e6,3e6,3.2e6),
                           historic.df=historic.df,length.df=length.df)

## -----------------------------------------------------------------------------
df <- sim(1)
head(df)
tail(df)

## -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
ggplot(df %>% filter(Year > 2005),aes(x=Year,y=N,colour=factor(Target)))+geom_line()

## -----------------------------------------------------------------------------
library(furrr)
plan(multiprocess)
system.time(df <- future_map_dfr(1:101,sim))

## -----------------------------------------------------------------------------
ggplot(df %>% filter(Year <= 2005),
       aes(x=Year,y=N,group=Run))+
  geom_line(colour=rgb(0,0,0,0.1))+ylim(0,5E07)

## -----------------------------------------------------------------------------
ggplot(df %>% filter(Year > 2005),
       aes(x=Year,y=N,group=Run))+
  geom_line(colour=rgb(0,0,0,0.1))+
  facet_wrap(~Target,ncol=1)+ylim(0,5E07)


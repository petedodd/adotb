library(here)
library(data.table)
library(ggplot2)
library(readxl)

## ----------------- BMI/thinness -------
## BMI distribution
## 18.5 lower quintile
## 95th obese = 30

## fit gamma distribution to these assumptions:
F <- function(x) (qgamma(0.2,shape=x[1],scale=x[2])-18.5)^2 + (qgamma(0.95,shape=x[1],scale=x[2])-30)^2 #error
(ans <- optim(par=c(2,20),fn=F,method='L-BFGS-B',lower=0.1,upper=100)) #ans:  25.4526201  0.8751584
curve(dgamma(x,shape=ans$par[1],scale=ans$par[2]),from=1,to=50) #check

## now what distribution gives a particular level of thinness?
(m <- ans$par[1]*ans$par[2])
(s <- ans$par[1]*ans$par[2]^2)
xr <- pnorm(-2) #ref point for -2z
bmi2z <- qgamma(xr,shape=ans$par[1],scale=ans$par[2])
pgamma(bmi2z,shape=ans$par[1],scale=ans$par[2]) #prob below in std dist
(pc <- pgamma(bmi2z,shape=ans$par[1],scale=ans$par[2]/1.3)) #EG prob below in std dist

## function to scale this so as to achive thinness proportion
getscale <- function(propthin)uniroot(function(x) propthin-pgamma(bmi2z,shape=ans$par[1],scale=ans$par[2]/x),interval = c(0.5,10))$root
getscale(0.2) #test


## https://academic.oup.com/ije/article/39/1/149/713956?login=false
## 13.8% per unit BMI 13.8% (95% CI 13.4−14.2)
## MGF = E[exp(tX)] = 1/(1-theta*t)^k
alph <- log(1-0.138) #per BMI

## un-normalized RR
unRR <- function(a,k,th) 1/(1-th*a)^k
unRR(alph,ans$par[1],ans$par[2]/1.2) / unRR(alph,ans$par[1],ans$par[2]) #test

getRR <- function(scale) unRR(alph,ans$par[1],ans$par[2]/scale) / unRR(alph,ans$par[1],ans$par[2])
getRR(1.2) #test

TH <- read_excel(here('rawdata/Adolescent Thiness Data.xlsx'))
TH <- as.data.table(TH)
TH <- TH[2:nrow(TH)]
names(TH)[2] <- 'pcthin'
names(TH)[3] <- 'go'
TH[,go:=NULL]

## calculate the RRs from thinness
TH[,RR:=1.0]
for(i in 1:nrow(TH)){
  scl <- getscale(TH$pcthin[i]/1e2)
  TH$RR[i] <- getRR(scl)
}

fwrite(TH,file=here('outdata/thinness.RRs.csv'))

## TODO incorporate uncertainty in alph

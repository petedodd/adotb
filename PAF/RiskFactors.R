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
getscale <- function(propthin) uniroot(function(x) propthin-pgamma(bmi2z,shape=ans$par[1],scale=ans$par[2]/x),
                                      interval = c(0.5,10))$root
getscale(0.2) #test

png(here('plots/BMIeg.png'))
curve(dgamma(x,shape=ans$par[1],scale=ans$par[2]/getscale(0.1)),
      from=1,to=50,col=2,
      xlab='BMI',ylab='density') #check
curve(dgamma(x,shape=ans$par[1],scale=ans$par[2]),
      from=1,to=50,add=TRUE) #check
grid()
dev.off()


## https://academic.oup.com/ije/article/39/1/149/713956?login=false
## 13.8% per unit BMI 13.8% (95% CI 13.4−14.2)
## MGF = E[exp(tX)] = 1/(1-theta*t)^k
alph <- log(1-0.138) #per BMI
alph.lo <- log(1-0.142) #per BMI
alph.hi <- log(1-0.134) #per BMI

## un-normalized RR
unRR <- function(a,k,th) 1/(1-th*a)^k
unRR(alph,ans$par[1],ans$par[2]/1.2) / unRR(alph,ans$par[1],ans$par[2]) #test

getRR <- function(scale,A) unRR(A,ans$par[1],ans$par[2]/scale) / unRR(A,ans$par[1],ans$par[2])
getRR(1.2,alph) #test

## NOTE getRR is monotonic decreasing in a (for scale>1):
## (1-B*a)^k/(1-C*a)^k
curve(((1-ans$par[2]*x)/(1-ans$par[2]*x/1.01))^ans$par[1],from=-0.1,to=-0.2)
## and increasing if scale<1

TH <- read_excel(here('rawdata/Adolescent Thiness Data.xlsx'))
TH <- as.data.table(TH)
TH <- TH[2:nrow(TH)]
names(TH)[2] <- 'pcthin'
names(TH)[3] <- 'go'
TH[,go:=NULL]

## calculate the RRs from thinness
TH[,c('RR','RR.lo','RR.hi'):=1.0]
for(i in 1:nrow(TH)){
  scl <- getscale(TH$pcthin[i]/1e2)
  TH$RR[i] <- getRR(scl,alph)
  if(scl>1){
    TH$RR.hi[i] <- getRR(scl,alph.lo)
    TH$RR.lo[i] <- getRR(scl,alph.hi)
  } else {
    TH$RR.hi[i] <- getRR(scl,alph.hi)
    TH$RR.lo[i] <- getRR(scl,alph.lo)
  }
}
TH[,RR.sd:=(RR.hi-RR.lo)/3.92]

TH

fwrite(TH,file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin

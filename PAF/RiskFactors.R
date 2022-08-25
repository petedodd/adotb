library(here)
library(data.table)
library(ggplot2)
library(readxl)
library(metafor)
library(HEdtree) #install via devtools::install_github('petedodd/HEdtree')
lo <- function(x) quantile(x,0.025)
hi <- function(x) quantile(x,1-0.025)
rds <- function(x) round(x,1)
rds(10.1); rds(1.01); rds(-1e-5)
fmt <- function(x,y,z) paste0(rds(x)," (",rds(y)," to ",rds(z),")")
fmtpc <- function(x,y,z) paste0(rds(1e2*x)," (",rds(1e2*y)," to ",rds(1e2*z),")")
gh <- function(x) glue(here(x))

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
getscale <- function(propthin) uniroot(
                                 function(x){
                                   propthin-pgamma(bmi2z,shape=ans$par[1],scale=ans$par[2]/x)
                                 },
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



## ====== HIV ===

## HIV data
load(here('rawdata/H.Rdata')) #NOTE no uncertainty
load(here('rawdata/ckey.Rdata'))
H <- H[iso3 %in% ckey$iso3] #restrict to relevant countries

## IRR estimates
PD <- read.csv(here('rawdata/HIVirrs.csv'))
PD[,c('NAME','DESCRIPTION')]
## hivp,hivpi,artp
P <- parse.parmtable(PD)
names(P)

## reshape
H[,country:=NULL]
H <- melt(H,id='iso3')
H[,c('qty','acat'):=tstrsplit(variable,"\\.")]
H <- dcast(H,iso3+acat~qty,value.var = 'value')

## --- PSAify

## extend H
nrep <- 1e3 #replicates
nn <- nrow(H)
H <- H[rep(1:nn,each = nrep)]
H[,replicate:=rep(1:nrep,nn)]

## make IRRs
IRR <- makePSA(nrep,P)
IRR[,replicate:=1:nrep]
IRR2 <- copy(IRR)
IRR[,acat:='10-14']
IRR2[,acat:='15-19']
IRR <- rbind(IRR,IRR2)
IRR <- merge(IRR,H,by=c('replicate','acat'))
IRR[,c('a','h'):=.(artpc/1e2,hivpc/1e2)]
IRR[,c('artpc','hivpc'):=NULL]
IRR[,irr:=hivpi*(artp*a+1-a)]
IRR[,irr2:=hivp*(artp*a+1-a)]
IRR[,HIVinTB:=h*irr / (h*irr+1-h)]
IRR[,HIVinTB2:=h*irr2 / (h*irr2+1-h)]
IRM <- IRR[,.(HIVinTB=mean(HIVinTB),HIVinTB2=mean(HIVinTB2),hiv=mean(h)),by=.(iso3,acat)] #mean
IRR

## add in IRRs for thinness
ckey[,newcountry:=UN]
setdiff(ckey$UN,TH$Country)
setdiff(TH$Country,ckey$UN)
ckey[UN=="Democratic People's Republic of Korea",newcountry:="DPR Korea"]
ckey[UN=="Democratic Republic of the Congo",newcountry:="DR Congo"]
ckey[UN=="United Republic of Tanzania",newcountry:="UR Tanzania"]
ckey[UN=="Papua New Guinea",newcountry:="PNG"]

tmp <- merge(TH,ckey,by.x = 'Country',by.y = 'newcountry')
tmp <- tmp[,.(iso3,thin=pcthin/1e2,RR,RR.sd)]
tmp[,g.theta:=RR.sd^2/RR]
tmp[,g.k:=RR/g.theta]
tmp <- tmp[rep(1:nrow(tmp),2)]
tmp[,acat:=rep(c('10-14','15-19'),each=nn/2)]
tmp <- tmp[rep(1:nn,nrep)]
tmp[,replicate:=rep(1:nrep,nn)]
tmp[,IRRthin:=rgamma(nrow(tmp),shape=g.k,scale=g.theta)]

## combine
IRR <- merge(IRR,tmp[,.(replicate,iso3,acat,thin,IRRthin)],
             by=c('replicate','iso3','acat'))



## === SHS ===
load(here('rawdata/whokey.Rdata'))

## exposure data
S <- fread(here('rawdata/SHS.csv'),skip=1)
nn <- names(S)
nn <- nn[c(1,2,5,8,11)]
SP <- S[,..nn]

## look
SPM <- melt(SP,id='Country')
ggplot(SPM,aes(variable,value,col=Country,group=Country))+
  geom_point()+geom_line()
ggsave(here('plots/SHS_inputs.pdf'),w=10,h=5,device = cairo_pdf)


## relative
SP[,RR1:=`≥1 day`/`≥1 day`]
SP[,RR3:=`≥3 days`/`≥1 day`]
SP[,RR5:=`≥5 days`/`≥1 day`]
SP[,RR7:=`7 days`/`≥1 day`]

SPM <- melt(SP[,.(Country,RR1,RR3,RR5,RR7)],id='Country')

ggplot(SPM,aes(variable,value,col=Country,group=Country))+
  geom_point()+geom_line()
ggsave(here('plots/SHS_inputsRR.pdf'),w=10,h=5)


## regional meta-analysis to fill in missing
S <- merge(S,ckey,by.x = 'Country',by.y = 'newcountry')
S <- merge(S,whokey,all.y=FALSE,by='iso3')
S[,median(`≥1 day`,na.rm=TRUE),by=g_whoregion]
S[is.na(`≥1 day`),table(g_whoregion)]
S[,median(`7 days`,na.rm=TRUE),by=g_whoregion]
S[is.na(`7 days`),table(g_whoregion)]

S7 <- S[!is.na(`7 days`),.(iso3,g_whoregion,`7 days`,se=(`upper CI_7`-`lower CI_7`)/3.92)]
S7$g_whoregion <- as.factor(S7$g_whoregion)
S7

ma <- rma(yi=S7$`7 days`,sei=S7$se,mods = ~1+S7$g_whoregion)
S7[,.(unique(g_whoregion),as.integer(unique(g_whoregion)))]

shsp <- unique(cbind(as.data.table(predict(ma)),g_whoregion=S7$g_whoregion))

S <- merge(S,shsp[,.(g_whoregion,pred,ci.lb,ci.ub)],by='g_whoregion',all.x = TRUE,all.y=FALSE)
S[,sde:=`7 days`]
S[,sde.hi:=`upper CI_7`]
S[,sde.lo:=`lower CI_7`]

S[is.na(sde),c('sde','sde.lo','sde.hi'):=.(pred,ci.lb,ci.ub)]

## complete version
S <- S[,.(iso3,sde,sde.lo,sde.hi)]
SL <- S[rep(1:nrow(S),each=nrep)]
SL[,replicate:=rep(1:nrep,nrow(S))]
SL <- SL[,.(iso3,replicate,
            sde=sde/1e2,
            E=sde*(100-sde)/((sde.hi-sde.lo)/3.92)^2-1)]
SL[,a:=E*sde]
SL[,b:=E*(1-sde)]
SL[,shs:=rbeta(nrow(SL),a,b)]


## merge into IRR
IRR <- merge(IRR,SL[,.(iso3,replicate,shs)],by=c('iso3','replicate'),all.x=TRUE)

## 1·59, 95% CI 1·11–2·27 from Dogar et al
shsp <- getLNparms(1.59,(2.27-1.11)^2/3.92^2) #parametrize as log normal
curve(dlnorm(x,shsp$mu,shsp$sig),from=0.1,to=5,n=500) #check looks OK

## add in IRR for shs
IRR[,IRRshs:=rlnorm(nrow(IRR),shsp$mu,shsp$sig)]


save(IRR,file=here('PAF/data/IRR.Rdata'))


## [(1-h+h*irr) - (1)] / (1-h+h*irr) = 1 - 1/(1-h+h*irr)
IRR[,PAF.hiv:=1-1/(1-h+h*irr)]
IRR[,PAF.thin:=1-1/(1-thin+thin*IRRthin)]
IRR[,PAF.shs:=1-1/(1-shs+shs*IRRshs)]

IRRS <- IRR[,.(hiv.mid=mean(PAF.hiv),hiv.lo=lo(PAF.hiv),hiv.hi=hi(PAF.hiv),
               thinness.mid=mean(PAF.thin),thinness.lo=lo(PAF.thin),thinness.hi=hi(PAF.thin),
               shs.mid=mean(PAF.shs),shs.lo=lo(PAF.shs),shs.hi=hi(PAF.shs)),
            by=.(iso3,acat)]

IRRS[,thinness:=fmtpc(thinness.mid,thinness.lo,thinness.hi)]
IRRS[,hiv:=fmtpc(hiv.mid,hiv.lo,hiv.hi)]
IRRS[,shs:=fmtpc(shs.mid,shs.lo,shs.hi)]


IRRS <- IRRS[order(iso3,acat),.(iso3,acat,thinness,hiv,shs)] 

fwrite(IRRS,file=here('outdata/PAF.csv'))

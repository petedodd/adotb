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
F <- function(x) (qgamma(0.2,shape=x[1],scale=x[2])-18.5)^2 +
                   (qgamma(0.95,shape=x[1],scale=x[2])-30)^2 #error
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

## ===== BMI ===
load(here('rawdata/BMIREF.Rdata'))
load(here('rawdata/DRA.Rdata'))

## https://academic.oup.com/ije/article/39/1/149/713956?login=false
## 13.8% per unit BMI 13.8% (95% CI 13.4−14.2)
## MGF = E[exp(tX)] = 1/(1-theta*t)^k
alph <- log(1-0.138) #per BMI
alph.lo <- log(1-0.142) #per BMI
alph.hi <- log(1-0.134) #per BMI

## un-normalized RR
unRR <- function(a,k,th) 1/(1-th*a)^k

getRR <- function(A,k,theta,kref,thetaref) unRR(A,k,theta) / unRR(A,kref,thetaref)

## test:
top <- DRA[Country=='India' & Sex=='Girls' & age==19,.(k,theta)]
ref <- BMIREF[Sex=='Girls' & age==19,.(k=kref,theta=thetaref)]
top[,k*theta]
ref[,k*theta]
getRR(alph,top$k,top$theta,ref$k,ref$theta) #test

## TODO new plot of ref and IND

## calculate IRRs
DRA <- merge(DRA,BMIREF,by=c('Sex','age'))
DRA[,c('RR','RR.lo','RR.hi'):=1.0]
DRA[,c('RR','RR.lo','RR.hi'):=.(
       getRR(alph,k,theta,kref,thetaref),
       getRR(alph.hi,k,theta,kref,thetaref),
       getRR(alph.lo,k,theta,kref,thetaref)
     )]

DRAM <- DRA[age>10,.(RRn=mean(RR)),by=Country]
TH <- fread(file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin
DRAM <- merge(DRAM,TH[,.(Country,RRo=RR)],by='Country')

ggplot(DRAM,aes(RRo,RRn,label=Country))+
  geom_point() + ggrepel::geom_text_repel()+
  geom_abline(intercept = 0,slope=1,col=2)+
  geom_smooth(method='lm',forumla = y~x-1)+
  xlab('Previous RR') + ylab('New RR (average over ages/sex)')+
  ggtitle('Comparison of undernutrition BMIs')

ggsave('plots/BMIcf.png',w=6,h=6)

lm(data=DRAM,RRn~RRo-1)

CFA <- merge(DRA[age>=10,.(RRn=mean(RR)),by=Country],
             DRA[age>=10 & age <15,.(RRy=mean(RR)),by=Country],by='Country')

CFA <- merge(CFA,DRA[age>=15,.(RRo=mean(RR)),by=Country],by='Country')

CFA <- DRA[age>=10,.(RR=mean(RR)),by=.(Sex,Country)]
CFA <- dcast(CFA,Country ~ Sex,value.var = 'RR')

ggplot(CFA,aes(Boys,Girls,label=Country))+ geom_point()+ggrepel::geom_text_repel()+geom_abline(intercept=0,slope=1,col=2)
ggsave('CFsex.png',w=7,h=7)

DRA[,bmi:=k*theta]
CFA <- dcast(DRA[age==15],Country ~ Sex,value.var = 'bmi' )


## average by LTBI?
## load LTBI:
library(glue)
load(file=gh('LTBI/data/rnra.Rdata'))
## rnra <- merge(rnra,pzz,by='acat',all.x=TRUE)
## rnra[,prog.recent:=rbeta(nrow(rnra),shape1=A,shape2=B)]
tmp <- rnra[mixing=='assortative'][,.(iso3,age,P)]


ckey

rnra <- merge(rnra,pzzb,by='acat',all.x=TRUE)
rnra <- rnra[order(mixing,replicate,iso3,age)]
## ensure paired over mixing:
rnra[mixing=='assortative',prog.recent1:=rbeta(sum(mixing=='assortative'),shape1=A1,shape2=B1)]
rnra[,prog.recent1:=rep(prog.recent1[1:sum(mixing=='assortative')],2)]
rnra[mixing=='assortative',prog.recent2:=rlnorm(sum(mixing=='assortative'),meanlog = m2,sdlog = s2)]
rnra[,prog.recent2:=rep(prog.recent2[1:sum(mixing=='assortative')],2)]
## rnra[,prog.recent2:=rlnorm(nrow(rnra),meanlog = m2,sdlog = s2)]
## rnra[,prog.recent2:=rbeta(nrow(rnra),shape1=A2,shape2=B2)]

tmp <- rep(rlnorm(nrow(rnra)/2,meanlog=eps$meanlog,sdlog=eps$sdlog),2)
rnra[,prog.slow:=tmp]
## rnra[,prog.slow:=rlnorm(nrow(rnra),meanlog=eps$meanlog,sdlog=eps$sdlog)]

## rnra[,inc0:=(P1) * prog.recent + (P-P1) * prog.slow] #baseline incidence
rnra[,inc0:=(P1) * prog.recent1 + (P2-P1) * prog.recent2 + (P-P2) * prog.slow] #baseline incidence



## fwrite(TH,file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin

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
setdiff(ckey$UN,DRA$Country)
setdiff(DRA$Country,ckey$UN)
ckey[UN=="Democratic People's Republic of Korea",newcountry:="DPR Korea"]
ckey[UN=="Democratic Republic of the Congo",newcountry:="DR Congo"]
ckey[UN=="United Republic of Tanzania",newcountry:="UR Tanzania"]
ckey[UN=="Papua New Guinea",newcountry:="PNG"]

save(ckey,file=here('progression/data/ckey.Rdata'))


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


save(IRR,file=here('PAF/data/IRR.Rdata'))


## [(1-h+h*irr) - (1)] / (1-h+h*irr) = 1 - 1/(1-h+h*irr)
IRR[,PAF.hiv:=1-1/(1-h+h*irr)]
IRR[,PAF.thin:=1-1/(1-thin+thin*IRRthin)]

## summary table
IRRS <- IRR[,.(h.mid=mean(h),h.lo=lo(h),h.hi=hi(h),
               a.mid=mean(a),a.lo=lo(a),a.hi=hi(a),
               th.mid=mean(thin),th.lo=lo(thin),th.hi=hi(thin),
               hiv.mid=mean(PAF.hiv),hiv.lo=lo(PAF.hiv),hiv.hi=hi(PAF.hiv),
               thinness.mid=mean(PAF.thin),thinness.lo=lo(PAF.thin),thinness.hi=hi(PAF.thin)),
            by=.(iso3,acat)]
IRRS[,thinness:=paste0(rds(th.mid*1e2))]
IRRS[,hiv:=fmtpc(h.mid,h.lo,h.hi)]
IRRS[,art:=fmtpc(a.mid,a.lo,a.hi)]
IRRS[,thinness.PAF:=fmtpc(thinness.mid,thinness.lo,thinness.hi)]
IRRS[,hiv.PAF:=fmtpc(hiv.mid,hiv.lo,hiv.hi)]
IRRS <- IRRS[order(iso3,acat),.(iso3,acat,thinness,hiv,art,thinness.PAF,hiv.PAF)]
IRRS <- dcast(IRRS,iso3 ~ acat, value.var = c('thinness','thinness.PAF','hiv','art','hiv.PAF'))
IRRS <- merge(IRRS,ckey[,.(iso3,country=newcountry)],by='iso3')
setcolorder(IRRS,neworder = c("iso3","country",
                              "thinness_10-14",
                              "thinness.PAF_10-14",
                              "hiv_10-14",
                              "art_10-14",
                              "hiv.PAF_10-14",
                              "thinness_15-19",
                              "thinness.PAF_15-19",
                              "hiv_15-19",
                              "art_15-19",
                              "hiv.PAF_15-19"))

fwrite(IRRS,file=here('outdata/PAF.csv'))



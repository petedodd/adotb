library(here)
library(glue)
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

load(here('rawdata/ckey.Rdata'))

## ===== BMI ===
## TODO document these
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

## plot of ref and IND for 19 yo F
png(here('plots/BMIeg.png'))
curve(dgamma(x,shape=ref$k,scale=ref$theta),
      from=1,to=50,col=1,n=1e3,
      xlab='BMI',ylab='density') #check
curve(dgamma(x,shape=top$k,scale=top$theta),
      from=1,to=50,col=2,n=1e3,add=TRUE) #check
text(x=15,y=0.12,,col=2,
     labels = paste0(round(getRR(alph,top$k,top$theta,ref$k,ref$theta),1),'x'))
grid()
dev.off()


## calculate IRRs
DRA <- merge(DRA,BMIREF,by=c('Sex','age'))
DRA[,c('RR','RR.lo','RR.hi'):=1.0]
DRA[,c('RR','RR.lo','RR.hi'):=.(
       getRR(alph,k,theta,kref,thetaref),
       getRR(alph.hi,k,theta,kref,thetaref),
       getRR(alph.lo,k,theta,kref,thetaref)
     )]
save(DRA,file=here('PAF/data/DRAgamma.Rdata'))

DRAM <- DRA[age>10,.(RRn=mean(RR)),by=Country]
TH <- fread(file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin
DRAM <- merge(DRAM,TH[,.(Country,RRo=RR)],by='Country')

ggplot(DRAM,aes(RRo,RRn,label=Country))+
  geom_point() + ggrepel::geom_text_repel()+
  geom_abline(intercept = 0,slope=1,col=2)+
  geom_smooth(method='lm',forumla = y~x-1)+
  xlab('Previous RR using UN summaries') + ylab('NCD-RiSC-based RR (average over ages/sex)')+
  ggtitle('Comparison of undernutrition BMIs')

ggsave('plots/BMIcf.png',w=6,h=6)

## lm(data=DRAM,RRn~RRo-1)

CFA <- merge(DRA[age>=10,.(RRn=mean(RR)),by=Country],
             DRA[age>=10 & age <15,.(RRy=mean(RR)),by=Country],by='Country')

CFA <- merge(CFA,DRA[age>=15,.(RRo=mean(RR)),by=Country],by='Country')

CFA <- DRA[age>=10,.(RR=mean(RR)),by=.(Sex,Country)]
CFA <- dcast(CFA,Country ~ Sex,value.var = 'RR')

ggplot(CFA,aes(Boys,Girls,label=Country))+ geom_point()+ggrepel::geom_text_repel()+geom_abline(intercept=0,slope=1,col=2)
ggsave(here('plots/CFsex.png'),w=7,h=7)

DRA[,bmi:=k*theta]
CFA <- dcast(DRA[age==15],Country ~ Sex,value.var = 'bmi' )

## TODO -- will remove below here

## ## average by LTBI?
## ## load LTBI:

## load(file=gh('LTBI/data/rnra.Rdata'))
## ## rnra <- merge(rnra,pzz,by='acat',all.x=TRUE)
## ## rnra[,prog.recent:=rbeta(nrow(rnra),shape1=A,shape2=B)]
## tmp <- rnra[mixing=='assortative'][,.(iso3,age,P)]


## ckey

## rnra <- merge(rnra,pzzb,by='acat',all.x=TRUE)
## rnra <- rnra[order(mixing,replicate,iso3,age)]
## ## ensure paired over mixing:
## rnra[mixing=='assortative',prog.recent1:=rbeta(sum(mixing=='assortative'),shape1=A1,shape2=B1)]
## rnra[,prog.recent1:=rep(prog.recent1[1:sum(mixing=='assortative')],2)]
## rnra[mixing=='assortative',prog.recent2:=rlnorm(sum(mixing=='assortative'),meanlog = m2,sdlog = s2)]
## rnra[,prog.recent2:=rep(prog.recent2[1:sum(mixing=='assortative')],2)]
## ## rnra[,prog.recent2:=rlnorm(nrow(rnra),meanlog = m2,sdlog = s2)]
## ## rnra[,prog.recent2:=rbeta(nrow(rnra),shape1=A2,shape2=B2)]

## tmp <- rep(rlnorm(nrow(rnra)/2,meanlog=eps$meanlog,sdlog=eps$sdlog),2)
## rnra[,prog.slow:=tmp]
## ## rnra[,prog.slow:=rlnorm(nrow(rnra),meanlog=eps$meanlog,sdlog=eps$sdlog)]

## ## rnra[,inc0:=(P1) * prog.recent + (P-P1) * prog.slow] #baseline incidence
## rnra[,inc0:=(P1) * prog.recent1 + (P2-P1) * prog.recent2 + (P-P2) * prog.slow] #baseline incidence



## fwrite(TH,file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin

## ====== HIV ===
## ART
ART <- fread(here('rawdata/Adolescent ART Coverage for TB Modeling Study 2.csv'))
ART[Value=='<100',Value:='50']
ART[Value=='<200',Value:='100']
ART[Value=='<500',Value:='250']
ART[,Value:=as.numeric(gsub(',','',Value))]
ART[,sex:=substr(Sex,start=1,stop=1)]
ART[,acat:=gsub('Age ','',Age)]
ART <- merge(ART,ckey[,.(iso3,UN)],by.x='Country.Region',by.y='UN')
## HIV
H <- fread(here('rawdata/IHME-GBD_2019_DATA-5a125199-1.csv'))
H <- merge(H,ckey,by.x='location_name',by.y='ihme',all.x=FALSE,all.y=TRUE)
H[,acat:=gsub(' years','',age_name)]
H[,sex:=substr(sex_name,start=1,stop=1)]
## popn
load(here('rawdata/N80MF.Rdata'))
N80 <- N80[iso3 %in% ckey$iso3 & AgeGrp %in% c('10-14','15-19')]
N80[,Year:=NULL]
N80 <- melt(N80,id=c('iso3','AgeGrp'))
N80[,sex:=gsub('T','B',substr(variable,start=4,stop=4))]

## merge
H <- merge(H,N80,by.x=c('iso3','acat','sex'),by.y=c('iso3','AgeGrp','sex'),all.x=TRUE,all.y=FALSE)
H <- merge(H,ART,by=c('iso3','acat','sex'))
H[,art:=Value/(val*value*1e3)]
H[art>1,art:=1.0]
H <- H[,.(iso3,sex,acat,hiv=val,hiv.lo=lower,hiv.hi=upper,pop=1e3*value,art)]

## ## HIV data
## load(here('rawdata/H.Rdata')) #NOTE no uncertainty TODO delete
## H <- H[iso3 %in% ckey$iso3] #restrict to relevant countries

## IRR estimates
PD <- read.csv(here('rawdata/HIVirrs.csv'))
PD[,c('NAME','DESCRIPTION')]
## hivp,hivpi,artp
P <- parse.parmtable(PD)
names(P)

## ## reshape
## H[,country:=NULL]
## H <- melt(H,id='iso3')
## H[,c('qty','acat'):=tstrsplit(variable,"\\.")]
## H <- dcast(H,iso3+acat~qty,value.var = 'value')

## --- PSAify

## extend H
nrep <- 1e3 #replicates
nn <- nrow(H)
H <- H[sex!='B']
H <- H[rep(1:nn,each = nrep)]
H[,replicate:=rep(1:nrep,nn)]
H[,hiv.sd:=(hiv.hi-hiv.lo)/3.92]
H[,hiv:=rgamma(nrow(H), #gamma distribution
               shape=(hiv/(hiv.sd+1e-6))^2,
               scale = hiv.sd^2/(hiv+1e-6))]

## make IRRs
IRR <- makePSA(nrep,P)
IRR[,replicate:=1:nrep]
sexage <- data.table(sex=rep(c('M','F'),2),acat=rep(c('10-14','15-19'),each=2))
IRR <- IRR[rep(1:nrep,each=nrow(sexage))]
sexage <- sexage[rep(1:4,nrep)]
IRR[,c('sex','acat'):=sexage]
IRR <- merge(IRR,H,by=c('replicate','acat','sex'))
IRR[,irr:=hivpi*(artp*art+1-art)]
IRR[,irr2:=hivp*(artp*art+1-art)]
IRR[,HIVinTB:=hiv*irr / (hiv*irr+1-hiv)]
IRR[,HIVinTB2:=hiv*irr2 / (hiv*irr2+1-hiv)]

IRM <- IRR[,.(HIVinTB=mean(HIVinTB),HIVinTB2=mean(HIVinTB2),hiv=mean(hiv)),by=.(iso3,acat)] #mean


## change below to use means over TODO
DRA[,c('RR','RR.lo','RR.hi'):=.(
       getRR(alph,k,theta,kref,thetaref),
       getRR(alph.hi,k,theta,kref,thetaref),
       getRR(alph.lo,k,theta,kref,thetaref)
     )]
save(DRA,file=here('PAF/data/DRA.Rdata'))

DRAM <- DRA[age>10,.(RRn=mean(RR)),by=Country]
TH <- fread(file=here('outdata/thinness.RRs.csv')) #NOTE no uncertainty given for pcthin
DRAM <- merge(DRAM,TH[,.(Country,RRo=RR)],by='Country')


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
load(file=here('PAF/data/IRR.Rdata'))


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



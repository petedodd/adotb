library(here)
library(glue)
library(data.table)
library(ggplot2)
library(scales)
library(ggrepel)
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
ssum <- function(x) sqrt(sum(x^2))

load(here('rawdata/ckey.Rdata'))

## ===== BMI ===
## NOTE see rawdata/BMIdataprep.R
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

## for uncertainty
getRRhilosd <- function(A,A.sd,k,theta,kref,thetaref,N=1e3){
  HLS <- matrix(nrow=length(k),ncol=3)
  az <- rnorm(N,A,A.sd)
  for(i in 1:length(k)){
    rrz <- getRR(az,k[i],theta[i],kref[i],thetaref[i])
    HLS[i,] <- c(quantile(rrz,1-0.025),quantile(rrz,0.025),sd(rrz))
  }
  HLS
}

getRRhilosd(alph,(alph.hi-alph.lo)/3.92,top$k,top$theta,ref$k,ref$theta) #test


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
DRA[,RR:=getRR(alph,k,theta,kref,thetaref)]
HLS <- DRA[,getRRhilosd(alph,(alph.hi-alph.lo)/3.92,k,theta,kref,thetaref)]
DRA[,c('RR.sd','RR.lo','RR.hi'):=.(HLS[,3],HLS[,2],HLS[,1])]
## check
DRA[,any(RR.sd<0)]
DRA[,any(RR.hi<RR.lo)]
DRA[,any(RR.hi<RR)]
DRA[,any(RR.lo>RR)]
## save
save(DRA,file=here('PAF/data/DRAgamma.Rdata'))

CFA <- merge(DRA[age>=10,.(RRn=mean(RR)),by=Country],
             DRA[age>=10 & age <15,.(RRy=mean(RR)),by=Country],by='Country')

CFA <- merge(CFA,DRA[age>=15,.(RRo=mean(RR)),by=Country],by='Country')

CFA <- DRA[age>=10,.(RR=mean(RR)),by=.(Sex,Country)]
CFA <- dcast(CFA,Country ~ Sex,value.var = 'RR')

ggplot(CFA,aes(Boys,Girls,label=Country))+ geom_point()+ggrepel::geom_text_repel()+geom_abline(intercept=0,slope=1,col=2)
ggsave(here('plots/CFsex.png'),w=7,h=7)

DRA[,bmi:=k*theta]
CFA <- dcast(DRA[age==15],Country ~ Sex,value.var = 'bmi' )

## ==== build in non-BMI risk factors
load(here('PAF/data/DRAgamma.Rdata')) #BMI data
ART <- fread(here('rawdata/Adolescent ART Coverage for TB Modeling Study 2.csv')) #ART data
H <- fread(here('rawdata/IHME-GBD_2019_DATA-5a125199-1.csv')) #HIV data
load(here('rawdata/N80MF.Rdata'))                             #demographic data


## ====== HIV ===
## ART
ART[Value=='<100',Value:='50']
ART[Value=='<200',Value:='100']
ART[Value=='<500',Value:='250']
ART[,Value:=as.numeric(gsub(',','',Value))]
ART[,sex:=substr(Sex,start=1,stop=1)]
ART[,acat:=gsub('Age ','',Age)]
ART <- merge(ART,ckey[,.(iso3,UN)],by.x='Country.Region',by.y='UN')
## HIV
H <- merge(H,ckey,by.x='location_name',by.y='ihme',all.x=FALSE,all.y=TRUE)
H[,acat:=gsub(' years','',age_name)]
H[,sex:=substr(sex_name,start=1,stop=1)]
## popn
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

## IRR estimates
PD <- read.csv(here('rawdata/HIVirrs.csv'))
PD[,c('NAME','DESCRIPTION')]
## hivp,hivpi,artp
P <- parse.parmtable(PD)
names(P)

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

## average
DRAM <- DRA[age>9,.(Country,Sex,age,RR,RR.sd)]
DRAM[,acat:=ifelse(age<=14,'10-14','15-19')]
DRAM[,sex:=ifelse(Sex=='Boys','M','F')]
setdiff(DRAM$Country,ckey$newcountry) #these get dropped
DRAML <- DRAM[Country %in% ckey$newcountry] #for plotting
DRAM <- merge(DRAM,ckey[,.(iso3,newcountry)],by.x = 'Country',by.y='newcountry',all.y=TRUE,all.x = FALSE)
DRAM <- DRAM[,.(RR=mean(RR),RR.sd=ssum(RR.sd)/sqrt(5)),by=.(iso3,sex,acat)]

## check haven't spuriously shrunk uncertainty:
DRA[,1e2*mean(RR.sd/RR)]
DRAM[,1e2*mean(RR.sd/RR)]
DRA[,1e2*median(RR.sd/RR)]
DRAM[,1e2*median(RR.sd/RR)]

## lengthen
nn <- nrow(DRAM)
DRAM <- DRAM[rep(1:nn,nrep)]
DRAM[,replicate:=rep(1:nrep,each = nn)]
DRAM[,RR:=rnorm(n=nrow(DRAM),mean=RR,sd=RR.sd)]

## generate/capture thinness measure
load(file=here('rawdata/UW.Rdata'))
UW <- merge(UW,ckey[,.(iso3,newcountry)],by.y='newcountry',by.x = 'Country')
names(UW) <- c('Country','Sex','Year','Age group',
               'p1SD','l1SD','h1SD','p2SD','l2SD','h2SD','iso3')
UW <- UW[`Age group`>9]
UW[,sex:=ifelse(Sex=='Boys','M','F')]
UW[,acat:=ifelse(`Age group`>14,'15-19','10-14')]
UW <- UW[,.(iso3,sex,acat,
            p1SD,p1SD.sd=(h1SD-l1SD)/3.92,
            p2SD,p2SD.sd=(h2SD-l2SD)/3.92)]
UWS <- UW[,.(p1SD=mean(p1SD),p1SD.sd=ssum(p1SD.sd)/sqrt(10),
             p2SD=mean(p2SD),p2SD.sd=ssum(p2SD.sd)/sqrt(10)),
          by=.(iso3,acat)]

## how similar?
GP2 <- ggplot(UWS,aes(p1SD,p2SD,label=iso3))+
  geom_point()+geom_text_repel(show.legend = FALSE)+
  facet_wrap(~acat)+
  scale_x_continuous(label=percent)+scale_y_continuous(label=percent)+
  theme_linedraw()+## theme(legend.position='top')+
  xlab('BMI < 1SD')+ylab('BMI < 2SD')

ggsave(GP2,file=here('plots/UWcf.png'),h=5,w=10)



## combine
IRR <- merge(IRR,DRAM[,.(replicate,iso3,sex,acat,RRbmi=RR,thin=RR)],
             by=c('replicate','iso3','sex','acat'))


save(IRR,file=here('PAF/data/IRR.Rdata'))

load(file=here('PAF/data/IRR.Rdata'))


## [(1-h+h*irr) - (1)] / (1-h+h*irr) = 1 - 1/(1-h+h*irr)
IRR[,PAF.hiv:=1-1/(1-hiv+hiv*irr)]
IRR[,PAF.bmi:=1-1/(1-1.0+1.0*RRbmi)] #coverage is 1 since average for all population

## TODO is HIV uncertainty included?

## summary table
IRRSS <- IRR[,.(h.mid=mean(hiv),h.lo=lo(hiv),h.hi=hi(hiv),
               a.mid=mean(art),a.lo=lo(art),a.hi=hi(art),
               th.mid=mean(thin),th.lo=lo(thin),th.hi=hi(thin),
               hiv.mid=mean(PAF.hiv),hiv.lo=lo(PAF.hiv),hiv.hi=hi(PAF.hiv),
               BMI.mid=mean(PAF.bmi),BMI.lo=lo(PAF.bmi),BMI.hi=hi(PAF.bmi)),
            by=.(iso3,sex,acat)]

tmp <- dcast(data=IRRSS,iso3+acat ~ sex,value.var = c('BMI.mid','hiv.mid'))
tmp[,sex:=NA]

GP <- ggplot(IRRSS,aes(hiv.mid,BMI.mid,col=sex,shape=sex,label=iso3))+
  scale_x_continuous(label=percent)+scale_y_continuous(label=percent)+
  geom_point()+
  geom_text_repel(show.legend = FALSE,max.overlaps=50)+
  facet_wrap(~acat)+
  geom_abline(intercept = 0,slope=1,col=2)+
  theme_linedraw()+theme(legend.position='top')+
  xlab('PAF for HIV/ART')+ylab('PAF for BMI')

GP + geom_segment(data=tmp,aes(x=hiv.mid_F,xend=hiv.mid_M,
                            y=BMI.mid_F,yend=BMI.mid_M,
                            group=paste0(iso3,acat)),col='grey')

ggsave(GP,file=here('plots/PAFsexCF.png'),h=8,w=15)

## TODO perhaps incude segments in the above to higlight sex pairings?

IRRS <- IRR[,.(h.mid=mean(hiv),h.lo=lo(hiv),h.hi=hi(hiv),
                a.mid=mean(art),a.lo=lo(art),a.hi=hi(art),
                hiv.mid=mean(PAF.hiv),hiv.lo=lo(PAF.hiv),hiv.hi=hi(PAF.hiv),
                BMI.mid=mean(PAF.bmi),BMI.lo=lo(PAF.bmi),BMI.hi=hi(PAF.bmi)),
             by=.(iso3,acat)]
IRRS <- merge(IRRS,
              UWS[,.(iso3,acat,th.mid=p2SD,th.lo=pmax(p2SD-1.96*p2SD.sd,0),th.hi=p2SD+1.96*p2SD.sd)],
              by=c('iso3','acat')) #merge underweight data

## below for table 1 TODO check unc
## TODO also include unceretainty in thinness estimates
## IRRS[,thinness:=paste0(rds(th.mid*1e2))]
IRRS[,thinness:=fmtpc(th.mid,th.lo,th.hi)]
IRRS[,hiv:=fmtpc(h.mid,h.lo,h.hi)]
IRRS[,art:=fmtpc(a.mid,a.lo,a.hi)]
IRRS[,BMI.PAF:=fmtpc(BMI.mid,BMI.lo,BMI.hi)]
IRRS[,hiv.PAF:=fmtpc(hiv.mid,hiv.lo,hiv.hi)]
IRRS <- IRRS[order(iso3,acat),.(iso3,acat,thinness,hiv,art,BMI.PAF,hiv.PAF)]
IRRS <- dcast(IRRS,iso3 ~ acat, value.var = c('thinness','BMI.PAF','hiv','art','hiv.PAF'))
IRRS <- merge(IRRS,ckey[,.(iso3,country=newcountry)],by='iso3')
setcolorder(IRRS,neworder = c("iso3","country",
                              "thinness_10-14",
                              "BMI.PAF_10-14",
                              "hiv_10-14",
                              "art_10-14",
                              "hiv.PAF_10-14",
                              "thinness_15-19",
                              "BMI.PAF_15-19",
                              "hiv_15-19",
                              "art_15-19",
                              "hiv.PAF_15-19"))

fwrite(IRRS,file=here('outdata/PAF.csv'))



## quick look at RRs for BMI
GP <- ggplot(DRAML,aes(age,RR,ymin=RR-RR.sd*1.96,ymax=RR+RR.sd*1.96,col=sex,group=sex,fill=sex))+
  geom_ribbon(col=NA,alpha=0.3)+
  geom_line()+
  facet_wrap(~Country)+
  theme_linedraw()+theme(legend.position='top')+
  xlab('Age in years')+ylab('Risk ratio for TB due to BMI distribution')

ggsave(GP,file=here('plots/bmiRRbyage.png'),h=10,w=10)


## MF plot
MF <- IRR[,.(replicate,iso3,sex,acat,RRB = RRbmi * (1-hiv+hiv*irr) )]
MF <- dcast(MF,replicate+iso3+acat ~ sex,value.var = 'RRB')
MF[,mfratio:=M/F]
MF <- MF[,.(mf=mean(mfratio),mf.lo=lo(mfratio),mf.hi=hi(mfratio)),by=.(iso3,acat)]

save(MF,file=here('PAF/data/MF.Rdata'))

load(file=here('PAF/data/MF.Rdata'))

tmp <- MF[acat=='10-14']
der <- order(tmp$mf)
lvls <- unique(tmp[der,iso3])
length(lvls)
MF$iso3 <- factor(MF$iso3,levels=lvls,ordered=TRUE)
MF[,age:=acat]
## text data
TXT <- dcast(data=MF[,.(iso3,age,mf)],iso3 ~ age,value.var='mf'  )
TXT[,txtr:=`15-19`/`10-14`]
TXT[,txt:=round(txtr,2)]
TXT[,c('mf','mf.lo','mf.hi'):=1.5]
TXT[,c('sex','age'):=NA]
TXT <- TXT[iso3 %in% MF$iso3]

pd <- position_dodge(0.25)
GPP <- ggplot(MF,aes(iso3,mf,ymin=mf.lo,ymax=mf.hi,col=age))+
  geom_pointrange(position=pd,shape=1)+
  coord_flip(clip='off')+
  theme_classic()+theme(legend.position='top')+ggpubr::grids()+
  xlab('Country')+ylab('M:F ratio of risk ratios due to HIV & BMI')+
  geom_hline(yintercept = 1,col='grey',lty=2) +
  geom_text(data=TXT,aes(label=txt),show.legend=FALSE)+
  annotate('text',y=1.485,x=31.5,label='Older/Younger\nMF ratios',size=3)
## GPP

ggsave(GPP,file=here('plots/MFcountry.png'),h=7,w=7)

## comparison with notifications
fn <- gh('progression/data/NR.Rdata')
if(file.exists(fn)){
  load(fn)
} else {
  N <- fread('~/Dropbox/Documents/WHO_TBreports/data2020/TB_notifications_2020-10-15.csv')
  NR <- N[year==2019,.(iso3,
                       n1014=newrel_m1014+newrel_f1014,
                       n1519=newrel_m1519+newrel_f1519,
                       n1014m=newrel_m1014,n1014f=newrel_f1014,
                       n1519m=newrel_m1519,n1519f=newrel_f1519)]
  NR <- melt(NR,id='iso3')
  NR[,acat:='10-14']
  NR[grepl('19',variable),acat:='15-19']
  NR[,sex:='B']
  NR[grepl('m',variable),sex:='M']; NR[grepl('f',variable),sex:='F']
  NR <- NR[,.(iso3,acat,sex,
              notified=value)]
  NR <- NR[iso3 %in% unique(rnra$iso3)]
  save(NR,file=fn)
}

NR <- dcast(NR[sex!='B'],iso3+acat~sex,value.var='notified')
NR[,nmf:=M/F]
NR <- merge(MF,NR[,.(iso3,acat,nmf)],by=c('iso3','acat'))
NR <- NR[!is.na(nmf)]


ggplot(NR,aes(nmf,mf,ymin=mf.lo,ymax=mf.hi,col=age,shape=age,label=iso3))+
  geom_pointrange()+
  geom_text_repel(show.legend=FALSE)+
  geom_abline(slope=1,intercept=0,col=2)+
  geom_abline(slope=0,intercept=1,col=2,lty=2)+
  geom_vline(xintercept=1,col=2,lty=2)+
  theme_classic()+ggpubr::grids()+
  xlab('MF ratio in TB notifications')+
  ylab('MF ratio due to BMI & HIV')+
  theme(legend.position='top')

ggsave(file=here('plots/MFcountryVNotes.png'),h=7,w=10)

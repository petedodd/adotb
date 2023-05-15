## TODO more reps?
## using previous outputs to compute incidence
library(here)
library(data.table)
library(scales)
library(ggthemes)
library(ggplot2)
library(ggrepel)
library(glue)
library(HEdtree)
library(ggarchery)
library(ungeviz)

## utility functions
ssum <- function(x) sqrt(sum(x^2))
hi <- function(x,p=0.05) quantile(x,probs=1-p/2)
lo <- function(x,p=0.05) quantile(x,probs=p/2)
rd <- function(x) formatC(round(x),big.mark = ",",format='d')
rd(1); rd(1234); rd(1e9)
rdb <- function(x) format(
                     signif(x,3),
                     digits = 3,
                     nsmall = 0L,
                     big.mark = " ",
                     justify = 'right',
                     drop0trailing = TRUE,
                     scientific = FALSE
                   )
rdb(1); rdb(123456); rdb(1e9)
fmt <- function(x,y,z) paste0(rd(x)," (",rd(y)," to ",rd(z),")")
fmtb <- function(x,y,z) paste0(rdb(x)," (",rdb(y)," to ",rdb(z),")")
gh <- function(x) glue(here(x))
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

load(here('progression/data/ckey.Rdata'))


## Martinez 2-year progression rates:
## 8.8 (3.7-19.7) 10-14
## 10.6 (4.4-23.3)
pmy <- 8.8/1e2; pvy <- (19.7-3.7)^2/392^2
pzy <- getAB(pmy,pvy)
pmo <- 10.6/1e2; pvo <- (23.3-4.4)^2/392^2
pzo <- getAB(pmo,pvo)
pzz <- data.table(acat=c('10-14','15-19'),A=c(pzy$a,pzo$a),B=c(pzy$b,pzo$b))

## 89% within first year, 98% within 2 years
f <- 0.89/0.98 #fraction of 2 year progression within 1 year
## young
pmy1 <- pmy * f; pmy2 <- pmy * (1-f);
pvy1 <- pvy * sqrt(f); pvy2 <- pvy * sqrt(1-f);
pzy1 <- getAB(pmy1,pvy1);
pzy2 <- getLNparms(pmy2,pvy2)
## curve(dbeta(x,pzy1$a,pzy1$b));curve(dbeta(x,pzy2$a,pzy2$b))
## curve(dlnorm(x,pzy2$mu,pzy2$sig))
## old
pmo1 <- pmo * f; pmo2 <- pmo * (1-f);
pvo1 <- pvo * sqrt(f); pvo2 <- pvo * sqrt(1-f);
pzo1 <- getAB(pmo1,pvo1);
## pzo2 <- getAB(pmo2,pvo2)
pzo2 <- getLNparms(pmo2,pvo2)
## curve(dbeta(x,pzo1$a,pzo1$b));curve(dbeta(x,pzo2$a,pzo2$b))
pzzb <- data.table(acat=c('10-14','15-19'),
                   A1=c(pzy1$a,pzo1$a),B1=c(pzy1$b,pzo1$b),
                   ## A2=c(pzy2$a,pzo2$a),B2=c(pzy2$b,pzo2$b)
                   m2=c(pzy2$mu,pzo2$mu),s2=c(pzy2$sig,pzo2$sig)
                   )

## Ragonnet slow progression for >2 years:
eps <- list(meanlog=-6.89,sdlog=0.58)         #nu: Ragonnet
1e2*mean(rlnorm(1e4,meanlog=-6.89,sdlog=0.58)) #~0.1%/y mean

## load LTBI:
load(file=gh('LTBI/data/rnra.Rdata'))
## rnra <- merge(rnra,pzz,by='acat',all.x=TRUE)
## rnra[,prog.recent:=rbeta(nrow(rnra),shape1=A,shape2=B)]
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



## notification data for comparison
fn <- gh('progression/data/NR.Rdata')
if(file.exists(fn)){
  load(fn)
} else {
  N <- fread('~/Dropbox/Documents/WHO_TBreports/data2020/TB_notifications_2020-10-15.csv')
  NR <- N[year==2019,.(iso3,n1014=newrel_m1014+newrel_f1014,n1519=newrel_m1519+newrel_f1519)]
  NR <- melt(NR,id='iso3')
  NR[,acat:='10-14']
  NR[grepl('19',variable),acat:='15-19']
  NR <- NR[,.(iso3,acat,
              notified=value)]
  NR <- NR[iso3 %in% unique(rnra$iso3)]
  save(NR,file=fn)
}

## aggregate over ages
rnr <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2),inc0=mean(inc0)),
               by=.(iso3,acat,replicate,mixing)] #mean over ages

## relative risk data
load(gh('PAF/data/IRR.Rdata'))

## include HIV & pop
rnr <- merge(rnr,
             IRR[,.(value=pop,iso3,acat,replicate,
                    h,irr,thin,IRRthin)],
             by=c('iso3','acat','replicate'))
rnr[,LTBI:=P*value]
rnr[,LTBI1:=P1*value]
rnr[,LTBI2:=P2*value]
rnr[,inc.num0:=inc0*value]
rnr[,inc:=inc0 * (1-h+h*irr)*(1-thin + thin*IRRthin)]
rnr[,inc.num:=inc*value]


## summarize
rnrss <- rnr[,.(P.mid=mean(P),inc0.mid=mean(inc0),inc.mid=mean(inc),
                LTBI.mid=mean(LTBI),
                inc.num0.mid=mean(inc.num0),inc.num.mid=mean(inc.num),
                LTBI1.mid=mean(LTBI1),P1.mid=mean(P1),
                LTBI1.lo=lo(LTBI1),P1.lo=mean(P1),
                LTBI1.hi=hi(LTBI1),P1.hi=mean(P1),
                LTBI2.mid=mean(LTBI2),P2.mid=mean(P2),
                LTBI2.lo=lo(LTBI2),P2.lo=mean(P2),
                LTBI2.hi=hi(LTBI2),P2.hi=mean(P2),
                P.lo=lo(P),P.hi=hi(P),
                inc0.hi=hi(inc0),inc0.lo=lo(inc0),
                inc.hi=hi(inc),inc.lo=lo(inc),
                LTBI.lo=lo(LTBI), LTBI.hi=hi(LTBI),
                inc.num0.hi=hi(inc.num0),inc.num0.lo=lo(inc.num0),
                inc.num.hi=hi(inc.num),inc.num.lo=lo(inc.num)
                ),
             by=.(iso3,acat,mixing)] #mean/hi/lo

## LTBI TOTAL
rnrtot <- rnr[,.(P=weighted.mean(P,value),
                 P1=weighted.mean(P1,value),
                 P2=weighted.mean(P2,value),
                 LTBI=sum(LTBI),
                 LTBI1=sum(LTBI1),
                 LTBI2=sum(LTBI2)
                 ), by=.(acat,replicate,mixing)]
rnrtot <- rnrtot[,.(P.mid=mean(P),P1.mid=mean(P1),P2.mid=mean(P2),
                    LTBI.mid=mean(LTBI),LTBI1.mid=mean(LTBI1),LTBI2.mid=mean(LTBI2),
                    P.lo=lo(P),P1.lo=lo(P1),P2.lo=lo(P2),
                    LTBI.lo=lo(LTBI),LTBI1.lo=lo(LTBI1),LTBI2.lo=lo(LTBI2),
                    P.hi=hi(P),P1.hi=hi(P1),P2.hi=hi(P2),
                    LTBI.hi=hi(LTBI),LTBI1.hi=hi(LTBI1),LTBI2.hi=hi(LTBI2)),
                 by=.(acat,mixing)]
## inc totals
rnrtoti <- rnrss[,.(inc.num.mid=sum(inc.num.mid),inc.num0.mid=sum(inc.num0.mid),
                    inc.num.sd=ssum(inc.num.hi-inc.num.lo)/3.92,
                    inc.num0.sd=ssum(inc.num0.hi-inc.num0.lo)/3.92
         ),by=.(acat,mixing)]
rnrtoti[,c('inc.num0.lo','inc.num0.hi','inc.num.lo','inc.num.hi'):=
           .(inc.num0.mid - 1.96*inc.num0.sd, inc.num0.mid + 1.96*inc.num0.sd,
             inc.num.mid - 1.96*inc.num.sd, inc.num.mid + 1.96*inc.num.sd)]
rnrtoti[,c('inc.num0.sd','inc.num.sd'):=NULL]
rnrtot <- merge(rnrtot,rnrtoti,by=c('acat','mixing'))
rnrtot[,iso3:='TOTAL']

(tnmz <- names(rnrtot))
smy <- rbind(rnrss[,..tnmz],rnrtot) #TODO
lvls <- rnr[,unique(iso3)]
lvls <- lvls[lvls!='TOTAL']
lvls <- sort(as.character(lvls))
lvls <- c(lvls,'TOTAL')
smy$iso3 <- factor(smy$iso3,levels=lvls,ordered = TRUE)

smy[,LTBI.fmt:=fmtb(LTBI.mid,LTBI.lo,LTBI.hi)]
smy[,inc.num0.fmt:=fmtb(inc.num0.mid,inc.num0.lo,inc.num0.hi)]
smy[,inc.num.fmt:=fmtb(inc.num.mid,inc.num.lo,inc.num.hi)]
smy[,.(iso3,acat,LTBI.fmt,inc.num0.fmt,inc.num.fmt)]
smys <- smy[,.(iso3,mixing,acat,LTBI.fmt,
               inc.num0.fmt,inc.num.fmt,
              LTBI.mid,LTBI.lo,LTBI.hi,
              inc.num0.mid,inc.num0.lo,inc.num0.hi,
              inc.num.mid,inc.num.lo,inc.num.hi)]

print(smy[iso3=='TOTAL' & mixing=='assortative',.(sum(inc.num0.mid)/1e6,sum(inc.num.mid)/1e6)])

smy <- merge(smy,NR,by=c('iso3','acat'),all.x=TRUE,all.y = FALSE)
smy[,CDR0:=1e2*notified/inc.num0.mid]
smy[,CDR:=1e2*notified/inc.num.mid]

(tmp <- smy[!is.na(CDR),.(iso3,mixing,acat,CDR0,CDR)][order(mixing)])

fwrite(tmp,file=gh('outdata/CDR.csv'))

## reformat for plotting
smy2 <- smy[,.(iso3,mixing,acat,notified,
               inc.num0.mid,inc.num0.lo,inc.num0.hi,
               inc.num.mid,inc.num.lo,inc.num.hi
               )]
smy2 <- melt(smy2,id=c('iso3','acat','notified','mixing'))
smy2[,method:='with risk factors']
smy2[grepl('0',variable),method:='without risk factors']
smy2[,variable:=gsub('0','',variable)]
smy2 <- dcast(smy2,iso3+acat+notified+method + mixing ~ variable,value.var = 'value')
smy2 <- merge(smy2,ckey,by = 'iso3',all.x=TRUE)


## plot
m <- 1e6*0.75
plt <- ggplot(smy2[method=='with risk factors' & mixing=='assortative'],
              aes(x=notified,y=inc.num.mid,
                  ymin=inc.num.lo,ymax=inc.num.hi,
                  label=newcountry)) +
  geom_abline(slope=1,intercept = 0,col=2)+
  geom_arrowsegment(data=smy2,
                    aes(x = notified-1e5*sqrt(notified/2e5)/4,
                        xend = notified, y = inc.num.mid, yend = inc.num.mid,
                        fill=paste0(method,', ',mixing),
                        col=paste0(method,', ',mixing)),
                    arrows = arrow(type = 'closed',length = unit(0.07, "inches")))+
  geom_point(size=2,shape=1) +
  geom_errorbar(width=10,alpha=0.75)+
  geom_text_repel(show.legend = FALSE,max.overlaps = Inf,nudge_x = 1e2,nudge_y = 20)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  scale_fill_colorblind(name=NULL)+
  scale_color_colorblind(name=NULL)+
  facet_wrap(~acat)+coord_fixed()+# + xlim(0,m)+ylim(0,m)+
  ylab('Estimated tuberculosis incidence 2019 (sqrt scale)')+
  xlab('Notified tuberculosis 2019 (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
## plt


ggsave(plt,file=here('plots/IvN.pdf'),h=7,w=14)
ggsave(plt,file=here('plots/IvN.png'),h=7,w=14)


## barplot version
smy2[,fmeth:=paste0(method,', ',mixing)]
smy2[is.na(newcountry),newcountry:='TOTAL']
cnys <- unique(smy2$newcountry)
cnys <- c(sort(cnys[cnys!='TOTAL']),'TOTAL')
smy2$newcountry <- factor(smy2$newcountry,levels=cnys,ordered = TRUE)
dog <- position_dodge()

## ymax=inc.num.hi,ymin=inc.num.lo
## geom_errorbar(width=0,col=2,position = dog)+ #NOTE tricky to dodge and unhelpful global
plt <- ggplot(smy2,aes(acat,y=inc.num.mid,fill=fmeth))+
  geom_bar(stat='identity',position = dog)+
  geom_point(aes(acat,notified),col=2,show.legend = FALSE)+
  geom_hpline(aes(y = notified, x = acat),col=2,width=1)+
  facet_wrap(~newcountry,scales='free_y')+
  scale_fill_colorblind(name=NULL)+
  scale_y_continuous(label = comma)+
  xlab('Age group (years)')+
  ylab('Tuberculosis incidence 2019')+
  theme(legend.position = c(0.5,0.1/2),legend.direction = 'horizontal')
plt

ggsave(plt,file=here('plots/Ibar.pdf'),h=9,w=12)
ggsave(plt,file=here('plots/Ibar.png'),h=9,w=12)


## comparison table for totals
II <- fread(gh('rawdata/IHME-GBD_2019_DATA-6d1c5bfb-1.csv'))
II <- merge(II,ckey,by.x='location_name',by.y='ihme',all.x=TRUE)
II[,isd:=(upper-lower)/3.92]
II[,acat:=gsub(' years','',age_name)]
II <- II[,.(ihme=sum(val),ihme.sd=ssum(isd)),by=.(iso3,newcountry,acat)]
IIT <- II[,.(ihme=sum(ihme),ihme.sd=ssum(ihme.sd)),by=acat]
IIT <- IIT[,.(acat,fmeth='IHME',inc.num.mid=ihme,
              inc.num.lo=ihme-1.96*ihme.sd,
              inc.num.hi=ihme+1.96*ihme.sd)]

## snow approach, same splits
## ========== active TB =======
fn <- gh('rawdata/TB_burden_age_sex_2020-10-15.csv')
snow <- fread(fn)
exa <- c('0-14','15plus','all')
ado <- c('0-4','5-14','15-24')
snow[,V:=((hi-lo)/3.92)^2]

## sanity checks
snow[sex!='a',]
snow[sex!='a' & !age_group %in% exa ,sum(best)*1e-6]    #yup

## extract
snowC <- snow[sex!='a' & age_group %in% ado[-1]  & iso3 %in% ckey$iso3,
          .(incidence=1.0*sum(best),V=1.0*sum(V)),by=.(age_group,iso3)]


## Kathryn Snow split
## kids
kidf <- 0.206
snowC[age_group %in% c('5-14'),c('incidence','V'):=.(incidence*kidf,V*kidf)]

## older
adof <- (0.300 + 0.454)/2## 0.3
snowC[age_group %in% c('15-24'),c('incidence','V'):=.(incidence*adof,V*adof)]

snowC[,acat:=ifelse(age_group=='5-14','10-14','15-19')]
snowC[,c('age_group'):=NULL]

## check
snowC[,.(sum(incidence)*1e-6),by=acat] #yup

## make total
snowCT <- snowC[,.(inc.num.mid=sum(incidence),V=sum(V)),by=acat]
snowCT <- snowCT[,.(acat,fmeth='Snow',inc.num.mid,
                    inc.num.lo=inc.num.mid-1.96*sqrt(V),
                    inc.num.hi=inc.num.mid+1.96*sqrt(V))]


## ==== ours
TCF <- smy2[newcountry=='TOTAL',.(acat,fmeth,inc.num.mid,inc.num.lo,inc.num.hi)]

## combine
TCF <- rbind(TCF,IIT)
TCF <- rbind(TCF,snowCT)

TCF[,txt:=fmtb(inc.num.mid, inc.num.lo, inc.num.hi)]
TCFe <- dcast(TCF,acat~fmeth,value.var = 'txt')

setcolorder(TCFe,c('acat',
                   'with risk factors, assortative','with risk factors, random',
                   'without risk factors, assortative','without risk factors, random',
                   'IHME','Snow'))

fn <- gh('outdata/cftab.csv')
fwrite(TCFe,file=fn)

print(TCFe)


## country-level comparisons
CC <- list()
CC[[1]] <- snowC[,.(iso3,acat,method='Snow',incidence,
                    incidence.lo=pmax(0,incidence-1.96*sqrt(V)),incidence.hi=incidence+1.96*sqrt(V))]
CC[[2]] <- II[,.(iso3,acat,method='IHME',incidence=ihme,
                 incidence.lo=pmax(0,ihme-1.96*ihme.sd),incidence.hi=ihme+1.96*ihme.sd)]
CC[[3]] <- smy2[newcountry!='TOTAL',.(iso3,acat,method=fmeth,incidence=inc.num.mid,
                           incidence.lo=inc.num.lo,incidence.hi=inc.num.hi)]
CC <- rbindlist(CC)
addon <- CC[method=='with risk factors, assortative',.(iso3,acat,ref=incidence)]
CC <- merge(CC,addon,by=c('iso3','acat'),all.x=TRUE)
CC[,iso3t:=ifelse(method=='IHME',iso3,'')]

M <- 5e5
plt <- ggplot(CC,aes(ref,incidence,ymin=incidence.lo,ymax=incidence.hi,
              col=method,shape=method,label=iso3t))+
  geom_point(size=2,shape=1) +
  geom_errorbar(width=10,alpha=0.75)+
  facet_wrap(~acat)+coord_fixed()+
  geom_text_repel(show.legend = FALSE,max.overlaps = Inf,nudge_x = -12,nudge_y = 50,alpha=0.5)+
  scale_y_sqrt(limits=c(0,M),label=comma)+scale_x_sqrt(limits=c(0,M),label=comma)+
  scale_color_colorblind(name=NULL)+
  geom_abline(intercept = 0,slope=1,col=2)+
  ylab('Estimated tuberculosis incidence 2019 (sqrt scale)')+
  xlab('Reference estimate (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
plt

fac <- 2
WW <- 10
ggsave(plt,file=here('plots/IvE.pdf'),h=WW,w=fac*WW)
ggsave(plt,file=here('plots/IvE.png'),h=WW,w=fac*WW)


## ratio by country: old to young; effect of mixing

rats <- smy2[newcountry!='TOTAL',.(iso3,acat,method=fmeth,incidence=inc.num.mid,
                                      incidence.lo=inc.num.lo,incidence.hi=inc.num.hi)]

rats1 <- dcast(rats,iso3+acat~method,value.var = 'incidence')
rats1 <- rats1[acat=='15-19',.(iso3,
 `with risk factors, assortative` = `with risk factors, assortative` / `without risk factors, random`,
 `without risk factors, assortative` = `without risk factors, assortative` / `without risk factors, random`,
 `with risk factors, random` = `with risk factors, random` / `without risk factors, random`)]

ratsm <- melt(rats1,id='iso3')

ggplot(ratsm,aes(x=iso3,y=value,col=variable,shape=variable)) +
  geom_hline(yintercept = 1,col=2)+
  geom_point(size=2)+
  coord_flip()+
  expand_limits(y=0)+
  xlab('Ratio relative to random mixing & no risk factors')+
  theme(legend.position =  'top',legend.title = element_blank() )## +
  ## guides(col=guide_legend(ncol=2),shape=guide_legend(ncol=2))


ggsave(file=here('plots/ratio_method.png'),h=5,w=7)


## age ratios
rats2 <- dcast(rats[method=='with risk factors, assortative'],iso3~acat,value.var = 'incidence')
rats2[,`ratio by age category`:=`15-19`/`10-14`]

ggplot(rats2,aes(x=iso3,y=`ratio by age category`)) +
  geom_hline(yintercept = 1,col=2)+
  geom_point(size=2)+
  geom_segment(aes(xend=iso3,y=0,yend=`ratio by age category`))+
  coord_flip()

ggsave(file=here('plots/ratio_age.png'),h=5,w=5)



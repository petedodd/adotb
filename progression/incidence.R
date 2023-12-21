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
library(ggh4x)
library(ggpubr)
library(grid)

## utility functions
ssum <- function(x) sqrt(sum(x^2))
hi <- function(x,p=0.05) quantile(x,probs=1-p/2)
lo <- function(x,p=0.05) quantile(x,probs=p/2)
rd <- function(x) formatC(round(x),big.mark = ",",format='d')
rd1 <- function(x) format(round(x,1),nsmall = 1)
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
fmt1 <- function(x,y,z) paste0(rd1(x)," (",rd1(y)," to ",rd1(z),")")
gh <- function(x) glue(here(x))
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

set.seed(1234)

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
## load(file=gh('LTBI/data/rnra.Rdata'))
load(file=gh('LTBI/tmpdata/rnra.full.Rdata'))

## add in acats
rnra[,acat:=ifelse(age>14,'15-19','10-14')]
table(rnra$acat) #OK
rnra$acat <- factor(rnra$acat,levels=c('10-14','15-19'),ordered = TRUE)
rnra <- merge(rnra,pzzb,by='acat',all.x=TRUE)
rnra <- rnra[order(mixing,replicate,iso3,age)]

## ensure paired over mixing:
rnra[mixing=='assortative',prog.recent1:=rbeta(sum(mixing=='assortative'),shape1=A1,shape2=B1)]
rnra[,prog.recent1:=rep(prog.recent1[1:sum(mixing=='assortative')],2)]
rnra[mixing=='assortative',prog.recent2:=rlnorm(sum(mixing=='assortative'),meanlog = m2,sdlog = s2)]
rnra[,prog.recent2:=rep(prog.recent2[1:sum(mixing=='assortative')],2)]

tmp <- rep(rlnorm(nrow(rnra)/2,meanlog=eps$meanlog,sdlog=eps$sdlog),2)
rnra[,prog.slow:=tmp]
## incidence
rnra[,inc0.a:=(P1) * prog.recent1 + (P2-P1) * prog.recent2 + (P-P2) * prog.slow] #baseline incidence
rnra[,inc0.m:=(m.P1) * prog.recent1 + (m.P2-m.P1) * prog.recent2 + (m.P-m.P2) * prog.slow] #baseline incidence, males
rnra[,inc0.f:=(f.P1) * prog.recent1 + (f.P2-f.P1) * prog.recent2 + (f.P-f.P2) * prog.slow] #baseline incidence, females


## notification data for comparison
fn <- gh('progression/data/NR.Rdata')
load(fn)

## relative risk data
load(gh('PAF/data/IRR.Rdata'))
load(here('PAF/data/DRAgamma.Rdata')) #BMI data by age & sex

## include HIV & pop
DRA <- merge(DRA[,.(Country,Sex,age,RR,RR.sd)],ckey[,.(iso3,Country=newcountry)],by='Country')
DRA[,sex:=ifelse(Sex=='Boys','M','F')]
DRA <- DRA[age>9]

## include sex in rnra
rnra <- melt(rnra[,.(acat,iso3,mixing,replicate,age,
                     m.P,m.P1,m.P2,m.inc0=inc0.m,
                     f.P,f.P1,f.P2,f.inc0=inc0.f)],
             id=c('acat','iso3','mixing','replicate','age'))
rnra[,c('sex','variable'):=tstrsplit(variable,split='\\.')]
rnra[,sex:=toupper(sex)]
rnra <- dcast(data=rnra,acat+iso3+mixing+replicate+age+sex ~ variable,value.var='value')

## checks for pattern: OK
rnra[,sum(inc0),by=.(sex,mixing)]
rnra[,sum(inc0),by=.(sex,mixing,acat)]

## merge in
rnra <- merge(rnra,DRA[,.(iso3,sex,age,RR,RR.sd)],
              by=c('iso3','sex','age'))
rnra[,RR:=rnorm(nrow(rnra),mean=RR,sd=RR.sd)]

rnr <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2),inc0=mean(inc0),IRRthin=mean(RR)),
               by=.(iso3,sex,acat,replicate,mixing)] #mean over ages
rnr <- merge(rnr,
             IRR[,.(value=pop,iso3,acat,sex,replicate,
                    hiv,irr)],
             by=c('iso3','sex','acat','replicate'))
rnr[,LTBI:=P*value]
rnr[,LTBI1:=P1*value]
rnr[,LTBI2:=P2*value]
rnr[,inc.num0:=inc0*value]
rnr[,inc:=inc0 * (1-hiv+hiv*irr)*(1-1.0 + 1.0*IRRthin)]
rnr[,inc.num:=inc*value]


## summarize
rnrss <- rnr[,.(P=mean(P),inc0=mean(inc0),inc=mean(inc),
                LTBI=sum(LTBI),
                inc.num0=sum(inc.num0),inc.num=sum(inc.num),
                LTBI1=sum(LTBI1),P1=mean(P1),
                LTBI2=sum(LTBI2),P2=mean(P2)
                ),
             by=.(iso3,acat,mixing,replicate)] #mean/hi/lo
rnrss <- rnrss[,.(P.mid=mean(P),inc0.mid=mean(inc0),inc.mid=mean(inc),
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
rnrss2 <- rnr[,.(P.mid=mean(P),inc0.mid=mean(inc0),inc.mid=mean(inc),
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
             by=.(iso3,sex,acat,mixing)] #mean/hi/lo


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
## by sex and age
rnrtot2 <- rnr[,.(P=weighted.mean(P,value),
                 P1=weighted.mean(P1,value),
                 P2=weighted.mean(P2,value),
                 LTBI=sum(LTBI),
                 LTBI1=sum(LTBI1),
                 LTBI2=sum(LTBI2)
                 ), by=.(sex,acat,replicate,mixing)]
rnrtot2 <- rnrtot2[,.(P.mid=mean(P),P1.mid=mean(P1),P2.mid=mean(P2),
                    LTBI.mid=mean(LTBI),LTBI1.mid=mean(LTBI1),LTBI2.mid=mean(LTBI2),
                    P.lo=lo(P),P1.lo=lo(P1),P2.lo=lo(P2),
                    LTBI.lo=lo(LTBI),LTBI1.lo=lo(LTBI1),LTBI2.lo=lo(LTBI2),
                    P.hi=hi(P),P1.hi=hi(P1),P2.hi=hi(P2),
                    LTBI.hi=hi(LTBI),LTBI1.hi=hi(LTBI1),LTBI2.hi=hi(LTBI2)),
                 by=.(sex,acat,mixing)]


## inc totals
rnrtoti <- rnr[,.(inc.num0=sum(inc.num0),inc.num=sum(inc.num)),by=.(acat,mixing,replicate)]
rnrtoti <- rnrtoti[,.(inc.num.mid=mean(inc.num),inc.num0.mid=mean(inc.num0),
                      inc.num.lo=lo(inc.num),inc.num0.lo=lo(inc.num0),
                      inc.num.hi=hi(inc.num),inc.num0.hi=hi(inc.num0)),
                   by=.(acat,mixing)]
rnrtot <- merge(rnrtot,rnrtoti,by=c('acat','mixing'))
rnrtot[,iso3:='TOTAL']

## by sex and age
rnrtoti2 <- rnr[,.(inc.num0=sum(inc.num0),inc.num=sum(inc.num)),by=.(sex,acat,mixing,replicate)]
rnrtoti2 <- rnrtoti2[,.(inc.num.mid=mean(inc.num),inc.num0.mid=mean(inc.num0),
                      inc.num.lo=lo(inc.num),inc.num0.lo=lo(inc.num0),
                      inc.num.hi=hi(inc.num),inc.num0.hi=hi(inc.num0)),
                   by=.(sex,acat,mixing)]
rnrtot2 <- merge(rnrtot2,rnrtoti2,by=c('sex','acat','mixing'))
rnrtot2[,iso3:='TOTAL']

## by sex only
rnrtoti3 <- rnr[,.(inc.num0=sum(inc.num0),inc.num=sum(inc.num)),by=.(sex,mixing,replicate)]
rnrtoti3 <- rnrtoti3[,.(inc.num.mid=mean(inc.num),inc.num0.mid=mean(inc.num0),
                        inc.num.lo=lo(inc.num),inc.num0.lo=lo(inc.num0),
                        inc.num.hi=hi(inc.num),inc.num0.hi=hi(inc.num0)),
                     by=.(sex,mixing)]


## for getting levels correct
lvls <- rnr[,unique(iso3)]
lvls <- lvls[lvls!='TOTAL']
lvls <- sort(as.character(lvls))
lvls <- c(lvls,'TOTAL')


(tnmz <- names(rnrtot))
smy <- rbind(rnrss[,..tnmz],rnrtot)
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

smy <- merge(smy,NR[sex=='B'],by=c('iso3','acat'),all.x=TRUE,all.y = FALSE)
smy[,CDR0:=1e2*notified/inc.num0.mid]
smy[,CDR:=1e2*notified/inc.num.mid]

(tmp <- smy[!is.na(CDR),.(iso3,mixing,acat,CDR0,CDR)][order(mixing)])

fwrite(tmp,file=gh('outdata/CDR.csv'))
fwrite(smy,file=gh('outdata/smy.csv'))

## CDR stats for text
cdry <- tmp[mixing=='assortative' & acat=='10-14']
cdry <- cdry[order(CDR)]

fwrite(cdry,file=gh('outdata/cdry.csv'))

cdro <- tmp[mixing=='assortative' & acat=='15-19']
cdro <- cdro[order(CDR)]

fwrite(cdro,file=gh('outdata/cdro.csv'))


## reformat
smy[,LTBI.fmt:=fmtb(LTBI.mid,LTBI.lo,LTBI.hi )]
smy[,LTBI1.fmt:=fmtb(LTBI1.mid,LTBI1.lo,LTBI1.hi )]
smy[,LTBI2.fmt:=fmtb(LTBI2.mid,LTBI2.lo,LTBI2.hi )]

## LTBI output by country
out.ltbi <- dcast(smy[mixing=='assortative'],
                  iso3~acat,
                  value.var = c('LTBI.fmt','LTBI1.fmt','LTBI2.fmt'))
out.ltbi$iso3 <- factor(out.ltbi$iso3,levels=lvls,ordered = TRUE)
setkey(out.ltbi,iso3)
setcolorder(out.ltbi,c("iso3","LTBI.fmt_10-14","LTBI2.fmt_10-14","LTBI1.fmt_10-14",
                                     "LTBI.fmt_15-19","LTBI2.fmt_15-19","LTBI1.fmt_15-19"))

fwrite(out.ltbi,file=gh('outdata/out.ltbi.csv'))

## incidence etc by country
out.inc <- smy[,.(iso3,acat,mixing,inc.num.fmt,inc.num0.fmt,notified,CDR,CDR0)]
out.inc <- melt(out.inc,id=c('iso3','acat','mixing'))
out.inc[,RF:=ifelse(grepl('0',variable),'Without risk factors','With risk factors')]
out.inc[,variable:=gsub('0','',variable)]
out.inc <- dcast(out.inc,iso3+mixing+RF ~ variable + acat)
out.inc$iso3 <- factor(out.inc$iso3,levels=lvls,ordered = TRUE)
setkey(out.inc,iso3,mixing,RF)
setcolorder(out.inc,c("iso3","mixing","RF",
                      "inc.num.fmt_10-14","notified_10-14","CDR_10-14",
                      "inc.num.fmt_15-19","notified_15-19","CDR_15-19"   ))

fwrite(out.inc[!is.na(`notified_10-14`)],file=gh('outdata/out.inc.csv'))



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
plt


ggsave(plt,file=here('plots/IvN.pdf'),h=7,w=14)
ggsave(plt,file=here('plots/IvN.png'),h=7,w=14)


## barplot version
smy2[,fmeth:=paste0(method,', ',mixing)]
smy2[is.na(newcountry),newcountry:='TOTAL']
cnys <- unique(smy2$newcountry)
cnys <- c(sort(cnys[cnys!='TOTAL']),'TOTAL')
smy2$newcountry <- factor(smy2$newcountry,levels=cnys,ordered = TRUE)
dog <- position_dodge()

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


## per capita version of Ibar
popss <- unique(IRR[,.(iso3,sex,acat,pop)])
pops <- popss[,.(pop=sum(pop)),by=.(iso3,acat)]
popst <- pops[,.(pop=sum(pop)),by=acat]
popst[,iso3:='TOTAL']
pops <- rbind(pops,popst)
smy3 <- merge(smy2,pops,by=c('iso3','acat'),all.x=TRUE)

plt <- ggplot(smy3,aes(acat,y=1e5*inc.num.mid/pop,fill=fmeth))+
  geom_bar(stat='identity',position = dog)+
  geom_point(aes(acat,1e5*notified/pop),col=2,show.legend = FALSE)+
  geom_hpline(aes(y = 1e5*notified/pop, x = acat),col=2,width=1)+
  facet_wrap(~newcountry,scales='free_y')+
  scale_fill_colorblind(name=NULL)+
  scale_y_continuous(label = comma)+
  xlab('Age group (years)')+
  ylab('Tuberculosis incidence 2019 (per 100,000)')+
  theme(legend.position = c(0.5,0.1/2),legend.direction = 'horizontal')
plt

ggsave(plt,file=here('plots/IbarPC.pdf'),h=9,w=12)
ggsave(plt,file=here('plots/IbarPC.png'),h=9,w=12)



## === sex-based versions for plotting
(tnmz <- names(rnrtot2))
smys <- rbind(rnrss2[,..tnmz],rnrtot2)
smys <- merge(smys,NR[sex!='B'],by=c('iso3','sex','acat'),all.x=TRUE,all.y = FALSE)
smys$iso3 <- factor(smys$iso3,levels=lvls,ordered = TRUE)
smys2 <- smys[,.(iso3,mixing,sex,acat,notified,
                 inc.num0.mid,inc.num0.lo,inc.num0.hi,
                 inc.num.mid,inc.num.lo,inc.num.hi
                 )]

smys2 <- melt(smys2,id=c('iso3','sex','acat','notified','mixing'))
smys2[,method:='with risk factors']
smys2[grepl('0',variable),method:='without risk factors']
smys2[,variable:=gsub('0','',variable)]
smys2 <- dcast(smys2,iso3+sex+acat+notified+method + mixing ~ variable,value.var = 'value')
smys2 <- merge(smys2,ckey,by = 'iso3',all.x=TRUE)
smys2[iso3=='TOTAL',newcountry:='TOTAL']

## plot
m <- 3.1e5
plt <- ggplot(smys2[method=='with risk factors' & mixing=='assortative'],
              aes(x=notified,y=inc.num.mid,
                  ymin=inc.num.lo,ymax=inc.num.hi,
                  label=newcountry)) +
  geom_abline(slope=1,intercept = 0,col=2)+
  geom_arrowsegment(data=smys2,
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
  facet_grid(sex~acat)+coord_fixed()+# + xlim(0,m)+ylim(0,m)+
  ylab('Estimated tuberculosis incidence 2019 (sqrt scale)')+
  xlab('Notified tuberculosis 2019 (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
plt


ggsave(plt,file=here('plots/IvNsex.pdf'),h=14,w=14)
ggsave(plt,file=here('plots/IvNsex.png'),h=14,w=14)


## barplot version
smys2[,fmeth:=paste0(method,', ',mixing)]
cnys <- unique(smys2$newcountry)
cnys <- c(sort(cnys[cnys!='TOTAL']),'TOTAL')
smys2$newcountry <- factor(smys2$newcountry,levels=cnys,ordered = TRUE)
smys2[,ymx:=pmax(inc.num.mid,notified,na.rm=TRUE),by=.(iso3,acat)] #y-max


dog <- position_dodge()
plist <- list()
for(cn in cnys){
  plist[[cn]] <- ggplot(smys2[newcountry==cn],aes(acat,y=inc.num.mid,fill=fmeth))+
    geom_bar(stat='identity',position = dog)+
    geom_point(aes(acat,notified),col=2,show.legend = FALSE)+
    geom_hpline(aes(y = notified, x = acat),col=2,width=1)+
    facet_nested_wrap(~newcountry+sex)+
    scale_y_continuous(label = comma)+
    scale_fill_colorblind(name=NULL)+
    xlab('')+
    ylab('')+
    theme(legend.position = 'top',legend.direction = 'horizontal')+
    theme(strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "grey"),
          plot.margin = unit(c(0,0.0,0,0), 'lines'))
}

GA <- ggarrange(plotlist=plist,nrow=8,ncol=4,common.legend=TRUE)
GA <- annotate_figure(GA,
                      left = textGrob("Tuberculosis incidence 2019",
                                      rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                      bottom = textGrob("Age group (years)", gp = gpar(cex = 1.3)))

ggsave(GA,file=here('plots/Ibarsex.pdf'),h=12,w=12)
ggsave(GA,file=here('plots/Ibarsex.png'),h=12,w=12)



## per capita version of Ibar
popss <- unique(IRR[,.(iso3,sex,acat,pop)])
popsst <- popss[,.(pop=sum(pop)),by=.(sex,acat)]
popsst[,iso3:='TOTAL']
popss <- rbind(popss,popsst)
smys3 <- merge(smys2,popss,by=c('iso3','sex','acat'),all.x=TRUE)
smys3[,ymx:=pmax(1e5*inc.num.mid/pop,1e5*notified/pop,na.rm=TRUE),by=.(iso3,acat)] #y-max


dog <- position_dodge()
plist <- list()
for(cn in cnys){
  plist[[cn]] <- ggplot(smys3[newcountry==cn],
                        aes(acat,y=1e5*inc.num.mid/pop,fill=fmeth))+
    geom_bar(stat='identity',position = dog)+
    geom_point(aes(acat,1e5*notified/pop),col=2,show.legend = FALSE)+
    geom_hpline(aes(y = 1e5*notified/pop, x = acat),col=2,width=1)+
    facet_nested_wrap(~newcountry+sex)+
    scale_fill_colorblind(name=NULL)+
    scale_y_continuous(label = comma)+
    xlab('')+
    ylab('')+
    theme(legend.position = 'top',legend.direction = 'horizontal')+
    theme(strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "grey"),
          plot.margin = unit(c(0,0.0,0,0), 'lines'))
}

## check
plist[['Brazil']]

GA <- ggarrange(plotlist=plist,nrow=8,ncol=4,common.legend=TRUE)
GA <- annotate_figure(GA,
                      left = textGrob("Tuberculosis incidence 2019 (per 100,000)",
                                      rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                      bottom = textGrob("Age group (years)", gp = gpar(cex = 1.3)))

ggsave(plt,file=here('plots/IbarPCsex.pdf'),h=12,w=12)
ggsave(plt,file=here('plots/IbarPCsex.png'),h=12,w=12)


## --- MF plot -------
MF <- rnr[,.(inc0,inc),
          by=.(iso3,sex,acat,mixing,replicate)]
MF <- melt(MF,id=c('iso3','sex','acat','mixing','replicate'))
MF[,variable:=ifelse(grepl(0,variable),'no risk factors','with risk factors')]
MF <- dcast(data=MF,iso3+mixing+variable+replicate+acat ~ sex)
MF[,mf:=M/F]
MF <- MF[,.(mf=mean(mf),mf.lo=lo(mf),mf.hi=hi(mf)),by=.(iso3,mixing,variable,acat)]
## MF <- MF[variable!='no risk factors']
tmp <- MF[acat=='10-14' & variable=='with risk factors' & mixing=='random']
der <- order(tmp$mf)
lvls <- unique(tmp[der,iso3])
MF$iso3 <- factor(MF$iso3,levels=lvls,ordered=TRUE)
MF[,age:=acat]
MF$mixing <- factor(MF$mixing,levels=rev(c('assortative','random')),ordered=TRUE)
## data for text
TXT <- dcast(data=MF,iso3 ~ mixing + variable + acat,value.var='mf')
TXT[,txtr:=`assortative_with risk factors_15-19`/`assortative_with risk factors_10-14`]
TXT[,txt:=round(txtr,2)]
TXT[,c('mf','mf.lo','mf.hi'):=1.7]
TXT[,variable:='with risk factors']
TXT[,sex:='M']
TXT[,acat:='10-14']
TXT[,mixing:='assortative']
TXT[,age:=acat]

## keep only relevant data
tmp <- MF[!(acat=='10-14' & mixing=='assortative')]
tmp <- tmp[!(acat=='15-19' & mixing=='assortative' & variable=='no risk factors')]
tmp <- tmp[!(mixing=='random' & variable=='no risk factors')]

## plot
pd <- position_dodge(0.5)
GPP <- ggplot(tmp,
              aes(iso3,mf,ymin=mf.lo,ymax=mf.hi,
                     col=age,shape=mixing))+
  geom_pointrange(position=pd)+
  scale_shape(solid=FALSE)+
  coord_flip(clip='off')+
  theme_classic()+theme(legend.position='top',)+ggpubr::grids()+
  xlab('Country')+ylab('M:F ratio of risk ratios')+
  geom_hline(yintercept = 1,col='grey',lty=2) +
  geom_text(data=TXT,aes(label=txt),show.legend=FALSE,col='black')+
  annotate('text',y=1.7,x=31.5,label='Older/Younger\nMF ratios',size=3)+
  ylim(c(0.85,1.72))
GPP

ggsave(GPP,file=here('plots/MFcountryFULL.png'),h=7,w=7)
ggsave(GPP,file=here('plots/MFcountryFULL.pdf'),h=7,w=7)


## percentages
## dQ/Q=dlog(A/B) = dA/A + dB/B
PC <- rnrss[mixing=='assortative',
            .(iso3,acat,inc.num.mid,dAoA=abs(inc.num.lo-inc.num.hi)/inc.num.mid/3.92)]
pct <- rnrtot[mixing=='assortative',
              .(acat,tot=inc.num.mid,dBoB=abs(inc.num.lo-inc.num.hi)/inc.num.mid/3.92)]
PC <- merge(PC,pct,by=c('acat'),all.x=TRUE,all.y=FALSE)
PC[,pc:=1e2*inc.num.mid/tot]
PC[,pc.sd:=pc*sqrt(dAoA^2+dBoB^2)]
PC[,pc.lo:=pc-1.96*pc.sd]
PC[,pc.hi:=pc+1.96*pc.sd]
PC[,percentage:=fmt1(pc,pc.lo,pc.hi)]
youngpc <- PC[acat=='10-14'][order(pc,decreasing=TRUE),.(iso3,percentage)]
oldpc <- PC[acat=='15-19'][order(pc,decreasing=TRUE),.(iso3,percentage)]

fwrite(youngpc,file=gh('outdata/PC_ranked_1014.csv'))
fwrite(oldpc,file=gh('outdata/PC_ranked_1519.csv'))

## comparison table for totals
II <- fread(gh('rawdata/IHME-GBD_2019_DATA-6d1c5bfb-1.csv'))
II <- merge(II,ckey,by.x='location_name',by.y='ihme',all.x=TRUE)
II[,isd:=(upper-lower)/3.92]
II[,acat:=gsub(' years','',age_name)]
II <- II[,.(ihme=sum(val),ihme.sd=ssum(isd)),by=.(iso3,newcountry,acat)]
IIT <- II[,.(ihme=sum(ihme),ihme.sd=ssum(ihme.sd)),by=acat]
IIT2 <- IIT[,.(acat='10-19',ihme=sum(ihme),ihme.sd=ssum(ihme.sd))]
IIT <- rbind(IIT,IIT2)
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
allage <- copy(snow)

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
snowCT2 <- snowCT[,.(acat='10-19',inc.num.mid=sum(inc.num.mid),V=sum(V))]
snowCT <- rbind(snowCT,snowCT2)
snowCT <- snowCT[,.(acat,fmeth='Snow',inc.num.mid,
                    inc.num.lo=inc.num.mid-1.96*sqrt(V),
                    inc.num.hi=inc.num.mid+1.96*sqrt(V))]


## ==== ours
TCF <- smy2[newcountry=='TOTAL',.(acat,fmeth,inc.num.mid,inc.num.lo,inc.num.hi)]
TCF2 <- TCF[,.(acat='10-19',inc=sum(inc.num.mid),inc.sdx=ssum(inc.num.hi-inc.num.lo)),by=fmeth]
TCF <- rbind(TCF,TCF2[,.(acat,fmeth,inc.num.mid=inc,inc.num.lo=inc-inc.sdx/2,inc.num.hi=inc+inc.sdx/2)])

## combine
TCF <- rbind(TCF,IIT)
TCF <- rbind(TCF,snowCT)

TCF[,txt:=fmtb(inc.num.mid, inc.num.lo, inc.num.hi)]
TCFe <- dcast(TCF,acat~fmeth,value.var = 'txt')

acata <- c('10-14','15-19','10-19')
TCFe$acat <- factor(TCFe$acat,levels=acata,ordered=TRUE)
setkey(TCFe,acat)

setcolorder(TCFe,c('acat',
                   'with risk factors, assortative','with risk factors, random',
                   'without risk factors, assortative','without risk factors, random',
                   'IHME','Snow'))

fn <- gh('outdata/cftab.csv')
fwrite(TCFe,file=fn)

print(TCFe)

## version with sex 
TWS <- copy(rnrtoti2)
tmp <- copy(rnrtoti3)
tmp[,acat:='10-19']
TWS <- rbind(TWS,tmp,use.names=TRUE)
TWS <- melt(TWS,id=c('sex','mixing','acat'))
TWS[,meth:='with risk factors']
TWS[grepl(0,variable),meth:='without risk factors']
TWS[,variable:=gsub('0','',variable)]
TWS[,fmeth:=paste0(meth,', ',mixing)]
TWS <- dcast(data=TWS,sex+fmeth+acat ~ variable,value.var = 'value')
TWS[,txt:=fmtb(inc.num.mid, inc.num.lo, inc.num.hi)]
TWS <- melt(TWS[,.(sex,acat,fmeth,txt)],id=c('sex','acat','fmeth'))
TWS <- dcast(TWS,acat+sex~fmeth,value.var = 'value')

TWS$acat <- factor(TWS$acat,levels=acata,ordered=TRUE)
setkey(TWS,acat)


setcolorder(TWS,c('acat','sex',
                   'with risk factors, assortative','with risk factors, random',
                   'without risk factors, assortative','without risk factors, random'))

TWS

fn <- gh('outdata/cfWS.csv')
fwrite(TWS,file=fn)


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
  ylab('Ratio relative to random mixing & no risk factors')+
  theme(legend.position =  'top',legend.title = element_blank() )


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

## against mean age of TB
allage <- allage[risk_factor=='all' & sex!='a' & 
                 !age_group %in% c('all','0-14','15plus','18plus'),
                 .(iso3,sex,age_group,best,V)]

allage <- allage[,.(best=sum(best)),by=.(iso3,age_group)]
allage[,agemid:=as.numeric(gsub('(.*)-(.*)','\\1',age_group))]
allage[,agetop:=as.numeric(gsub('(.*)-(.*)','\\2',age_group))]
allage[is.na(agemid),agemid:=65]
allage[is.na(agetop),agetop:=70]
allage[,agemid:=(agemid+agetop+1)/2]
meanage <- allage[,.(`mean age of TB`=weighted.mean(agemid,best)),by=iso3]

rats2 <- merge(rats2,meanage,by='iso3',all.x = TRUE,all.y=FALSE)

ggplot(rats2,aes(label=iso3,x=`mean age of TB`,y=`ratio by age category`)) +
  ## geom_smooth(method = 'lm')+
  geom_point(size=2)+
  geom_text_repel()


ggsave(file=here('plots/ratio_age_scatter.png'),h=5,w=5)

popl <- unique(IRR[,.(iso3,sex,acat,pop)])
popl <- popl[,.(pop=sum(pop)),by=.(iso3,acat)]
nrts <- merge(rats[method=='with risk factors, assortative',.(iso3,acat,incidence)],
              popl,by=c('iso3','acat'))

nrts[,pci:=1e5*incidence/pop]

rrdata <- fread(gh('LTBI/data/RR.csv'))

nrts2 <- dcast(nrts,iso3~acat,value.var = 'pci')
nrts2[,pciratio:=`15-19`/`10-14`]
nrts2 <- merge(nrts2,rrdata,by='iso3')

ggplot(nrts2,aes(rr,pciratio,label=iso3))+
  geom_point(size=2)+
  geom_text_repel()+
  xlab('Estimated infection risk ratio for 15-19 vs 10-14 year olds')+
  ylab('Ratio in estimated per capita TB incidence for 15-19 vs 10-14 year olds')+
  theme_classic()+ggpubr::grids()+
  geom_abline(slope=1,intercept=0,col=2)


ggsave(file=here('plots/ratio_pcivRR_scatter.png'),h=7,w=7)



## HIV/TB by region
load(here('rawdata/whokey.Rdata'))
PAF <- fread(here('outdata/PAF.csv'))
PAF <- PAF[,.(iso3,yng=`hiv.PAF_10-14`,old=`hiv.PAF_15-19`)]
getnum <- function(X) #get mid from string as num
  as.numeric(unlist(lapply(X,function(x) strsplit(x,split=' ')[[1]][1])))
PAF[,midy:=getnum(yng)]
PAF[,mido:=getnum(old)]
PAF <- melt(PAF[,.(iso3,`10-14`=midy/1e2,`15-19`=mido/1e2)],id='iso3')
PAF <- PAF[,.(iso3,acat=variable,phiv=value)]
nhiv <- merge(smy2[method=='with risk factors' & mixing=='assortative',
                   .(iso3,acat,inc.num.mid,inc.num.sdw=inc.num.hi-inc.num.lo)],
              PAF,
              by=c('iso3','acat'))
nhiv <- merge(nhiv,whokey,by='iso3',all.x=TRUE,all.y=FALSE)
nhiva <- nhiv[,.(hivtb=sum(inc.num.mid*phiv),hivtb.w=ssum(inc.num.sdw*phiv)),by=.(acat,g_whoregion)]
nhiv <- nhiv[,.(hivtb=sum(inc.num.mid*phiv),hivtb.w=ssum(inc.num.sdw*phiv)),by=.(g_whoregion)]
nhiv <- rbind(nhiva,nhiv[,.(g_whoregion,acat='10-19',hivtb,hivtb.w)])
nhiv[,c('hivtb.lo','hivtb.hi'):=.(pmax(0,hivtb-hivtb.w/2),hivtb+hivtb.w/2)]
nhiv[,txt:=fmtb(as.integer(hivtb), as.integer(hivtb.lo), as.integer(hivtb.hi))]
nhiv <- nhiv[order(acat,g_whoregion)]

fwrite(nhiv[,.(acat,g_whoregion,txt)],file=here('outdata/nhiv.csv'))

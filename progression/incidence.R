## using previous outputs to compute incidence
library(here)
library(data.table)
library(scales)
library(ggthemes)
library(ggplot2)
library(ggrepel)
library(glue)
library(HEdtree)

## utility functions
hi <- function(x,p=0.05) quantile(x,probs=1-p/2)
lo <- function(x,p=0.05) quantile(x,probs=p/2)
rd <- function(x) formatC(round(x),big.mark = ",",format='d')
rd(1); rd(1234); rd(1e9)
fmt <- function(x,y,z) paste0(rd(x)," (",rd(y)," to ",rd(z),")")
gh <- function(x) glue(here(x))


## Martinez progression rates:
## 8.8 (3.7-19.7) 10-14
## 10.6 (4.4-23.3)
pzy <- getAB(8.8/1e2,(19.7-3.7)^2/392^2)
pzo <- getAB(10.6/1e2,(23.3-4.4)^2/392^2)
pzz <- data.table(acat=c('10-14','15-19'),A=c(pzy$a,pzo$a),B=c(pzy$b,pzo$b))

## Ragonnet slow progression:
eps <- list(meanlog=-6.89,sdlog=0.58)         #nu: Ragonnet
1e2*mean(rlnorm(1e4,meanlog=-6.89,sdlog=0.58)) #~0.1%/y mean

## load LTBI:
load(file=gh('LTBI/data/rnra.Rdata'))
rnra <- merge(rnra,pzz,by='acat',all.x=TRUE)
rnra[,prog.recent:=rbeta(nrow(rnra),shape1=A,shape2=B)]
rnra[,prog.slow:=rlnorm(nrow(rnra),meanlog=eps$meanlog,sdlog=eps$sdlog)]
rnra[,inc0:=(P2) * prog.recent + (P-P2) * prog.slow] #baseline incidence

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
rnr <- rnra[,.(P=mean(P),P2=mean(P2),inc0=mean(inc0)),
               by=.(iso3,acat,replicate)] #mean over ages

## relative risk data
load(gh('PAF/data/IRR.Rdata'))

## include HIV & pop
rnr <- merge(rnr,
             IRR[,.(value=pop,iso3,acat,replicate,
                    h,irr,thin,IRRthin)],
             by=c('iso3','acat','replicate'))
rnr[,LTBI:=P*value]
rnr[,LTBI2:=P2*value]
rnr[,inc.num0:=inc0*value]
rnr[,inc:=inc0 * (1-h+h*irr)*(1-thin + thin*IRRthin)]
rnr[,inc.num:=inc*value]


## TODO check uncertainty
## TOTAL
rnrtot <- rnr[,.(P=weighted.mean(P,value),
                 P2=weighted.mean(P2,value),
                 inc0=weighted.mean(inc0,value),
                 inc=weighted.mean(inc,value),
                 LTBI=sum(LTBI),
                 LTBI2=sum(LTBI2),
                 inc.num0=sum(inc.num0),
                 inc.num=sum(inc.num),
                 value=sum(value)
                 ), by=.(acat,replicate)]
rnrtot[,iso3:='TOTAL']
names(rnr)
(tnmz <- names(rnrtot))
rnr <- rbind(rnr[,..tnmz],rnrtot) #TODO
lvls <- rnr[,unique(iso3)]
lvls <- lvls[lvls!='TOTAL']
lvls <- sort(as.character(lvls))
lvls <- c(lvls,'TOTAL')
rnr$iso3 <- factor(rnr$iso3,levels=lvls,ordered = TRUE)


## summarize
rnrss <- rnr[,.(P.mid=mean(P),inc0.mid=mean(inc0),inc.mid=mean(inc),
                LTBI.mid=mean(LTBI),
                inc.num0.mid=mean(inc.num0),inc.num.mid=mean(inc.num),
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
             by=.(iso3,acat)] #mean/hi/lo

smy <- rnrss[,.(iso3,acat,
                LTBI.mid,LTBI.lo,LTBI.hi,
                inc.num0.mid,inc.num0.lo,inc.num0.hi,
                inc.num.mid,inc.num.lo,inc.num.hi,
                LTBI2.mid,LTBI2.lo,LTBI2.hi,
                inc.mid,inc.lo,inc.hi,
                inc0.mid,inc0.lo,inc0.hi
                )]
smy[,LTBI.fmt:=fmt(LTBI.mid,LTBI.lo,LTBI.hi)]
smy[,inc.num0.fmt:=fmt(inc.num0.mid,inc.num0.lo,inc.num0.hi)]
smy[,inc.num.fmt:=fmt(inc.num.mid,inc.num.lo,inc.num.hi)]
smy[,.(iso3,acat,LTBI.fmt,inc.num0.fmt,inc.num.fmt)]
smys <- smy[,.(iso3,acat,LTBI.fmt,
               inc.num0.fmt,inc.num.fmt,
              LTBI.mid,LTBI.lo,LTBI.hi,
              inc.num0.mid,inc.num0.lo,inc.num0.hi,
              inc.num.mid,inc.num.lo,inc.num.hi)]
print(smy[iso3=='TOTAL',.(sum(inc.num0.mid)/1e6,sum(inc.num.mid)/1e6)])

fwrite(smy,file=gh('outdata/smy.csv')) #~1.4M or 1.7M
save(smy,file=gh('progression/data/smy.Rdata'))

smy <- merge(smy,NR,by=c('iso3','acat'),all.x=TRUE,all.y = FALSE)
smy[,CDR0:=1e2*notified/inc.num0.mid]
smy[,CDR:=1e2*notified/inc.num.mid]

tmp <- smy[!is.na(CDR),.(iso3,acat,CDR0,CDR)]
fwrite(tmp,file=gh('outdata/CDR.csv'))

## reformat for plotting
smy2 <- smy[,.(iso3,acat,notified,
               inc.num0.mid,inc.num0.lo,inc.num0.hi,
               inc.num.mid,inc.num.lo,inc.num.hi
               )]
smy2 <- melt(smy2,id=c('iso3','acat','notified'))
smy2[,method:='with risk factors']
smy2[grepl('0',variable),method:='without risk factors']
smy2[,variable:=gsub('0','',variable)]
smy2 <- dcast(smy2,iso3+acat+notified+method~variable,value.var = 'value')

## plot
m <- 1e6*0.75
plt <- ggplot(smy2,aes(x=notified,y=inc.num.mid,
                      ymin=inc.num.lo,ymax=inc.num.hi,
                      label=iso3,col=method)) +
  geom_point(shape=1,size=2) +
  geom_errorbar(width=10,alpha=0.75)+
  geom_text_repel(show.legend = FALSE)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  geom_abline(slope=1,intercept = 0,col=2)+
  facet_wrap(~acat)+coord_fixed()+# + xlim(0,m)+ylim(0,m)+
  ylab('Estimated incidence 2019 (sqrt scale)')+
  xlab('Notified TB 2019 (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
plt

## TODO smarten
ggsave(plt,file=here('plots/IvN.pdf'),h=7,w=14)
ggsave(plt,file=here('plots/IvN.png'),h=7,w=14)


## quick compare with KS data
load('~/Dropbox/Documents/comms/KatharinaKranzer/ado2/graphs/TBC.Rdata')

TBC
smy2

smy2 <- merge(smy2,TBC[,.(iso3,KSI=incidence,acat)],
              by=c('iso3','acat'),all.x=TRUE,all.y=FALSE)

plt <- ggplot(smy2,aes(x=KSI,y=inc.num.mid,
                ymin=inc.num.lo,ymax=inc.num.hi,
                label=iso3,col=method)) +
  geom_point(shape=1,size=2) +
  geom_errorbar(width=10,alpha=0.75)+
  geom_text_repel(show.legend = FALSE)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  geom_abline(slope=1,intercept = 0,col=2)+
  facet_wrap(~acat)+coord_fixed()+# + xlim(0,m)+ylim(0,m)+
  ylab('Estimated incidence 2019 (sqrt scale)')+
  xlab('KS estimate TB 2017 (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
plt

ggsave(plt,file='~/Dropbox/Documents/comms/KatharinaKranzer/ado2/graphs/TBC.png')

tmp <- smy2[iso3!='TOTAL',.(Snow=sum(KSI)/1e6,prog=sum(inc.num.mid)/1e6),by=.(acat,method)]
tmp
fwrite(tmp,file='~/Dropbox/Documents/comms/KatharinaKranzer/ado2/graphs/tbc.csv')

## comparison of totals
tot <- smy2[iso3!='TOTAL',.(snowtot=sum(KSI),inc=sum(inc.num.mid)),by=.(acat,method)]
tot[,inc/snowtot]
tot <- smy2[iso3!='TOTAL',.(snowtot=sum(KSI),inc=sum(inc.num.mid)),by=.(method)]
tot[,inc/snowtot]


m <- 1e6*0.75
plt <- ggplot(smy2,aes(x=notified,y=KSI,
                       label=iso3,col=method)) +
  geom_point(shape=1,size=2) +
  geom_text_repel(show.legend = FALSE)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  geom_abline(slope=1,intercept = 0,col=2)+
  facet_wrap(~acat)+coord_fixed()+# + xlim(0,m)+ylim(0,m)+
  ylab('KSI 2017 (sqrt scale)')+
  xlab('Notified TB 2019 (sqrt scale)')+
  theme_light()+theme(legend.position = 'top')
plt

ggsave(plt,file='~/Dropbox/Documents/comms/KatharinaKranzer/ado2/graphs/TBC2.png')


## TODO list
## check prog factor
## compare snow, compare IHME?
## include w/ and w/o mixing
## go back over checks in LTBI code

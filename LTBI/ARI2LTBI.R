rm(list=ls())
library(here)
library(data.table)
library(ggplot2)
library(scales)
library(ggthemes)
library(glue)
library(ggrepel)

## utility functions
hi <- function(x,p=0.05) quantile(x,probs=1-p/2)
lo <- function(x,p=0.05) quantile(x,probs=p/2)
rd <- function(x) formatC(round(x),big.mark = ",",format='d')
rd(1); rd(1234); rd(1e9)
fmt <- function(x,y,z) paste0(rd(x)," (",rd(y)," to ",rd(z),")")
gh <- function(x) glue(here(x))


## read in ARI run data from LTBI/tmpdata
fz <- dir(path=here('LTBI/tmpdata'),pattern = "z")
rnr <- list()
for(f in fz){
  print(f)
  fn <- gh('LTBI/tmpdata/{f}')
  load(fn)
  rnr[[f]] <- as.data.table(runsdf)
}
rnr <- rbindlist(rnr)

## relative ARIs
load(gh('LTBI/data/RR.Rdata'))
rnr <- merge(rnr,RR,by='iso3',all.x = TRUE, all.y = FALSE)

## work on converting this to LTBI
rnr[,ari:=exp(lari)]          #true ARI
rnr[,year:= 2019 - year]       #age
rnr <- rnr[order(replicate,iso3,year)] #order
rnr <- rnr[year<20]

## rrmask for ages
newrr <- paste0('rr',0:19)
rnewrr <- paste0('R_',newrr)
rnr[,c(newrr):=rr]
rnr[,c(rnewrr):=rr]
nmnewrr <- paste0('nm.rr',0:19) #no mix versions
nmrnewrr <- paste0('nm.R_',newrr)
rnr[,c(nmnewrr):=rr]
rnr[,c(nmrnewrr):=rr]
nmnewrr <- paste0('nm.rr',0:19) #no mix versions
y1rnewrr <- paste0('y1R_',newrr) #yrs 1
nmy1rnewrr <- paste0('nm.y1R_',newrr)

## NOTE can't acquire risk for ages not yet lived
for(k in 0:19){
  ## all
  ## -- mixing
  nm <- paste0('rr',k)
  rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year>k,c(nm):=0.0]
  ## -- no mixing
  nm <- paste0('nm.rr',k)
  rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year>k,c(nm):=0.0]
  ## recent
  ## -- mixing
  nm <- paste0('R_rr',k)
  rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year>min(1,k),c(nm):=0.0]
  ## -- no mixing
  nm <- paste0('nm.R_rr',k)
  rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year>min(1,k),c(nm):=0.0]
  ## really recent = 1 year
  ## -- mixing
  nm <- paste0('y1R_rr',k)
  ## rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year==0,c(nm):=ifelse(k>=14,rr,1.0)] #step @ 15
  rnr[year>0,c(nm):=0.0]
  ## -- no mixing
  nm <- paste0('nm.y1R_rr',k)
  ## rnr[year>=max(0,k+1-14),c(nm):=1.0] #step @ 15
  rnr[year==0,c(nm):=1.0] #step @ 15
  rnr[year>0,c(nm):=0.0]
}

rnr[replicate==1 & iso3=='AGO'] #check


## sums
## -- mix versions
rnr[,c(newrr):=lapply(.SD,function(x) x*ari),.SDcols=newrr] #multiply by ARI: all
rnr[,c(rnewrr):=lapply(.SD,function(x) x*ari),.SDcols=rnewrr]   #multiply by ARI: 2y recent
rnr[,c(y1rnewrr):=lapply(.SD,function(x) x*ari),.SDcols=y1rnewrr] #multiply by ARI: 1y recent
## -- no mix versions
rnr[,c(nmnewrr):=lapply(.SD,function(x) x*ari),.SDcols=nmnewrr] #multiply by ARI: all
rnr[,c(nmrnewrr):=lapply(.SD,function(x) x*ari),.SDcols=nmrnewrr]   #multiply by ARI: 2y recent
rnr[,c(nmy1rnewrr):=lapply(.SD,function(x) x*ari),.SDcols=nmy1rnewrr] #multiply by ARI: 1y recent

## cumulatives
## -- mix versions
rnr[,c(newrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=newrr] #cumulative sum
rnr[,c(rnewrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=rnewrr] #cumulative sum
rnr[,c(y1rnewrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=y1rnewrr] #cumulative sum
## -- no mix versions
rnr[,c(nmnewrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=nmnewrr] #cumulative sum
rnr[,c(nmrnewrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=nmrnewrr] #cumulative sum
rnr[,c(nmy1rnewrr):=lapply(.SD,sum),by=.(iso3,replicate),.SDcols=nmy1rnewrr] #cumulative sum

## reshape results
keep <- c('iso3','replicate','year',
          newrr,rnewrr,y1rnewrr,
          nmnewrr,nmrnewrr,nmy1rnewrr)
tmp <- melt(rnr[,..keep],id=c('iso3','replicate','year'))

## safety
tmp[,year2:=gsub('y1','',variable)] #remove y1
tmp[,year2:=gsub("\\.", "", year2)]
tmp[,year2:=as.numeric(gsub("[a-zA-z]", "", year2))]
tmp <- tmp[year==year2] #diagonal only (wasteful)
tmp[,year2:=NULL]       #drop
tmp[,mixing:='random']
tmp[!grepl('nm',variable),mixing:='assortative']
tmp[,variable:=gsub('nm\\.','',variable)]
tmp[grepl('y1R_',variable),variable:='dH1']
tmp[grepl('R_',variable),variable:='dH2']
tmp[grepl('rr',variable),variable:='H']

## checks
tmp[iso3=='AGO' & replicate==1]
tmp[iso3=='AGO' & replicate==1,unique(variable)]
tmp[iso3=='AGO' & replicate==1,table(variable)]
tmp[iso3=='AGO' & replicate==1,table(mixing)]

## reshape
tmp <- dcast(tmp,iso3+mixing+replicate+year~variable,value.var = 'value')

## merge back & drop
drop <- c(newrr,rnewrr,y1rnewrr,nmnewrr,nmrnewrr,nmy1rnewrr)
rnr[,c(drop):=NULL]
rnr <- merge(rnr,tmp,by=c('iso3','replicate','year'))

rnr[,dH2:=H-dH2] #NOTE dH is actually cumulative EXCEPT last 2 years in calx
rnr[,dH1:=H-dH1] #NOTE dH is actually cumulative EXCEPT last 1 years in calx


## classic approach (w/o mixing)
rnr <- rnr[order(mixing,replicate,iso3,year)] #order

## rnr[,H0:=cumsum(ari),by=.(iso3,replicate)] #cumulative ARI, not including RR
## ## distinguish past 2 years
## mask <- rep(1,length(unique(rnr$year)))
## mask[1:2] <- 0                          #all except last 2 years
## rnr[,dH0:=cumsum(ari*mask),by=.(iso3,replicate)] #cumhaz!2y no mix

## convert to prevalence
rnr[,P:=1-exp(-H)]                  #ever
rnr[,Pf1:=-exp(-H)+exp(-dH1)]         #1st recent=prob ever - prob not<1
rnr[,Pf2:=-exp(-H)+exp(-dH2)]         #1st recent=prob ever - prob not<2

## infection protection from LTBI - Andrews: 0.79 .7-.86
pm <- 0.79
pv <- (0.86-0.7)^2/3.92^2
apb <- pm*(1-pm)/pv-1
pa <- pm*apb                            #77.88
pb <- (1-pm)*apb                        #20.70
## curve(dbeta(x,shape1 = pa,shape2=pb),from=0,to=1)
## abline(v=pm,col=2);abline(v=.86,col=2,lty=2);abline(v=.7,col=2,lty=2);
## swap
alph <- rbeta(nrow(rnr),shape1=pb,shape2=pa)

rnr[,P1:=alph*(H-dH1) + (1-alph)*(exp(-dH1)-exp(-H))]      #anyrecent <1
rnr[,P2:=alph*(H-dH2) + (1-alph)*(exp(-dH2)-exp(-H))]      #anyrecent <2
## rnr[,P20:=alph*(H0-dH0) + (1-alph)*(exp(-dH0)-exp(-H0))]      #anyrecent

## check
rnr[iso3=='AGO' & replicate==1]
rnr[iso3=='AGO' & replicate==1 & year==18]

## restrict to relevant age groups
rnra <- rnr[year>=10 & year<20] #ADO only

save(rnra,file=here('LTBI/tmpdata/rnra.full.Rdata'))

rnras <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2)),
              by=.(iso3,age=year,mixing)]

rnras[,table(mixing)]
## P = ever
## P1 = any 1 year
## P2 = any infection within 2 years


## plot
rnras2 <- rnra[,.(`infection >=2 years`=mean(P-P2),
                  `infection 1-2 years`=mean(P2-P1),
                  `infection  <1 years`=mean(P1)),by=.(iso3,age=year,mixing)]
rnrasm <- melt(rnras2,id=c('iso3','age','mixing'))


plt <- ggplot(rnrasm[mixing=='assortative'],
              aes(as.factor(age),value,fill=variable))+
  geom_bar(stat='identity')+
  scale_y_continuous(label=percent)+
  scale_fill_colorblind()+
  ylab('TBI prevalence')+xlab('Age (years)')+
  facet_wrap(~iso3)+
  theme_light()+
  theme(legend.position = 'top',legend.title = element_blank())
## plt

ggsave(plt,file=here('plots/LTBI.jpg'),w=20,h=15)
ggsave(plt,file=here('plots/LTBI.pdf'),w=20,h=15)

## check
rnrc <- dcast(rnrasm[variable=='infection  <1 years'],iso3 + age ~ mixing,value.var = 'value')
rnrc[age>=15,assortative/random] #TODO check country-spec RRs?

## include age categories & save
rnra[,acat:='10-14']
rnra[year>14,acat:='15-19']
table(rnra$acat) #OK
rnra$acat <- factor(rnra$acat,levels=c('10-14','15-19'),ordered = TRUE)

save(rnra,file=gh('LTBI/data/rnra.Rdata'))


## ## comparing the change in LTBI due to mixing
## tmp <- rnra[acat=='15-19',.(Iold_mix=mean(P-P2),Inew_mix=mean(P2),
##                             Iold_nomix=mean(P0-P20),Inew_nomix=mean(P20)),
##                by=.(iso3,replicate)]

## tmp <- tmp[,.(Iold_mix_mid=mean(Iold_mix),Inew_mix_mid=mean(Inew_mix),
##               Iold_nomix_mid=mean(Iold_nomix),Inew_nomix_mid=mean(Inew_nomix),
##               Iold_mix_sd=sd(Iold_mix),Inew_mix_sd=sd(Inew_mix),
##               Iold_nomix_sd=sd(Iold_nomix),Inew_nomix_sd=sd(Inew_nomix)),by=iso3]

## tmp <- melt(tmp,id='iso3')
## tmp[,quantity:=ifelse(grepl('mid',variable),'mid','sd')]
## tmp[,mixing:=ifelse(grepl('nomix',variable),'without mixing RR','with mixing RR')]
## tmp[,variable:=ifelse(grepl('old',variable),'infection >2 years','infection <= 2 years')]
## tmp <- dcast(tmp,iso3+variable~quantity+mixing)
## names(tmp) <- gsub('mid_','',names(tmp))

## GP <- ggplot(data=tmp,aes(`with mixing RR`,`without mixing RR`,label=iso3))+
##   geom_point()+
##   facet_wrap(~variable,scales='free')+
##   geom_abline(intercept=0,slope=1,col=2)+
##   geom_errorbar(aes(ymin=`without mixing RR`-`sd_without mixing RR`,
##                     ymax=`without mixing RR`+`sd_without mixing RR`),width=0)+
##   geom_errorbarh(aes(xmin=`with mixing RR`-`sd_with mixing RR`,
##                      xmax=`with mixing RR`+`sd_with mixing RR`),height=0)+
##   geom_text_repel()

## ggsave(GP,file=here('plots/LTBIcf.jpg'),w=20,h=15)
## ggsave(GP,file=here('plots/LTBIcf.pdf'),w=20,h=15)


## tmp[,median(`with mixing RR`/`without mixing RR`),by=variable]
## ## variable       V1
## ## <char>    <num>
## ## 1: infection <= 2 years 1.211151
## ## 2:   infection >2 years 1.019446

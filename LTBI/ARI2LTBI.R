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
rd1 <- function(x) round(x,digits = 1)
rd(1); rd(1234); rd(1e9)
fmt <- function(x,y,z) paste0(rd(x)," (",rd(y)," to ",rd(z),")")
fmt1 <- function(x,y,z) paste0(rd1(x)," (",rd1(y)," to ",rd1(z),")")
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
fmtb <- function(x,y,z) paste0(rdb(x)," (",rdb(y)," to ",rdb(z),")")
gh <- function(x) glue(here(x))
cbPalette <- c("#000000", "#E69F00", "#56B4E9","#009E73",
               "#F0E442", "#0072B2","#D55E00", "#CC79A7")
load(here('progression/data/ckey.Rdata'))

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

## sex ratios?
SR <- fread(here('rawdata/TB_burden_age_sex_2020-10-15.csv'))
SR <- SR[age_group=='15plus' & risk_factor=='all' & sex!='a',.(iso3,sex,best)]
SR <- dcast(data=SR,iso3~sex,value.var='best')
SR[,sr:=(m+1)/(f+1)]
av <- SR[,sum(m)/sum(f)]
SR[m+f<1e2,sr:=av] #average for countries where this is small
SR <- SR[,.(iso3,sr)]

## work on converting this to LTBI
rnr[,ari:=exp(lari)]          #true ARI
rnr[,year:= 2019 - year]       #age
rnr <- rnr[order(replicate,iso3,year)] #order
rnr <- rnr[year<20]

## merge in sex stuff
## KH review: M  56% (IQR 54%–58%), F 59 (IQR 57%–63%)
## m=56.5 IQR 54-63 -> SD=(63-54)/1.35
## HEdtree::getAB(0.565,((63-54)/135)^2)
rnr <- merge(rnr,SR,by='iso3',all.x=TRUE)
rnr[,assort:=rbeta(nrow(rnr),30.67915,23.62023)]
rnr[,Rf:=(assort/(1+sr) + (1-assort)/(1+1/sr))*2]
rnr[,Rm:=((1-assort)/(1+sr) + assort/(1+1/sr))*2]

## better way to do below calx?
## year = years ago
rnr <- rnr[order(replicate,iso3,year)] #order

## loop over ages
tmpl <- list()
for(a in 10:19){
  rnr[,rrt:=1.0]
  rnr[,c('rrtf','rrtm'):=1.0]
  rnr[year<=a-15,rrt:=rr,by=.(iso3,replicate)] #RR applies only to years lived over 15
  rnr[year<=a-15,rrtf:=Rf,by=.(iso3,replicate)] #RR applies only to years lived over 15
  rnr[year<=a-15,rrtm:=Rm,by=.(iso3,replicate)] #RR applies only to years lived over 15
  tmp <- rnr[year<a,.(random.H=sum(ari),
                      random.dH2=sum(ari[1:2]),
                      random.dH1=sum(ari[1]),
                      assortative.H=sum(ari*rrt),
                      assortative.dH2=sum(ari[1:2]*rrt[1:2]),
                      assortative.dH1=sum(ari[1]*rrt[1]),
                      ## males
                      assortative.H_m=sum(ari*rrt*rrtm),
                      assortative.dH2_m=sum(ari[1:2]*rrt[1:2]*rrtm[1:2]),
                      assortative.dH1_m=sum(ari[1]*rrt[1]*rrtm[1]),
                      ## females
                      assortative.H_f=sum(ari*rrt*rrtf),
                      assortative.dH2_f=sum(ari[1:2]*rrt[1:2]*rrtf[1:2]),
                      assortative.dH1_f=sum(ari[1]*rrt[1]*rrtf[1])),
             by=.(iso3,replicate)]
  tmp[,age:=a]
  tmpl[[a]] <- tmp
}
tmpl <- rbindlist(tmpl)


## reshape
tmp <- melt(tmpl,id=c('iso3','replicate','age'))
tmp[,c('mixing','variable'):=tstrsplit(variable,'\\.')]
tmp <- dcast(tmp,iso3+mixing+replicate+age~variable,value.var = 'value')
tmp <- tmp[order(mixing,replicate,iso3,age)]
tmp[is.na(dH1_f),c('dH1_f','dH1_m'):=dH1] #random same by sex
tmp[is.na(dH2_m),c('dH2_f','dH2_m'):=dH2]
tmp[is.na(H_m),c('H_f','H_m'):=H]
tmp[,dH2:=H-dH2] #NOTE dH is actually cumulative EXCEPT last 2 years in calx
tmp[,c('dH2_f','dH2_m'):=.(H_f-dH2_f,H_m-dH2_m)] #as above by sex
tmp[,dH1:=H-dH1] #NOTE dH is actually cumulative EXCEPT last 1 years in calx
tmp[,c('dH1_f','dH1_m'):=.(H_f-dH1_f,H_m-dH1_m)] #as above by sex
rnr <- tmp

## convert to prevalence
rnr[,P:=1-exp(-H)]                  #ever
rnr[,c('f.P','m.P'):=.(1-exp(-H_f),1-exp(-H_m))] #as above by sex
rnr[,Pf1:=-exp(-H)+exp(-dH1)]         #1st recent=prob ever - prob not<1
rnr[,Pf2:=-exp(-H)+exp(-dH2)]         #1st recent=prob ever - prob not<2
rnr[,c('f.Pf1','m.Pf1'):=.(-exp(-H_f)+exp(-dH1_f),-exp(-H_m)+exp(-dH1_m))] #by sex
rnr[,Pf2:=-exp(-H)+exp(-dH2)]         #1st recent=prob ever - prob not<2
rnr[,c('f.Pf2','m.Pf2'):=.(-exp(-H_f)+exp(-dH2_f),-exp(-H_m)+exp(-dH2_m))] #by sex

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
rnr[,c('f.P1','m.P1'):=.(alph*(H_f-dH1_f) + (1-alph)*(exp(-dH1_f)-exp(-H_f)),
                         alph*(H_m-dH1_m) + (1-alph)*(exp(-dH1_m)-exp(-H_m)))]
rnr[,P2:=alph*(H-dH2) + (1-alph)*(exp(-dH2)-exp(-H))]      #anyrecent <2
rnr[,c('f.P2','m.P2'):=.(alph*(H_f-dH2_f) + (1-alph)*(exp(-dH2_f)-exp(-H_f)),
                         alph*(H_m-dH2_m) + (1-alph)*(exp(-dH2_m)-exp(-H_m)))]
## rnr[,P20:=alph*(H0-dH0) + (1-alph)*(exp(-dH0)-exp(-H0))]      #anyrecent

## check
rnr[iso3=='AGO' & replicate==1]

## ## restrict to relevant age groups
## rnra <- rnr[year>=10 & year<20] #ADO only

rnra <- rnr

save(rnra,file=here('LTBI/tmpdata/rnra.full.Rdata'))

rnras <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2)),
              by=.(iso3,age,mixing)]

rnras[,table(mixing)]
## P = ever
## P1 = any 1 year
## P2 = any infection within 2 years

## some LTBI outputs by sex?
## plot
rnras2 <- rnra[,.(`infection >=2 years`=mean(P-P2),
                  `infection 1-2 years`=mean(P2-P1),
                  `infection  <1 years`=mean(P1)),by=.(iso3,age,mixing)]
rnrasm <- melt(rnras2,id=c('iso3','age','mixing'))
mxl <- levels(rnrasm$variable)
rnrasm$variable <- factor(rnrasm$variable,levels=mxl,ordered = TRUE)
rnrasm <- merge(rnrasm,ckey,by = 'iso3',all.x=TRUE)

plt <- ggplot(rnrasm[mixing=='assortative'],
              aes(as.factor(age),value,fill=variable))+
  geom_bar(stat='identity')+
  scale_y_continuous(label=percent)+
  scale_fill_manual(values=rev(cbPalette[1:3]))+
  ylab('TBI prevalence')+xlab('Age (years)')+
  facet_wrap(~newcountry)+
  theme_light()+
  theme(legend.position = 'top',legend.title = element_blank())
## plt

ggsave(plt,file=here('plots/LTBI.jpg'),w=20,h=15)
ggsave(plt,file=here('plots/LTBI.pdf'),w=20,h=15)

## check
rnrc <- dcast(rnrasm[variable=='infection  <1 years'],iso3 + age ~ mixing,value.var = 'value')
rnrc[age>=15,assortative/random] #check country-spec RRs
rnrc[age<15,assortative/random]

## include age categories & save
rnra[,acat:='10-14']
rnra[age>14,acat:='15-19']
table(rnra$acat) #OK
rnra$acat <- factor(rnra$acat,levels=c('10-14','15-19'),ordered = TRUE)

save(rnra,file=gh('LTBI/data/rnra.Rdata'))

## NOTE some of this is repeated in progression/incidence.R
## additional outputs
load(file=gh('LTBI/data/rnra.Rdata'))
load(gh('PAF/data/IRR.Rdata'))
popn <- unique(IRR[,.(iso3,acat,pop)])

rnra <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2)),by=.(iso3,mixing,replicate,acat)] #over ages
rnra <- merge(rnra,popn,by=c('iso3','acat'))
rnr <- rnra[,.(P=weighted.mean(P,pop),
               P1=weighted.mean(P1,pop),
               P2=weighted.mean(P2,pop),
               pop=sum(pop)),by=.(iso3,mixing,replicate)] #no age
rnrc <- rnr[,.(P=1e2*mean(P),
               P1=1e2*mean(P1),
               P2=1e2*mean(P2)),by=.(iso3,mixing)] #for country variation
rnrtota <- rnra[,.(P=1e2*weighted.mean(P,pop),
                   P1=1e2*weighted.mean(P1,pop),
                   P2=1e2*weighted.mean(P2,pop),
                   TBI=sum(P*pop),
                   TBI1=sum(P1*pop),
                   TBI2=sum(P2*pop)),by=.(mixing,replicate,acat)] #with age
rnrtot <- rnra[,.(P=1e2*weighted.mean(P,pop),
                   P1=1e2*weighted.mean(P1,pop),
                   P2=1e2*weighted.mean(P2,pop),
                   TBI=sum(P*pop),
                   TBI1=sum(P1*pop),
                   TBI2=sum(P2*pop)),by=.(mixing,replicate)] #w/o age

## outputs
rnrtotasPC <- rnrtota[,lapply(.SD,function(x) fmt1(mean(x),lo(x),hi(x))),
                    by=.(mixing,acat),.SDcols=c('P','P1','P2')]
rnrtotasN <- rnrtota[,lapply(.SD,function(x) fmtb(mean(x),lo(x),hi(x))),
                      by=.(mixing,acat),.SDcols=c('TBI','TBI1','TBI2')]
rnrtotsPC <- rnrtot[,lapply(.SD,function(x) fmt1(mean(x),lo(x),hi(x))),
                      by=.(mixing),.SDcols=c('P','P1','P2')]
rnrtotsN <- rnrtot[,lapply(.SD,function(x) fmtb(mean(x),lo(x),hi(x))),
                     by=.(mixing),.SDcols=c('TBI','TBI1','TBI2')]
ltbi.crange <- rnrc[,.(Pmin=min(P),Pmax=max(P),
                       Cmin=iso3[which.min(P)],Cmax=iso3[which.max(P)]),
                    by=mixing]

fwrite(rnrtotasPC,file=gh('outdata/ltbi.PC.age.csv'))
fwrite(rnrtotasN,file=gh('outdata/ltbi.tot.age.csv'))
fwrite(rnrtotsPC,file=gh('outdata/ltbi.PC.noage.csv'))
fwrite(rnrtotsN,file=gh('outdata/ltbi.tot.noage.csv'))
fwrite(ltbi.crange,file=gh('outdata/ltbi.crange.csv'))

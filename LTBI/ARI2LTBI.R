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
fmt <- function(x,y,z) paste0(rd(x)," (",rd(y)," - ",rd(z),")")
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
rnr[,rr:=ifelse(year>14,rr,1.0)]       #RR step at 15
rnr[,H:=cumsum(rr*ari),by=.(iso3,replicate)] #cumulative ARI, including RR
rnr[,H0:=cumsum(ari),by=.(iso3,replicate)] #cumulative ARI, not including RR

## distinguish past 2 years
mask <- rep(1,length(unique(rnr$year)))
mask[1:2] <- 0                          #all except last 2 years
rnr[,dH:=cumsum(rr*ari*mask),by=.(iso3,replicate)] #cumhaz!2y
rnr[,dH0:=cumsum(ari*mask),by=.(iso3,replicate)] #cumhaz!2y no mix
rnr[,P:=1-exp(-H)]                  #ever
rnr[,P1:=-exp(-H)+exp(-dH)]         #1st recent=prob ever - prob not<2
rnr[,P0:=1-exp(-H0)]                  #ever
rnr[,P10:=-exp(-H0)+exp(-dH0)]         #1st recent=prob ever - prob not<2

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
rnr[,P2:=alph*(H-dH) + (1-alph)*(exp(-dH)-exp(-H))]      #anyrecent
rnr[,P20:=alph*(H0-dH0) + (1-alph)*(exp(-dH0)-exp(-H0))]      #anyrecent

## restrict to relevant age groups
rnra <- rnr[year>=10 & year<20] #ADO only
rnras <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2),
                 P0=mean(P0),P10=mean(P10),P2=mean(P20)),
              by=.(iso3,age=year)]
## P = ever
## P1 = 1st recent=prob ever - prob not<2
## P2 = any infection within 2 years
## 0 post-pended -> without mixing RR

## plot
rnras2 <- rnra[,.(`infection >=2 years`=mean(P-P2),
                  `infection  <2 years`=mean(P2)),by=.(iso3,age=year)]
rnrasm <- melt(rnras2,id=c('iso3','age'))
plt <- ggplot(rnrasm,aes(as.factor(age),value,fill=variable))+
  geom_bar(stat='identity')+
  scale_y_continuous(label=percent)+
  scale_fill_colorblind()+
  ylab('LTBI prevalence')+xlab('Age (years)')+
  facet_wrap(~iso3)+
  theme_light()+
  theme(legend.position = 'top',legend.title = element_blank())
## plt

ggsave(plt,file=here('plots/LTBI.jpg'),w=20,h=15)
ggsave(plt,file=here('plots/LTBI.pdf'),w=20,h=15)


## include age categories & save
rnra[,acat:='10-14']
rnra[year>14,acat:='15-19']
table(rnra$acat) #OK
rnra$acat <- factor(rnra$acat,levels=c('10-14','15-19'),ordered = TRUE)

save(rnra,file=gh('LTBI/data/rnra.Rdata'))


## comparing the change in LTBI due to mixing
tmp <- rnra[acat=='15-19',.(Iold_mix=mean(P-P2),Inew_mix=mean(P2),
                            Iold_nomix=mean(P0-P20),Inew_nomix=mean(P20)),
               by=.(iso3,replicate)]

tmp <- tmp[,.(Iold_mix_mid=mean(Iold_mix),Inew_mix_mid=mean(Inew_mix),
              Iold_nomix_mid=mean(Iold_nomix),Inew_nomix_mid=mean(Inew_nomix),
              Iold_mix_sd=sd(Iold_mix),Inew_mix_sd=sd(Inew_mix),
              Iold_nomix_sd=sd(Iold_nomix),Inew_nomix_sd=sd(Inew_nomix)),by=iso3]

tmp <- melt(tmp,id='iso3')
tmp[,quantity:=ifelse(grepl('mid',variable),'mid','sd')]
tmp[,mixing:=ifelse(grepl('nomix',variable),'without mixing RR','with mixing RR')]
tmp[,variable:=ifelse(grepl('old',variable),'infection >2 years','infection <= 2 years')]
tmp <- dcast(tmp,iso3+variable~quantity+mixing)
names(tmp) <- gsub('mid_','',names(tmp))

GP <- ggplot(data=tmp,aes(`with mixing RR`,`without mixing RR`,label=iso3))+
  geom_point()+
  facet_wrap(~variable,scales='free')+
  geom_abline(intercept=0,slope=1,col=2)+
  geom_errorbar(aes(ymin=`without mixing RR`-`sd_without mixing RR`,
                    ymax=`without mixing RR`+`sd_without mixing RR`),width=0)+
  geom_errorbarh(aes(xmin=`with mixing RR`-`sd_with mixing RR`,
                     xmax=`with mixing RR`+`sd_with mixing RR`),height=0)+
  geom_text_repel()

ggsave(GP,file=here('plots/LTBIcf.jpg'),w=20,h=15)
ggsave(GP,file=here('plots/LTBIcf.pdf'),w=20,h=15)


tmp[,median(`with mixing RR`/`without mixing RR`),by=variable]
## variable        V1
## 1: infection <= 2 years 0.9943804
## 2:   infection >2 years 1.0403761

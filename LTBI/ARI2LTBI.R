library(here)
library(data.table)
library(ggplot2)
library(scales)
library(ggthemes)
library(glue)

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

## work on converting this to LTBI
rnr[,ari:=exp(lari)]          #true ARI
rnr[,year:= 2019 - year]       #age
rnr <- rnr[order(replicate,iso3,year)] #order
rnr[,H:=cumsum(ari),by=.(iso3,replicate)] #cumulative ARI

## distinguish past 2 years
mask <- rep(1,length(unique(rnr$year)))
mask[1:2] <- 0                          #all except last 2 years
rnr[,dH:=cumsum(ari*mask),by=.(iso3,replicate)] #cumhaz!2y
rnr[,P:=1-exp(-H)]                  #ever
rnr[,P1:=-exp(-H)+exp(-dH)]         #1st recent=prob ever - prob not<2

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

## restrict to relevant age groups
rnra <- rnr[year>=10 & year<20] #ADO only
rnras <- rnra[,.(P=mean(P),P1=mean(P1),P2=mean(P2)),by=.(iso3,age=year)]
## P = ever
## P1 = 1st recent=prob ever - prob not<2
## P2 = any infection within 2 years


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


## TODO include mixing in above

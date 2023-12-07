library(here)
library(data.table)
library(ggplot2)
library(ggrepel)


## https://ncdrisc.org/data-downloads-adiposity-ado.html
D <- fread(here('NCD_RisC_Lancet_2020_BMI_child_adolescent_country.csv')) #NOTE not in repo as large


D[,table(Year,is.na(`Prevalence of BMI<minus2SD lower 95% uncertainty interval`))]

DR <- D[Year==2019] #so the distribution info is available before 2017

## questions:
## 1. How does the mean BMI compare to the data we have?
## 2. How big are the shifts in mean BMI from 2017 to 2019?

## question 1: comparison (BTW our data is 2016)
DO <- readxl::read_excel(here('rawdata/Adolescent Thiness Data.xlsx'),skip=0)
DO <- as.data.table(DO)
DO <- DO[-1,]
names(DO)[2] <- 'thin'
tmp <- D[Year %in% c(2016),.(Year,Sex,age=`Age group`,Country,
                             thin2=100*`Prevalence of BMI<minus2SD (moderate & severe underweight)`)]
tmp <- tmp[,.(thin2=mean(thin2)),by=.(Country)]
tmp <- merge(tmp,DO,by='Country')

ggplot(tmp,aes(thin,thin2,label=Country)) +
  geom_point() + geom_text_repel()+
  geom_abline(slope = 1,intercept = 0,col=2)+
  xlab('UNICEF 2016 data') + ylab('Average of 2016 NCD RisC estimates')+
  ggtitle('Comparison of percentage >= 2SD low-BMI')

ggsave('plots/BMIcompareChange.png',w=7,h=7)

## question 2: BMI shifts?
tmp <- D[Year %in% c(2017,2019),.(Year,Sex,age=`Age group`,Country,mBMI=`Mean BMI`)]
tmp <- dcast(tmp,Country+Sex+age ~ Year,value.var='mBMI')
tmp[,pcnt:=100*(`2019`/`2017`-1)]
tmp[,summary(pcnt)] #mean 0.4%, max 2.1%

tmp[,qplot(pcnt)] + xlab('Percentage difference') + ggtitle('Comparison of mean BMI 2019 vs 2016 (all ages, sexes, countries)')
names(DR)
DR[,length(unique(Country))]

## fit BMI distributiosn to this data (gamma distribution)
nnmz <- c('m.BMI','h2.BMI','h1.BMI','l1.BMI','l2.BMI')
anmz <- c(c('Country','Year','Sex','Age group'),nnmz)
tmp <- D[Country=='Afghanistan' & Year == 2016]
names(tmp)[c(5,9,12,15,18)] <- nnmz
tmp <- tmp[,..anmz]


## want to know: BMIs for z scores -2,-1,+1,+2 in ref pop for sex and age
## from: https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators/bmi-for-age
B <- readxl::read_excel(here('rawdata/bmi-boys-z-who-2007-exp.xlsx')) #
G <- readxl::read_excel(here('rawdata/bmi-girls-z-who-2007-exp.xlsx')) #
B <- as.data.table(B)
G <- as.data.table(G)
G[,Sex:='Girls']
B[,Sex:='Boys']
BG <- rbind(B,G)

gammaStats <- function(k,theta,
                       b2112 #BMI refs for +2,+1,-1,-2 SD
                       ){
  M <- k * theta
  c(M,                                        #mean
    1-pgamma(b2112[1:2],shape=k,scale=theta), #high tails
    pgamma(b2112[3:4],shape=k,scale=theta)    #low tails
    )
}


ba5 <- unlist(B[Month==5*12+6,.(SD2,SD1,SD1neg,SD2neg)])
gammaStats(2,9,ba5)

## loss function
gamErr <- function(x,tgt,ref){
  x <- exp(x)
  res <- gammaStats(x[1],x[2],ref)
  mean((res/(tgt+1e-10)-1)^2) #~SSE relative, equal weight
}

## test
TGT <- unlist(tmp[1,..nnmz])
gamErr(c(0,0),TGT,ba5)
## optimize
(out <- optim(par=c(0,0),fn=function(x) gamErr(x,TGT,ba5)))
exp(out$par)
gammaStats(exp(out$par[1]),exp(out$par[2]),ba5)
TGT
## check
curve(dgamma(x,shape=exp(out$par[1]),scale=exp(out$par[2])),from=10,to=30,n=1e3)


##  apply this fitting across the data:
D

## harmonize names
setdiff(DO$Country,unique(D$Country))
grep('North Korea',unique(D$Country),value=TRUE)
grep('Tanzania',unique(D$Country),value=TRUE)
grep('Papua',unique(D$Country),value=TRUE)
## rename
D[grepl('North Korea',Country),Country:="DPR Korea"]
D[grepl('Tanzania',Country),Country:="UR Tanzania"]
D[grepl('Papua',Country),Country:="PNG"]

## merge
DR <- D[Country %in% DO$Country & Year == 2016]
unique(DR$Country); nrow(DO) #check
## rename & restrict
names(DR)[c(5,9,12,15,18)] <- nnmz
DR <- DR[,..anmz]
DR[,Month:=12*`Age group`] #low-point (otherwise miss 19-20)
DR[Month==60,Month:=61]    #otherwise missing from BG
DR[,range(Month)]; BG[,range(Month)]
akey <- unique(DR[,.(age=`Age group`,Month)])
save(akey,file=here('rawdata/akey.Rdata'))

## merge:
DR <- merge(DR,BG,by=c('Sex','Month'),all.x = TRUE)

## Loop above work flow
DRA <- DR[,{
  print(Country);print(Sex);print(Month);
  TGT <- c(m.BMI,  h2.BMI, h1.BMI, l1.BMI, l2.BMI); #targets in NCD RisC data
  saref <- c(SD2, SD1, SD1neg, SD2neg);             #BMIs in WHO reference data
  out <- optim(par=c(0,0),fn=function(x) gamErr(x,TGT,saref));
  if(abs(out$convergence)>0) print('*** has not converged! ***')
  list(k=exp(out$par[1]),theta=exp(out$par[2]))
},by=.(Country,Month,Sex)]

DRA[,age:=as.integer(Month/12)]
DRA[,table(age,Month)]

## checks
DRA[,mn:=k*theta]
DRA[,summary(mn)]
DRA[,qplot(mn)]

## convergence problem for:
## Angola Girls 61
## SLE Girls 156
## Congo Girls
## RUS girls 204
## THA G 204
## Congo G 216
## Libera girls 228


DRA[mn>25] #PNG boys 120

## inspect problems
bad <- paste(c('Angola','Sierra Leone','Russian Federation','Thailand','Congo','Libera'),
             rep('Girls',7),
             c(61,156,204,204,216,228),
             sep='_')
bad <- c(bad,'PNG_Boys_120')

DRA[,id:=paste(Country,Sex,Month,sep='_')]
DR[,id:=paste(Country,Sex,Month,sep='_')]

DRA[id %in% bad]
DR[id %in% bad]

## replace these units with age/sex averages
DRAM <- DRA[,.(km=median(k),thetam=median(theta)),by=.(Sex,Month)]
DRAM[,km*thetam,by=.(Sex,Month)] #check
DRA <- merge(DRA,DRAM,by=c('Sex','Month'))
DRA[id %in% bad,c('k','theta'):=.(km,thetam)] #replace

## check
DRA[,mn:=k*theta]
DRA[mn>25] #PNG boys 120
DRA[,summary(mn)] #OK

## drop extra vars
DRA[,c('mn','Month','km','thetam','id'):=NULL]

## save
save(DRA,file=here('rawdata/DRA.Rdata'))


## BG refs
load(file=here('rawdata/DRA.Rdata'))
load(file=here('rawdata/akey.Rdata'))

BG <- BG[Month %in% akey$Month]
BG <- merge(BG,akey,by = 'Month')


pnorm(5) #mapping z to cumulative
qgamma(0.5,2,scale=9) #mapping cumuilative to x

gz <- function(z,k,theta)qgamma(pnorm(z),shape=k,scale=theta)
gz(c(-2,-1,0,1,2),2,9)


gamErr2 <- function(x,tgt){
  x <- exp(x)
  res <- gz(c(-2,-1,0,1,2),k=x[1],theta=x[2])
  mean((res/(tgt+1e-10)-1)^2) #~SSE relative, equal weight
}


BMIREF <- BG[,{
  print(Sex);print(Month);
  WHOref <- c(SD2neg, SD1neg, SD0,SD1, SD2);             #BMIs in WHO reference data
  out <- optim(par=c(0,0),fn=function(x) gamErr2(x,WHOref));
  if(abs(out$convergence)>0) print('*** has not converged! ***')
  list(kref=exp(out$par[1]),thetaref=exp(out$par[2]))
},by=.(Sex,age)]


## check
BMIREF[,kref*thetaref]

## save
save(BMIREF,file=here('rawdata/BMIREF.Rdata'))

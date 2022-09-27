library(here)
library(glue)
library(ggplot2)
library(data.table)
library(ggtext)
gh <- function(x) glue(here(x))
setDTthreads(1) #if parallelizing above
## gh('bla{cn}')

## setwd('../LTBIest')

## Analysis code for Houben & Dodd 2016, distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
args <- commandArgs(TRUE)
jj <- as.numeric(args[1])
print(jj)
library(MASS)
library(Matrix)
lin <- 0    #linear or constant  --- CHANGE HERE!


## ------------- dataprep
## now prevalence in same fashion
## log-normal: mu=1.678 sig=0.371
mu <- 1.678
sig <- 0.371
F <- exp(mu + .5*sig^2)                                  #styblo factor!! (lnorm mn)
vF <- (exp(sig^2)-1)*exp(2*mu+sig^2)                     #log normal variance

sqrt(vF)/F

## fix country
load(here('rawdata/ckey.Rdata'))
cn <- ckey[jj,iso3]
cnl <- ckey[jj,ihme]
print(c(cn,cnl))

## data input
load(gh('LTBI/data/K.Rdata'))

## NOTE see data prep
HM <- fread(gh('rawdata/IHME-GBD_2019_DATA-7d7d66ed-1.csv'))
HM <- HM[grepl('susceptible',cause_name) & grepl(cnl,location_name),.(year,val,upper,lower)]

HM[,prevSD:=(upper-lower)/3.92e5]
HM[,iso3:=cn]


HM <- merge(HM,K,by='iso3',all.x = TRUE,all.y=FALSE)
HM[,c('S','vlS'):=.(1,0)]

HM[,ari:=val*1e-5 * mFr * S * F]
HM[,lari:=log(ari)]
HM[,E1:=(prevSD)/(val/1e5)] #sd/mean from prev
HM[,E:=sqrt(E1^2+vF/F^2 + vlS + vFr/mFr^2)]
## HM[,ariSD:=ari*sqrt(vF/F^2+(val/1e5)^2/(prevSD)^2)]

load('~/Dropbox/Documents/Rwork/GPLTBI/LTBIest/data/All_3.Rdata') #old version for comparison in India & for TST svys
All <- All[All$lari!=-Inf,]
All <- All[All$iso3 %in% cn,]
HM[,type:='Prevalence estimate']
HM[,src:='IHME']
All$src <- 'WHO (old)'
CF <- rbind(HM[,.(iso3,year,lari,E,type,src)],
            All[,c('iso3','year','lari','E','type','src')])
CF[grep('TST',type),src:='Survey']


## NOTE just looking at India
if(cn=='IND'){

  ## adding real prevalence survey NOTE india only
  RPS <- data.table(year=2019,val=316,lower=296,upper=337)

  HMC <- ggplot(data=HM,
                aes(year,y=val,ymin=lower,ymax=upper))+
    geom_pointrange()+expand_limits(y=0)+
    geom_point(data=RPS,col="red",size=3)+
    geom_errorbar(data=RPS,aes(ymin=lower,ymax=upper),col="red")+
    xlab('Year') + ylab('TB prevalence per 100,000')+
    labs(title =  "<span style = 'color: red;'>Prevalence survey</span> vs IHME estimates")+
    theme_classic()+ggpubr::grids() + theme(plot.title = element_markdown())
  ggsave(HMC,file=gh('plots/IND_TBPScf.pdf'),w=9,h=5)
  ggsave(HMC,file=gh('plots/IND_TBPScf.png'),w=9,h=5)

  GGP <- ggplot(data=CF,
                aes(year,y=lari,col=src,shape=type,
                    ymin=lari-E,ymax=lari+E))+
    geom_pointrange(position = position_dodge(width=0.5))+
    theme_classic()+ggpubr::grids()+xlab('Year')+ylab('Log(ARI)')

  ggsave(GGP,file=gh('plots/IND_data.pdf'),w=9,h=5)
  ggsave(GGP,file=gh('plots/IND_data.png'),w=9,h=5)


} #NOTE end looking at IND only




## new version of All that uses IHME + survey
All <- CF[!grep('old',src)]
All <- All[order(year)]
All <- as.data.frame(All)

## =========== regression ==============================

getKKtonly <- function(t1,t2,k=function(x,y) exp(-abs(x-y)),Wm){
    K <- outer(t1,t2,FUN=k)
    K
}                                       #make K(X,Y) matrices

## function to do double-backsolve inversion given cholesky decomp
bsi <- function(CH,v) backsolve(CH,backsolve(CH, v, transpose = TRUE))

## this is basically 1,t x I(country==j)
getHtonly <- function(t,n=1){             #n is highest power
    H <- Matrix(0,nrow=(n+1),ncol=length(t)) #3 for 1,t,t^2
    for(k in 0:n)
        H[1 + k,] <- t^k
    as(H,'sparseMatrix')
}

tr <- function(x) sum(diag(x))

## gets s_k^2, L, tscale
PRMconvert2 <- function(x) c( exp(x[1]), exp(x[2]/2))

getMLnGradT <- function(x,grad=TRUE){                    # see Eqn 2.45 in R&W
    ## preliminaries
    ## -- covariances --
    a <- x[1]; b <- x[2]
    K <- outer(tdz,tdz,FUN=function(x,y) exp(a-exp(-b)*(x-y)^2) )
    K <- Matrix(K)                      #kxx
    K2 <- outer(tdz,tdz,FUN=function(x,y) (x-y)^2*exp(a-exp(-b)*(x-y)^2) )
    K2 <- Matrix(K2)                     #dK/da
    ## -- derived matrices (as used above) --
    ## new version
    cvy <- K + sigz
    U <- chol(cvy)
    Uy <- backsolve(U, y, transpose = TRUE)
    ky <- backsolve(U, Uy) 
    hky <- H %*% ky
    AM <- symmpart(H %*% bsi(U,t(H)))
    V <- chol(AM)
    Vy <- backsolve(V, hky, transpose = TRUE)
    ## -- marginal loglikelihood --
    LML <- -sum(Uy^2) + sum(Vy^2) #data likelihood
    LML <- LML - 2*sum(log(diag(U)))
    LML <- LML - 2*sum(log(diag(V)))
    LML <- LML/2
    ## extras for dLML
    VVy <- backsolve(V, Vy)
    ympy <- bsi(U, t(H) %*% VVy)       #K^{-1} x ...
    dLML <- NULL
    if(grad){
        ## -- gradient --
        ## --- gradient helper function ---
        dHelp <- function(dK){              #takes the local vars from up-level
            dML <- t(ky) %*% dK %*% ky + t(ympy) %*% dK %*% ympy
            dML <- dML - 2*t(ympy) %*% dK %*% ky
            tmp <- bsi(U,dK)        #K^{-1}dK
            dML <- dML - tr(tmp)
            tmp <- bsi(V,H%*%tmp)
            tmp <- bsi(U, t(H)%*%tmp)
            dML <- dML + tr(tmp)  
            return(as.numeric(dML/2))
        }
        ## --- get gradient ---
        dLML <- c( dHelp(K), dHelp(K2) )
    }
    if(LML>0){LML <- -1e3; dLML <- -1e3*rep(1,2)}
    ## return
    return(list(LML=LML,dLML=dLML))
}

getPredztonly <- function(x,tdz,tez,y,Vz){
    ## from here
    usek <- function(i,j)x[1]*exp(-abs(i-j)^2/x[2]^2)
    ## H
    H <- getHtonly(tdz,n=lin)
    Hs <- getHtonly(tez,n=lin)
    ## make matrices
    kxx <- getKKtonly(tdz,tdz,k=usek)
    kxxs <- getKKtonly(tdz,tez,k=usek)
    kxsx <- t(kxxs)
    kxsxs <- getKKtonly(tez,tez,k=usek)
    sigz <- Diagonal(x=Vz)                         #noise
    covy <- kxx + sigz
    U <- chol(covy)
    ## regress
    reg <- 1                                #flag
    HKH <- H %*% bsi(U,t(H))
    V <- chol(symmpart(HKH))
    R <- Hs - H %*% bsi(U,kxxs)
    mn <- 0
    y <- y-mn
    ## mean/covar
    mf <- kxsx %*% bsi(U,y)          #mean prediction
    cf <- kxsxs  - kxsx %*% bsi(U,kxxs)
    ## mean stuff
    bbar <- bsi(V,(H %*% bsi(U,y)))
    mg <- mf + reg * t(R) %*% bbar
    if(nrow(V)>1)
        cg <- cf + reg * t(R) %*% bsi(V,R)
    else
        cg <- cf + reg * t(R) %*% (R/V[1,1])
    ## return
    return(list(mg=mg,cg=cg))
}


## ============== work ====================

## time/country vectors for data/extrapolation
fyear <- 1933
tdz <- All$year-fyear
tez <- 1934:2019 - fyear    #extrapolation times: now 2019 not 2014

## H as a global
## dcnz, tdz, sigz, y as globals
## H
H <- getHtonly(tdz,n=lin)
Hs <- getHtonly(tez,n=lin)
y <- All$lari
vz <- (All$E)^2
sigz <- Diagonal(x=vz)                         #noise

mz <- c(log(.5),2*1.5*log(2));sz <- c(1,1)*100
LMLfun2 <- function(x) -(getMLnGradT(x,grad=FALSE)$LML
                         - sum(.5*(x-mz)^2/sz^2))
dLMLfun2 <- function(x) -(getMLnGradT(x)$dLML-(x-mz)/sz^2)
x02 <- mz
## optimize
system.time({                           #<1s
    testo2 <- optim(par=x02,fn = LMLfun2,gr = dLMLfun2)
})
pab <- testo2$par
ab <- PRMconvert2(pab)
## ab
print(ab);
## NB trade-off flexibility vs growing noise going back

xx <- ab
tot <- getPredztonly(xx,tdz,tez,y,vz)
scf <- as.numeric(sqrt(diag(tot$cg)))
erw <- data.frame(year=tez+fyear,iso3=as.character(unique(All$iso3)),
                  lari=as.numeric(tot$mg),
                  upper=as.numeric(tot$mg) + 1.96*scf,
                  lower=as.numeric(tot$mg) - 1.96*scf)

save(erw,file=gh('LTBI/tmpdata/erw_{cn}.Rdata'))

runs <- mvrnorm(n=2e2,mu=as.numeric(tot$mg),
                Sigma=as.matrix(symmpart(tot$cg)))
runsdf <- data.frame(year=tez+fyear,
                     iso3=as.character(unique(All$iso3)),
                     lari=c(t(runs)),
                     replicate=rep(1:nrow(runs),each=ncol(runs)))

save(runsdf,file=gh('LTBI/tmpdata/zz_{cn}.Rdata'))

erw$src <- NA
erw$type <- NA
erw$E <- NA


if(cn=='IND'){
  GGP2 <- GGP +
    geom_line(data=erw,aes(x=year,y=lari))+
    geom_line(data=erw,aes(x=year,y=lower),lty=2)+
    geom_line(data=erw,aes(x=year,y=upper),lty=2)

  ggsave(GGP2,file=gh('plots/IND_data_fit.pdf'),w=9,h=5)
  ggsave(GGP2,file=gh('plots/IND_data_fit.png'),w=9,h=5)
}

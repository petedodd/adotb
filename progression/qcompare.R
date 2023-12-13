library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
library(here)
load(here('progression/data/smy.Rdata'))
load(here('progression/data/ckey.Rdata'))

W <- fread('/Users/petedodd/Dropbox/gtb2020/Pete/MCEE/outdata/db_incidence_disaggregated_MCEEages_allyear0.csv')
H <- fread('/Users/petedodd/Downloads/IHME-GBD_2019_DATA-644cb085-1/IHME-GBD_2019_DATA-644cb085-1.csv')

H <- merge(H,ckey[,.(iso3,ihme)],by.x = 'location_name',by.y='ihme',all.x=TRUE)
H[,acat:=ifelse(age_name=='10-14 years','10-14','15-19')]

HS <- H[cause_name=='Drug-susceptible tuberculosis'] #DS only
W <- W[year==2019]
W <- W[iso3 %in% smy$iso3]
W <- W[age_group %in% c('10_14','15_19')]
W[,acat:=ifelse(age_group=='10_14','10-14','15-19')]
W <- W[,.(mcinc=sum(inc.nh.num)),by=.(iso3,acat)]
W

## comparison object
cf <- merge(smy[,.(iso3,acat,mixing,
                   inc.num.mid,inc.num.lo,inc.num.hi,
                   inc.num0.mid,inc.num0.lo,inc.num0.hi)],
            W,by=c('iso3','acat'),all.x=TRUE)
tmp <- HS[,.(iso3,acat,ihme=val,ihme.lo=lower,ihme.hi=upper)]
cf <- merge(cf,tmp,by=c('iso3','acat'),all.x=TRUE)


m <- 1e6*0.75
ggplot(cf[mixing=='assortative'],aes(inc.num.mid,mcinc,label=iso3))+
  geom_point()+
  geom_text_repel()+
  facet_wrap(~acat)+
  geom_abline(intercept = 0,slope=1,col=2)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  xlab('Infection->progression')+ylab('Snow/MCEE')

ggsave(file=here('plots/CF1.png'),h=7,w=14)

summary(cf)
cf[,mean(inc.num.mid/mcinc,na.rm = TRUE),by=.(mixing,acat)] #4x in 10-14; 1.1-1.2 x in 15-19



m <- 1e6*0.75
ggplot(cf[mixing=='assortative'],aes(inc.num.mid,ihme,label=iso3,
                                     ymin=ihme.lo,ymax=ihme.hi,
                                     xmin=inc.num.lo,xmax=inc.num.hi))+
  geom_errorbar(width=0)+geom_errorbarh(height=0)+
  geom_point()+
  geom_text_repel()+
  facet_wrap(~acat)+
  geom_abline(intercept = 0,slope=1,col=2)+
  scale_x_sqrt(limits=c(0,m),label=comma)+
  scale_y_sqrt(limits=c(0,m),label=comma)+
  xlab('Infection->progression')+ylab('IHME')

ggsave(file=here('plots/CF2.png'),h=7,w=14)


cf[,mean(inc.num.mid/ihme,na.rm = TRUE),by=.(mixing,acat)] #1.5x in 10-14; 0.75 x in 15-19
cf[,median(inc.num.mid/ihme,na.rm = TRUE),by=.(mixing,acat)] #1.5x in 10-14; 0.75 x in 15-19


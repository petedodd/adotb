library(here)
library(data.table)
library(ggplot2)
library(viridis)

## key to WHO regions
load(file=here('rawdata/whokey.Rdata'))


## --- mixing data
## contact data
load(here('rawdata/synthetic_contacts_2020.csv.Rdata'))
CD <- synthetic_contacts_2020[location_contact=="all",
                              .(iso3=iso3c,age_contactor,
                                age_contactee=age_cotactee,
                                ctx=mean_number_of_contacts)]

## --- TB estimates
## read in WHO age-specific incidence
E <- fread('~/Dropbox/Documents/WHO_TBreports/data2020/TB_burden_age_sex_2020-10-15.csv')
E[,unique(age_group)]
exa <- c('all','0-14',"15plus","18plus") #exlude age groups
E <- E[!age_group %in% exa & risk_factor=='all',
       .(iso3,sex,age_group,TB=best)] #choose right age groups
E <- E[,.(TB=sum(TB)),by=.(iso3,age_group)] #aggregate over sex

## --- WPP19 demography for 2020
load(file=here('rawdata/N80.Rdata'))

## age key
akey <- data.table(AgeGrp=unique(N80$AgeGrp))
akey[,age_group:=AgeGrp]
akey[age_group %in% c('65-69','70-74','75-79','80+'),age_group:='65plus']
akey[age_group %in% c('5-9','10-14'),age_group:='5-14']
akey[age_group %in% c('15-19','20-24'),age_group:='15-24']
akey[age_group %in% c('25-29','30-34'),age_group:='25-34']
akey[age_group %in% c('35-39','40-44'),age_group:='35-44']
akey[age_group %in% c('45-49','50-54'),age_group:='45-54']
akey[age_group %in% c('55-59','60-64'),age_group:='55-64']
akey[,acat:=gsub("plus","+",age_group)]
akey[,cage:=AgeGrp]
akey[AgeGrp %in% c('75-79','80+'),cage:='75+'] #for contacts

agz <- akey[,unique(acat)]


NS <- merge(N80,akey)
NS <- NS[,.(pop=sum(pop)),by=.(iso3,age_group)]


## join demo + TB
E <- merge(E,NS,by=c('iso3','age_group'),all.x=TRUE,all.y=FALSE)
E[,pcTB:=1e5*TB/pop] #per capita tb

E <- merge(E,unique(akey[,.(age_group,acat)]),by='age_group')
E[is.na(pcTB),pcTB:=0]
E$acat <- factor(E$acat,levels=agz,ordered=TRUE)
refs <- E[acat=='15-24',.(iso3,refpcTB=pcTB)]
E <- merge(E,refs,by=c('iso3'))
E[,relpcTB:=pcTB/refpcTB]
E[!is.finite(relpcTB),relpcTB:=NA]

## examine:
E <- merge(E,whokey,by='iso3') #merge region
ER <- E[,.(relpcTB=mean(relpcTB,na.rm = TRUE)),by=.(g_whoregion,acat)]
ER[,iso3:=g_whoregion]


GP <- ggplot(E[!acat %in% c('0-4','5-14')],
       aes(acat,relpcTB,group=iso3))+
  geom_line()+
  facet_wrap(~g_whoregion)+
  geom_hline(yintercept=1,col=2)+
  scale_y_sqrt()+
  ylab('Square-root of relative per capita TB incidence')+
  xlab('Age category')+
  ggtitle('Step 1: Relative per capita TB incidence (WHO estimates)')+
  theme_light()
## GP

ggsave(GP,file=here('plots/ARI_step1_percapTB.pdf'),w=12,h=7)
ggsave(GP,file=here('plots/ARI_step1_percapTB.png'),w=12,h=7)


## aggregating contacts
CD[,AO:=gsub(" to ","-",age_contactor)]
CD[,AI:=gsub(" to ","-",age_contactee)]
CD <- merge(CD,unique(akey[,.(cage,acato=acat)]),
            by.x = 'AO',by.y = 'cage',
            all.x=TRUE,all.y=FALSE
            )
CD <- merge(CD,unique(akey[,.(cage,acati=acat)]),
            by.x = 'AI',by.y = 'cage',
            all.x=TRUE,all.y=FALSE
            )
CD <- CD[,.(ctx=sum(ctx)),by=.(iso3,acato,acati)]
CD <- merge(CD,whokey,by='iso3') #regions

## summary(CD)


## inspect
CDR <- CD[,.(contacts=mean(ctx)),
          by=.(g_whoregion,acato,acati)] #regional ave
CDR[,acato:=factor(acato,levels=agz,ordered=TRUE)]
CDR[,acati:=factor(acati,levels=agz,ordered=TRUE)]


GP <- ggplot(data=CDR,aes(x=acato,y=acati,fill=contacts))+
  geom_tile()+
  scale_fill_viridis()+
  theme(legend.position='bottom')+
  facet_wrap(~g_whoregion)+
  xlab('Age of contactor')+
  ylab('Age of contactee')+
  ggtitle('Step 2: Regional average contact patterns')
## GP

ggsave(GP,file=here('plots/ARI_step2_contacts.pdf'),w=12,h=8)
ggsave(GP,file=here('plots/ARI_step2_contacts.png'),w=12,h=8)

## merge in contact data

EC <- merge(E[,.(iso3,acat,TB,pop,relpcTB,g_whoregion)],
            CD,
            by.x = c('iso3','g_whoregion','acat'),
            by.y = c('iso3','g_whoregion','acati'),
            all.x = TRUE, all.y=FALSE)
EC <- EC[!is.na(ctx)]

length(unique(EC$iso3)) #177 countries


## under 15 not infectious
EC[acat %in% c('0-4','5-14'),relpcTB:=0]
EC[,ari:= relpcTB * ctx]

EC <- EC[,.(ari=sum(ari)),by=.(g_whoregion,iso3,acat=acato)]
refs <- EC[acat=='5-14',.(iso3,refari=ari)]
EC <- merge(EC,refs,by='iso3',all.x = TRUE)
EC[,relari:=ari/refari]
EC <- EC[,.(iso3,g_whoregion,acat,relari)]
EC$acat <- factor(EC$acat,levels=agz,ordered=TRUE)

save(EC,file=here('LTBI/data/EC.Rdata'))


GP <- ggplot(EC,aes(acat,relari,group=iso3))+
  geom_line()+
  facet_wrap(~g_whoregion)+
  geom_hline(yintercept=1,col=2)+
  xlab('Age')+
  ylab('Relative ARI')+
  ggtitle('Step 3: Relative ARI implied by mixing')
## GP

ggsave(GP,file=here('plots/ARI_step3_relari.pdf'),w=12,h=8)
ggsave(GP,file=here('plots/ARI_step3_relari.png'),w=12,h=8)


## 12.5 vs 17.5
EC[,c('la','ua'):=tstrsplit(acat,"-")]
EC[,mida:=as.numeric(la) + (as.numeric(ua)-as.numeric(la)+1)/2]
EC <- EC[!is.na(relari)]

RR <- list()
for(cn in EC[,unique(iso3)]){
  tmp <- EC[iso3==cn & !is.na(mida)]
  rar <- approx(x=tmp$mida,y=tmp$relari,xout = c(12.5,17.5))$y
  RR[[cn]] <- data.table(iso3=cn,rr=rar[2]/rar[1])
}
RR <- rbindlist(RR)

save(RR,file=here('LTBI/data/RR.Rdata'))

## look
RR <- merge(RR,whokey,by='iso3',all.x = TRUE)




GP <- ggplot(RR,aes(iso3,rr))+
  geom_point()+
  facet_wrap(~g_whoregion,scales='free')+
  coord_flip()+
  geom_hline(yintercept=1,col=2)+
  xlab('Country')+
  ylab('Relative ARI')+
  ggtitle('RR age 17.5 vs 12.5')
## GP

ggsave(GP,file=here('plots/ARI_step4_RR.pdf'),w=12,h=10)
ggsave(GP,file=here('plots/ARI_step4_RR.png'),w=12,h=10)

##### mod5 only: conspecific v conmycorrhizal heterspecific:
##### mod2 also: conspecific rel to all heterospecific

#0. rich
#1. lat
#2. myc
#3. abund

# load packages
library(tidyverse)
library(lsmeans)
library(ggeffects)
library(lme4)
library(lmerTest)
library(gridExtra)
library(taxonlookup)
library(spatstat)
library(MCMCglmm)
library(phytools)
library(mgcv)  
library(pez)

# read in lat:
lat<- read.table("data/FGeo_LatLon.txt",header =TRUE) %>%
  mutate(site=toupper(site))%>%
  select(-site)%>%
  rename(site=full_site) %>%
  mutate(abslat=abs(lat))

#########################
###### GLM OR GAM #######
###### COEF OR PROP #####
######## Net A ##########
#########################

## CHOOSE ONE
#GLM coef net A
dat<-read.csv("data/output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.coef.A) %>% mutate(resp=round(resp,digits=3))
#GLM coef mod 2
dat<-read.csv("data/output.study.net.m1-2.A_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.coef.A) %>% mutate(resp=round(resp,digits=3))
#GLM prop net A
dat<-read.csv("data/output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))
#GLM prop mod 2
dat<-read.csv("data/output.study.net.m1-2.A_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))
#GAM prop net A (1-5)*****
dat<-read.csv("data/GAM.output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))
#GAM prop mod 2 (1-5)
dat<-read.csv("data/GAM.output.study.net.m1-5.A_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))
#GAM prop net A (1-2)
dat<-read.csv("data/GAM.output.study.net.m1-2.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))
#GAM prop mod 2 (1-2)
dat<-read.csv("data/GAM.output.study.net.m1-2.A_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A) %>% mutate(resp=round(resp,digits=3))

dat<-dat %>%left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(site=as.factor(site)) %>%  #for GAMs
  mutate(myc = fct_relevel(myc,"AM","EM"))

dat.remove<- dat %>% 
  group_by(site, myc,.drop=FALSE) %>% 
  distinct(latin) %>% 
  summarize(N = n())%>%
  mutate(freq = N / sum(N)) %>%
  filter(freq==1)

dat<- dat %>% filter(!site %in% dat.remove$site)   

#######################
### CONIFER SECTION ###
#######################

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

dat.sp<- unique(dat$latin)
dat.conif<- lookup_table(c(dat.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

dat<- dat %>% left_join(dat.conif,by="latin") %>%
  mutate(myc.conif= paste(myc,Coniferae)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM.NC", myc.conif)) %>% # replace with AM for NC
  mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM.C", myc.conif)) %>% # replace with AM for C
  mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM NA", "AM.NC",myc.conif))%>%
  select(-myc) %>%
  rename(myc=myc.conif) 

#######################
#######################

richdat5<-
  dat %>% 
  group_by(site) %>%                         
  summarize(sprich = length(latin)) 

latdat5<-dat %>% 
  group_by(site) %>%          
  summarise(medcndd = weighted.median(resp, lgabund)) %>%     
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  dat %>% 
  group_by(site,myc) %>%                       
  summarize(sprich = length(latin)) 

latdat5.myc<-dat %>% 
  group_by(site) %>%                     
  summarize(medcndd = weighted.median(resp, lgabund)) %>%     
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

latdat5.EM<- latdat5.myc %>%
  filter(myc=="EM")%>%
  rename(sprichEM=sprich) %>%
  left_join(richdat5, by="site") %>%
  mutate(propEM= sprichEM/sprich)

#0 C: GLM
lm.rich<-glm(medcndd~sprich, dat=latdat5)               
summary(lm.rich)
saveRDS(lm.rich,"cs_lmrich.rds")

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
#Fig1.A<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Species richness")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-1.25,0.50)

#0 C: GLM myc
lm.rich.myc<-glm(medcndd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"), raw=TRUE)
Fig1.B<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = lat, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Absolute latitude")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-1.25,0.50)

#1 C: GLM
lm.lat.myc<-glm(medcndd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#1 C: GLM version prop EM
lm.pmyc<-glm(medcndd~propEM, dat=latdat5.EM)               
summary(lm.pmyc)

datplot <- ggpredict(lm.pmyc, terms=c("propEM"))
Fig1.C<-ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5.EM, mapping = aes(x = propEM, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Proportion EM")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-1.25,0.50)

#2 C: GLM all sites
lm.myc<-lmer(resp~myc + (1|site), dat=dat)      
summary(lm.myc)
anova(lm.myc)

#2 C: GLM each site
lm.myc.s<-glm(resp~myc*site,  dat=dat)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= dat)
ref.table<-as.data.frame(ref$lsmeans) %>%
  left_join(latdat5, by = "site") %>%
  arrange(abslat) 

ggplot(ref.table, aes(x=reorder(myc, abslat), lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-10,5)

#2 C: MCMCGLM all sites
#load tree
load("data/GBOTB.extended.rda")
phy<-GBOTB.extended

#merge congenerics
phy.ex<-congeneric.merge(phy, dat$latin, split = "_", cite = FALSE)

#subset tree to observed species in data
alltreesp<- c(phy.ex$tip.label) %>% as.data.frame() 
colnames(alltreesp)<-"latin"
nodrop<- dat %>% select(latin)                                         
todrop<- setdiff(alltreesp, nodrop)%>% as.vector()                        
todrop<-c(t(todrop))                                                 
p.tree<- drop.tip(phy.ex, todrop)                              

#clean up the tree
p.tree$node.label <- NULL
phy$node.label<-NULL
rnd <- inverseA(p.tree)$Ainv # inverse of matrix

# subset dat to include only latin in tips of tree:
dat<-dat %>% filter(latin %in% p.tree$tip.label)

#define priors (phylo and site):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002), 
                                             G2 = list(V=1, nu=0.002))) 

lm.myc.mcmc<-MCMCglmm(resp~myc,random= ~latin+site,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.mcmc)

ref<-lsmeans(lm.myc.mcmc,pairwise~myc, data= dat)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lower.HPD, ymax=upper.HPD), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15), alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-5,1)

#2 C: MCMCGLM each site
#define priors (phylo only):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002)))

lm.myc.mcmc<-MCMCglmm(resp~myc*site,random= ~latin,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.mcmc)

ref<-lsmeans(lm.myc.mcmc,pairwise~myc*site, data= dat)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lower.HPD, ymax=upper.HPD), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15), alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-5,1)+
  facet_wrap(~site, ncol=5)

#3 C: GLM all sites
lm.myc.abund<-lmer(resp~myc*lgabund+(1|site),  dat=dat)
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
Fig3.A<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = resp),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-2.5,1)

#3 C: GAM all sites
lm.myc.abund.gam<-gam(resp~s(lgabund)+s(site, bs = 're'),  dat=dat) 
summary(lm.myc.abund.gam)

pred <- data.frame(resp= dat$resp,
                   lgabund = dat$lgabund,
                   predict(lm.myc.abund.gam,newdata = dat,exclude = 's(site)',se=TRUE))

ggplot(data = pred, mapping = aes(x = lgabund, y = fit))+
  geom_line() +
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = resp),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-2.5,1)

#3 C: GLM each site
lm.myc.abund.s<-glm(resp~myc*lgabund*site, dat=dat)
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

#3 C: MCMCGLM all sites
#define priors (phylo and site):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002), 
                                             G2 = list(V=1, nu=0.002))) 

lm.myc.abund.mcmc<-MCMCglmm(resp~myc*lgabund,random= ~latin+site,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.abund.mcmc)

#3 C: MCMCGLM each site
#define priors (phylo only):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002)))

lm.myc.abund.mcmc<-MCMCglmm(resp~myc*lgabund*site,random= ~latin,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.abund.mcmc)

#4 Visualize
# Density plot
ggplot(data=dat, mapping =aes(x=resp,color=myc, fill=myc))+
  geom_density(alpha=0.3) +
  xlab("Conspecific Recruitment")+
  ylab("Density") +
  xlim(-2,2) 

# Fig1: lat and propEM
grid.arrange(Fig1.A, Fig1.B, Fig1.C, ncol=3)

#########################
####### GLM OR GAM ######
###### COEF OR PROP #####
######## Net CM #########
#########################

#GLM coef net CM
dat<-read.csv("data/output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = resp) %>% mutate(resp=round(resp,digits=3))
#GAM prop net A (1-5)*****
dat<-read.csv("data/GAM.output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.CMH) %>% mutate(resp=round(resp,digits=3))

dat<-
  dat%>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(site=as.factor(site)) %>%  #for GAMs
  mutate(myc = fct_relevel(myc,"AM","EM"))

dat.remove<- dat %>% 
  group_by(site, myc,.drop=FALSE) %>% 
  distinct(latin) %>% 
  summarize(N = n())%>%
  mutate(freq = N / sum(N)) %>%
  filter(freq==1)

dat<- dat %>% filter(!site %in% dat.remove$site) 

#######################
### CONIFER SECTION ###
#######################

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

dat.sp<- unique(dat$latin)
dat.conif<- lookup_table(c(dat.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

dat<- dat %>% left_join(dat.conif,by="latin") %>%
  mutate(myc.conif= paste(myc,Coniferae)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM.NC", myc.conif)) %>% # replace with AM for NC
  mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM.C", myc.conif)) %>% # replace with AM for C
  mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM NA", "AM.NC",myc.conif))%>%
  select(-myc) %>%
  rename(myc=myc.conif) 

#######################
#######################

richdat5<-
  dat %>% 
  group_by(site) %>%                         
  summarize(sprich = length(latin)) 

latdat5<-dat %>% 
  group_by(site) %>% 
  summarise(medcmdd = weighted.median(resp, lgabund)) %>%     
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  dat %>% 
  group_by(site,myc) %>%                          
  summarize(sprich = length(latin)) 

latdat5.myc<-dat %>% 
  group_by(site) %>%                    
  summarize(medcmdd = weighted.median(resp, lgabund)) %>%    
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

latdat5.EM<- latdat5.myc %>%
  filter(myc=="EM")%>%
  rename(sprichEM=sprich) %>%
  left_join(richdat5, by="site") %>%
  mutate(propEM= sprichEM/sprich)

#0 C: GLM
lm.rich<-glm(medcmdd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcmdd),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Species Richness")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-1.25,0.50)

#0 C: GLM myc
lm.rich.myc<-glm(medcmdd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 CM: GLM: 
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"), raw=TRUE)
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = lat, y = medcmdd),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Absolute latitude")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-0.5,0.5)

#1 C: GLM
lm.lat.myc<-glm(medcmdd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#1 C: GLM version prop EM
lm.pmyc<-glm(medcmdd~propEM, dat=latdat5.EM)               
summary(lm.pmyc)

datplot <- ggpredict(lm.pmyc, terms=c("propEM"))
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5.EM, mapping = aes(x = propEM, y = medcmdd),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Proportion EM")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-0.5,0.5)

#2 CM: GLM all sites:NS
lm.myc<-lmer(resp~myc + (1|site),  dat=dat)      
summary(lm.myc)
anova(lm.myc,test="LRT")

ref<-lsmeans(lm.myc,pairwise~myc, data= dat)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 CM: GLM each site
lm.myc.s<-glm(resp~myc*site,  dat=dat)      
summary(lm.myc.s)
anova(lm.myc.s, test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= dat)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: MCMCGLM all sites

#load tree
load("data/GBOTB.extended.rda")
phy<-GBOTB.extended

#merge congenerics
phy.ex<-congeneric.merge(phy,dat$latin, split = "_", cite = FALSE)

#subset tree to observed species in data
alltreesp<- c(phy.ex$tip.label) %>% as.data.frame() 
colnames(alltreesp)<-"latin"
nodrop<- dat %>% select(latin)                                         
todrop<- setdiff(alltreesp, nodrop)%>% as.vector()                        
todrop<-c(t(todrop))                                                 
p.tree<- drop.tip(phy.ex, todrop)                              

#clean up the tree
p.tree$node.label <- NULL
phy$node.label<-NULL
rnd <- inverseA(p.tree)$Ainv # inverse of matrix

# subset dat to include only latin in tips of tree:
dat<-dat %>% filter(latin %in% p.tree$tip.label)

#define priors (phylo and site):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002), 
                                             G2 = list(V=1, nu=0.002))) 

lm.myc.mcmc<-MCMCglmm(resp~myc,random= ~latin+site,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.mcmc)

ref<-lsmeans(lm.myc.mcmc,pairwise~myc, data= dat)
ref.table<-as.data.frame(ref$lsmeans) %>% mutate(site="Overall") %>% mutate(abslat=0)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lower.HPD, ymax=upper.HPD), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=resp,color=myc),position=position_jitter(width =.15), alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-2.5,1)

#2 C: MCMCGLM each site
#define priors (phylo only):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002)))

lm.myc.mcmc<-MCMCglmm(resp~myc*site,random= ~latin,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.mcmc)

#3 C: GLM all sites
lm.myc.abund<-lmer(resp~myc*lgabund+(1|site), dat=dat)
summary(lm.myc.abund)
anova(lm.myc.abund,test="LRT")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
Fig3.B<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = resp),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("log(abundance)")+
  ylim(-2.5,1)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=dat,aes(x=lgabund,y=resp,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("log(abundance)")+
  ylim(-5,5)

#3 C: GAM all sites
lm.myc.abund.gam<-gam(resp~s(lgabund)+s(site, bs = 're'),  dat=dat) 
summary(lm.myc.abund.gam)

pred <- data.frame(resp= dat$resp,
                   lgabund = dat$lgabund,
                   predict(lm.myc.abund.gam,newdata = dat,exclude = 's(site)',se=TRUE))

ggplot(data = pred, mapping = aes(x = lgabund, y = fit))+
  geom_line() +
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = resp),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed") +
  ylim(-2.5,1)


dat.AM<-dat%>%filter(myc=="AM")
lm.AM.abund<-gam(resp~s(lgabund)+s(site, bs = 're'),  dat=dat.AM)
summary(lm.AM.abund)

pred <- data.frame(resp= dat.AM$resp,
                   lgabund = dat.AM$lgabund,
                   predict(lm.AM.abund, newdata = dat.AM,exclude = 's(site)', se=TRUE))

ggplot(pred, aes(x = lgabund)) +
  geom_point(aes(y = resp), size = 1, alpha = 0.5) + 
  geom_line(aes(y = fit))+
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha = .5) + 
  ylab("CMDD")+
  xlab("log (abundance)")

dat.EM<-dat%>%filter(myc=="EM")
lm.EM.abund<-gam(net.prop.change.A~s(lgabund),  dat=dat.EM)
summary(lm.EM.abund)
plot(lm.EM.abund)

pred <- data.frame(resp= dat.EM$resp,
                   lgabund = dat.EM$lgabund,
                   predict(lm.EM.abund, newdata = dat.EM,exclude = 's(site)', se=TRUE))

ggplot(pred, aes(x = lgabund)) +
  geom_point(aes(y = resp), size = 1, alpha = 0.5) + 
  geom_line(aes(y = fit))+
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha = .5) + 
  ylab("CMDD")+
  xlab("log (abundance)")

#3 C: GLM each site
lm.myc.abund.s<-glm(resp~myc*lgabund*site,dat=dat)
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#3 C: MCMCGLM all sites
#define priors (phylo and site):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002), 
                                             G2 = list(V=1, nu=0.002))) 

lm.myc.abund.mcmc<-MCMCglmm(resp~myc*lgabund,random= ~latin+site,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.abund.mcmc)

ex.grid<- expand.grid(myc=dat$myc)
lgabund.range <- seq(from=quantile(dat$lgabund, 0.001),to=quantile(dat$lgabund, 0.999),length.out=5126) 

pred.dat<-cbind(lgabund.range,ex.grid)
pred<- predict.MCMCglmm(lm.myc.abund.mcmc,type="response",newdata = pred.dat) # predict results for these values


datplot <- ggpredict(lm.myc.mcmc, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = resp),alpha=0.25)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  ylim(-20,10)

#3 C: MCMCGLM each site
#define priors (phylo only):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002)))

lm.myc.abund.mcmc<-MCMCglmm(resp~myc*lgabund*site,random= ~latin,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.abund.mcmc)

#4 Visualize 

# density plot AM v EM
ggplot(data=dat, mapping =aes(x=resp,color=myc, fill=myc))+
  geom_density(alpha=0.3) +
  xlab("Conmycorrhizal Recruitment")+
  ylab("Density")+
  xlim(-1,2)

grid.arrange(Fig3.A,Fig3.B,ncol=2)

#########################
###### GLM OR GAM #######
###### Net A & CM #######
#########################

# GLM net A & CM coef
dat<-read.csv("data/output.study.net.A.CM_10.15.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(site=as.factor(site)) %>%  #for GAMs
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  select(site,latin,myc,net.coef.A, net.coef.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.coef.A", "net.coef.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

# GLM net A & CM prop
dat<-read.csv("data/output.study.net.A.CM_10.15.21.csv", sep=",",header=TRUE) %>% 
  mutate(lgabund=log(abund)) %>%
  mutate(site=as.factor(site)) %>%  #for GAMs
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  select(site,latin,myc,net.prop.change.A, net.prop.change.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.prop.change.A", "net.prop.change.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

# GAM net A & CM prop (1-5)***
dat<-read.csv("data/GAM.output.study.net.m1-5.A.CM_10.15.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(site=as.factor(site)) %>%  #for GAMs
  mutate(myc = fct_relevel(myc,"AM","EM")) %>%
  select(site,latin,myc,net.prop.change.A, net.prop.change.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.prop.change.A", "net.prop.change.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

dat.remove<- dat %>% 
  group_by(site, myc,.drop=FALSE) %>% 
  distinct(latin) %>% 
  summarize(N = n())%>%
  mutate(freq = N / sum(N)) %>%
  filter(freq==1)

dat<- dat %>% filter(!site %in% dat.remove$site) 

#######################
### CONIFER SECTION ###
#######################

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

dat.sp<- unique(dat$latin)
dat.conif<- lookup_table(c(dat.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

dat<- dat %>% left_join(dat.conif,by="latin") %>%
  mutate(myc.conif= paste(myc,Coniferae)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM.NC", myc.conif)) %>% # replace with AM for NC
  mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM.C", myc.conif)) %>% # replace with AM for C
  mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM NA", "AM.NC",myc.conif))%>%
  select(-myc) %>%
  rename(myc=myc.conif) 

#######################
#######################

#2 CM: GLM all sites
lm.myc<-lmer(value~myc*type + (1|site),   dat=dat)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc*type, data= dat, pbkrtest.limit = 10000)
ref.table<-as.data.frame(ref$lsmeans)

ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)+
  ylim(-2,1)

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site,  dat=dat)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)

#overall values
ref<-lsmeans(lm.myc,pairwise~myc*type, data= dat)
ref.tableA<-as.data.frame(ref$lsmeans) %>% mutate(site="Overall")

#site values
ref<-lsmeans(lm.myc,pairwise~myc*type*site, data= dat)
ref.table<-as.data.frame(ref$lsmeans) 
ref.table<-rbind(ref.table,ref.tableA) %>%
  mutate(site = fct_relevel(site, 
                            "Overall", "BCI", "Danum_Valley", 
                            "Heishiding", "Indian_Cave", "Michigan_Big_Woods", "SCBI","SERC","Sinharaja",
                            "Wabikon","Wind_River","Yosemite")) 

ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=dat,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Recruitment")+
  xlab("")+
  coord_flip()+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#2 C: MCMCGLM all sites
#load tree
load("data/GBOTB.extended.rda")
phy<-GBOTB.extended

#merge congenerics
phy.ex<-congeneric.merge(phy,dat$latin, split = "_", cite = FALSE)

#subset tree to observed species in data
alltreesp<- c(phy.ex$tip.label) %>% as.data.frame() 
colnames(alltreesp)<-"latin"
nodrop<- dat %>% select(latin)                                         
todrop<- setdiff(alltreesp, nodrop)%>% as.vector()                        
todrop<-c(t(todrop))                                                 
p.tree<- drop.tip(phy.ex, todrop)                              

#clean up the tree
p.tree$node.label <- NULL
phy$node.label<-NULL
rnd <- inverseA(p.tree)$Ainv # inverse of matrix

# subset dat to include only latin in tips of tree:
dat<-dat %>% filter(latin %in% p.tree$tip.label)

#define priors (phylo and site):
priors <- list(R=list(V=1, nu=0.002), G=list(G1 = list(V=1, nu=0.002), 
                                             G2 = list(V=1, nu=0.002))) 


lm.myc.mcmc<-MCMCglmm(value~myc*type,random= ~latin+site,ginverse=list(latin=rnd), dat=dat, prior =priors)      
summary(lm.myc.mcmc)

ref<-lsmeans(lm.myc.mcmc,pairwise~myc*type, data= dat, pbkrtest.limit = 10000)
ref.table<-as.data.frame(ref$lsmeans)
Fig2.A<-
  ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lower.HPD, ymax=upper.HPD), width=.5,size=1, position=position_dodge(1),alpha=0.3) + 
  geom_point(data=dat,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+ylab("")+
  xlab("")+
  theme_minimal()+
  theme(legend.position = "none") +
  facet_grid(~type)+
  ylim(-3,3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()
  )+
  theme(text = element_text(size = 20),
        axis.text.y=element_text(size = 15))

#2 C: MCMCGLM each site

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site), dat=dat)
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = dat, mapping = aes(x = lgabund, y = value))+
  ylab("Recruitment")+
  xlab("log(abundance)")+
  facet_grid(~group)+
  ylim(-5,2)
  

#3 C: GLM each site
lm.myc.abund.s<-glm(value~myc*lgabund*type*site,dat=dat)
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#3 C: MCMCGLM all sites

#3 C: MCMCGLM each site

#4 Visualize density plot AM v EM
Fig2.B<-
  ggplot(data=dat, mapping =aes(x=value,color=myc, fill=myc))+
  geom_density(alpha=0.3) +
  geom_vline(xintercept = 0, linetype="dashed")+
  coord_flip()+
  xlab("")+
  ylab("") +
  xlim(-3,3)  +
  theme_minimal()+
  facet_grid(cols = vars(type),scales="free")+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text = element_text(size = 20),
        axis.text.y=element_text(size = 15))

grid.arrange(Fig2.B, Fig2.A,nrow=2)

##########################
#### MAP OF LOCATIONS ####
##########################

library(maps)
library(sf)

world <- map_data("world")

# make lat of plots we have samples for:
lat.samp<- lat %>% filter(site %in% c("BCI","San_Lorenzo","Ordway_Swisher","Luquillo","Lambir","Fushan",
                             "Heishiding","Kenting","Lienhuachih", "Wytham_Woods", "Zofin",
                             "Indian_Cave","University_of_California_Santa_Cruz",
                             "University_of_Maryland_Baltimore_County","Utah","Wabikon", "Wind_River",
                             "Yosemite","Laupaheoheo","Speulderbos","Tyson_Research_Center","SERC","SCBI",
                             "Little_Dickie_Woods","Harvard_Forest","Palamanui",
                             "FOREG-31_East","FOREG-31_West","FOREG-27",                               
                             "FOREG-02","FOREG-38","FOREG-28","FOREG-22","FOREG-21"))
lat.nosamp<-lat %>% filter(!site %in% lat.samp)

png("figures/MAPS_FGEO/mapsitessamp_10.15.21.jpg", width=15, height= 7, units='in', res=300)
ggplot(data = world) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group))+
  geom_point(data=lat.dat, 
             aes(x=lon, y=lat), colour="deepskyblue", 
             fill="deepskyblue",pch=21, size=5, alpha=I(0.6))+
  geom_point(data=lat.nodat, 
             aes(x=lon, y=lat), colour="darkgreen", 
             fill="forestgreen",pch=21, size=5, alpha=I(0.35))+
  xlab("Latitude")+ylab ("Longitude")
dev.off()

# make lat of plots we have data for:
#GAM prop mod 2 (1-2)
dat<-read.csv("data/GAM.output.study.net.m1-2.A_10.15.21.csv", sep=",",header=TRUE) %>% mutate(resp = net.prop.change.A)

lat.dat<- lat %>% filter(site %in% dat$site)
lat.nodat<-lat %>% filter(!site %in% dat$site)

png("figures/MAPS_FGEO/mapsitesdat_10.15.21.jpg", width=15, height= 7, units='in', res=300)
ggplot(data = world) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group))+
  geom_point(data=lat.dat, 
             aes(x=lon, y=lat), colour="deepskyblue", 
             fill="deepskyblue",pch=21, size=5, alpha=I(0.6))+
  geom_point(data=lat.nodat, 
             aes(x=lon, y=lat), colour="darkgreen", 
             fill="forestgreen",pch=21, size=5, alpha=I(0.35))+
  xlab("Latitude")+ylab ("Longitude")
dev.off()

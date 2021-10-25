##### mod5 only: conspecific v conmycorrhizal heterspecific:
##### mod2 also: conspecific rel to all heterospecific

#0. rich
#1. Lat
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

# read in lat:
lat<- read.table("data/FGeo_LatLon.txt",header =TRUE) %>%
  mutate(site=toupper(site))%>%
  select(-site)%>%
  rename(site=full_site) %>%
  mutate(abslat=abs(lat))

# set plot theme:
theme_set(theme(text = element_text(size = 10, family = "Tahoma"),
                axis.title = element_text(face="bold"),
                axis.text.x=element_text(size = 7),
                axis.text.y=element_text(size = 7),
                legend.position = "bottom"))

#########################
##### COEF data GLM #####
######## Net A ##########
#########################

#net A
output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE)
#mod 2
output5<-read.csv("data/output.study.net.A_9.29.21.csv", sep=",",header=TRUE) 

output5<-output5 %>% left_join(lat, by="site") %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM")) %>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  filter(!myc=="NM")

richdat5<-
  output5 %>% 
  group_by(site) %>%                         
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                    
  summarize(medcndd= median(net.coef.A)) %>%     
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                       
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                     
  summarize(medcndd= median(net.coef.A)) %>%     
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

#0 C: GLM
lm.rich<-glm(medcndd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcndd),alpha=0.25)+
  ylab("CNDD")+
  xlab("Species Richness")

#0 C: GLM myc
lm.rich.myc<-glm(medcndd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = lat, y = medcndd),alpha=0.25)+
  ylab("CNDD")+
  xlab("Absolute latitude")

#1 C: GLM
lm.lat.myc<-glm(medcndd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#2 C: GLM all sites
lm.myc<-lmer(net.coef.A~myc + (1|site),weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.coef.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CNDD")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: GLM each site
lm.myc.s<-glm(net.coef.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.coef.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CNDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-20,5)

#3 C: GLM all sites
lm.myc.abund<-lmer(net.coef.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.coef.A),alpha=0.25)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  ylim(-20,3)

#3 C: GLM each site
lm.myc.abund.s<-glm(net.coef.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

#########################
##### COEF data GLM #####
######## Net CM #########
#########################

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  mutate(net.coef.CMH= coef_CMH-coef_HMH) %>%
  filter(!myc=="NM")

richdat5<-
  output5 %>% 
  group_by(site) %>%                         
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                   
  summarize(medcmdd= median(net.coef.CMH)) %>%     
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                          
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                    
  summarize(medcmdd= median(net.coef.CMH)) %>%    
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

#0 C: GLM
lm.rich<-glm(medcmdd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcmdd),alpha=0.25)+
  ylab("CMDD")+
  xlab("Species Richness")

#0 C: GLM myc
lm.rich.myc<-glm(medcmdd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

datplot <- ggpredict(lm.rich.myc, terms=c("sprich","myc"), raw=TRUE)
plot(datplot)
ggplot(data = datplot, mapping = aes(x = x, y = predicted,color=myc))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5.myc, mapping = aes(x = sprich, y = medcmdd),alpha=0.25)+
  ylab("CMDD")+
  xlab("Species Richness")

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

#1 C: GLM
lm.lat.myc<-glm(medcmdd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#2 CM: GLM all sites:NS
lm.myc<-lmer(net.coef.CMH~myc + (1|site), weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.coef.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 CM: GLM each site
lm.myc.s<-glm(net.coef.CMH~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s, test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.coef.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.coef.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.coef.CMH),alpha=0.25)+
  ylab("CMDD")+
  xlab("log(abundance)")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=output5,aes(x=lgabund,y=net.coef.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CMDD")+
  xlab("log(abundance)")

#3 C: GLM each site
lm.myc.abund.s<-glm(net.coef.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#########################
##### COEF data GLM #####
###### Net A & CM #######
#########################

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  mutate(net.coef.CMH= coef_CMH-coef_HMH) %>%
  select(site,latin,myc,net.coef.A, net.coef.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.coef.A", "net.coef.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

#2 CM: GLM all sites
lm.myc<-lmer(value~myc*type + (1|site), weights=relabund,  dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CNDD")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)

#overall values
ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.tableA<-as.data.frame(ref$lsmeans) %>% mutate(site="Overall")

#site values
ref<-lsmeans(lm.myc,pairwise~myc*type*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans) 
ref.table<-rbind(ref.table,ref.tableA) %>%
  mutate(site = fct_relevel(site, 
                            "Overall", "BCI", "Danum_Valley", 
                            "Heishiding", "Indian_Cave", "Michigan_Big_Woods", "SCBI","SERC","Sinharaja",
                            "Wabikon","Wind_River","Yosemite")) 

ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("CNDD")+
  xlab("")+
  coord_flip()+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = value))+
  ylab("Recruitment")+
  xlab("log(abundance)")

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  mutate(net.coef.CMH= coef_CMH-coef_HMH) %>%
  filter(!myc=="NM")

#5 most abundant species by site for each myc group:
abund<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.coef.A","net.coef.CMH") %>%
  gather(type,value,5:6)
rare<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = -lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.coef.A","net.coef.CMH") %>%
  gather(type,value,5:6)

ggplot(data = abund, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))

ggplot(data = rare, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))

#########################
##### PROP data GLM #####
######### Net A #########
#########################

#net A
output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE)
#mod 2
output5<-read.csv("data/output.study.net.A_9.29.21.csv", sep=",",header=TRUE) 

output5<-output5 %>%left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM")) %>%
  filter(!myc=="NM")

richdat5<-
  output5 %>% 
  group_by(site) %>%                         
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                   
  summarize(medcndd= median(net.prop.change.A)) %>%      
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                        
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                    
  summarize(medcndd= median(net.prop.change.A)) %>%     
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

#0 C: GLM
lm.rich<-glm(medcndd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Species Richness")

#0 C: GLM myc
lm.rich.myc<-glm(medcndd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = lat, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Absolute latitude")

#1 C: GLM
lm.lat.myc<-glm(medcndd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#2 C: GLM all sites
lm.myc<-lmer(net.prop.change.A~myc + (1|site),weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: GLM each site
lm.myc.s<-glm(net.prop.change.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.prop.change.A),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  ylim(-2,1)

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

#########################
##### PROP data GLM #####
######### Net CM ########
#########################

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  filter(!myc=="NM")

richdat5<-
  output5 %>% 
  group_by(site) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

#0 C: GLM
lm.rich<-glm(medcmdd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcmdd),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Species Richness")

#0 C: GLM myc
lm.rich.myc<-glm(medcmdd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

#1 C: GLM
lm.lat.myc<-glm(medcmdd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#2 CM: GLM all sites:NS
lm.myc<-lmer(net.prop.change.CMH~myc + (1|site), weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

#2 CM: GLM each site
lm.myc.s<-glm(net.prop.change.CMH~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s, test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.CMH,color=myc),position=position_jitter(width =.15), alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.prop.change.CMH),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("log(abundance)")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("red","cyan3"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=output5,aes(x=lgabund,y=net.coef.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("log(abundance)")

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#########################
##### PROP data GLM #####
###### Net A & CM #######
#########################

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  select(site,latin,myc,net.prop.change.A, net.prop.change.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.prop.change.A", "net.prop.change.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

#2 CM: GLM all sites
lm.myc<-lmer(value~myc*type + (1|site), weights=relabund,  dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)

ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Recruitment")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)

ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.tableA<-as.data.frame(ref$lsmeans) %>% mutate(site="Overall")

ref<-lsmeans(lm.myc,pairwise~myc*type*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans) 
ref.table<-rbind(ref.table,ref.tableA) %>%
  mutate(site = fct_relevel(site, 
                            "Overall", "BCI", "Danum_Valley", 
                            "Heishiding", "Indian_Cave", "Michigan_Big_Woods", "SCBI","SERC","Sinharaja",
                            "Wabikon","Wind_River","Yosemite")) 

ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Recruitment")+
  xlab("")+
  coord_flip()+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = value))+
  ylab("Recruitment")+
  xlab("log(abundance)")

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  filter(!myc=="NM")

#5 most abundant species by site for each myc group:
abund<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)
rare<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = -lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)

png("figures/p_quantchange_8.20.21/M5.A.CM.mycabundrecruit_9.29.21.jpg", width=7, height= 15, units='in', res=300)
ggplot(data = abund, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))
dev.off()

png("figures/p_quantchange_8.20.21/M5.A.CM.mycrarerecruit_9.29.21.jpg", width=7, height= 15, units='in', res=300)
ggplot(data = rare, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))
dev.off()

#########################
##### PROP data GAM #####
######### Net A #########
#########################

##### PLUS CONIFER ######
####### DIVISION ########

#net A
output5<-read.csv("data/GAM.output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE)
#mod 2
output5<-read.csv("data/GAM.output.study.net.A_9.29.21.csv", sep=",",header=TRUE)

output5<-output5 %>% left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))
  

###################
# conifer section #
###################

# if using this section, don't do lat/propEM section

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

output5.sp<- unique(output5$latin)
output5.conif<- lookup_table(c(output5.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

output5<- output5 %>% left_join(output5.conif,by="latin") %>%
            mutate(myc.conif= paste(myc,Coniferae)) %>%
            mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM.NC", myc.conif)) %>% # replace with AM for NC
            mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM.C", myc.conif)) %>% # replace with AM for C
            mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
            mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
            mutate(myc.conif= ifelse(myc.conif=="AM NA", "AM.NC",myc.conif))%>%
            select(-myc) %>%
            rename(myc=myc.conif) 

###################
# conifer section #
###################

richdat5<-
  output5 %>% 
  group_by(site) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcndd= median(net.prop.change.A)) %>%      # summarize median cndd
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcndd= median(net.prop.change.A)) %>%      # summarize median cndd
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site")           # join with lat data

latdat5.EM<- latdat5.myc %>%
  filter(myc=="EM")%>%
  rename(sprichEM=sprich) %>%
  left_join(richdat5, by="site") %>%
  mutate(propEM= sprichEM/sprich) 

#0 C: GLM
lm.rich<-glm(medcndd~sprich, dat=latdat5)               
summary(lm.rich)

png("figures/GAM/CS_sprich_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Species Richness")
dev.off()

#0 C: GLM myc
lm.rich.myc<-glm(medcndd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

png("figures/GAM/CS_lat_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.lat, terms=c("abslat"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = lat, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Absolute latitude")
dev.off()

#1 C: GLM
lm.lat.myc<-glm(medcndd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#1 C: GLM version prop EM
lm.pmyc<-glm(medcndd~propEM, dat=latdat5.EM)               
summary(lm.pmyc)

png("figures/GAM/CS_propEM_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.pmyc, terms=c("propEM"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5.EM, mapping = aes(x = propEM, y = medcndd),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Proportion EM")
dev.off()

#2 C: GLM all sites
lm.myc<-lmer(net.prop.change.A~myc + (1|site), dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans) %>% mutate(site="Overall") %>% mutate(abslat=0)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.A,color=myc),position=position_jitter(width =.15), alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-2.5,0.25)

#2 C: GLM each site
lm.myc.s<-glm(net.prop.change.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans) %>%
      left_join(latdat5, by = "site") %>%
      arrange(abslat) 

ggplot(ref.table, aes(x=reorder(myc, abslat), lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-2,1)

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM")))
summary(lm.myc.abund)

png("figures/GAM/CS_abund_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.prop.change.A),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  ylim(-2,1)
dev.off()

png("figures/GAM/CS_abund*myc_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("red","cyan3"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=output5,aes(x=lgabund,y=net.prop.change.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conspecific Recruitment")+
  xlab("log(abundance)")+
  ylim(-2,1)
dev.off()

#3 C: GLM each site: ugly
lm.myc.abund.s<-glm(net.prop.change.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

#########################
##### PROP data GAM #####
######### Net CM ########
#########################

output5<-read.csv("data/GAM.output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))

###################
# conifer section #
###################

# if using this section, don't do lat/propEM section

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

output5.sp<- unique(output5$latin)
output5.conif<- lookup_table(c(output5.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

output5<- output5 %>% left_join(output5.conif,by="latin") %>%
  mutate(myc.conif= paste(myc,Coniferae)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM.NC", myc.conif)) %>% # replace with AM for NC
  mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM.C", myc.conif)) %>% # replace with AM for C
  mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="AM NA", "AM.NC",myc.conif))%>%
  select(-myc) %>%
  rename(myc=myc.conif) 

###################
# conifer section #
###################

richdat5<-
  output5 %>% 
  group_by(site) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(richdat5, by="site") %>%
  left_join(lat, by ="site") 

richdat5.myc<-
  output5 %>% 
  group_by(site,myc) %>%                          # group by site
  summarize(sprich = length(latin)) 

latdat5.myc<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(richdat5.myc, by="site") %>%
  left_join(lat, by ="site") 

#0 C: GLM
lm.rich<-glm(medcmdd~sprich, dat=latdat5)               
summary(lm.rich)

datplot <- ggpredict(lm.rich, terms=c("sprich"), raw=TRUE)
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = latdat5, mapping = aes(x = sprich, y = medcmdd),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Species Richness")

#0 C: GLM myc
lm.rich.myc<-glm(medcmdd~sprich*myc, dat=latdat5.myc)               
summary(lm.rich.myc)

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

#1 C: GLM
lm.lat.myc<-glm(medcmdd~abslat*myc, dat=latdat5.myc)               
summary(lm.lat.myc)

#2 CM: GLM all sites:NS
lm.myc<-lmer(net.prop.change.CMH~myc,dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Mycorrhizal type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylim(-2,.25)

#2 CM: GLM each site
lm.myc.s<-glm(net.prop.change.CMH~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s, test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=net.prop.change.A,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = net.prop.change.CMH),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("log(abundance)")

png("figures/GAM/CM_abund*myc_10.1.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("red","cyan3"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data=output5,aes(x=lgabund,y=net.prop.change.CMH,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Conmycorrhizal Recruitment")+
  xlab("log(abundance)")+
  ylim(-2,1)
dev.off()

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#########################
##### PROP data GAM #####
###### Net A & CM #######
#########################

##### PLUS CONIFER ######
####### DIVISION ########

output5<-read.csv("data/GAM.output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM")) %>%
  select(site,latin,myc,net.prop.change.A, net.prop.change.CMH,relabund,lgabund)%>%
  gather(type,value,4:5) %>%
  mutate(type = factor(type, levels = c("net.prop.change.A", "net.prop.change.CMH"),
                       labels = c("Conspecific", "Conmycorrhizal")
  ))

###################
# conifer section #
###################

# if using this section, don't do lat/propEM section

Coniferae<-add_higher_order()%>% filter(Coniferae=="Coniferae") %>% select(c("family","genus","Coniferae"))

output5.sp<- unique(output5$latin)
output5.conif<- lookup_table(c(output5.sp), by_species=TRUE, version="1.1.5") %>% 
  mutate(latin=rownames(.)) %>%
  left_join(Coniferae, by= c("family","genus")) %>%
  replace_na(list(Coniferae = "Non-Coniferae")) %>%
  select(c("latin","Coniferae"))

output5<- output5 %>% left_join(output5.conif,by="latin") %>%
  mutate(myc.conif= paste(myc,Coniferae)) %>%
  #mutate(myc.conif= ifelse(myc.conif=="AM Non-Coniferae", "AM", myc.conif)) %>% # replace with AM for NC
  #mutate(myc.conif= ifelse(myc.conif=="AM Coniferae", "AM", myc.conif)) %>% # replace with AM for C
  mutate(myc.conif= ifelse(myc.conif=="EM Non-Coniferae", "EM.NC", myc.conif)) %>%
  mutate(myc.conif= ifelse(myc.conif=="EM Coniferae", "EM.C", myc.conif)) %>%
  select(-myc) %>%
  rename(myc=myc.conif) 

###################
# conifer section #
###################

#2 CM: GLM all sites
lm.myc<-lmer(value~myc*type + (1|site), weights=relabund,  dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)

png("figures/GAM/CSCM*myc_10.1.21.jpg", width=6, height= 6, units='in', res=300)
ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5, pbkrtest.limit = 10000)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Recruitment")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)+
  ylim(-1,0.1)
dev.off()

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)

ref<-lsmeans(lm.myc,pairwise~myc*type*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans) %>% 
  left_join(latdat5, by = "site") %>%
  arrange(abslat) 

png("figures/GAM/CSCM*myc*site_10.1.21.jpg", width=6, height= 10, units='in', res=300)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  geom_point(data=output5,aes(x=myc,y=value,color=myc),position=position_jitter(width =.15),alpha=0.25)+
  ylab("Recruitment")+
  xlab("")+
  coord_flip()+
  ylim(-1.5,1.5)+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))
dev.off()

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("red","cyan3","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  geom_point(data = output5, mapping = aes(x = lgabund, y = value),alpha=0.25)+
  ylab("Recruitment")+
  xlab("log(abundance)")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=facet))+
  scale_colour_manual(values = c("red","cyan3","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("log(abundance)")+
  facet_grid(~group)

#3 C: GLM each site:
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/GAM.output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"AM","EM"))%>%
  filter(!myc=="NM")

#5 most abundant species by site for each myc group:
abund<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)
rare<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = -lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)

ggplot(data = abund, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))

ggplot(data = rare, mapping = aes(x = reorder(latin, value), y = value,color=myc))+
  geom_point()+
  coord_flip()+
  ylab("Recruitment")+
  xlab("Species")+
  ylim(-2,1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(myc),
    cols=vars(type))

#########################
##### PROP data GAM #####
##### CONIFER MYC #######
######### Net A #########
#########################

output5<-read.csv("data/GAM.conif.output.study.net.A_9.29.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  select(site,myc.con,net.prop.change.A,relabund,lgabund)

#2 C: GLM all sites
lm.myc<-lm(net.prop.change.A~myc.con*site,weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

##########################
#### MAP OF LOCATIONS ####
##########################

library(maps)
library(sf)

world <- map_data("world")

# make lat of plots we have samples for:

lat$site

lat.samp<- lat %>% filter(site %in% c("BCI","San_Lorenzo","Ordway_Swisher","Luquillo","Lambir","Fushan",
                             "Heishiding","Kenting","Lienhuachih", "Wytham_Woods", "Zofin",
                             "Indian_Cave","University_of_California_Santa_Cruz",
                             "University_of_Maryland_Baltimore_County","Utah","Wabikon", "Wind_River",
                             "Yosemite","Laupaheoheo","Speulderbos","Tyson_Research_Center","SERC","SCBI",
                             "Little_Dickie_Woods","Harvard_Forest","Palamanui",
                             "FOREG-31_East","FOREG-31_West","FOREG-27",                               
                             "FOREG-02","FOREG-38","FOREG-28","FOREG-22","FOREG-21"))
lat.nosamp<-lat %>% filter(!site %in% c("BCI","San_Lorenzo","Ordway_Swisher","Luquillo","Lambir","Fushan",
                                        "Heishiding","Kenting","Lienhuachih", "Wytham_Woods", "Zofin",
                                        "Indian_Cave","University_of_California_Santa_Cruz",
                                        "University_of_Maryland_Baltimore_County","Utah","Wabikon", "Wind_River",
                                        "Yosemite","Laupaheoheo","Speulderbos","Tyson_Research_Center","SERC","SCBI",
                                        "Little_Dickie_Woods","Harvard_Forest","Palamanui",
                                        "FOREG-31_East","FOREG-31_West","FOREG-27",                               
                                        "FOREG-02","FOREG-38","FOREG-28","FOREG-22","FOREG-21"))

png("figures/mapsitessamp_9.29.21.jpg", width=15, height= 7, units='in', res=300)
ggplot(data = world) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group))+
  #coord_sf(xlim = c(-100, -50), ylim = c(-5, 50), expand = FALSE)+
  geom_point(data=lat.samp, 
             aes(x=lon, y=lat), colour="deepskyblue", 
             fill="deepskyblue",pch=21, size=5, alpha=I(0.6))+
  geom_point(data=lat.nosamp, 
             aes(x=lon, y=lat), colour="darkgreen", 
             fill="forestgreen",pch=21, size=5, alpha=I(0.35))+
  
  xlab("Latitude")+ylab ("Longitude")
dev.off()

# make lat of plots we have data for:

lat$site

lat.dat<- lat %>% filter(site %in% c("Cocoli", "Wind_River", "Lienhuachih","Speulderbos","SERC",
                                     "SCBI","BCI","San_Lorenzo","Yosemite","Huai_Kha_Khaeng",
                                     "Wanang","Sinharaja","Wabikon","Khoa_Chong","Palanan",
                                     "Michigan_Big_Woods","Heishiding","Indian_Cave","Laupahoehoe","Rabi",
                                      "FOREG-31_East","FOREG-31_West","FOREG-27",                               
                                      "FOREG-02","FOREG-38","FOREG-28","FOREG-22","FOREG-21"))
lat.nodat<-lat %>% filter(!site %in% c("Cocoli", "Wind_River", "Lienhuachih","Speulderbos","SERC",
                                       "SCBI","BCI","San_Lorenzo","Yosemite","Huai_Kha_Khaeng",
                                       "Wanang","Sinharaja","Wabikon","Khoa_Chong","Palanan",
                                       "Michigan_Big_Woods","Heishiding","Indian_Cave","Laupahoehoe","Rabi",
                                       "FOREG-31_East","FOREG-31_West","FOREG-27",                               
                                       "FOREG-02","FOREG-38","FOREG-28","FOREG-22","FOREG-21"))

png("figures/mapsitesdat_9.29.21.jpg", width=15, height= 7, units='in', res=300)
ggplot(data = world) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group))+
  #coord_sf(xlim = c(-100, -50), ylim = c(-5, 50), expand = FALSE)+
  geom_point(data=lat.dat, 
             aes(x=lon, y=lat), colour="deepskyblue", 
             fill="deepskyblue",pch=21, size=5, alpha=I(0.6))+
  geom_point(data=lat.nodat, 
             aes(x=lon, y=lat), colour="darkgreen", 
             fill="forestgreen",pch=21, size=5, alpha=I(0.35))+
  
  xlab("Latitude")+ylab ("Longitude")
dev.off()


##### mod5 only: conspecific v conmycorrhizal heterspecific:

# load packages
library(tidyverse)
library(lsmeans)
library(ggeffects)
library(lme4)
library(lmerTest)
library(gridExtra)

# read in lat:
lat<- read.table("data/FGeo_LatLon.csv",header =TRUE) %>%
  mutate_if(is.numeric, round, digits = 4) %>%
  mutate(site=toupper(site))%>%
  select(-site)%>%
  rename(site=full_site) %>%
  mutate(abslat=abs(lat))

# set plot theme:
theme_set(theme(text = element_text(size = 20, family = "Tahoma"),
                axis.title = element_text(face="bold"),
                axis.text.x=element_text(size = 10),
                axis.text.y=element_text(size = 10),
                legend.position = "bottom"))

#########################
##### PROP data GLM #####
###### Net A: COEF ######
#########################

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM")) %>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcndd= median(net.coef.A)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("Absolute latitude")

#2 C: GLM all sites
lm.myc<-lmer(net.coef.A~myc + (1|site),weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.B<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("Recruitment")+
  xlab("Mycorrhizal type")+
  ylim(-1,0.01)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: GLM each site
lm.myc.s<-glm(net.coef.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.C<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("CNDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  ylim(-1.50,0.10)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

png("figures/GLM.propFig1.9.6.21.jpg", width=20, height= 6, units='in', res=300)
grid.arrange(Fig1.A,Fig1.B,Fig1.C,ncol=3)
dev.off()

#3 C: GLM all sites
lm.myc.abund<-lmer(net.coef.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

png("figures/GLM.propCNDD.lgabund.9.6.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: ugly
lm.myc.abund.s<-glm(net.coef.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","site"))
ggplot(data = datplot%>% rename(site =group), mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)

#########################
##### coef data GLM #####
######### Net CM ########
#########################

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  mutate(net.coef.A= coef_A-coef_HMH) %>%
  mutate(net.coef.CMH= coef_CMH-coef_HMH) %>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.coef.CMH)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

#2 CM: GLM all sites:NS
lm.myc<-lmer(net.coef.CMH~myc + (1|site), weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

#2 CM: GLM each site
lm.myc.s<-glm(net.coef.CMH~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s, test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.coef.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
Fig3.A<-ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
Fig3.B<-ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")

Fig3.B
png("figures/GLM.propCMDD.Fig3.9.6.21.jpg", width=12, height= 6, units='in', res=300)
grid.arrange(Fig3.A,Fig3.B, ncol=2)
dev.off()

#3 C: GLM each site
lm.myc.abund.s<-glm(net.coef.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","myc","site"))
ggplot(data = datplot%>% rename(site =facet), mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)
dev.off()

##### consp v conmyc: mod 5

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM") %>%
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
library(lsmeans)

png("figures/GLM.prop.myc*type.9.6.21.jpg", width=6, height= 6, units='in', res=300)
ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
FigA<- ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("CNDD")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)
dev.off()

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)
png("figures/GLM.propCMDD.myc*type*site.9.6.21.jpg", width=6, height= 12, units='in', res=300)

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
  ylab("CNDD")+
  xlab("")+
  coord_flip()+
  ylim(-2,1)+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

png("figures/GLM.propCMDD.lgabund*type.9.6.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: crap.
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
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

png("figures/p_quantchange_8.20.21/M5.A.CM.mycabundrecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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

png("figures/p_quantchange_8.20.21/M5.A.CM.mycrarerecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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
##### PROP data GLM #####
######### Net A #########
#########################

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM")) %>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcndd= median(net.prop.change.A)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"))
Fig1.A<-ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("Absolute latitude")

#2 C: GLM all sites
lm.myc<-lmer(net.prop.change.A~myc + (1|site),weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.B<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("Recruitment")+
  xlab("Mycorrhizal type")+
  ylim(-0.30,0.01)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: GLM each site
lm.myc.s<-glm(net.prop.change.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.C<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("CNDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  ylim(-1.50,0.10)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

png("figures/GLM.propFig1.9.6.21.jpg", width=20, height= 6, units='in', res=300)
grid.arrange(Fig1.A,Fig1.B,Fig1.C,ncol=3)
dev.off()

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

png("figures/GLM.propCNDD.lgabund.9.6.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: ugly
lm.myc.abund.s<-glm(net.prop.change.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","site"))
ggplot(data = datplot%>% rename(site =group), mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)

#########################
##### PROP data GLM #####
######### Net CM ########
#########################

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

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
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
Fig3.A<-ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
Fig3.B<-ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")

png("figures/GLM.propCMDD.Fig3.9.6.21.jpg", width=12, height= 6, units='in', res=300)
grid.arrange(Fig3.A,Fig3.B, ncol=2)
dev.off()

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","myc","site"))
ggplot(data = datplot%>% rename(site =facet), mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)
dev.off()

##### consp v conmyc: mod 5

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM") %>%
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

png("figures/GLM.prop.myc*type.9.6.21.jpg", width=6, height= 6, units='in', res=300)
ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("Recruitment")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)
dev.off()

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)
png("figures/GLM.propCMDD.myc*type*site.9.6.21.jpg", width=6, height= 12, units='in', res=300)

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
  ylab("Recruitment")+
  xlab("")+
  coord_flip()+
  ylim(-2,1)+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

png("figures/GLM.propCMDD.lgabund*type.9.6.21.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: crap.
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/output.study.net.A.CM_8.16.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM")

#5 most abundant species by site for each myc group:
abund<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)
rare<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = -lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)

png("figures/p_quantchange_8.20.21/M5.A.CM.mycabundrecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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

png("figures/p_quantchange_8.20.21/M5.A.CM.mycrarerecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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

output5<-read.csv("data/GAM.output.study.net.A.CM_9.9.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM")) %>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcndd= median(net.prop.change.A)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 C: GLM
lm.lat<-glm(medcndd~abslat, dat=latdat5)               
summary(lm.lat)

datplot <- ggpredict(lm.lat, terms=c("abslat"))
Fig1.A<-ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("Absolute latitude")

#2 C: GLM all sites
lm.myc<-lmer(net.prop.change.A~myc + (1|site),weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc)

ref<-lsmeans(lm.myc,pairwise~myc, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.B<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("Recruitment")+
  xlab("Mycorrhizal type")+
  ylim(-0.5,0.1)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#2 C: GLM each site
lm.myc.s<-glm(net.prop.change.A~myc*site, weights=relabund, dat=output5)      
summary(lm.myc.s)
anova(lm.myc.s,test="LRT")

ref<-lsmeans(lm.myc.s,pairwise~myc*site, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
Fig1.C<-ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("CNDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  ylim(-2,1.5)+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

png("figures/GLM.propFig1.9.921.jpg", width=20, height= 6, units='in', res=300)
grid.arrange(Fig1.A,Fig1.B,Fig1.C,ncol=3)
dev.off()

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.A~myc*lgabund+(1|site), weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

png("figures/GLM.propCNDD.lgabund.9.921.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted,color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: ugly
lm.myc.abund.s<-glm(net.prop.change.A~myc*lgabund*site, weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s,test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","myc","site"))
ggplot(data = datplot%>% rename(site =facet), mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)

#########################
##### PROP data GAM #####
######### Net CM ########
#########################

output5<-read.csv("data/GAM.output.study.net.A.CM_9.9.21.csv", sep=",",header=TRUE) %>%
  left_join(lat, by="site")%>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM")

latdat5<-output5 %>% 
  group_by(site) %>%                     # group by site
  summarize(medcmdd= median(net.prop.change.CMH)) %>%      # summarize median cndd
  left_join(lat, by ="site")             # join with lat data

#1 CM: GLM: NS
lm.lat<-glm(medcmdd~abslat, dat=latdat5)              
summary(lm.lat)

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
  ylab("CMDD")+
  xlab("Mycorrhizal type")+
  facet_grid(~site,scales="free")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")

#3 C: GLM all sites
lm.myc.abund<-lmer(net.prop.change.CMH~myc*lgabund+(1|site), weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)

datplot <- ggpredict(lm.myc.abund, terms=c("lgabund"))
Fig3.A<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")
Fig3.A
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","myc"))
#Fig3.B<-
  ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CMDD")+
  xlab("log(abundance)")

png("figures/GLM.propCMDD.Fig3.9.921.jpg", width=12, height= 6, units='in', res=300)
grid.arrange(Fig3.A,Fig3.B, ncol=2)
dev.off()

#3 C: GLM each site
lm.myc.abund.s<-glm(net.prop.change.CMH~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

datplot <- ggpredict(lm.myc.abund.s, terms=c("lgabund","myc","site"))
ggplot(data = datplot%>% rename(site =facet), mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("CNDD")+
  xlab("log(abundance)")+
  facet_grid(~site)
dev.off()

##### consp v conmyc: mod 5

output5<-read.csv("data/GAM.output.study.net.A.CM_9.9.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM") %>%
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

png("figures/GLM.prop.myc*type.9.921.jpg", width=6, height= 6, units='in', res=300)
ref<-lsmeans(lm.myc,pairwise~myc*type, data= output5)
ref.table<-as.data.frame(ref$lsmeans)
ggplot(ref.table, aes(myc, lsmean,color=myc)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  ylab("Recruitment")+
  xlab("")+
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(~type)
dev.off()

#2 CM: GLM each site
lm.myc<-lm(value~myc*type*site, weights=relabund, dat=output5)      
summary(lm.myc)
anova(lm.myc,test="LRT")

# check assumptions
par(mfrow=c(2,2))
plot(lm.myc)
library(lsmeans)
png("figures/GLM.propCMDD.myc*type*site.9.921.jpg", width=6, height= 12, units='in', res=300)

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
  ylab("Recruitment")+
  xlab("")+
  coord_flip()+
  #ylim(-2,1)+
  labs(color = "Mycorrhizal Type")+
  geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  facet_grid(
    rows=vars(site),
    cols=vars(type))

#3 C: GLM all sites
lm.myc.abund<-lmer(value~myc*lgabund*type+(1|site),weights=relabund, dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund)
anova(lm.myc.abund, test="LRT")

png("figures/GLM.propCMDD.lgabund*type.9.921.jpg", width=6, height= 6, units='in', res=300)
datplot <- ggpredict(lm.myc.abund, terms=c("lgabund","type"))
ggplot(data = datplot, mapping = aes(x = x, y = predicted, color=group))+
  scale_colour_manual(values = c("cyan3", "red","darkgreen"))+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("Recruitment")+
  xlab("log(abundance)")
dev.off()

#3 C: GLM each site: crap.
lm.myc.abund.s<-glm(net.prop.change.CM~myc*lgabund*site,weights=relabund,dat=output5 %>% mutate(myc = fct_relevel(myc,"AM","EM","NM")) %>% filter(!myc=="NM"))
summary(lm.myc.abund.s)
anova(lm.myc.abund.s, test="LRT")

#####PLOT MOST ABUNDANT SPECIES AM V EM:

output5<-read.csv("data/GAM.output.study.net.A.CM_9.9.21.csv", sep=",",header=TRUE) %>%
  mutate(lgabund=log(abund)) %>%
  mutate(myc = fct_relevel(myc,"NM","AM","EM"))%>%
  filter(!myc=="NM")

#5 most abundant species by site for each myc group:
abund<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)
rare<-output5 %>% distinct(latin,.keep_all=TRUE) %>% group_by(site,myc) %>% slice_max(order_by = -lgabund, n = 5) %>%
  select("site","latin","myc","lgabund","net.prop.change.A","net.prop.change.CMH") %>%
  gather(type,value,5:6)

png("figures/p_quantchange_8.20.21/M5.A.CM.mycabundrecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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

png("figures/p_quantchange_8.20.21/M5.A.CM.mycrarerecruit_8.16.21.jpg", width=7, height= 15, units='in', res=300)
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
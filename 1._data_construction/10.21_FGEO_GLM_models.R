## This script creates all initial models for all sites and species 

###################################################
###################################################
################# SET UP PACKAGES #################
###################################################
###################################################

library(tidyverse)  # tidyverse
library(mgcv)       # gam models
library(doParallel) # parallelization
library(foreach)    # parallelization


###################################################
###################################################
################## SET UP PARS ####################
###################################################
###################################################

abundance = 'stem'			            # Option to use stem abundance or basal area to calculate distance-weighted abundances ('stem' = stem abundances; 'ba' = basal area)
null.type = 'allometric.disp.null'	# Selection of null-model type: 'allometric.disp.null' = dispersal-kernel model with dispersal differing among species as a function of their max. height; 'constant.disp.null' = dispersal-kernel model with constant dispersal among species (set in function code above

###################################################
###################################################
################## BRING IN DAT ###################
###################################################
###################################################
###################################################

dat<- readRDS(dat.myc, 'data/allcensusdata_withmyc.rds')

site = unique(dat$site)                                                       # unique sites

dist.weighted.list<- readRDS('data/density_dependence.input.data.rds')

###################################################
###################################################
################## RUN ANALYSES ###################
###################################################
###################################################
###################################################

set.seed(609)

mod.list <- list() # full list site, sp, model
for (k in site) { 
  site.list<-list()
  # bring in model inputs
  lp.dat = dist.weighted.list[[k]]
  # Run model per species:
  for(i in 1:SP) {
    if(lp.dat$N[i] > tr & lp.dat$adultquads[i] > 9 & lp.dat$sapquads[i] > 9) 
    {
      tbl = data.frame(A = lp.dat$adultDistWeightedAbund[,i], ln_A =log(lp.dat$adultDistWeightedAbund[,i]),Ha = lp.dat$HadultDistWeightedAbund[,i], Hs = lp.dat$HsapDistWeightedAbund[,i],
                       CMa = lp.dat$CMadultDistWeightedAbund[,i], CMs = lp.dat$CMsapDistWeightedAbund[,i], CMHa = lp.dat$CMHadultDistWeightedAbund[,i],
                       CMHs = lp.dat$CMHsapDistWeightedAbund[,i],HMHa = lp.dat$HMHadultDistWeightedAbund[,i], HMHs = lp.dat$HMHsapDistWeightedAbund[,i], S = lp.dat$saplings[,i])
      # model 1
      mod1 = glm(S ~ A ,offset=(ln_A), family = 'poisson', data = tbl) 
      # model 2
      mod2 = glm(S ~ A+ Ha ,offset=(ln_A), family = 'poisson', data = tbl) 
      # model 3
      mod3 = glm(S ~ A+ CMa ,offset=(ln_A), family = 'poisson', data = tbl) 
      # model 4    
      mod4 = glm(S ~ A+ CMa +Ha,offset=(ln_A), family = 'poisson', data = tbl) 
      # model 5
      mod5 = glm(S ~ A+ CMHa +HMHa,offset=(ln_A), family = 'poisson', data = tbl) 
      # save site.list
      site.list[[as.character(species[i])]] <- list(mod1,mod2,mod3,mod4,mod5) # put model outputs in list
      names(site.list[[as.character(species[i])]])<-c("mod1","mod2","mod3","mod4","mod5") # name models for this species
    }
    mod.list[[k]] <- site.list 
  }
}

#save model list as .rds object.
saveRDS(mod.list, 'data/GLMmodeloutputs_m1-5_10.21.rds')

print('GLM models complete')

##################################
########## GLM MOD 1-5 ###########
###### CALC PROP CHANGE 1-2 ######
####### CALC COEFFICIENTS ########
##################################

mod.list<-readRDS('data/GLMmodeloutputs_m1-5_10.21.rds')
mod.list['Luquillo'] <- NULL

# simplify dat for myc status later
dat.myc<-dat %>% select(latin,myc) %>%
  distinct(latin,.keep_all=TRUE)

# count abund and rel.abund per species per site where dbh>0
dat.abund<- dat %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(abund = n()) %>%
  mutate(relabund = abund / sum(abund)) 

# count ba and rel.ba per species per site where dbh>0
dat.ba<- dat %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(ba = sum(dbh)) %>%
  mutate(relba = ba / sum(ba)) 

### calculate A effects
tree.range <- c(0.01,0.61) 
A <- tree.range 
output.study.A <- data.frame() 

for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]] 
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(A, t(colMeans(species.dat[[1]]$data[,c('Ha','CMa','CMHa','HMHa','ln_A')])))
      pred<-predict.glm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (recruitment_per_capita2/ recruitment_per_capita1) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      coef_A<-species.dat[[m]]$coefficients[['A']]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2,coef_A) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_A"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_A"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.A <- rbind(output.site, output.study.A)
}  

colnames(output.study.A) <- c('site','species','model','prop.change.A','recruit_pc1','recruit_pc2','se1','se2','coef_A')

output.study.A <- output.study.A %>%
  rename(latin=species)%>%
  left_join(dat.myc, by=c("latin")) %>%
  left_join(dat.abund, by= c("site","latin")) %>%
  left_join(dat.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate CM effects
CMa <- tree.range 
output.study.CM <- data.frame()    

for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod3","mod4")] # subset model 3 and 4 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(CMa, t(colMeans(species.dat[[1]]$data[,c('A','Ha','CMHa','HMHa','ln_A')])))
      pred<-predict.lm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      coef_CM<-species.dat[[m]]$coefficients[['CMa']]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2,coef_CM) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_CM"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_CM"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.CM <- rbind(output.site, output.study.CM)
}  

colnames(output.study.CM) <- c('site','species','model','prop.change.CM','recruit_pc1','recruit_pc2','se1','se2','coef_CM')

output.study.CM <- output.study.CM %>%
  rename(latin=species)%>%
  left_join(dat.myc, by=c("latin")) %>%
  left_join(dat.abund, by= c("site","latin")) %>%
  left_join(dat.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate CMH effects
CMHa <- tree.range 
output.study.CMH <- data.frame()    

for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod5")] # subset model 5 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(CMHa, t(colMeans(species.dat[[1]]$data[,c('A','Ha','CMa','HMHa','ln_A')])))
      pred<-predict.lm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      coef_CMH<-species.dat[[m]]$coefficients[['CMHa']]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2,coef_CMH) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_CMH"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_CMH"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.CMH <- rbind(output.site, output.study.CMH)
}  

colnames(output.study.CMH) <- c('site','species','model','prop.change.CMH','recruit_pc1','recruit_pc2','se1','se2','coef_CMH')

output.study.CMH <- output.study.CMH %>%
  rename(latin=species)%>%
  left_join(dat.myc, by=c("latin")) %>%
  left_join(dat.abund, by= c("site","latin")) %>%
  left_join(dat.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate H effects
Ha <- tree.range 
output.study.H <- data.frame()    

for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod2","mod4")] # subset model 2 and 4 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(Ha, t(colMeans(species.dat[[1]]$data[,c('A','CMHa','CMa','HMHa','ln_A')])))
      pred<-predict.lm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      coef_H<-species.dat[[m]]$coefficients[['Ha']]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2,coef_H) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_H"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_H"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.H <- rbind(output.site, output.study.H)
} 

colnames(output.study.H) <- c('site','species','model','prop.change.H','recruit_pc1','recruit_pc2','se1','se2','coef_H')

output.study.H <- output.study.H %>%
  rename(latin=species)%>%
  left_join(dat.myc, by=c("latin")) %>%
  left_join(dat.abund, by= c("site","latin")) %>%
  left_join(dat.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate HMHa effects
HMHa <- tree.range 
output.study.HMH <- data.frame()    

for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod5")] # subset model 5 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(HMHa, t(colMeans(species.dat[[1]]$data[,c('A','Ha','CMa','CMHa','ln_A')])))
      pred<-predict.lm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      coef_HMH<-species.dat[[m]]$coefficients[['HMHa']]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2,coef_HMH) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_HMH"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2","coef_HMH"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.HMH <- rbind(output.site, output.study.HMH)
} 

colnames(output.study.HMH) <- c('site','species','model','prop.change.HMH','recruit_pc1','recruit_pc2','se1','se2','coef_HMH')

output.study.HMH <- output.study.HMH %>%
  rename(latin=species)%>%
  left_join(dat.myc, by=c("latin")) %>%
  left_join(dat.abund, by= c("site","latin")) %>%
  left_join(dat.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# Conspecific v conspecific:heterospecific
output.study.net.A<-merge(output.study.A, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod2") %>% # only keep model that has both A and H
  mutate(net.prop.change.A=prop.change.A-prop.change.H,
         net.coef.A= coef_A-coef_H)

# Conmyc v heterospecific
output.study.net.CM<-merge(output.study.CM, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod4") %>% # only keep model that has both CM and H
  mutate(net.prop.change.CM=prop.change.CM-prop.change.H)

# A  v  HMH; CMH v HMH
output.study.net.A.CM<- merge(output.study.A, output.study.CMH, by=c("site","latin","model","myc","relabund","relba","abund","ba")) #
output.study.net.A.CM<- merge(output.study.net.A.CM, output.study.HMH, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  mutate(net.prop.change.A=prop.change.A-prop.change.HMH,
         net.prop.change.CMH=prop.change.CMH-prop.change.HMH,
         net.coef.A= coef_A-coef_HMH,
         net.coef.CMH= coef_CMH-coef_HMH) 

write_csv(output.study.net.A,"data/GLMoutput.study.net.m1-5.A_10.15.21.csv")
write_csv(output.study.net.CM,"data/GLMoutput.study.net.m1-5.CM_10.15.21.csv")
write_csv(output.study.net.A.CM,"data/GLMoutput.study.net.m1-5.A.CM_10.15.21.csv")



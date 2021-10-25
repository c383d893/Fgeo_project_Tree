######## ANALYSES FOR INITIAL AM-EM TEST ##########
########## UPDATED: FORESTGEO 5.5.21 #############

# weights = 1/SE = aka inversely to its' variance (lower SE, higher weight)
# need to scale abundance if used.

####################
##### SET UP #######
####################

#load packages
library(lme4);library(lmerTest);library(ggeffects);library(ggplot2);library(tidyverse);library(lsmeans)

##################################
##################################
#### COMPARING COEFICCIENTS ######
##################################
##################################

##################################
####### ADULT SPECIES ONLY #######
##################################

# Compare cndd coef from glm including Hs and Ha or not:
scbi<-read.table("data/scbi.dat_A.csv",header=TRUE)  # read dat
scbi<-scbi %>% filter(!sp == "66")
# 1. Open jpeg file
png("figures/scbi.coef.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi, mapping = aes(x = cndd, y = cndd.noh))+
  geom_point()+
  #geom_errorbar(aes(xmin=cndd-cndd.se, xmax=cndd+cndd.se), width=.2)+
  #geom_errorbar(aes(ymin=cndd.noh-cndd.se.noh, ymax=cndd.noh+cndd.se.noh), width=.2)+
  labs(title = "cndd comparison species: heterospecific effects ",
       x = "cndd with hetero effects",
       y = "cndd without hetero effects") +
  xlim (-4,4) +
  ylim (-4,4) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()
#colour = myc, not very different

##################################
####### ADULT GENUS ONLY #########
##################################

# Compare cndd coef from glm including Hs and Ha or not:
scbi.g<-read.table("data/scbi.dat_A_genus.csv",header=TRUE)  # read dat
# 1. Open jpeg file
png("figures/scbi.g.coef.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.g, mapping = aes(x = cndd, y = cndd.noh))+
  geom_point()+
  #geom_errorbar(aes(xmin=cndd-cndd.se, xmax=cndd+cndd.se), width=.2)+
  #geom_errorbar(aes(ymin=cndd.noh-cndd.se.noh, ymax=cndd.noh+cndd.se.noh), width=.2)+
  labs(title = "cndd comparison genus: heterospecific effects ",
       x = "cndd with hetero effects",
       y = "cndd without hetero effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()
#colour = myc, not very different

#############################################
######## ADULT & CONMYC SPECIES ONLY ########
#############################################

scbi.cm<-read.table("data/scbi.dat_A_conmyc.csv",header=TRUE)  # read dat
#compare cndd for each model (pairs)
#only cndd
pairs(scbi.cm[ ,c(3,7,9,13)]) # all very close to 1:1 line

# cndd between conmyc and not
# 1. Open jpeg file
png("figures/scbi.cm.cnddcoef.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cm, mapping = aes(x = cndd.cm, y = cndd.A))+
  geom_point()+
  labs(title = "cndd comparison: conmyc effects",
       x = "cndd without conmyc effect",
       y = "cndd with conmyc effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#cmndd
# 1. Open jpeg file
png("figures/scbi.cm.cmnddcoef.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cm, mapping = aes(x = cmndd.cm, y = cmndd.all))+
  geom_point()+
  labs(title = "cmndd comparison: hetero effects",
       x = "cmndd without hetero effect",
       y = "cmndd with hetero effects") +
  xlim (-0.3,0.3) +
  ylim (-0.3,0.3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#compare magnitude of cmndd and cndd (within model S ~A +CmA+CmS)
# grouped boxplot
# 1. Open jpeg file
png("figures/scbi.cm.cnddvcmndd.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
scbi.cm.long <- scbi.cm %>% gather(type, value, c(9,11), factor_key=TRUE) #make data long format
ggplot(scbi.cm.long, aes(x=type, y=value, fill=myc)) + 
  geom_boxplot()
# 3. Close the file
dev.off()

#############################################
######## ADULT & CONMYC SPECIES ONLY ########
################### R2 ######################
#############################################

scbi.cmaa<-read.table("data/scbi.dat_A_conmycallad.csv",header=TRUE)  # read dat
# S~ A +cma to all adults
# S ~A +all adults

# cndd between conmyc and not
# 1. Open jpeg file
png("figures/scbi.cm.cnddcoefallad.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cmaa, mapping = aes(x = cndd, y = cnddnocm))+
  geom_point()+
  labs(title = "cndd comparison all adults: conmyc effects",
       x = "cndd without conmyc effect",
       y = "cndd with conmyc effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#compare magnitude of cmndd and cndd (within model S~ A+cmA +all adults)
# grouped boxplot
# 1. Open jpeg file
png("figures/scbi.cm.cnddvcmnddvallad.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
scbi.cmaa.long <- scbi.cmaa %>% gather(type, value, c(2,4,6), factor_key=TRUE) #make data long format
ggplot(scbi.cmaa.long, aes(x=type, y=value, fill=myc)) + 
  geom_boxplot()
# 3. Close the file
dev.off()

#############################################
######## ADULT & CONMYC SPECIES ONLY ########
############## R3: NO OFFSET ################
#############################################

scbi.cmaa.no<-read.table("data/scbi.dat_A_conmycallad_nooff.csv",header=TRUE)  # read dat
# S~ A +cma to all adults
# S ~A +all adults

# cndd between conmyc and not
# 1. Open jpeg file
png("figures/scbi.cm.cnddcoefallad.no.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cmaa.no, mapping = aes(x = cndd, y = cnddnocm))+
  geom_point()+
  labs(title = "cndd comparison all adults: conmyc effects",
       x = "cndd without conmyc effect",
       y = "cndd with conmyc effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#compare magnitude of cmndd and cndd (within model S~ A+cmA +all adults)
# grouped boxplot
# 1. Open jpeg file
png("figures/scbi.cm.cnddvcmnddvallad.no.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
scbi.cmaa.no.long <- scbi.cmaa.no %>% gather(type, value, c(2,4,6), factor_key=TRUE) #make data long format
ggplot(scbi.cmaa.no.long, aes(x=type, y=value, fill=myc)) + 
  geom_boxplot()
# 3. Close the file
dev.off()

#############################################
######## ADULT & CONMYC SPECIES ONLY ########
########### R4: ALL ADULT OFFSET ############
#############################################

scbi.cmaa.Aao<-read.table("data/scbi.dat_A_conmycallad_Aaoffset.csv",header=TRUE)  # read dat
# S~ A +cma to all adults
# S ~A +all adults

# cndd between conmyc and not
# 1. Open jpeg file
png("figures/scbi.cm.cnddcoefallad.Aao.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cmaa.Aao, mapping = aes(x = cndd, y = cnddnocm))+
  geom_point()+
  labs(title = "cndd comparison all adults: conmyc effects",
       x = "cndd without conmyc effect",
       y = "cndd with conmyc effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#compare magnitude of cmndd and cndd (within model S~ A+cmA +all adults)
# grouped boxplot
# 1. Open jpeg file
png("figures/scbi.cm.cnddvcmnddvallad.Aao.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
scbi.cmaa.Aao.long <- scbi.cmaa.Aao %>% gather(type, value, c(2,4,6), factor_key=TRUE) #make data long format
ggplot(scbi.cmaa.Aao.long, aes(x=type, y=value, fill=myc)) + 
  geom_boxplot()
# 3. Close the file
dev.off()

#############################################
######## ADULT & CONMYC SPECIES ONLY ########
############## R4: ALL ADULT OFFSET #########
#############################################

scbi.cmaa.Cmo<-read.table("data/scbi.dat_A_conmycallad_Cmaoffset.csv",header=TRUE)  # read dat
# S~ A +cma to all adults
# S ~A +all adults

# cndd between conmyc and not
# 1. Open jpeg file
png("figures/scbi.cm.cnddcoefallad.Cmo.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(data = scbi.cmaa.Cmo, mapping = aes(x = cndd, y = cnddnocm))+
  geom_point()+
  labs(title = "cndd comparison all adults: conmyc effects",
       x = "cndd without conmyc effect",
       y = "cndd with conmyc effects") +
  xlim (-3,3) +
  ylim (-3,3) +
  geom_abline(intercept = 0, slope = 1, colour = "blue")    # Add 1:1 line
# 3. Close the file
dev.off()

#compare magnitude of cmndd and cndd (within model S~ A+cmA +all adults)
# grouped boxplot
# 1. Open jpeg file
png("figures/scbi.cm.cnddvcmnddvallad.Cmo.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
scbi.cmaa.Cmo.long <- scbi.cmaa.Cmo %>% gather(type, value, c(2,4,6), factor_key=TRUE) #make data long format
ggplot(scbi.cmaa.Cmo.long, aes(x=type, y=value, fill=myc)) + 
  geom_boxplot()
# 3. Close the file
dev.off()


##################################
##################################
########## GLM MODELS ############
##################################
##################################

##################################
####### ADULT SPECIES ONLY #######
##################################

scbi<-read.table("data/scbi.dat_A.csv",header=TRUE)  # read dat
scbi<-scbi %>% filter(!sp == "66")

scbi.lm.myc<-lm(cndd~myc, weights = 1/cndd.se,dat=scbi)                   # lm
saveRDS(scbi.lm.myc,"2._data_analysis/scbi.lm.myc.rds")
summary(scbi.lm.myc)   
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) -0.13021    0.08065  -1.614   0.1185  
#mycEM       -0.06968    0.11794  -0.591   0.5598  
#mycNM       -1.40130    0.67419  -2.078   0.0477 *

ref<-lsmeans(scbi.lm.myc,pairwise~myc, data= scbi)
ref.table<-as.data.frame(ref$lsmeans)

# 1. Open jpeg file
png("figures/scbi.cndd.myctype.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(ref.table, aes(myc, lsmean)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  theme_bw()+
  scale_y_continuous(name = "CNDD")+
  xlab("Mycorrhizal type")+
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") 
# 3. Close the file
dev.off()

scbi.lm.myc.N<-lm(cndd~myc*abund, weights = 1/cndd.se,dat=scbi)               # model 2 with rarity
saveRDS(scbi.lm.myc.N,"2._data_analysis/scbi.lm.myc.N.rds")
summary(scbi.lm.myc.N)
#(Intercept) -2.380e-01  1.150e-01  -2.071   0.0493 *
#mycEM       -9.555e-03  1.809e-01  -0.053   0.9583  
#mycNM       -1.297e+00  6.798e-01  -1.907   0.0685 .
#abund           2.015e-05  1.529e-05   1.318   0.1999  
#mycEM:abund      2.902e-05  1.143e-04   0.254   0.8017  
#mycNM:abund            NA         NA      NA       NA 

datplot4 <- ggpredict(scbi.lm.myc.N, terms=c("abund","myc"))
plot(datplot4,add.data = TRUE)

par(mfrow=c(1,1), las=1) 

AM.dat4<-datplot4[datplot4$group == "AM",]
EM.dat4<-datplot4[datplot4$group == "EM",]
NM.dat4<-datplot4[datplot4$group == "NM",]

plot(AM.dat4$x,AM.dat4$predicted,type="l",col="#5F9EA04D",cex.lab=1, lwd =2, xlim=c(-2000,2000),ylim=c(-2000,2000),
     xlab="N",ylab="cndd")
points(EM.dat4$x,EM.dat4$predicted, type="l",col="#FFD7004D",cex.lab=1, lwd =2)
points(NM.dat4$x,NM.dat4$predicted, type="l",col="#B222224D",cex.lab=1, lwd =2)

polygon(c(AM.dat4$x,rev(AM.dat4$x)),c(AM.dat4$conf.low,rev(AM.dat4$conf.high)),col="#5F9EA04D", lty=0, lwd=0.05)
polygon(c(EM.dat4$x,rev(EM.dat4$x)),c(EM.dat4$conf.low,rev(EM.dat4$conf.high)),col="#FFD7004D", lty=0, lwd= 0.05)
polygon(c(NM.dat4$x,rev(NM.dat4$x)),c(NM.dat4$conf.low,rev(NM.dat4$conf.high)),col="#B222224D",lty=0, lwd= 0.05)

AM.scbi<-scbi[scbi$myc == "AM",]
EM.scbi<-scbi[scbi$my == "EM",]
NM.scbi<-scbi[scbi$my == "NM",]

points(AM.scbi$abund,AM.scbi$cndd, type="l",col="#FFD7004D",cex.lab=1, lwd =2)
points(EM.scbi$abund,EM.scbi$cndd, type="l",col="#FFD7004D",cex.lab=1, lwd =2)
points(NM.scbi$abund,NM.scbi$cndd, type="l",col="#B222224D",cex.lab=1, lwd =2)

legend(x=0.7, y =7.5, legend= c("AM","EM","NM"), col = c("cadetblue","gold","firebrick"),border = "black", pch=19)

##################################
######## ADULT GENUS ONLY ########
##################################

#CNDD: AMNM >AM; with N:AM>EM; AMNM>AM
#Rarity: NS

scbi.g<-read.table("data/scbi.dat_A_genus.csv",header=TRUE)   # read dat
scbi.g.lm.myc<-lm(cndd~myc, weights = 1/cndd.se,dat=scbi.g)                        # lm
saveRDS(scbi.g.lm.myc,"2._data_analysis/scbi.g.lm.myc.rds")
summary(scbi.g.lm.myc)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.15891    0.08466  -1.877   0.0778 .
#mycAMNM     -1.37363    0.72179  -1.903   0.0741 .                                # AMNM stronger CNDD than AM (NM influence)
#mycEM        0.14850    0.12464   1.192   0.2498                                  # EM NS weaker CNDD than AM

ref<-lsmeans(scbi.g.lm.myc,pairwise~myc, data= scbi.g)
ref.table<-as.data.frame(ref$lsmeans)

# 1. Open jpeg file
png("figures/scbi.cndd.g.myctype.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(ref.table, aes(myc, lsmean)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  theme_bw()+
  scale_y_continuous(name = "CNDD")+
  xlab("Mycorrhizal type")+
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") 
# 3. Close the file
dev.off()

scbi.g.lm.myc.N<-lm(cndd~myc*abund, weights = 1/cndd.se,dat=scbi.g)                    # model 2 with rarity
saveRDS(scbi.g.lm.myc.N,"2._data_analysis/scbi.g.lm.myc.N.rds")
summary(scbi.g.lm.myc.N)
#              Estimate Std. Error t value Pr(>|t|) 
#(Intercept) -2.851e-01  1.145e-01  -2.489   0.0250 *
#mycAMNM     -1.251e+00  7.004e-01  -1.786   0.0943 .                              # AMNM stronger CNDD than AM
#mycEM        4.182e-01  2.150e-01   1.945   0.0707 .                              # EM weaker CNDD than AM
#abund            2.441e-05  1.554e-05   1.571   0.1371  
#mycAMNM:abund           NA         NA      NA       NA  
#mycEM:abund     -7.832e-05  6.178e-05  -1.268   0.2242  


#############################################
######## ADULT & CONMYC SPECIES ONLY ########
#############################################

#just looking at S~A model == first models above at the species level, less datapoints
scbi.cm<-read.table("data/scbi.dat_A_conmyc.csv",header=TRUE)   # read dat
scbi.cma.lm.myc<-lm(cndd.A~myc, weights = 1/cndd.A.se ,dat=scbi.cm)                       # lm
summary(scbi.cma.lm.myc)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) -0.01842    0.02852  -0.646    0.527
#mycEM       -0.43130    0.26784  -1.610    0.126
saveRDS(scbi.cma.lm.myc,"2._data_analysis/scbi.cma.lm.myc.rds")

ref<-lsmeans(scbi.cma.lm.myc,pairwise~myc, data= scbi.cm)
ref.table<-as.data.frame(ref$lsmeans)

# 1. Open jpeg file
png("figures/scbi.cndd.conmycmyctype.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(ref.table, aes(myc, lsmean)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  theme_bw()+
  scale_y_continuous(name = "CNDD")+
  xlab("Mycorrhizal type")+
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") 
# 3. Close the file
dev.off()

##################################
######## MUTLI SITE GLM ##########
##################################

scbi<-read.table("data/scbi.dat_A.csv",header=TRUE)  # read dat
scbi<-scbi %>% filter(!sp == "66")

trc<-read.table("data/trc.dat_A.csv",header=TRUE)  # read dat

fgeo<-rbind(scbi,trc)
fgeo<-fgeo[!fgeo$myc=="NM",]
fgeo$abund<-scale(fgeo$abund)
fgeo.lm.myc<-lmer(cndd~myc + (1|site), weights = 1/cndd.se, dat=fgeo)                   # lm
summary(fgeo.lm.myc)

#Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)  -0.1081     0.1435 47.0000  -0.753    0.455
#mycEM        -0.0755     0.2213 47.0000  -0.341    0.734
#mycNM        -0.2058     0.6483 47.0000  -0.318    0.752


ref<-lsmeans(fgeo.lm.myc,pairwise~myc, data= fgeo)
ref.table<-as.data.frame(ref$lsmeans)

# 1. Open jpeg file
png("figures/scbi.cndd.scbitrc.jpg", width=6, height= 6, units='in', res=300)
# 2. Create the plot
ggplot(ref.table, aes(myc, lsmean)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.5,size=1, position=position_dodge(1)) + 
  theme_bw()+
  scale_y_continuous(name = "CNDD")+
  xlab("Mycorrhizal type")+
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") 
# 3. Close the file
dev.off()

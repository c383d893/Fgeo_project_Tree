## This script creates all null models for all sites and species 

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
############### SET UP FUNCTIONS ##################
###################################################
###################################################

source("functions/dispersal1.R")
source("functions/dispersal2.R")
source("functions/dispersal3.R")
source("functions/DistWeighted.R")

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

dat2<- readRDS('data/allcensusdata_withmyc.rds')

site = unique(dat2$site)                                                       # unique sites

###################################################
###################################################
############### NULL ITTERATIONS ##################
###################################################
###################################################

null.itt = 200	 # set null reps	

registerDoParallel(3) 
#registerDoParallel(cores=(Sys.getenv("SLURM_CPUS_PER_TASK"))) # for use on the cluster
system.time({master.list <- foreach(z = 1:null.itt, .combine='c', .multicombine=TRUE)%dopar%{    
  itt.list<-list()
  for(k in site) {
      site.list <- list()
      lp.dat = filter(dat2, site == k)                                       # Make loop dataframe
      print(k)                                                               # Print site
      # Keep only main stems
      lp.dat<- lp.dat %>% group_by(treeID) %>% arrange(dbh) %>% slice(1) %>% ungroup() %>% as.data.frame()
      # record values    
      sp = lp.dat$latin 				              # Species ID for each stem
      gx = lp.dat$gx					                # X-coordinate of each stem
      gy = lp.dat$gy				         	        # Y-coordinate of each stem
      dbh =lp.dat$dbh                         # DBH
      species = unique(lp.dat$latin)		      # Species in the plot
      SP = length(species)				            # Number of species in the plot
      Lx = max(lp.dat$gx) ; Ly = max(lp.dat$gy); Area = (Lx*Ly)/10000	  # East-west width (Lx) and north-south height (Ly) of plot along with number of total ha of plot
      myc = lp.dat$myc                        # Myc category
      # Make model inputs
      rm(lp.dat.ran)
      lp.dat.ran = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type)
      # Run model per species:
      species_names<-vector()
      for(i in 1:SP) {
        if(lp.dat.ran$N[i] > tr & lp.dat.ran$adultquads[i] > 9 & lp.dat.ran$sapquads[i] > 9) {
          tbl = data.frame(A = lp.dat.ran$adultDistWeightedAbund[,i],ln_A =log(lp.dat.ran$adultDistWeightedAbund[,i]),Ha = lp.dat.ran$HadultDistWeightedAbund[,i], Hs = lp.dat.ran$HsapDistWeightedAbund[,i],
                           CMa = lp.dat.ran$CMadultDistWeightedAbund[,i], CMs = lp.dat.ran$CMsapDistWeightedAbund[,i], CMHa = lp.dat.ran$CMHadultDistWeightedAbund[,i],
                           CMHs = lp.dat.ran$CMHsapDistWeightedAbund[,i],HMHa = lp.dat.ran$HMHadultDistWeightedAbund[,i], HMHs = lp.dat.ran$HMHsapDistWeightedAbund[,i], S = lp.dat.ran$saplings[,i])
          # model 1
          mod1 = gam(S ~ s(A,k=3), offset=(ln_A), family= "poisson", data = tbl) 
          # model 2
          mod2 = gam(S ~ s(A,k=3)+ s(Ha,k=3),offset=(ln_A), family= "poisson",data = tbl) 
          # model 3
          mod3 = gam(S ~ s(A,k=3)+ s(CMa,k=3),offset=(ln_A),family= "poisson",data = tbl) 
          # model 4    
          mod4 = gam(S ~ s(A,k=3)+ s(CMa,k=3)+ s(Ha,k=3), offset=(ln_A),family= "poisson", data = tbl) 
          # model 5
          mod5 = gam(S ~ s(A,k=3)+ s(CMHa,k=3)+ s(HMHa,k=3), offset=(ln_A),family= "poisson", data = tbl) 
          # save site.list
          site.list[[as.character(species[i])]] <- list(mod1,mod2,mod3,mod4,mod5) # put model outputs in list
          names(site.list[[as.character(species[i])]])<-c("mod1","mod2","mod3","mod4","mod5") # name models for this species
        }
      }
      itt.list[[k]] <- site.list
      return(itt.list)
    }
   }
  })

#save model list as .rds object.
saveRDS(master.list, 'data/Rand_GAMmodeloutputs_m1-5_10.21.rds')

print('GAM Null models complete')

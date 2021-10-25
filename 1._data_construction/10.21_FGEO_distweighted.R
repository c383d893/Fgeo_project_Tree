###################################################
###################################################
################# SET UP PACKAGES #################
###################################################
###################################################

library(tidyverse) # tidyverse

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
################## BRING IN DAT ###################
###################################################
###################################################
###################################################

dat<- readRDS('data/allcensusdata_withmyc.rds')   # load data
dat<- readRDS('data/allcensusdata_withmyc_Amacayacu.rds')   # load data

site = unique(dat$site)                           # unique sites

###################################################
###################################################
################ RUN DISTWEIGHTED #################
###################################################
###################################################
###################################################

dist.weighted.list <- list()
for (k in site) { 
  lp.dat = filter(dat, site == k)                                       # Make loop dataframe
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
  lp.dat = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
  dist.weighted.list[[k]]<- lp.dat
}

#save site list as .rds object.
saveRDS(dist.weighted.list, 'data/density_dependence.input.data.rds')


  
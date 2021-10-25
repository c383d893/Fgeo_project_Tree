# FOR LOOP VERSION OF J.LAMANNA'S ADAPTED CODE
# UDPATED 10.15.21 BY C.DELAVAUX

####################
### !!random to sort:
####################

# set difftest for SE between two values
difftest_sd <- function(x1, x2, model){
  vardiff <- ((pred[[2]][1])^2 + 
                (pred[[2]][2])^2) - (2*(vcov(model)[x1, x2])) 
  diffse <- sqrt(vardiff)
  return(diffse=round(diffse, digits =4))
}

#DIVIDE BY SQRT Neff

# models to run:
#1.S ~ consp
#2.S ~ consp + heterosp (all)
#3.S ~ consp + conmyc(all)
#4.S ~ consp + conmyc(all) + heterosp(all)
#5.S ~ consp + hetersp-conmyc + heterosp-heteromyc

# Code:
# Conspecific = C
# Heterospecific = H
# Conmyc_all = CM
# Conmyc_het = CMH
# Hetmyc_het = HMH

###################################################
###################################################
################# SET UP PACKAGES #################
###################################################
###################################################

library(tidyverse) # tidyverse
library(mgcv)      # gam package
library(tidymv)    # predict_gam()

###################################################
###################################################
############### SET UP FUNCTIONS ##################
###################################################
###################################################

# Clark's 2Dt seed dispersal kernel ; the null model dispersal kernel/ seed rain shape from crown out
# Alpha = mean dispersal distance
dispersal.fun = function(x, y, alpha, n) {
  seedp = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))*2*pi*r}	# Clark's 2Dt p = 1 (integrated in 2 polar dimensions)
  disp.dists = seq(0, 10000, by = 1)
  rho <- sample(disp.dists, size = n, replace = TRUE, prob = seedp(r = disp.dists, alpha = alpha))
  theta <- runif(n, 0, 2*pi)
  x1 <- (rho * cos(theta)) + x
  y1 <- (rho * sin(theta)) + y
  result = data.frame(x1, y1)
  return(result)
}

# If seed disperses our of edge of focal area, it comes back through the opposite side
dispersal.fun2 = function(loc, xlim = plotwidth, ylim = plotheight, alpha) {	# Recruits disperse across plot edges in a torus
  
  test = dispersal.fun(loc[1], loc[2], alpha = alpha, n = 1)
  test[,1] = test[,1] %% xlim		# Torus
  test[,2] = test[,2] %% ylim
  return(data.frame(x1 = test[,1], y1 = test[,2]))
}

# Thomson null: Relationship between max tree height and mean dispersal distance (plus observed error around regression fit)
# Determining distance of spread based on height
mean.disp.dist = function(max.height) {10 ^ ((0.1975 + (0.936 * log10(max.height))) + rnorm(length(max.height), mean = 0, sd = 0.245))}	   	# Relationship between maximum tree height and mean dispersal distance with observed error from Thomson et al 2011 J. of Ecol.
clark.alpha = function(mean.disp, shape = 1) {exp(log(mean.disp) - (0.5*log(pi) + log(gamma(shape - 1/2)) - log(2) - log(gamma(shape))))^2}	# Calculate alpha from Clark's 2Dt model given mean dispersal distance

# DBH-height allometry function from ForestGEO CTFS Package
# Determining height based on dbh
predht.asym=function(dbh,param)
{
  if(is.null(dim(param)))
  {
    ymax=param[1]
    a=param[2]
    b=param[3]
  }
  else
  {
    ymax=param[,1]
    a=param[,2]
    b=param[,3]
  }
  
  return(ymax*(1-exp(-a*dbh^b)))
}
htparam=c(41.7, 0.057, 0.748) # default CTFS height parameters from BCI data

# Distance weighted function
# Takes forestgeo df and calculates distance weighted abundance/BA around focal seedling
# Generates for conspecific, conmycorrhizal (2 v), heterospecific (2 v)
DistWeighted = function(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type) {

  # quadrat grid
  x = seq(0, Lx, by = DX)                             # x loc
  y = seq(0, Ly, by = DX)                             # y loc
  n = length(x) - 1                                   # Number of x location minus 1
  m = length(y) - 1                                   # Number of y location minus 1
  
  # inizialize variables
  species = unique(sp)                                # Unique species
  S = length(species)                                 # Number in species
  saplings = matrix(0, n*m, S)                        # Conspecific saps
  adults = matrix(0, n*m, S)                          # Conspecific adults
  Hsaps = matrix(0, n*m, S)                           # Het saps
  Hadults = matrix(0, n*m, S)                         # Het adults
  CMsaps = matrix(0, n*m, S)                          # Conmycorrhizal saps all
  CMadults = matrix(0, n*m, S)                        # Conmycorrhizal adults all
  CMHsaps = matrix(0, n*m, S)                         # Conmycorrhizal hetero saps
  CMHadults = matrix(0, n*m, S)                       # Conmycorrhizal hetero adults
  HMHsaps = matrix(0, n*m, S)                         # Heteromycorrhizal hetero saps
  HMHadults = matrix(0, n*m, S)                       # Heteromycorrhizal hetero adults
  adultDistWeightedAbund = matrix(0, n*m, S)          # Adult dist weighted abundance (conspecific)
  HsapDistWeightedAbund = matrix(0, n*m, S)           # Het sap dist weighted abundance
  HadultDistWeightedAbund = matrix(0, n*m, S)         # Het adult dist weighted abundance
  CMsapDistWeightedAbund = matrix(0, n*m, S)          # Conmyc sap dist weighted abundance
  CMadultDistWeightedAbund = matrix(0, n*m, S)        # Conmyc adult dist weighted abundance
  CMHsapDistWeightedAbund = matrix(0, n*m, S)         # Conmyc het sap dist weighted abundance
  CMHadultDistWeightedAbund = matrix(0, n*m, S)       # Conmyc het adult dist weighted abundance
  HMHsapDistWeightedAbund = matrix(0, n*m, S)         # Non Conmyc het sap dist weighted abundance: no conspecific
  HMHadultDistWeightedAbund = matrix(0, n*m, S)       # Non Conmyc het adult dist weighted abundance: no conspecific
  BasalArea = rep(NA, times = S)                      # BA (total) per species
  N = rep(NA, times = S)                              # N of each species in plot
  status = rep(NA, times = length(dbh))               # adult or sapling
  Dcutoff = rep(NA, times = S)                        # DBH cutoff (determined below)
  sapquads = rep(NA, times = S)                       # Number of quadrats in which saplings of that species occur
  adultquads = rep(NA, times = S)                     # Number of quadrats in which adults of that species occur
  ba = (pi / 4) * ((dbh / 100) ^ 2)	                	# Basal area (m^2)
  
  # For each species: determine cutoff for sap/adult; lowest is 2 cm, since we need 1 cm for saplings as a minimum
  for(s in 1:S) {                                                           # For each unique species,
    use = which(sp == species[s] & dbh > 0)                                 # Find where sp is species s and dbh > 0
    N[s] = length(use)                                                      # Sample size for species
    BasalArea[s] = pi / 4 * sum((dbh[use] / 100) ^ 2) / (Lx * Ly / 10000) 	# Species-level basal area sum (m^2 per ha)
    D = dbh[use]
    if(sum(D > 5)/length(D) < 0.2) {Dcutoff[s] = 2}		  # Record DBH cutoff point at 2 cm (fewer than 20% of individuals >= 5 cm DBH)
    if(sum(D > 5)/length(D) >= 0.2) {Dcutoff[s] = 5}		# Record DBH cutoff point at 5 cm (greater than 20% of individuals >= 5 cm DBH)
    if(sum(D > 10)/length(D) >= 0.2) {Dcutoff[s] = 10}	# Record DBH cutoff point at 10 cm (greater than 20% of individuals >= 10 cm DBH)
    status[use] = c("adult")
    status[use][which(dbh[use] <= Dcutoff[s])] = "sap"  # Call ind where dbh is less than or equal to cutoff saplings
  }
  
  for(s in 1:S) {
    if(N[s] > tr) {                                              # If greater than tr (total number of quadrats in which species occurs); only use species that have min # obs
      use = which(sp == species[s] & dbh > 0)                    # Conspecific species == conspecific and conmycorrhizal (by default)
      usehet = which(sp != species[s] & dbh > 0)                 # Heterspecific species == non conspecific
      usecm = which(dbh > 0 & myc == myc[s])                     # Conmyc species: all conmycorrhizal, including species
      usecmh = which(sp != species[s] & dbh > 0 & myc == myc[s])  # Non conspecific conmyc species: conmycorrhizal, excluding species
      usehmh = which(sp != species[s] & dbh > 0 & myc != myc[s])  # Heteromyc species == heterospecific but not conmycorrhizal
      x0 = gx[use]; y0 = gy[use]				                         # Conspecific locations 
      xH = gx[usehet]; yH = gy[usehet]			                     # Heterospecific locations 
      xCM = gx[usecm]; yCM = gy[usecm]			                     # Conmycorrhizal locations
      xCMH = gx[usecmh]; yCMH = gy[usecmh]			                 # Conmycorrhizal het locations 
      xHMH = gx[usehmh]; yHMH = gy[usehmh]			                 # Heteromyc het locations
      D = dbh[use]						                                   # List of DBH measurements for focal species
      # get adults and saplings position
      use1 = which(D <= Dcutoff[s])		           	# Conspecific saplings
      use2 = which(D > Dcutoff[s])		          	# Conspecific adults
      usehet1 = which(status[usehet] == "sap")	  # Heterospecific saplings, cutoff based on each species' cutoff
      usehet2 = which(status[usehet] == "adult")	# Heterospecific adults, cutoff based on each species' cutoff
      usecm1 = which(status[usecm] == "sap")	    # Conmycorrhizal saplings, cutoff based on each species' cutoff
      usecm2 = which(status[usecm] == "adult")   	# Conmycorrhizal adults, cutoff based on each species' cutoff
      usecmh1 = which(status[usecmh] == "sap")	  # Conmycorrhizal saplings, cutoff based on each species' cutoff: no conspecific
      usecmh2 = which(status[usecmh] == "adult")  # Conmycorrhizal adults, cutoff based on each species' cutoff: no conspecific
      usehmh1 = which(status[usehmh] == "sap")	    # Heterospecific saplings, cutoff based on each species' cutoff
      usehmh2 = which(status[usehmh] == "adult")	  # Heterospecific adults, cutoff based on each species' cutoff
      x1 = x0[use1]; y1 = y0[use1]			          # Conspecific sap locations
      x2 = x0[use2]; y2 = y0[use2]			          # Conspecific adults locations
      xH1 = xH[usehet1]; yH1 = yH[usehet1]		    # Heterospecific sap locations
      xH2 = xH[usehet2]; yH2 = yH[usehet2]		    # Heterospecific adults locations
      xCM1 = xCM[usecm1]; yCM1 = yCM[usecm1]		  # Conmycorrhizal sap locations
      xCM2 = xCM[usecm2]; yCM2 = yCM[usecm2]		  # Conmycorrhizal adults locations
      xCMH1 = xCMH[usecmh1]; yCMH1 = yCMH[usecmh1]		      # Conmycorrhizal sap locations: no conspecific
      xCMH2 = xCMH[usecmh2]; yCMH2 = yCMH[usecmh2]		      # Conmycorrhizal adults locations: no conspecific
      xHMH1 = xHMH[usehmh1]; yHMH1 = yHMH[usehmh1]		  # All sap locations
      xHMH2 = xHMH[usehmh2]; yHMH2 = yHMH[usehmh2]		  # All adults locations
      
      # Wiegand-Moloney Thomas Process null model (fixes adult locations to preserve habitat heterogeneity/preferences, then disperses young away from adults given Clark's 2dT dispersal kernel; allometric uses DBH-height allometry function from ForestGEO CTFS Package)
      if(type == 'constant.disp.null') {
        alpha = clark.alpha(30)			# Sets the mean dispersal distance for species, meaning mean dispersal distance will be set to 30 m across species, but individual recruits can disperse a wide range of dispersal distances (only the mean of the 2dT distribution is set to 30 m here)
        num.live = ceiling(length(use2)*(1-adult.mort)); num.dead = length(use2) - num.live   # Randomly select reproductive adults
        repro.adults <- sample(1:length(use2), size = (length(use1) + num.dead), replace = T) # Assigns saplings per adults
        saplocs = do.call('rbind', apply(data.frame(x2, y2)[repro.adults,], 1, dispersal.fun2, xlim = Lx, ylim = Ly, alpha = alpha)) # Determine sapling locations from reproducing adults
        live.adults = sample(1:length(use2), size = num.live, replace = F)                    # Adults to adults; no replacement
        new.adults = sample(1:nrow(saplocs), size = num.dead, replace = F)                    # Sapling to adults; no replacement
        x2 = c(x2[live.adults], saplocs[new.adults,1]); y2 = c(y2[live.adults], saplocs[new.adults,2])                            # Combine all locations new adults
        x1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),1]; y1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),2]	# All sapling locations						
      }
      if(type == 'allometric.disp.null') {
        max.height = predht.asym(dbh = max(D), param = htparam) # Dispersal set based on tree height
        mean.disp.distance = mean.disp.dist(max.height)	        # Find max height of tree
        alpha = clark.alpha(mean.disp.distance)                 # Get mean disp distance; changes every itteration
        num.live = ceiling(length(use2)*(1-adult.mort)); num.dead = length(use2) - num.live
        repro.adults <- sample(1:length(use2), size = (length(use1) + num.dead), replace = T)
        saplocs = do.call('rbind', apply(data.frame(x2, y2)[repro.adults,], 1, dispersal.fun2, xlim = Lx, ylim = Ly, alpha = alpha))
        live.adults = sample(1:length(use2), size = num.live, replace = F)
        new.adults = sample(1:nrow(saplocs), size = num.dead, replace = F)
        x2 = c(x2[live.adults], saplocs[new.adults,1]); y2 = c(y2[live.adults], saplocs[new.adults,2])
        x1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),1]; y1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),2]
      }
      
      for(i in 1:n) {
        use1 = which(x1 >= x[i] & x1 < (x[i] + DX))		        	# All conspecific saplings within a given plot column
        use2 = which(x2 >= x[i] & x2 < (x[i] + DX))			        # All conspecific adults within a given plot column
        usehet1 = which(xH1 >= x[i] & xH1 < (x[i] + DX))	    	# All heterospecific saplings within a given plot column
        usehet2 = which(xH2 >= x[i] & xH2 < (x[i] + DX))	    	# All heterospecific adults within a given plot column
        usecm1 = which(xCM1 >= x[i] & xCM1 < (x[i] + DX))	    	# All conmycorrhizal saplings within a given plot column
        usecm2 = which(xCM2 >= x[i] & xCM2 < (x[i] + DX))	    	# All conmycorrhizal adults within a given plot column
        usecmh1 = which(xCMH1 >= x[i] & xCMH1 < (x[i] + DX))	  # All conmycorrhizal saplings within a given plot column: no conspecific
        usecmh2 = which(xCMH2 >= x[i] & xCMH2 < (x[i] + DX))	 	# All conmycorrhizal adults within a given plot column: no conspecific
        usehmh1 = which(xHMH1 >= x[i] & xHMH1 < (x[i] + DX))	    	# All heteromycorrhizal saplings within a given plot column
        usehmh2 = which(xHMH2 >= x[i] & xHMH2 < (x[i] + DX))	    	# All heteromycorrhizal adults within a given plot column
        if(i == n) {							                            	# Correction for eastern-most plot boundary (so trees at border are not excluded)
          use1 = which(x1 >= x[i] & x1 <= (x[i] + DX))	    	  # All conspecific saplings within a given plot column
          use2 = which(x2 >= x[i] & x2 <= (x[i] + DX))		      # All conspecific adults within a given plot column
          usehet1 = which(xH1 >= x[i] & xH1 <= (x[i] + DX))		    # All heterospecific saplings within a given plot column
          usehet2 = which(xH2 >= x[i] & xH2 <= (x[i] + DX))		    # All heterospecific adults within a given plot column
          usecm1 = which(xCM1 >= x[i] & xCM1 <= (x[i] + DX))		  # All conmycorrhizal saplings within a given plot column
          usecm2 = which(xCM2 >= x[i] & xCM2 <= (x[i] + DX))		  # All conmycorrhizal adults within a given plot column
          usecmh1 = which(xCMH1 >= x[i] & xCMH1 <= (x[i] + DX))		# All conmycorrhizal saplings within a given plot column: no conspecific
          usecmh2 = which(xCMH2 >= x[i] & xCMH2 <= (x[i] + DX))	  # All conmycorrhizal adults within a given plot column: no conspecific
          usehmh1 = which(xHMH1 >= x[i] & xHMH1 <= (x[i] + DX))		  # All saplings within a given plot column
          usehmh2 = which(xHMH2 >= x[i] & xHMH2 <= (x[i] + DX))		  # All adults within a given plot column
          
        }
        dx = x[i] - x2 + (DX / 2)						        # Distance along x-axis from each conspecific adult to focal quadrat center
        dHx1 = x[i] - xH1 + (DX / 2)					    	# Distance along x-axis from each heterospecific sapling to focal quadrat center
        dHx2 = x[i] - xH2 + (DX / 2)					    	# Distance along x-axis from each heterospecific adult to focal quadrat center
        dCMx1 = x[i] - xCM1 + (DX / 2)						  # Distance along x-axis from each conmycorrhizal sapling to focal quadrat center
        dCMx2 = x[i] - xCM2 + (DX / 2)						  # Distance along x-axis from each conmycorrhizal adult to focal quadrat center
        dCMHx1 = x[i] - xCMH1 + (DX / 2)						# Distance along x-axis from each conmycorrhizal sapling to focal quadrat center: no conspecific
        dCMHx2 = x[i] - xCMH2 + (DX / 2)						# Distance along x-axis from each conmycorrhizal adult to focal quadrat center: no conspecific
        dHMHx1 = x[i] - xHMH1 + (DX / 2)						  # Distance along x-axis from each sapling to focal quadrat center
        dHMHx2 = x[i] - xHMH2 + (DX / 2)						  # Distance along x-axis from each adult to focal quadrat center
        for(j in 1:m) {
          LinInd = (i - 1) * m + j
          saplings[LinInd,s] = sum(y1[use1] >= y[j] & y1[use1] < (y[j] + DX))	      	      # Number of conspecific saplings within focal quadrat
          adults[LinInd,s] = sum(y2[use2] >= y[j] & y2[use2] < (y[j] + DX))		              # Number of conspecific adults within focal quadrat
          Hsaps[LinInd,s] = sum(yH1[usehet1] >= y[j] & yH1[usehet1] < y[j] + DX)	          # Number of heterospecific saplings within focal quadrat
          Hadults[LinInd,s] = sum(yH2[usehet2] >= y[j] & yH2[usehet2] < y[j] + DX)   	      # Number of heterospecific adults within focal quadrat
          CMsaps[LinInd,s] = sum(yCM1[usecm1] >= y[j] & yCM1[usecm1] < y[j] + DX)          	# Number of conmycorrhizal saplings within focal quadrat
          CMadults[LinInd,s] = sum(yCM2[usecm2] >= y[j] & yCM2[usecm2] < y[j] + DX)   	    # Number of conmycorrhizal adults within focal quadrat
          CMHsaps[LinInd,s] = sum(yCMH1[usecmh1] >= y[j] & yCMH1[usecmh1] < y[j] + DX)     	# Number of conmycorrhizal saplings within focal quadrat: no conspecific
          CMHadults[LinInd,s] = sum(yCMH2[usecmh2] >= y[j] & yCMH2[usecmh2] < y[j] + DX)   	# Number of conmycorrhizal adults within focal quadrat: no conspecific
          HMHsaps[LinInd,s] = sum(yHMH1[usehmh1] >= y[j] & yHMH1[usehmh1] < y[j] + DX)	          # Number of saplings within focal quadrat
          HMHadults[LinInd,s] = sum(yHMH2[usehmh2] >= y[j] & yHMH2[usehmh2] < y[j] + DX)   	    # Number of adults within focal quadrat
          if(j == m) {												                                              # Correction for northern-most plot boundary (so trees at border are not excluded)
            saplings[LinInd,s] = sum(y1[use1] >= y[j] & y1[use1] <= (y[j] + DX))		        # Number of conspecific saplings within focal quadrat
            adults[LinInd,s] = sum(y2[use2] >= y[j] & y2[use2] <= (y[j] + DX))	    	      # Number of conspecific adults within focal quadrat
            Hsaps[LinInd,s] = sum(yH1[usehet1] >= y[j] & yH1[usehet1] <= y[j] + DX)	        # Number of heterospecific saplings within focal quadrat
            Hadults[LinInd,s] = sum(yH2[usehet2] >= y[j] & yH2[usehet2] <= y[j] + DX)	      # Number of heterospecific adults within focal quadrat
            CMsaps[LinInd,s] = sum(yCM1[usecm1] >= y[j] & yCM1[usecm1] <= y[j] + DX)	      # Number of conmycorrhizal saplings within focal quadrat
            CMadults[LinInd,s] = sum(yCM2[usecm2] >= y[j] & yCM2[usecm2] <= y[j] + DX)	    # Number of conmycorrhizal adults within focal quadrat
            CMHsaps[LinInd,s] = sum(yCMH1[usecmh1] >= y[j] & yCMH1[usecmh1] <= y[j] + DX)	  # Number of conmycorrhizal saplings within focal quadrat: no conspecific
            CMHadults[LinInd,s] = sum(yCMH2[usecmh2] >= y[j] & yCMH2[usecmh2] <= y[j] + DX)	# Number of conmycorrhizal adults within focal quadrat: no conspecific
            HMHsaps[LinInd,s] = sum(yHMH1[usehmh1] >= y[j] & yHMH1[usehmh1] <= y[j] + DX)	      # Number of saplings within focal quadrat
            HMHadults[LinInd,s] = sum(yHMH2[usehmh2] >= y[j] & yHMH2[usehmh2] <= y[j] + DX)	    # Number of adults within focal quadrat
          }
          dy = y[j] - y2 + (DX / 2)					                  # Distance along y-axis from each conspecific adult to focal quadrat center
          dHy1 = y[j] - yH1 + (DX / 2)				                # Distance along y-axis from each heterospecific sapling to focal quadrat center
          dHy2 = y[j] - yH2 + (DX / 2)				                # Distance along y-axis from each heterospecific adult to focal quadrat center
          dCMy1 = y[j] - yCM1 + (DX / 2)			               	# Distance along y-axis from each conmycorrhizal sapling to focal quadrat center
          dCMy2 = y[j] - yCM2 + (DX / 2)				              # Distance along y-axis from each conmycorrhizal adult to focal quadrat center
          dCMHy1 = y[j] - yCMH1 + (DX / 2)				            # Distance along y-axis from each conmycorrhizal sapling to focal quadrat center: no conspecific
          dCMHy2 = y[j] - yCMH2 + (DX / 2)			             	# Distance along y-axis from each conmycorrhizal adult to focal quadrat center: no conspecific
          dHMHy1 = y[j] - yHMH1 + (DX / 2)			        	      # Distance along y-axis from each sapling to focal quadrat center
          dHMHy2 = y[j] - yHMH2 + (DX / 2)			                # Distance along y-axis from each adult to focal quadrat center
          r2 = (dx * dx) + (dy * dy)				                  # Squared distance from each conspecific adult to focal quadrat center
          r2H1 = (dHx1 * dHx1) + (dHy1 * dHy1)			          # Squared distance from each heterospecific sapling to focal quadrat center
          r2H2 = (dHx2 * dHx2) + (dHy2 * dHy2)			          # Squared distance from each heterospecific adult to focal quadrat center
          r2CM1 = (dCMx1 * dCMx1) + (dCMy1 * dCMy1)			      # Squared distance from each conmycorrhizal sapling to focal quadrat center
          r2CM2 = (dCMx2 * dCMx2) + (dCMy2 * dCMy2)			      # Squared distance from each conmycorrhizal adult to focal quadrat center
          r2CMH1 = (dCMHx1 * dCMHx1) + (dCMHy1 * dCMHy1)			# Squared distance from each conmycorrhizal sapling to focal quadrat center: no conspecific
          r2CMH2 = (dCMHx2 * dCMHx2) + (dCMHy2 * dCMHy2)			# Squared distance from each conmycorrhizal adult to focal quadrat center: no conspecific
          r2HMH1 = (dHMHx1 * dHMHx1) + (dHMHy1 * dHMHy1)			      # Squared distance from each sapling to focal quadrat center
          r2HMH2 = (dHMHx2 * dHMHx2) + (dHMHy2 * dHMHy2)			      # Squared distance from each adult to focal quadrat center
          adultDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2 + L ^ 2)) ^ 2) * 500    	      # Distance weighted formula for conspecific adults (using stem abundance) *Note: number at end of line (i.e. 500) is a scaling factor 
          HsapDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2H1 + L ^ 2)) ^ 2) * 500          # Distance weighted formula for heterospecific saplings (using stem abundance) 
          HadultDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2H2 + L ^ 2)) ^ 2) * 500        # Distance weighted formula for heterospecific adults (using stem abundance) 
          CMsapDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2CM1 + L ^ 2)) ^ 2) * 500        # Distance weighted formula for conmycorrhizal adults (using stem abundance) *Note: number at end of line (i.e. 500) is a scaling factor 
          CMadultDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2CM2 + L ^ 2)) ^ 2) * 500      # Distance weighted formula for conmycorrhizal adults (using stem abundance) *Note: number at end of line (i.e. 500) is a scaling factor 
          CMHsapDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2CMH1 + L ^ 2)) ^ 2) * 500      # Distance weighted formula for conmycorrhizal adults (using stem abundance) *Note: number at end of line (i.e. 500) is a scaling factor: non conspecific
          CMHadultDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2CMH2 + L ^ 2)) ^ 2) * 500    # Distance weighted formula for conmycorrhizal adults (using stem abundance) *Note: number at end of line (i.e. 500) is a scaling factor: non conspecific 
          HMHsapDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2HMH1 + L ^ 2)) ^ 2) * 500        # Distance weighted formula for saplings (using stem abundance) 
          HMHadultDistWeightedAbund[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2HMH2 + L ^ 2)) ^ 2) * 500      # Distance weighted formula for adults (using stem abundance) 
            }
      } 
      sapquads[s] = sum(saplings[,s] > 0) 
      adultquads[s] = sum(adults[,s] > 0)
    }
  }
  result = list(saplings, adultDistWeightedAbund, HsapDistWeightedAbund,  HadultDistWeightedAbund, CMsapDistWeightedAbund,  CMadultDistWeightedAbund,CMHsapDistWeightedAbund,  CMHadultDistWeightedAbund,  HMHsapDistWeightedAbund,HMHadultDistWeightedAbund, adults, species, Hsaps, Hadults, N, BasalArea, Dcutoff, sapquads, adultquads)
  names(result) = c("saplings", "adultDistWeightedAbund", "HsapDistWeightedAbund",  "HadultDistWeightedAbund", "CMsapDistWeightedAbund",  "CMadultDistWeightedAbund", "CMHsapDistWeightedAbund", "CMHadultDistWeightedAbund", "HMHsapDistWeightedAbund",  "HMHadultDistWeightedAbund","adults", "species", "Hsaps", "Hadults", "N", "BasalArea", "Dcutoff", "sapquads", "adultquads")
  return(result)
}

###################################################
###################################################
################## SET UP PARS ####################
###################################################
###################################################

quant = 0.5					                # Cutoff point for adults and saplings, this puts the cutoff point at the median DBH for each species (0.5 = 50th percentile). This was done to be comparable with Detto et al. 2019, but we recommend using a cutoff point specific to each species that is based on biological information
L = 2 * 20 / pi; 			            	# Mean dispersal = pi/2*L -- this sets the shape of the distance-weighted adult abundance function, not dispersal for the recruits in the diserpsal-kernel model
DX = 20					                    # Quadrat size (m), this equates to a 20x20 m quadrat size
adult.mort = 0.10				            # Adult mortality (proportion of adults that die and are replaced by saplings in the dispersal-kernel model)
abundance = 'stem'			            # Option to use stem abundance or basal area to calculate distance-weighted abundances ('stem' = stem abundances; 'ba' = basal area)
null.type = 'allometric.disp.null'	# Selection of null-model type: 'allometric.disp.null' = dispersal-kernel model with dispersal differing among species as a function of their max. height; 'constant.disp.null' = dispersal-kernel model with constant dispersal among species (set in function code above
tr=20						                    # Minimum abundance theshold for calculations (20 b/c CNDD only calculated for species with at least 20 individuals, 10 saplings in unique quadrats, and 10 adults in unique quadrats)

###################################################
###################################################
################## BRING IN DAT ###################
############# MYC SPECIES AND GENUS ###############
###################################################
###################################################

mycstat<-read.table("data/Funroot_species_08092021.csv",header=TRUE,sep=",")     # Species myc list

mycstat.g<-read.table("data/Funroot_genus_08092021.csv",header=TRUE,sep=",")     # Genus myc list

###################################################
###################################################
################## BRING IN DAT ###################
#################### PREP DAT #####################
###################################################
###################################################

# Read data and extract census 1 only and alive only.
dat<- read.csv("data/allcensusdata_9.20.2021.csv") %>%
  group_by(site) %>%
  slice_max(census) %>% # get latest census
  ungroup() %>%
  filter(status=="alive") %>%
  filter(!site %in% c("Traunstein", "Amacayacu"))  # remove Traunstein=managed forest; Amacayacu =TBD

# Add mycorrhizal category sp, then genus:
dat.myc<-dat %>%                           
  left_join(mycstat,by="latin")                                                  # Merge serc.dat with mycstat (species)

dat.myc.g<-dat %>% 
  left_join(mycstat.g,by="genus")%>%                                             # Merge sp list (alive) with by genus
  select(c("latin","myc")) %>%                                                   # Subset just latin and myc
  rename(myc.g=myc) %>%                                                          # Rename cols
  distinct(latin,.keep_all=TRUE)                                                 # Keep only on representative from each species

dat2<-dat.myc %>% left_join(dat.myc.g, by="latin")%>%                            # Merge sp myc list with genus myc list
  mutate(myc=as.character(myc),myc.g=as.character(myc.g)) %>%                    # Convert to characters
  mutate(consensus_myc = ifelse(is.na(myc), myc.g, myc)) %>%                     # if myc is na, assign myc.g, otherwise assigned, myc
  select(-c(myc,myc.g)) %>%
  rename(myc= consensus_myc) %>%
  filter(myc=="AM"| myc=="EM") %>% # only keep two cat
  as.data.frame() %>%
  select(-c("quadrat","census","status","genus","species"))

# Find sites that don't have both myc types: can't include for conmyc heteromyc
dat2.myc<- dat2 %>% 
    group_by(site,myc,.drop=FALSE) %>% 
    distinct(latin) %>% 
    summarize(N = n())%>%
    mutate(freq = N / sum(N)) %>%
    filter(freq==1)

##### CHOOSE ONE #####
# GLM 1-5: all
dat2 <- dat2 %>% filter(!site %in% dat2.myc$site)

# GLM 1-2: no extra info from GLM 1-5

# GAM 1-5: subset out those that can't run CM
dat2 <- dat2 %>% filter(!site %in% dat2.myc$site)

# GAM 1-2: all

site = unique(dat2$site)                                                       # unique sites

###################################################
###################################################
################## RUN ANALYSES ###################
###################################################
###################################################
###################################################

##################################
########## GLM MOD 1-5 ###########
#### STORE MODEL OUTPUTS LIST ####
##################################

mod.list <- list() # full list site, sp, model
for (k in site) { 
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
  lp.dat2 = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
  # Run model per species:
  species_names<-vector()
  for(i in 1:SP) {
    if(lp.dat2$N[i] > tr & lp.dat2$adultquads[i] > 9 & lp.dat2$sapquads[i] > 9) 
      {
      tbl = data.frame(A = lp.dat2$adultDistWeightedAbund[,i],ln_A =log(lp.dat2$adultDistWeightedAbund[,i]),Ha = lp.dat2$HadultDistWeightedAbund[,i], Hs = lp.dat2$HsapDistWeightedAbund[,i],
                       CMa = lp.dat2$CMadultDistWeightedAbund[,i], CMs = lp.dat2$CMsapDistWeightedAbund[,i], CMHa = lp.dat2$CMHadultDistWeightedAbund[,i],
                       CMHs = lp.dat2$CMHsapDistWeightedAbund[,i],HMHa = lp.dat2$HMHadultDistWeightedAbund[,i], HMHs = lp.dat2$HMHsapDistWeightedAbund[,i], S = lp.dat2$saplings[,i])
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
saveRDS(mod.list, 'data/GLMmodeloutputs_m1-5_10.15.21.rds')

##################################
########## GLM MOD 1-2 ###########
#### STORE MODEL OUTPUTS LIST ####
##################################

mod.list <- list() # full list site, sp, model
for (k in site) { 
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
  lp.dat2 = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
  # Run model per species:
  species_names<-vector()
  for(i in 1:SP) {
    if(lp.dat2$N[i] > tr & lp.dat2$adultquads[i] > 9 & lp.dat2$sapquads[i] > 9) 
    {
      tbl = data.frame(A = lp.dat2$adultDistWeightedAbund[,i],ln_A =log(lp.dat2$adultDistWeightedAbund[,i]),Ha = lp.dat2$HadultDistWeightedAbund[,i], Hs = lp.dat2$HsapDistWeightedAbund[,i],
                       CMa = lp.dat2$CMadultDistWeightedAbund[,i], CMs = lp.dat2$CMsapDistWeightedAbund[,i], CMHa = lp.dat2$CMHadultDistWeightedAbund[,i],
                       CMHs = lp.dat2$CMHsapDistWeightedAbund[,i],HMHa = lp.dat2$HMHadultDistWeightedAbund[,i], HMHs = lp.dat2$HMHsapDistWeightedAbund[,i], S = lp.dat2$saplings[,i])
      # model 1
      mod1 = glm(S ~ A ,offset=(ln_A), family = 'poisson', data = tbl) 
      # model 2
      mod2 = glm(S ~ A+ Ha ,offset=(ln_A), family = 'poisson', data = tbl) 
      # save site.list
      site.list[[as.character(species[i])]] <- list(mod1,mod2) # put model outputs in list
      names(site.list[[as.character(species[i])]])<-c("mod1","mod2") # name models for this species
    }
    mod.list[[k]] <- site.list 
  }
}

#save model list as .rds object.
saveRDS(mod.list, 'data/GLMmodeloutputs_m1-2_10.15.21.rds')

##################################
########## GAM MOD 1-5 ###########
#### STORE MODEL OUTPUTS LIST ####
##################################

mod.list <- list() # full list site, sp, model
for (k in site) { 
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
  lp.dat2 = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
  # Run model per species:
  species_names<-vector()
  for(i in 1:SP) {
    if(lp.dat2$N[i] > tr & lp.dat2$adultquads[i] > 9 & lp.dat2$sapquads[i] > 9) 
    {
      tbl = data.frame(A = lp.dat2$adultDistWeightedAbund[,i], ln_A =log(lp.dat2$adultDistWeightedAbund[,i]),Ha = lp.dat2$HadultDistWeightedAbund[,i], Hs = lp.dat2$HsapDistWeightedAbund[,i],
                       CMa = lp.dat2$CMadultDistWeightedAbund[,i], CMs = lp.dat2$CMsapDistWeightedAbund[,i], CMHa = lp.dat2$CMHadultDistWeightedAbund[,i],
                       CMHs = lp.dat2$CMHsapDistWeightedAbund[,i],HMHa = lp.dat2$HMHadultDistWeightedAbund[,i], HMHs = lp.dat2$HMHsapDistWeightedAbund[,i], S = lp.dat2$saplings[,i])
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
    mod.list[[k]] <- site.list 
  }
}

#save model list as .rds object.
saveRDS(mod.list, 'data/GAMmodeloutputs_m1-5_10.15.21.rds')

##################################
########## GAM MOD 1-2 ###########
#### STORE MODEL OUTPUTS LIST ####
##################################

mod.list <- list() # full list site, sp, model
for (k in site) { 
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
  lp.dat2 = DistWeighted(sp, dbh, myc, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
  # Run model per species:
  species_names<-vector()
  for(i in 1:SP) {
    if(lp.dat2$N[i] > tr & lp.dat2$adultquads[i] > 9 & lp.dat2$sapquads[i] > 9) 
    {
      tbl = data.frame(A = lp.dat2$adultDistWeightedAbund[,i], ln_A =log(lp.dat2$adultDistWeightedAbund[,i]),Ha = lp.dat2$HadultDistWeightedAbund[,i], Hs = lp.dat2$HsapDistWeightedAbund[,i],
                       CMa = lp.dat2$CMadultDistWeightedAbund[,i], CMs = lp.dat2$CMsapDistWeightedAbund[,i], CMHa = lp.dat2$CMHadultDistWeightedAbund[,i],
                       CMHs = lp.dat2$CMHsapDistWeightedAbund[,i],HMHa = lp.dat2$HMHadultDistWeightedAbund[,i], HMHs = lp.dat2$HMHsapDistWeightedAbund[,i], S = lp.dat2$saplings[,i])
      # model 1
      mod1 = gam(S ~ s(A,k=3), offset=(ln_A), family= "poisson", data = tbl) 
      # model 2
      mod2 = gam(S ~ s(A,k=3)+ s(Ha,k=3),offset=(ln_A), family= "poisson",data = tbl) 
      # save site.list
      site.list[[as.character(species[i])]] <- list(mod1,mod2) # put model outputs in list
      names(site.list[[as.character(species[i])]])<-c("mod1","mod2") # name models for this species
    }
    mod.list[[k]] <- site.list 
  }
}

#save model list as .rds object.
saveRDS(mod.list, 'data/GAMmodeloutputs_m1-2_10.15.21.rds')

##################################
########## GLM MOD 1-5 ###########
###### CALC PROP CHANGE 1-2 ######
####### CALC COEFFICIENTS ########
##################################

mod.list<-readRDS('data/GLMmodeloutputs_m1-5_10.15.21.rds')
mod.list['Luquillo'] <- NULL

# simplify dat2 for myc status later
dat2.myc<-dat2 %>% select(latin,myc) %>%
  distinct(latin,.keep_all=TRUE)

# count abund and rel.abund per species per site where dbh>0
dat2.abund<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(abund = n()) %>%
  mutate(relabund = abund / sum(abund)) 

# count ba and rel.ba per species per site where dbh>0
dat2.ba<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(ba = sum(dbh)) %>%
  mutate(relba = ba / sum(ba)) 

# determine CS quartile average (per quadrat)
output.study <- data.frame() 
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
    #Generate prediction dataframe.
    pred.dat <- as.data.frame(species.dat[[1]]$data[,c('A')]) # only take A
    colnames(pred.dat)<-"A"
    out.dat<-cbind(names(mod.list)[i],names(site.dat)[s],pred.dat$A)
    colnames(out.dat)<-c("site","species","A") # call col 
    output.species <- rbind(output.species, out.dat)
    #Add species level dataframe to site level dataframe.
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study <- rbind(output.site, output.study)
}  

output.study$A <- as.numeric(as.character(output.study$A))

str(output.study)
output.study.A<-output.study %>% group_by(site,species) %>% summarize(minA=min(A),medA=median(A),maxA=max(A))

range(output.study.A$maxA) # min of max is 0.62

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
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
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
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
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
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
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
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
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
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
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

##################################
########## GLM MOD 1-2 ###########
###### CALC PROP CHANGE 1-2 ######
##################################

mod.list<-readRDS('data/GLMmodeloutputs_m1-2_10.15.21.rds')

# simplify dat2 for myc status later
dat2.myc<-dat2 %>% select(latin,myc) %>%
  distinct(latin,.keep_all=TRUE)

# count abund and rel.abund per species per site where dbh>0
dat2.abund<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(abund = n()) %>%
  mutate(relabund = abund / sum(abund)) 

# count ba and rel.ba per species per site where dbh>0
dat2.ba<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(ba = sum(dbh)) %>%
  mutate(relba = ba / sum(ba)) 

tree.range<-c(0.01,0.62)

### calculate A effects
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
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-A)
      pred.dat <- cbind(A,pred.dat)
      pred<-predict.glm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.A <- rbind(output.site, output.study.A)
}  

colnames(output.study.A) <- c('site','species','model','prop.change.A','recruit_pc1','recruit_pc2','se1','se2')

output.study.A<- output.study.A %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate H effects
Ha <-  tree.range
output.study.H <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod2")] # subset model 2 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-Ha)
      pred.dat <- cbind(Ha,pred.dat)
      pred<-predict.glm(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.H <- rbind(output.site, output.study.H)
}  

colnames(output.study.H) <- c('site','species','model','prop.change.H','recruit_pc1','recruit_pc2','se1','se2')

output.study.H<- output.study.H %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# Conspecific v conspecific:heterospecific
output.study.net.A<-merge(output.study.A, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod2") %>% # only keep model that has both A and H
  mutate(net.prop.change.A=prop.change.A-prop.change.H,
         net.coef.A= coef_A-coef_H) 

# write out.
write_csv(output.study.net.A,"data/GLM.output.study.net.m1-2.A_10.15.21.csv")

##################################
########## GAM MOD 1-5 ###########
###### CALC PROP CHANGE 1-2 ######
##################################
#predict.gam(): Note that, in common with other prediction functions, 
#any offset supplied to gam as an argument is always ignored when predicting,
#unlike offsets specified in the gam model formula.
#but if the pred df gives A when offset=log(A), offset is effectively changed.

mod.list<-readRDS('data/GAMmodeloutputs_m1-5_10.15.21.rds')
mod.list[['Luquillo']]<-NULL # remove Luquillo; no EMs made it through

# simplify dat2 for myc status later
dat2.myc<-dat2 %>% select(latin,myc) %>%
  distinct(latin,.keep_all=TRUE)

# count abund and rel.abund per species per site where dbh>0
dat2.abund<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(abund = n()) %>%
  mutate(relabund = abund / sum(abund)) 

# count ba and rel.ba per species per site where dbh>0
dat2.ba<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(ba = sum(dbh)) %>%
  mutate(relba = ba / sum(ba)) 

tree.range<-c(0.01,0.61)

### calculate A effects
A <- tree.range
output.study.A <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat = mod.list[[i]] 
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat = site.dat[[s]]
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-A)
      pred.dat <- cbind(A,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
                                         mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.A <- rbind(output.site, output.study.A)
}  

colnames(output.study.A) <- c('site','species','model','prop.change.A','recruit_pc1','recruit_pc2','se1','se2')

output.study.A<- output.study.A %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate CM effects
CMa <-  tree.range
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
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-CMa)
      pred.dat <- cbind(CMa,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.CM <- rbind(output.site, output.study.CM)
}  

colnames(output.study.CM) <- c('site','species','model','prop.change.CM','recruit_pc1','recruit_pc2','se1','se2')

output.study.CM<- output.study.CM %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate CMH effects
CMHa <-  tree.range
output.study.CMH <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod5")] # subset model 3 and 4 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-CMHa)
      pred.dat <- cbind(CMHa,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.CMH <- rbind(output.site, output.study.CMH)
}  

colnames(output.study.CMH) <- c('site','species','model','prop.change.CMH','recruit_pc1','recruit_pc2','se1','se2')

output.study.CMH<- output.study.CMH %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate H effects
Ha <-  tree.range
output.study.H <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod2","mod4")] # subset model 3 and 4 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-Ha)
      pred.dat <- cbind(Ha,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.H <- rbind(output.site, output.study.H)
}  

colnames(output.study.H) <- c('site','species','model','prop.change.H','recruit_pc1','recruit_pc2','se1','se2')

output.study.H<- output.study.H %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate HMHa effects
HMHa <-  tree.range 
output.study.HMH <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat = site.dat[[s]]
    species.dat<-species.dat[c("mod5")] # subset model 3 and 4 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-HMHa)
      pred.dat <- cbind(HMHa,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.HMH <- rbind(output.site, output.study.HMH)
}  

colnames(output.study.HMH) <- c('site','species','model','prop.change.HMH','recruit_pc1','recruit_pc2','se1','se2')

output.study.HMH<- output.study.HMH %>%
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# Conspecific v conspecific:heterospecific
output.study.net.A<-merge(output.study.A, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod2") %>% # only keep model that has both A and H
  mutate(net.prop.change.A=prop.change.A-prop.change.H)

# Conmyc v heterospecific
output.study.net.CM<-merge(output.study.CM, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod4") %>% # only keep model that has both CM and H
  mutate(net.prop.change.CM=prop.change.CM-prop.change.H)

# A  v  HMH; CMH v HMH
output.study.net.A.CM<- merge(output.study.A, output.study.CMH, by=c("site","latin","model","myc","relabund","relba","abund","ba")) #
output.study.net.A.CM<- merge(output.study.net.A.CM, output.study.HMH, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
   mutate(net.prop.change.A=prop.change.A-prop.change.HMH,
          net.prop.change.CMH=prop.change.CMH-prop.change.HMH)

# write out.
write_csv(output.study.net.A,"data/GAM.output.study.net.m1-5.A_10.15.21.csv")
write_csv(output.study.net.CM,"data/GAM.output.study.net.m1-5.CM_10.15.21.csv")
write_csv(output.study.net.A.CM,"data/GAM.output.study.net.m1-5.A.CM_10.15.21.csv")

##################################
########## GAM MOD 1-2 ###########
###### CALC PROP CHANGE 1-2 ######
####### CALC COEFFICIENTS ########
##################################

mod.list<-readRDS('data/GAMmodeloutputs_m1-2_10.15.21.rds')
mod.list['Luquillo']<-NULL

# simplify dat2 for myc status later
dat2.myc<-dat2 %>% select(latin,myc) %>%
  distinct(latin,.keep_all=TRUE)

# count abund and rel.abund per species per site where dbh>0
dat2.abund<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(abund = n()) %>%
  mutate(relabund = abund / sum(abund)) 

# count ba and rel.ba per species per site where dbh>0
dat2.ba<- dat2 %>% filter(dbh>0) %>%
  group_by(site, latin) %>%
  summarise(ba = sum(dbh)) %>%
  mutate(relba = ba / sum(ba)) 

tree.range<-c(0.01,0.76)

### calculate A effects
A <- tree.range
output.study.A <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  print(i)
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-A)
      pred.dat <- cbind(A,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.A <- rbind(output.site, output.study.A)
}  

colnames(output.study.A) <- c('site','species','model','prop.change.A','recruit_pc1','recruit_pc2','se1','se2')

output.study.A<- output.study.A %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# calculate H effects
Ha <-  tree.range
output.study.H <- data.frame()    
for(i in 1:length(mod.list)){
  #grab species fit list for a given site.
  site.dat= mod.list[[i]]  
  #For each species in that site list....
  output.site <- data.frame()
  for(s in 1:length(site.dat)) {
    #grab the list of model fits for that given species, within that given site.
    species.dat= site.dat[[s]]
    species.dat<-species.dat[c("mod2")] # subset model 2 only
    #for each of those model fits....
    output.species <- data.frame() #empty dataframe to store results.
    for(m in 1:length(species.dat)) {
      #Generate prediction dataframe.
      pred.dat <- data.frame(t(colMeans(species.dat[[m]]$model))) %>% select(-Ha)
      pred.dat <- cbind(Ha,pred.dat)
      pred<-predict.gam(species.dat[[m]],newdata = pred.dat,se=TRUE)
      #calculate proportion change
      recruitment_per_capita1<-exp(pred[[1]][1])
      recruitment_per_capita2<-exp(pred[[1]][2])
      prop.change <- (exp(pred[[1]][2])/ exp(pred[[1]][1])) - 1
      se1<-pred[[2]][1]
      se2<-pred[[2]][2]
      #wrap output, store in species level dataframe.
      out.dat<- cbind(names(mod.list)[i],names(site.dat)[s],names(species.dat)[m],prop.change,recruitment_per_capita1,recruitment_per_capita2,se1,se2) # add site, species, and model to prop
      output.species <- rbind(output.species, out.dat)
    }
    #Add species level dataframe to site level dataframe.
    output.species <- output.species %>% mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.character) %>%
      mutate_at(c("prop.change","recruitment_per_capita1","recruitment_per_capita2","se1","se2"), as.numeric)
    output.site <- rbind(output.site, output.species)
  }
  #add site level output to study level output.
  output.study.H <- rbind(output.site, output.study.H)
}  

colnames(output.study.H) <- c('site','species','model','prop.change.H','recruit_pc1','recruit_pc2','se1','se2')

output.study.H<- output.study.H %>% 
  rename(latin=species)%>%
  left_join(dat2.myc, by=c("latin")) %>%
  left_join(dat2.abund, by= c("site","latin")) %>%
  left_join(dat2.ba, by= c("site","latin")) %>%
  drop_na() %>%
  select(-c('se1','se2','recruit_pc1','recruit_pc2'))

# Conspecific v conspecific:heterospecific
output.study.net.A<-merge(output.study.A, output.study.H, by=c("site","latin","model","myc","relabund","relba","abund","ba")) %>%
  filter(model=="mod2") %>% # only keep model that has both A and H
  mutate(net.prop.change.A=prop.change.A-prop.change.H) 

# write out.
write_csv(output.study.net.A,"data/GAM.output.study.net.m1-2.A_10.15.21.csv")

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


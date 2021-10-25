# If seed disperses our of edge of focal area, it comes back through the opposite side
dispersal.fun2 = function(loc, xlim = plotwidth, ylim = plotheight, alpha) {	# Recruits disperse across plot edges in a torus
  
  test = dispersal.fun(loc[1], loc[2], alpha = alpha, n = 1)
  test[,1] = test[,1] %% xlim		# Torus
  test[,2] = test[,2] %% ylim
  return(data.frame(x1 = test[,1], y1 = test[,2]))
}

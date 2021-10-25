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

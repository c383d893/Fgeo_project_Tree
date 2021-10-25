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
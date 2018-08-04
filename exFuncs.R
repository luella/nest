# NEST example functions


# Simulates data (s unif, mu normal)
simulateData.func = function(n, mu0, tau, smin, smax){
	mu=rnorm(n, mu0, tau)
  	s=runif(n, smin, smax)
  	eps=rnorm(n, 0, s)
  	x=mu+eps
  	y=list(x=x, mu=mu, s=s)
  	return(y)
}

# Finds 
findBest.func = function(metric.mat, xg, sg){
  minMetric = min(metric.mat)
  idx = which(metric.mat == minMetric, arr.ind = TRUE)
  y = list(hx = xg[idx[1]], hs = sg[idx[2]])
  return(y)
}
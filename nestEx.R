## Example code for using NEST
source("nestFuncs.R") 
source("exFuncs.R")

# simulated data 
n = 5000
mu0 = 3
tau = 1
smin = .1 
smax = 1.7
d = simulateData.func(n, mu0, tau, smin, smax)
x = d$x
s = d$s



## 1) Example of finding muhat.nest ----

# Parameters 
hx = .6
hs = .37

# NEST as described in paper, eqs 3.1 - 3.3
muhat.nest = nest.func(x, s, x, s, hx, hs)



## 2) Example of using cv+SURE to find hx, hs ----
# Described in paper, section 3.2, eqs 3.4, 3.5
# caution: very slow on large grids 

# grids
x.grid = seq(from = .5, to = .7, by = .1)
s.grid = seq(from = .35, to = .40, by = .05)
Nx = length(x.grid)
Ns = length(s.grid)

sure.mat = matrix(nrow = Nx, ncol = Ns)

Nk = 5 
testIdx = matrix(1:n, nrow=Nk)


# loop for x.grid
for (j in 1:Nx) {

	hx = x.grid[j]

	# loop for s.grid
	for (l in 1:Ns) {

		hs = s.grid[l]

    fhat.cv.nest = rep(NA, nrow = n)
  	dfhat.cv.nest = rep(NA, nrow = n)
  	ddfhat.cv.nest = rep(NA, nrow = n)

		# loop for kfolds
		for(k in 1:Nk){

			# split data
			test = testIdx[k,]
			x.train = x[-test]
			x.test = x[test]
			s.train = s[-test]
			s.test = s[test]

  		fhat.cv.nest[test] = dens.fhat.func2(x.test, s.test, x.train, s.train, hx, hs)
      dfhat.cv.nest[test] = dens.dfhat.func2(x.test, s.test, x.train, s.train, hx, hs)
      ddfhat.cv.nest[test] = dens.ddfhat.func2(x.test, s.test, x.train, s.train, hx, hs)

		} # end loop for kfolds

		sure.mat[j, l] = sure.func(s, fhat.cv.nest, dfhat.cv.nest, ddfhat.cv.nest)

	} # end loop for x.grid

} # end loop for s.grid


chosenBw = findBest.func(sure.mat, x.grid, s.grid)
chosenBw$hx
chosenBw$hs




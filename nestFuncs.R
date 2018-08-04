# Nest functions from the technical report: NONPARAMETRIC EMPIRICAL BAYES
# ESTIMATION ON HETEROGENEOUS DATA (2018, Fu, James & Sun)

# NEST ----
nest.func = function(g, gs, xjs, sjs, hx, hs){
  ## NEST estimator
  # NEST as described in paper, eqs 3.1 - 3.3
  ## Arguments:
  # g: grid where densities are evaluated
  # gs: sigma values corresponding to the chosen grid
  # xjs: observed values
  # sjs: known sds
  # hx: bandwidth for x, can be a value or vector
  # hs: bandwidth for sigma, can be a value or vector
  ## Values:
  # y: estimated effect sizes
  
  fhat.nest = dens.fhat.func2(x, s, x, s, hx, hs)
  dfhat.nest = dens.dfhat.func2(x, s, x, s, hx, hs)
  muhat.nest = tf.func(x, s, fhat.nest, dfhat.nest) 
  return(muhat.nest)
}


# density ----
dens.fhat.func2 = function(g, gs, xjs, sjs, hx, hs)
{
  ## kernel estimator for density
  ## Arguments:
  # g: grid where densities are evaluated
  # gs: sigma values corresponding to the chosen grid
  # xjs: observed values
  # sjs: known sds
  # hx: bandwidth for x, can be a value or vector
  # hs: bandwidth for sigma, can be a value or vector
  ## Values:
  # y: the densities evaluated at g
  
  n=length(xjs)
  m=length(g)
  y=rep(0, m)
  hxj=hx*sjs
  
  for (i in 1:m)
  {
    phigi=dnorm(g[i]-xjs, 0, hxj)
    phisi=dnorm(gs[i]-sjs, 0, hs)
    y[i]=sum(phisi*phigi)/sum(phisi)
  }
  return(y)
}


# derivative ---- 
dens.dfhat.func2 = function(g, gs, xjs, sjs, hx, hs)
{
  ## kernel estimator for density derivative
  ## Arguments:
  # g: grid where densities are evaluated
  # gs: sigma values corresponding to the chosen grid
  # xjs: observed values
  # sjs: known sds
  # hx: bandwidth for x, can be a value or vector
  # hs: bandwidth for sigma, can be a value or vector
  ## Values:
  # y: the densities evaluated at g
  
  n=length(xjs)
  m=length(g)
  y=rep(0, m)
  hxj=hx*sjs
  
  for (i in 1:m)
  {
    phigi=dnorm(g[i]-xjs, 0, hxj)
    phisi=dnorm(gs[i]-sjs, 0, hs)
    y[i]=sum(phisi*phigi*((xjs-g[i])/hxj^2))/sum(phisi)
  }
  return(y)
}

# second derivative kernel ----
dens.ddfhat.func2 = function(g, gs, xjs, sjs, hx, hs)
{
  ## kernel estimator for density derivative
  ## Arguments:
  # g: grid where densities are evaluated
  # gs: sigma values corresponding to the chosen grid
  # xjs: observed values
  # sjs: known sds
  # hx: bandwidth for x, can be a value or vector
  # hs: bandwidth for sigma, can be a value or vector
  ## Values:
  # y: the densities evaluated at g
  
  n=length(xjs)
  m=length(g)
  y=rep(0, m)
  hxj=hx*sjs
  
  for (i in 1:m)
  {
    phigi=dnorm(g[i]-xjs, 0, hxj)
    phisi=dnorm(gs[i]-sjs, 0, hs)
    y[i]=sum((1/hxj^2)*phisi*phigi*(((xjs-g[i])/hxj)^2 - 1))/sum(phisi)
  }
  return(y)
}



# Tweedie's formula ----
tf.func = function(x, s, fv, dfv)
{
  ## Tweedie's formula
  ## Arguments:
  # x: observations
  # s: SDs
  # fv: estimated density
  # dfv: estimated derivative
  # y: estimated effect sizes
  dl=dfv/fv
  y=x+s^2*dl
  return(y)
}


# Positive part TF ----
tf.pos.func = function(x, s, fv, dfv)
{
  ## Positive part Tweedie's formula
  # Use for sparse data. Described in paper, section 5.1, page 18 
  ## Arguments:
  # x: observations
  # s: SDs
  # fv: estimated density
  # dfv: estimated derivative
  # y: estimated effect sizes
  dl=dfv/fv
  temp = 1+s^2*dl/x
  y=(temp)*(temp>0)*x
  return(y)
}

# SURE function ----
sure.func = function(s, fv, dfv, ddfv){
  ## returns Stein's unbiased risk estimate
  # In paper, eq 3.4 and top of page 10
  ## Arguments
  # s: singleton or vector of sigma_i
  # fv: estimated density
  # dfv: estimated derivative
  # ddfv: estimated second derivative
  # y: estimated SURE
  y = mean(s^4*(2*fv*ddfv - dfv^2)/(fv^2))
  return(y)
}


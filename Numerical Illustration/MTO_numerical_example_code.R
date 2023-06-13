# MTO numerical example R code
# Last edited 06/06/2023


# *************The packages mvtnorm and nloptr must be installed.********************

library(mvtnorm)
library(nloptr)

# simcon does approximate calculation of containment probabilities (the probabilities on the left hand side of the inequalities (3.2)-(3.4) in the paper
# by simulation. This is to check on calculations done using the R package mvtnorm.
# The interval considered is [s,t]. The function returns probabilities for the intervals [-Inf,t] (prob_it), [s,t] (prob_st), [s,Inf] (prob_si).
# seed - random number generator seed
# nn - number of simulations
# mpar - values of model parameters
# spar values of structure parameters
# z - value of exogenous variables
simcon = function(seed,nn,mpar,spar,z,s,t, browse = F) {
  set.seed(seed)
  beta = mpar$beta;   alpha = mpar$alpha;   gamma1 = mpar$gamma1; gamma2 = mpar$gamma2
  a = spar$a; b = spar$b; c1 = spar$c1; c2 = spar$c2; d0 = spar$d0; d1 = spar$d1; s12 = spar$s12; s22 = spar$s22
  
  uvec = rmvnorm(nn,mean = c(0,0), sigma = rbind(c(1,s12),c(s12,s22)))
  
  y2 = rep(d0 + d1%*%z,nn) + uvec[,2]
  y1i = rep(b%*%z,nn) + a*y2 + uvec[,1]
  
  y1=rep(0.1,nn)
  y1[y1i <= c1] = 0
  y1[y1i > c1 & y1i <= c2] = 1
  y1[y1i > c2] = 2
  
  hold = rep(z%*%beta,nn) + alpha*y2
  
  prob_it = sum((y1 == 0)*(gamma1 - hold <= t))/nn + sum((y1 == 1)*(gamma2 - hold <= t))/nn      # (-Inf, t]
  prob_st = sum((y1 == 1)*(s < gamma1 - hold)*(gamma2 - hold <= t))/nn                         # [s, t]
  prob_si = sum((y1 == 1)*(s < gamma1 - hold))/nn + sum((y1 == 2)*(s < gamma2 - hold))/nn      # [s, Inf)
  
  py0 = sum(y1 == 0)/nn
  py1 = sum(y1 == 1)/nn
  py2 = sum(y1 == 2)/nn
  
  # average value of y2
  meany2 = mean(y2)
  
  return(list(c(prob_it,prob_st,prob_si),c(pnorm(t),pnorm(t)-pnorm(s),1-pnorm(s)), c(py0,py1,py2),meany2))
}

# conprob does calculation of containment probabilities (the probabilities on the left hand side of the inequalities (3.2)-(3.4) in the paper
# using bivariate normal distribution probabilities provided by the package mvtnorm
# The interval considered is [s,tt] returned as prob
# mpar - values of model parameters
# spar values of structure parameters
# z - value of exogenous variables
conprob = function(mpar,spar,z,s,tt) {
  beta = mpar$beta;   alpha = mpar$alpha;   gamma1 = mpar$gamma1; gamma2 = mpar$gamma2
  a = spar$a; b = spar$b; c1 = spar$c1; c2 = spar$c2; d0 = spar$d0; d1 = spar$d1; s12 = spar$s12; s22 = spar$s22
  
  meanv = rep(0,2)
  sigmav = rbind(c(a^2*s22+2*a*s12+1, a*alpha*s22+alpha*s12),c(a*alpha*s22+alpha*s12, alpha^2*s22))
  
  ahold = a*d0 + z%*%(b+a*d1)
  alphahold = alpha*d0 + z%*%(beta + alpha*d1)
  
  prob = 1
  
  if (abs(alpha) >= 0.00001) {
    
  if (is.infinite(s) & is.finite(tt)) {
    prob = pmvnorm(lower = c(-Inf,gamma1 - tt - alphahold), upper = c(c1 - ahold, Inf), mean = meanv, sigma = sigmav) + 
      pmvnorm(lower = c(c1 - ahold,gamma2 - tt - alphahold), upper = c(c2 - ahold, Inf), mean = meanv, sigma = sigmav)
  }
  
  if (is.finite(s) & is.finite(tt)) {
    if (tt-s <= gamma2 - gamma1) { prob = 0 } else {
      prob = pmvnorm(lower = c(c1 - ahold, gamma2 - tt - alphahold), upper = c(c2 - ahold, gamma1 - s -alphahold), mean = meanv, sigma = sigmav)
    }
  }

  if (is.finite(s) & is.infinite(tt)) {
    prob = pmvnorm(lower = c(c1 - ahold,-Inf), upper = c(c2 - ahold, gamma1 - s - alphahold), mean = meanv, sigma = sigmav) + 
      pmvnorm(lower = c(c2 - ahold,-Inf), upper = c(Inf, gamma2 - s - alphahold), mean = meanv, sigma = sigmav)
  }
    py=NULL
  } else { # alpha close to zero
    
    sd = sqrt(1+2*a*s12+a^2*s22)
    zbeta = z%*%beta
    
    py0 = pnorm((c1 - ahold)/sd)
    py1 = pnorm((c2 - ahold)/sd) - py0
    py2 = 1 - py0 - py1
    
    if (is.infinite(s) & is.finite(tt)) {
      prob = py0*(gamma1 - zbeta <= tt) + py1*(gamma2 - zbeta <= tt)
    }
    
    if (is.finite(s) & is.finite(tt)) {
      if (tt-s <= gamma2 - gamma1) { prob = 0 } else {
        prob = py1*(gamma1 - zbeta > s)*(gamma2 - zbeta <=tt)
      }
    }
    
    if (is.finite(s) & is.infinite(tt)) {
      prob = py1*(s < gamma1 - zbeta) + py2*(s < gamma2 - zbeta)
    }
    py = c(py0,py1,py2)
  }

  attr(prob,"msg") <- NULL;   attr(prob,"error") <- NULL
  
  return(list(prob=prob,py=py))
}

# ************** Uncomment and run commands below to compare probabilities calculated by simulation and using the mvtnorm package ******************

# mpar = list(beta = c(1,0),alpha = 1.5, gamma1 = -0.5, gamma2 = 0.5)
# spar = list(a = 1.5, b = c(1,0), c1 = -0.5, c2 = 0.5, d0 = 0, d1 = c(1,2), s12 = 1, s22 = 2)
# 
# simcon(1234,10000000,mpar,spar,c(-1,1),-0.8,1.0)
# 
# conprob(mpar,spar,c(-1,1),-Inf,1.0)
# conprob(mpar,spar,c(-1,1),-0.8,1.0)
# conprob(mpar,spar,c(-1,1),-0.8,Inf)
# 
# simcon(1234,1000000,mpar,spar,c(-1,1),-0.4,0.4)
# conprob(mpar,spar,c(-1,1),-0.4,0.4)
# 
# simcon(1234,1000000,mpar,spar,c(-1,1),-0.6,0.4)
# conprob(mpar,spar,c(-1,1),-0.6,0.4)
# 
# simcon(1234,10000000,mpar,spar,c(-1,1),-Inf,-0.4)
# conprob(mpar,spar,c(-1,1),Inf,-0.4)
# all agrees
# 
# mpar2 = list(beta = c(1,0),alpha = 1.5,gamma1 = -1.5,gamma2 = -1.4)
# 
# simcon(5678,10000000,mpar2,spar,c(-1,1),-0.8,1.0)
# conprob(mpar2,spar,c(-1,1),-Inf,1.0)
# conprob(mpar2,spar,c(-1,1),-0.8,1.0)
# conprob(mpar2,spar,c(-1,1),-0.8,Inf)
# # all agrees
# 
# simcon(1234,10000000,mpar,spar,c(-1,1),-0.49,0.55)
# 
# conprob(mpar,spar,c(-1,1),-0.49,0.55)
# 
# mpar3 = mpar; mpar3$alpha = 0
# simcon(1234,100000000,mpar3,spar,c(-1,1),0.4,Inf)
# conprob(mpar3,spar,c(-1,1),0.4,Inf)
# 
# mpar3 = mpar; mpar3$alpha = 0
# simcon(1234,100000000,mpar3,spar,c(-1,1),0.4,1.6)
# conprob(mpar3,spar,c(-1,1),0.4,1.6)
# 
# mpar3 = mpar; mpar3$alpha = 0
# simcon(1234,100000000,mpar3,spar,c(-1,1),-Inf,1.7)
# conprob(mpar3,spar,c(-1,1),-Inf,1.7)
# 
# mpar3 = mpar; mpar3$alpha = 0
# simcon(1234,100000000,mpar3,spar,c(-1,1),-Inf,-0.2)
# conprob(mpar3,spar,c(-1,1),-Inf,-0.2)
# 
# mpar3 = mpar; mpar3$alpha = 0.0001
# simcon(1234,100000000,mpar3,spar,c(-1,1),0.4,Inf)
# conprob(mpar3,spar,c(-1,1),0.4,Inf)
# 
# mpar = list(beta = c(1,0),alpha = 1.5, gamma1 = -0.5, gamma2 = 0.5)
# spar = list(a = 1.5, b = c(0,0), c1 = -1, c2 = 1, d0 = 0, d1 = c(1,2), s12 = 0, s22 = 0.00000000000001)
# simcon(1234,100000000,mpar,spar,c(-1,1),-0.4,1.4)
# conprob(mpar,spar,c(-1,1),-Inf,1.4)
# conprob(mpar,spar,c(-1,1),-0.4,1.4)
# conprob(mpar,spar,c(-1,1),-0.4,Inf)
# # all agree.
# 

# getsets: Creates a list of test sets which are all connected unions of nt-1 intervals that partition the real line, excluding the real line itself.
# sets are unions of intervales that partition [0,1] transformed onto he real line using the standard Gaussian quantile function.
getsets = function(nt) {
  testsets = matrix(NA,nrow = nt*(nt-1)/2,ncol=2)
  count = 0
  cells = seq(0,1,l = nt)
  for (i in (1:(nt-1))) {
    ss = cells[i]
    for (j in (i+1):(nt)) {
      count = count + 1
      tt = cells[j]
      testsets[count,] = c(ss,tt)
    }
  }
  testsets = qnorm(testsets[-(nt-1),]) # remove test set [0,1], map onto real line
  return(testsets)
}


# inset2 determines if te parameter value mpar in the identified set
#    when probabilities are delivered by a complete structure with parameter values spar.
# Test sets are all connected unions of nt-1 intervals that partition the real line as delivered by the function getsets
# The interval [-Inf,Inf] is excluded.
inset2 = function(mpar,spar,zvalues,nt) {
  testsets = getsets(nt)
  nz = nrow(zvalues)
  ins = TRUE
  mindisc = Inf
  for (it in (1:nrow(testsets))) {
    ss = testsets[it,1]
    tt = testsets[it,2]
    for (iz in 1:nz) {
      cprob = conprob(mpar,spar,zvalues[iz,],ss,tt)$prob
      gprob = pnorm(tt)-pnorm(ss)
      disc = gprob-cprob
      mindisc = min(mindisc, disc)
      if (disc<0) {
        ins = FALSE
        break()
      }
    }
    if (!ins) break()
  }
  return(list(ins,mindisc,mpar,spar,zvalues))
}

# Projections of identified sets obtained by optimization.

# ochoice_projection delivers a lower (lower = 1) or upper (lower = -1) bound on the interval which comprises the projection
# of an identified set in the ordered choice model for a particular parameter. To obtain a lower bound the algorithm COBYLA
# is used to minimize the value of the parameter under consideration by choice of all parameter values subject to the constraints that characterize the identified set.

# startvalue - start value for all parameters
# parindex - index (in the complete parameter list) of the parameter whose projection is be sought
# lower - plus 1 if a lower bound of a projection, minus 1 if upper bound
# spar - values of parameters of the structure
# zvalues - a matrix giving in rows the values of vector z
# constrained_index - locations of constrained parameters
# constrained_value - values of constrained parameters
# nt - number of cell boundaries used in determining test sets
# ub, lb - upper and lower bounds on all parameters, including constrained

# the parameters wrt which we optimize are arranged (gamma1, eta, alpha, beta) where eta = log(gamma2 - gamma1)
#      so gamma2 = gamma1 + exp(eta) enforcing the requirement gamma2 > gamma1.

# cfp true delivers bounds on counterfactual probabilities
# mcfp true delivers bounds on marginal effects of y2 on counterfactual probabilities
# if both are false then bounds on a parameter specified by parindex are calculated - results under this option tend to be unreliable.
# y1value, y2value - values of Y1 , Y2 at which counterfactual calculations are to be done.
# epsilon - a small value used to make constraints more demanding.
# algorithm, xtol_rel, printlevel, maxeval options supplied to the COBYLA algorithm

ochoice_projection = function(startvalue, parindex, lower, constrained_index, constrained_value, spar, zvalues, nt, ub=NULL, lb=NULL, 
                                algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1.0e-04, printlevel = 0, maxeval = 1, y2value = NULL,
                              cfp = FALSE, mcfp = FALSE,epsilon = 0, y1value = 0) {
  
  npar = length(startvalue) + length(constrained_index)
  testsets = getsets(nt)
  
  ineq_tol = rep(0, nrow(zvalues)*nrow(testsets)) # only effective using the global optimizing algorithm ISRES, not effective using COBYLA

  if (cfp) {
    qvalue = nloptr(x0 = startvalue,
                    eval_f = oc_obj_cfp,
                    eval_g_ineq = oc_ineq,
                    opts = list("algorithm" = algorithm,"xtol_rel" = xtol_rel,
                                "maxeval" = maxeval, "print_level" = printlevel, "ranseed" = "67230",
                                "tol_constraints_ineq" = ineq_tol),
                    parindex = parindex,
                    lower = lower,
                    npar = npar,
                    spar = spar,
                    zvalues = zvalues,
                    constrained_index = constrained_index,
                    constrained_value = constrained_value,
                    testsets = testsets,
                    ub = ub[-constrained_index],
                    lb = lb[-constrained_index],
                    y2value = y2value,
                    epsilon = epsilon,
                    y1value = y1value)
  }
  if (mcfp) {
    qvalue = nloptr(x0 = startvalue,
                    eval_f = oc_obj_mcfp,
                    eval_g_ineq = oc_ineq,
                    opts = list("algorithm" = algorithm,"xtol_rel" = xtol_rel,
                                "maxeval" = maxeval, "print_level" = printlevel, "ranseed" = "67230",
                                "tol_constraints_ineq" = ineq_tol),
                    parindex = parindex,
                    lower = lower,
                    npar = npar,
                    spar = spar,
                    zvalues = zvalues,
                    constrained_index = constrained_index,
                    constrained_value = constrained_value,
                    testsets = testsets,
                    ub = ub[-constrained_index],
                    lb = lb[-constrained_index],
                    y2value = y2value,
                    epsilon = epsilon,
                    y1value = y1value)
  }
  if(!cfp & !mcfp) {
    qvalue = nloptr(x0 = startvalue,
                    eval_f = oc_obj,
                    eval_g_ineq = oc_ineq,
                    opts = list("algorithm" = algorithm,"xtol_rel" = xtol_rel,
                                "maxeval" = maxeval, "print_level" = printlevel, "ranseed" = "67230",
                                "tol_constraints_ineq" = ineq_tol),
                    parindex = parindex,
                    lower = lower,
                    npar = npar,
                    spar = spar,
                    zvalues = zvalues,
                    constrained_index = constrained_index,
                    constrained_value = constrained_value,
                    testsets = testsets,
                    ub = ub[-constrained_index],
                    lb = lb[-constrained_index],
                    epsilon = epsilon)
  }
   

  return(qvalue)
}


# oc_obj returns the value of an objective function, namely the value of the parindex term in the parameter vector with constrained elements included.
# this for the case in which we calculate bounds on parameters. (not implemented in the example)
oc_obj = function(x, parindex, lower, constrained_index, constrained_value, spar, zvalues, testsets, npar, epsilon) {
 
  thispar = rep(NA,npar)
  thispar[constrained_index] = constrained_value
  thispar[setdiff(1:npar,constrained_index)] = x
  thispar[parindex]*lower
  #browser()
}

# returns the value of an objective function, namely the counterfactual probability of Y=y1value at a particular value of Y2 (and Z if included).
# used if cfp is TRUE
oc_obj_cfp = function(x, parindex, lower, constrained_index, constrained_value, spar, zvalues, testsets, npar, y2value, epsilon, y1value) {
  
  thispar = rep(NA,npar)
  thispar[constrained_index] = constrained_value
  thispar[setdiff(1:npar,constrained_index)] = x
  
  if (y1value == 0) {
    pnorm(thispar[1] - thispar[3]*y2value)*lower
  } else {
    if (y1value == 2) {
      (1 - pnorm(thispar[1] + exp(thispar[2]) - thispar[3]*y2value))*lower
    } else {
      if (y1value == 1) {
        (pnorm(thispar[1] + exp(thispar[2]) - thispar[3]*y2value) - pnorm(thispar[1] - thispar[3]*y2value))*lower
      }
    }
  }
}

# returns the value of an objective function, namely the marginal effect on the counterfactual probability of Y=y1value at a particular value of Y2 (and Z if included).
# used if mcfp is TRUE
oc_obj_mcfp = function(x, parindex, lower, constrained_index, constrained_value, spar, zvalues, testsets, npar, y2value, epsilon, y1value) {
  
  thispar = rep(NA,npar)
  thispar[constrained_index] = constrained_value
  thispar[setdiff(1:npar,constrained_index)] = x
  if (y1value == 0) {
    -thispar[3]*dnorm(thispar[1] - thispar[3]*y2value)*lower
  } else {
    if (y1value == 2) {
      thispar[3]*dnorm(thispar[1] + exp(thispar[2]) - thispar[3]*y2value)*lower
    } else {
      if (y1value == 1) {
        (-thispar[3]*dnorm(thispar[1] + exp(thispar[2]) - thispar[3]*y2value) + thispar[3]*dnorm(thispar[1] - thispar[3]*y2value))*lower
      }
    }
  }
}

# returns a list of values of m_j(theta) j running from 1 to J where the inequality constraints are m_j(theta) <= 0 for all j
oc_ineq = function(x, parindex, lower, constrained_index, constrained_value, spar, zvalues, testsets, npar, y2value = NULL,epsilon = 0,y1value) {

  nz = nrow(zvalues)
  nsets = nrow(testsets)
  ineq_values = rep(0,nz*nsets)
  count = 0
  
  thispar = rep(NA,npar)
  thispar[constrained_index] = constrained_value
  thispar[setdiff(1:npar,constrained_index)] = x
  
  mpar = list(alpha = thispar[3], gamma1 = thispar[1], gamma2 = thispar[1] + exp(thispar[2]), beta = thispar[-(1:3)])
  
  for (isets in 1:nsets) {
    ss = testsets[isets,1]
    tt = testsets[isets,2]
    for (iz in 1:nz) {
      count = count + 1
      z = zvalues[iz,]
      ineq_values[count] = conprob(mpar,spar,z,ss,tt)$prob - (pnorm(tt) - pnorm(ss)) + epsilon
    }
  }
  return(ineq_values)
}


# ineq_values delivers values of inequalities at a solution to an nlopt run. G(S) - Con(s) should all be nonnegative for values in the identified set.
ineq_values = function(solution,spar,zvalues,nt,constrained_index = c(4,5),constrained_value = c(0,0)) {
  
  npar = length(solution) + length(constrained_value)
  thispar = rep(0,npar)
  thispar[constrained_index] = constrained_value
  thispar[setdiff(1:npar,constrained_index)] = solution
  
  mpar = list(alpha = thispar[3], gamma1 = thispar[1], gamma2 = thispar[1] + exp(thispar[2]), beta = thispar[-(1:3)])
  
  testsets = getsets(nt)
  nz = nrow(zvalues)
  Jvalue = nz*nrow(testsets)
  ineq = rep(0,Jvalue)
  count = 0
  for (it in (1:nrow(testsets))) {
    ss = testsets[it,1]
    tt = testsets[it,2]
    for (iz in 1:nz) {
      count = count + 1
      cprob = conprob(mpar,spar,zvalues[iz,],ss,tt)$prob
      gprob = pnorm(tt)-pnorm(ss)
      disc = gprob-cprob
      ineq[count] = disc
    }
  }
  return(range(ineq))
}

# gettruth calculates values of CCP and MCCP in a structure for Y1 equal to 0, 1 and 2
# spar parameters of the structure delivering probabilities
# y2 - the value of Y2 employed
gettruth = function(spar,y2) {
  c1 = spar$c1; c2 = spar$c2; a = spar$a
  cfp0 = pnorm(c1 - a*y2)
  cfp2 = 1-pnorm(c2 - a*y2)
  cfp1 = pnorm(c2 - a*y2) - pnorm(c1 - a*y2)
  
  mfp0 = -a*dnorm(c1 - a*y2)
  mfp2 = a*dnorm(c2 - a*y2)
  mfp1 = -a*dnorm(c2 - a*y2) + a*dnorm(c1 - a*y2)
  return(list(cfp = c(cfp0,cfp1,cfp2),mcfp = c(mfp0,mfp1,mfp2)))
}


######################### ------------ Commands below deliver results in the numerical example

################################
################# asymmetric c's

# Y1 = 0
mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 1)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 CFP
holdasy_eps1l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 0)

ineq_values(holdasy_eps1l$solution,spar,zvalues,51) #

population_values_05 = gettruth(spar,0.5) # values of CCP and MCCP 

holdasy_eps1u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps1u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 CFP
holdasy_eps2l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps2l$solution,spar,zvalues,51) #

population_values_00 = gettruth(spar,0) # values of CCP and MCCP 

holdasy_eps2u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps2u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 CFP
holdasy_eps3l = ochoice_projection(startvalue = c(-0.5,0.262,0.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps3l$solution,spar,zvalues,51) #

population_values_negative_05 = gettruth(spar,-0.5) # values of CCP and MCCP 

holdasy_eps3u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps3u$solution,spar,zvalues,51) #


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 CFP
holdasy_eps4l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps4l$solution,spar,zvalues,51) #


holdasy_eps4u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps4u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps5l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps5l$solution,spar,zvalues,51) #

holdasy_eps5u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps5u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 CFP
holdasy_eps6l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps6l$solution,spar,zvalues,51) #

holdasy_eps6u = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, cfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps6u$solution,spar,zvalues,51) #


# MCFP


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 1)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 MCFP
holdasy_eps7l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 2000, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps7l$solution,spar,zvalues,51) #

holdasy_eps7u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps7u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 MCFP
holdasy_eps8l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps8l$solution,spar,zvalues,51) #

holdasy_eps8u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps8u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 MCFP
holdasy_eps9l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps9l$solution,spar,zvalues,51) #

holdasy_eps9u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-5)

ineq_values(holdasy_eps9u$solution,spar,zvalues,51) #

# s22 = 0.01

mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 MCFP
holdasy_eps10l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =12000, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps10l$solution,spar,zvalues,51) #

holdasy_eps10u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6)

ineq_values(holdasy_eps10u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps11l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6)

ineq_values(holdasy_eps11l$solution,spar,zvalues,51) #

holdasy_eps11u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6, y1value = 0)

ineq_values(holdasy_eps11u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 0.01 CFP
holdasy_eps12l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-5)

ineq_values(holdasy_eps12l$solution,spar,zvalues,51) #

holdasy_eps12u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6)

ineq_values(holdasy_eps12u$solution,spar,zvalues,51) #

####################### Y1 = 2
mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 1)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 CFP
holdasy_eps13l = ochoice_projection(startvalue = c(-0.5,0.262,1.5), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps13l$solution,spar,zvalues,51) #

holdasy_eps13u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps13u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 CFP
holdasy_eps14l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps14l$solution,spar,zvalues,51) #

holdasy_eps14u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps14u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 CFP
holdasy_eps15l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps15l$solution,spar,zvalues,51) #

holdasy_eps15u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps15u$solution,spar,zvalues,51) #


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 CFP
holdasy_eps16l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps16l$solution,spar,zvalues,51) #

holdasy_eps16u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps16u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps17l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps17l$solution,spar,zvalues,51) #

holdasy_eps17u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps17u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 0.01 CFP
holdasy_eps18l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps18l$solution,spar,zvalues,51) #

holdasy_eps18u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps18u$solution,spar,zvalues,51) #


# MCFP


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 1)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 MCFP
holdasy_eps19l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-5,y1value = 2)

ineq_values(holdasy_eps19l$solution,spar,zvalues,51) #

holdasy_eps19u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps19u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 MCFP
holdasy_eps20l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps20l$solution,spar,zvalues,51) #

holdasy_eps20u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps20u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 MCFP
holdasy_eps21l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps21l$solution,spar,zvalues,51) #

holdasy_eps21u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps21u$solution,spar,zvalues,51) #

# s22 = 0.01

mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 MCFP
holdasy_eps22l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =12000, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps22l$solution,spar,zvalues,51) #

holdasy_eps22u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps22u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps23l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps23l$solution,spar,zvalues,51) #

holdasy_eps23u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps23u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 0.01 CFP
holdasy_eps24l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps24l$solution,spar,zvalues,51) #

holdasy_eps24u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 2)

ineq_values(holdasy_eps24u$solution,spar,zvalues,51) #



####################### Y1 = 1
mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.05)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 CFP
holdasy_eps25l = ochoice_projection(startvalue = c(0,0,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 500, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps25l$solution,spar,zvalues,51) #

holdasy_eps25u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps25u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 CFP
holdasy_eps26l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 5e-5,y1value = 1)

ineq_values(holdasy_eps26l$solution,spar,zvalues,51) #

holdasy_eps26u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps26u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 CFP
holdasy_eps27l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 50, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps27l$solution,spar,zvalues,51) #

holdasy_eps27u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps27u$solution,spar,zvalues,51) #


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 CFP
holdasy_eps28l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-5,y1value = 1)

ineq_values(holdasy_eps28l$solution,spar,zvalues,51) #

holdasy_eps28u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.5,epsilon = 1e-5,y1value = 1)

ineq_values(holdasy_eps28u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps29l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps29l$solution,spar,zvalues,51) #

holdasy_eps29u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = 0.0,epsilon = 1e-5,y1value = 1)

ineq_values(holdasy_eps29u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 0.01 CFP
holdasy_eps30l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps30l$solution,spar,zvalues,51) #

holdasy_eps30u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, cfp = TRUE, y2value = -0.5,epsilon = 1e-5,y1value = 1)

ineq_values(holdasy_eps30u$solution,spar,zvalues,51) #


# MCFP


mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 1)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 1 MCFP
holdasy_eps31l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps31l$solution,spar,zvalues,51) #

holdasy_eps31u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps31u$solution,spar,zvalues,51) #


# y2 = 0.0 s22 = 1 MCFP
holdasy_eps32l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps32l$solution,spar,zvalues,51) #

holdasy_eps32u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps32u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 1 MCFP ?????????????
holdasy_eps33l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps33l$solution,spar,zvalues,51) #

holdasy_eps33u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps33u$solution,spar,zvalues,51) #

# s22 = 0.01

mpar = list(beta = c(0,0),alpha = 1, gamma1 = -0.5, gamma2 = 0.8)
spar = list(a = 1, b = c(0,0), c1 = -0.5, c2 = 0.8, d0 = 0, d1 = c(0,0.5), s12 = 0, s22 = 0.01)
zvalues = rbind(c(-1,1),c(1,-1))

# y2 = 0.5 s22 = 0.01 MCFP
holdasy_eps34l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =12000, mcfp = TRUE, y2value = 0.5,epsilon = 1e-5,y1value = 1)

ineq_values(holdasy_eps34l$solution,spar,zvalues,51) #

holdasy_eps34u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps34u$solution,spar,zvalues,51) #

# y2 = 0.0 s22 = 0.01 CFP
holdasy_eps35l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps35l$solution,spar,zvalues,51) #

holdasy_eps35u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval =1200, mcfp = TRUE, y2value = 0.0,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps35u$solution,spar,zvalues,51) #

# y2 = -0.5 s22 = 0.01 CFP
holdasy_eps36l = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = 1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps36l$solution,spar,zvalues,51) #

holdasy_eps36u = ochoice_projection(startvalue = c(-0.5,0.262,1), parindex = 3, lower = -1, constrained_index = c(4,5), constrained_value = c(0,0),
                                 spar = spar, zvalues = zvalues, nt = 51, printlevel = 0,maxeval = 1200, mcfp = TRUE, y2value = -0.5,epsilon = 1e-6,y1value = 1)

ineq_values(holdasy_eps36u$solution,spar,zvalues,51) #

print("Column 3 True Population Counterfactual Probabilities.")
print(round(matrix(t(cbind(population_values_05$cfp,population_values_00$cfp,population_values_negative_05$cfp)),nrow=9,ncol=1),2))

print("Column 4 Bounds on Counterfactual Probabilities. s22 = 1.\n")
cat("[",round(holdasy_eps1l$objective,2),", ",-round(holdasy_eps1u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps2l$objective,2),", ",-round(holdasy_eps2u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps3l$objective,2),", ",-round(holdasy_eps3u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps25l$objective,2),", ",-round(holdasy_eps25u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps26l$objective,2),", ",-round(holdasy_eps26u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps27l$objective,2),", ",-round(holdasy_eps27u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps13l$objective,2),", ",-round(holdasy_eps13u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps14l$objective,2),", ",-round(holdasy_eps14u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps15l$objective,2),", ",-round(holdasy_eps15u$objective,2),"]\n",sep="")

print("Column 5 Bounds on Counterfactual Probabilities.  s22 = 0.01.")
cat("[",round(holdasy_eps4l$objective,2),", ",-round(holdasy_eps4u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps5l$objective,2),", ",-round(holdasy_eps5u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps6l$objective,2),", ",-round(holdasy_eps6u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps28l$objective,2),", ",-round(holdasy_eps28u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps29l$objective,2),", ",-round(holdasy_eps29u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps30l$objective,2),", ",-round(holdasy_eps30u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps16l$objective,2),", ",-round(holdasy_eps16u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps17l$objective,2),", ",-round(holdasy_eps17u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps18l$objective,2),", ",-round(holdasy_eps18u$objective,2),"]\n",sep="")

print("Column 6 True Population Marginal Effects.")
print(round(matrix(t(cbind(population_values_05$mcfp,population_values_00$mcfp,population_values_negative_05$mcfp)),nrow=9,ncol=1),2))

print("Column 7 Bounds on Marginal Effects.  s22 = 1.")
cat("[",round(holdasy_eps7l$objective,2),", ",-round(holdasy_eps7u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps8l$objective,2),", ",-round(holdasy_eps8u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps9l$objective,2),", ",-round(holdasy_eps9u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps31l$objective,2),", ",-round(holdasy_eps31u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps32l$objective,2),", ",-round(holdasy_eps32u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps33l$objective,2),", ",-round(holdasy_eps33u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps19l$objective,2),", ",-round(holdasy_eps19u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps20l$objective,2),", ",-round(holdasy_eps20u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps21l$objective,2),", ",-round(holdasy_eps21u$objective,2),"]\n",sep="")

print("Column 8 Bounds on Marginal Effects.  s22 = 0.01.")
cat("[",round(holdasy_eps10l$objective,2),", ",-round(holdasy_eps10u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps11l$objective,2),", ",-round(holdasy_eps11u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps12l$objective,2),", ",-round(holdasy_eps12u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps34l$objective,2),", ",-round(holdasy_eps34u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps35l$objective,2),", ",-round(holdasy_eps35u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps36l$objective,2),", ",-round(holdasy_eps36u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps22l$objective,2),", ",-round(holdasy_eps22u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps23l$objective,2),", ",-round(holdasy_eps23u$objective,2),"]\n",sep="")
cat("[",round(holdasy_eps24l$objective,2),", ",-round(holdasy_eps24u$objective,2),"]\n",sep="")



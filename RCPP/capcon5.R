#######################################################################################################################################
# FUNCTION sncv computes the self-normal critical value from Chernozhukov, Chetverikov, and Kato (2019), as in Belloni, Bugni, and Chernozhukov (2018)
# INPUT ARGUMENTS
#   p: number of moments.
#   alpha: size of test.
#   n: sample size.
# RETURNS
#  critical value for a level alpha test.
#######################################################################################################################################*/
sncv <- function(p, alpha, n) {
  return(qnorm(1 - alpha/p)/sqrt(1-((qnorm(1 - alpha/p))^2)/n))
}

#######################################################################################################################################
# FUNCTION make_t_pairs computes a collection of (t1,t2) pairs with t1 < t2 for testing inequalities of the form (3.2) in CRS23 
# INPUT ARGUMENTS
#   n: sample size.
# RETURNS
#   Matrix of dimension n*(n+1)/2 by 2 in which each row represents an interval [t1,t2] for which to use inequality (3.2) in the discrepancy function.
#   Each row represents an intervals from one value of t to another, in which t values correspond to evently placed quantiles of the standard normal distribution.
#######################################################################################################################################*/
make_t_pairs <- function(n) {
  t1_t2_pairs <- matrix(0,nrow = n*(n+1)/2,ncol=2)
  count <- 0
  for (i in 0:(n-1)) {
    for (j in (i+1):n) {
      count <- count+1
      t1_t2_pairs[count,] <- c(i/n,j/n)
    }
  }
  tpairs = qnorm(t1_t2_pairs)
  return(tpairs[-n,]) # Remove (-Inf,Inf)
}

#######################################################################################################################################
# FUNCTION draw_random_start_theta draws a vector theta from the uniform distribution over the parameter space defined by the arguments.
#   Implemented by drawing uniformly over the multivariate rectange defined by upper and lower bounds on each components, and then
#   discarding any draws in violation of any additional constraints, and re-drawing until a vector is drawn not in violation.
# INPUT ARGUMENTS
#   theta_lb: lower bounds on parameter vector.
#   theta_ub: upper bounds on parameter vector.
#   pfun: function that defined target parameter, i.e. the target parameter is pfun(theta).
#   param_lb: lower bound constraint on the target parameter.
#   param_ub: upper bound constraint on the target parameter.
#   arglist: a list of additional arguments needed for pfun.
# RETURNS
#   Matrix of dimension n*(n+1)/2 by 2 in which each row represents an interval [t1,t2] for which to use inequality (3.2) in the discrepancy function.
#   Each row represents an intervals from one value of t to another, in which t values correspond to evently placed quantiles of the standard normal distribution.
#######################################################################################################################################*/
draw_random_start_theta <- function(theta_lb, theta_ub, pfun, param_lb, param_ub, arglist) {
  k <- length(theta_lb)
  theta_start <- runif(k,theta_lb,theta_ub)
  if (is.null(pfun)) {
    theta_start <- rep(0,k)
    for (i in 1:k) {
      theta_start[i] <- runif(1,theta_lb[i],theta_ub[i])
    }
  } else {
    while( (pfun(theta_start, arglist) < param_lb) || (pfun(theta_start, arglist) > param_ub) ) {
      # If the constraints on the profile function are in violation, then pick a new random starting value
      theta_start <- runif(k,theta_lb,theta_ub)
    }
  }
  return(theta_start)
}

#######################################################################################################################################
# FUNCTION LoadMTOData loads the Moving to Opportunity data.
# INPUT ARGUMENTS
#   datapath: the path and filename for the data file.
# RETURNS
#   A list containing the different data variables used in the analysis.
#######################################################################################################################################*/
LoadMTOData <- function(datapath) {
  mto.data <- read.csv(paste(datapath,"/34860-0001-Data-REST.csv", sep=""))
  Y <- (as.numeric(factor(mto.data$happy_scale123_ad)) - 1)
  W <- matrix(0,nrow = length(Y),ncol=2)
  W[,1] <- mto.data$f_c9010t_perpov_dw_z
  W[,2] <- mto.data$f_c9010t_pminorty_dw_z
  X <- matrix(0,nrow = length(Y),ncol=4)
  X[,1] <- mto.data$x_f_site_balt
  X[,2] <- mto.data$x_f_site_bos
  X[,3] <- mto.data$x_f_site_chi
  X[,4] <- mto.data$x_f_site_la
  Z <- mto.data$ra_group
  return(list(Y=Y,W=W,X=X,Z=Z))
}

##########################################################################################################
# FUNCTION ConMarginalEffect computes the derivative of Pr(Y(w0)=y0|X=x0) with respect to w0 (theta[3])
#   for a given value of y0, x0, and w0.
# INPUT ARGUMENTS:
#   theta: parameter vector (c1,c2,beta,gamma) at which to evaluate the conditional marginal effect
#   moreargs$...
#     y0: fixed value of y at which to compute the conditional marginal effect.
#     w0: fixed value of w at which to compute the conditional marginal effect.
#     x0: fixed value of x at which to compute the conditional marginal effect.
# RETURNS:
#   the conditional marginal effect.
################################################################################################################
ConMarginalEffect <- function(theta,moreargs) {
  w0 <- moreargs$w0
  wbeta_plus_xgamma <- sum(w0 * theta[seq(3,2+length(w0),1)]) + sum(moreargs$x0 * theta[seq(3+length(w0),length(theta),1)])
  cutoffs <- c(-Inf,theta[1],theta[2],Inf)
  return( ( dnorm(cutoffs[moreargs$y0+1]-wbeta_plus_xgamma) - dnorm(cutoffs[moreargs$y0+2]-wbeta_plus_xgamma) ) * theta[3] )
}

##########################################################################################################
# FUNCTION ConRP computes the counterfactual response probability Pr(Y(w0)=y0|X=x0)
#   for a given value of y0, x0, and w0.
# INPUT ARGUMENTS:
#   theta: parameter vector (c1,c2,beta,gamma) at which to evaluate the conditional marginal effect
#   moreargs$...
#     y0: fixed value of y at which to compute the counterfactual response probability.
#     w0: fixed value of w at which to compute the counterfactual response probability.
#     x0: fixed value of x at which to compute the counterfactual response probability.
# RETURNS:
#   the counterfactual response probability.
################################################################################################################
ConRP <- function(theta,moreargs) {
  w0 <- moreargs$w0
  wbeta_plus_xgamma <- sum(w0 * theta[seq(3,2+length(w0),1)]) + sum(moreargs$x0 * theta[seq(3+length(w0),length(theta),1)])
  cutoffs <- c(-Inf,theta[1],theta[2],Inf)
  return( pnorm(cutoffs[moreargs$y0+2]-wbeta_plus_xgamma) - pnorm(cutoffs[moreargs$y0+1]-wbeta_plus_xgamma) )
}

##########################################################################################################
# FUNCTION Simulate_Data creates simulated data from the complete triangular models described in the MTO draft at the specified parameter values.
# INPUT ARGUMENTS:
#   thresholds: J*1 vector of latent variables determining thresholds for 1,...,J 
#   beta: Coefficients on endogenous variables W.
#   gamma: Coefficients on exogenous variables X.
#   delta_x: Coefficients on X in the "first stage". The first column are coefficients to use when Z==0, second column when Z==1, third column when Z==2.
#   R: Vector of covariances of U with V.
#   SigmaV: Variance of (V1,V2).
#   XZprobs_file: Name of file from which to load the support and probability mass function for (x,z)
# RETURNS:
#   A list comprising [Y,W,X,Z] of data generated using the given arguments.
Simulate_Data <- function(thresholds,beta,gamma,delta_x,R,SigmaV,XZprobs_file) {
  
  load(XZprobs_file) 
  XZindices <- sample(x=seq(1,15), size=3273, replace=TRUE, prob=XZsupport_means)
  XZmatrix <- XZsupport[XZindices,]

  dimV <- 1
  dimX <- length(gamma)
  J <- length(thresholds) # J = maximum value of Y
  n <- dim(XZmatrix)[1]
  
  SigmaUV <- matrix(nrow = dimV+1, ncol = dimV+1)
  SigmaUV[1,1] <- 1 # Variance of U
  SigmaUV[2:(dimV+1),2:(dimV+1)] <- SigmaV
  SigmaUV[2:(dimV+1),1] <- R
  SigmaUV[1,2:(dimV+1)] <- R
  
  set.seed(82789384)
  UVdraws <- rmvnorm(n,matrix(0L,nrow=1,ncol=(dimV+1)),SigmaUV)
  U <- UVdraws[,1]
  V <- UVdraws[,2:(dimV+1)]
  X <- XZmatrix[1:n,1:dimX,drop=FALSE] # drop=FALSE to suppress matrix reduction to a vector when dimX==1
  X_sat <- cbind(X,1-rowSums(X))
  Z <- XZmatrix[1:n,dim(XZmatrix)[2]] # last column of XZmatrix
  W <- (X_sat %*% delta_x[,1]) * (Z==1) + (X_sat %*% delta_x[,2]) * (Z==2) + (X_sat %*% delta_x[,3]) * (Z==3) + V
  Ystar <- W %*% beta + X%*%gamma + U
  
  Y <- matrix(0,nrow=n,ncol=1) + 1*((Ystar>thresholds[1])&(Ystar<=thresholds[2])) + 2*(Ystar>thresholds[2])
  Z <- XZmatrix[1:n,(dimX+1):(dim(XZmatrix)[2]),drop=FALSE]; # drop=FALSE to suppress matrix reduction to a vector when dimX==1
  
  return(list(Y=Y,W=W,X=X,Z=Z))
}

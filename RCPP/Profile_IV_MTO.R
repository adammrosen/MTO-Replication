# Profile_IV_MTO
# Replication Code for Chesher, Rosen, and Siddique (2022): "Estimating Endogenous Effects on Ordinal Outcomes."
# In the single equation Gaussian error IV model search over the specified region
# for the specified target parameter and conduct estimation and inference.
# This file can be used to compute set estimates and confidence sets for Conditional Marginal Effects
# and Counterfactual Response Probabilities reported in Tables 3--5.

# The MTO ICPSR data is required for replication. The data is not public, but may be acquired from ICPSR.
# See https://www.icpsr.umich.edu/web/ICPSR/studies/34860

##### TO PRODUCE RESULTS REPORTED IN TABLE 4 #####################################################################
# SET include_NBH_MIN = F (line 99).
# SET target_parameter = CME (line 92) .
# SOURCE file 4 times with all configurations of isPosRegion = {T,F} (lines 95,96) and y0 = {0,2} (line 107).
# This will produce 4 output files, each with name beginning "PROFILE_SEIV_".
# Set estimates and confidence sets for each configuration can be found at the bottom of each of these files.
##################################################################################################################

##### TO PRODUCE RESULTS REPORTED IN TABLE 5 #####################################################################
# SET include_NBH_MIN = T (line 99).
# SET target_parameter = CME (line 92).
# SOURCE file 4 times with all configurations of isPosRegion = {T,F} (lines 95,96) and y0 = {0,2} (line 107).
# This will produce 4 output files, each with name beginning "PROFILE_SEIV_".
# Set estimates and confidence sets for each configuration can be found at the bottom of each of these files.
##################################################################################################################

##### TO PRODUCE RESULTS REPORTED IN TABLE 6 #####################################################################
# SET include_NBH_MIN = F (line 99).
# SET isPosRegion = F on (line 96).
# SET target_parameter = CRP (line 93).
# SOURCE file 6 times with all configurations of isMedianPov = {T,F} (lines 103, 104) and y0 = {0,1,2} (line 107).
# This will produce 6 output files, each with name beginning "PROFILE_SEIV_".
# Set estimates and confidence sets for each configuration can be found at the bottom of each of these files.
##################################################################################################################

##!!!!! Without access to the MTO data, this code may be run on data simulated from a complete triangular model using
##!!!!! the function Simulate_Data is capcon5.R.
##!!!!! To do this set use_MTO_Data = F (line 64) and use_simulated_Data = T (line 65).

######## Load libraries as needed. ##########################
if (!require("mvtnorm")) install.packages("mvtnorm")
library(mvtnorm)
if (!require("ordinal")) install.packages("ordinal")
library(ordinal)
if (!require("stargazer")) install.packages("stargazer")
library(stargazer)
if (!require("nloptr")) install.packages("nloptr")
library(nloptr)
if (!require("prodlim")) install.packages("prodlim")
library(prodlim)
if (!require("Matrix")) install.packages("Matrix")
library(Matrix)
if (!require("coop")) install.packages("coop")
library(coop)
if (!require("inline")) install.packages("inline")
library(inline)
if (!require("Rcpp")) install.packages("Rcpp")
library(Rcpp)
if (!require("RcppArmadillo")) install.packages("RcppArmadillo")
library(RcppArmadillo)
rm(list=ls())
start_time <- proc.time()

use_MTO_Data <- F # T for real MTO data, F for simulated or pseudo data
use_simulated_Data <- T
DRB_singleton_theta_set <- T # Set to T to only use theta_hat sample estimate when evaluating bootstrap rather than the set of all theta satisfying the profile constraint and achieving the optimal profile discrepancy
R_boot <- 99
echo <- T # Set to T for more detailed output during execution
search_grid_gap <- 0.005 # set grid point distance for subvector/function search

if ((use_MTO_Data + use_simulated_Data) != 1) stop("Exactly one of use_MTO_Data and use_simulated_Data must be true.")

if (use_MTO_Data) {
    sourcepath <- "D:/data/12942365/ICPSR_34860/DS0001/R"
} else {
    sourcepath <- "/Users/amr331/dropbox/tex documents/MTO/mtoempirics/Replication Files/RCPP" # "/Users/amr331/Git/MTO/RCPP" # For laptop
}

######## Set working directory and instantiate functions from source files ######################
setwd(sourcepath)
source("capcon5.R")
sourceCpp("DiscrepancyFunctions.cpp")
#################################################################################################

cv_type <- "DRB" # for Discard Resampling Bootstrap critical values
computeDRB <- T
alpha <- 0.05 # 1-alpha is nominal coverage level of CIs
c_omega_bar <- 0.1

#target_parameter <- "parameter_x" # for structural parameter component x
target_parameter <- "CME" # for conditional marginal effect
#target_parameter <- "CRP" # for conditional response probability

isPosRegion <- T     # Set to T to search for positive value of the target parameter
#isPosRegion <- F     # Set to F to search for negative value of the target parameter

#include_NBH_MIN <- T # Set to T to include neighborhood minority as an additional included endogenous variable
include_NBH_MIN <- F # Set to F if not.

print_val <- 0       # changing to 1, 2, or 3 will result in nloptr giving some output at each iteration

isMedianPov <- T # Set to T to use median poverty level for New York
#isMedianPov <- F # Set to F to use median poverty minus 1 sd for New York

# Set values at which to evaluate the conditional marginal effect or counterfactual response probability
y0 <- 0 # value of y0 for which to compute the marginal effect
x0 <- c(0,0,0,0) # value of x0 at which to compute the marginal effect

sinkfile <- gsub(":","-",paste("PROFILE_SEIV_", target_parameter, "_", Sys.time(), sep =""))
sink(file = sinkfile, split = T) # Write output both to file and console.

cat("Running the SEIV specification with neighborhood minority ",
    ifelse(include_NBH_MIN,"included","excluded"),
    ", y0 = ", y0, " and x0 = (", sep="")
cat(x0,sep = ", ")
cat(") with target parameter ", target_parameter, ".\n", sep = "")
cat("Profile grid search step size is ", search_grid_gap, "\n")

################ LOAD DATA ########################################
if (use_MTO_Data) {
  cat("** USING MTO DATA ***\n");
  YWXZ_obs <- LoadMTOData(sourcepath)
}
if (use_simulated_Data) {
  cat("** USING SIMULATED DATA ***\n");
  gamma_baltimore <- 0.2778
  gamma_boston <- 0.1435
  gamma_chicago <- 0.2989
  gamma_la <- 0.0893 
  gamma <- c(gamma_baltimore,gamma_boston,gamma_chicago,gamma_la)
  delta_x2 <- c(-1.0912,-1.2798,-0.3068,-0.8787,-0.8052) # coefficients in equation for W that multiply city dummy for Z=2, experimental voucher
  delta_x1 <- c(-1.0427,-1.0880,-0.1905,-0.8139,-0.3742) # coefficients in equation for W that multiply city dummy for Z=1, Section 8 voucher 
  delta_x0 <- c(-0.5220,-0.7145,0.2299,0.1584,0.1371) # coefficients in equation for W that multiply city dummy for Z=0, no voucher 
  delta_x <- cbind(delta_x0,delta_x1,delta_x2)
  YWXZ_obs <- Simulate_Data(c(-0.4750,0.9119),-0.1487,gamma,delta_x,0.0668,1,paste(sourcepath,"/XZprobs.RData",sep=""))
}

XZmatrix <- cbind(YWXZ_obs$X,YWXZ_obs$Z)
Xsupport <- unique(YWXZ_obs$X)
Zsupport <- unique(YWXZ_obs$Z)
XZsupport <- unique(XZmatrix)
nSupportX <- dim(Xsupport)[1]
nSupportZ <- length(Zsupport)
nSupportXZ <- dim(XZsupport)[1]
n_obs <- dim(XZmatrix)[1]
XZsupport_indices <- matrix(row.match(data.frame(XZmatrix),data.frame(XZsupport)),ncol=1)
Xsupport_indices <- matrix(row.match(data.frame(YWXZ_obs$X),data.frame(Xsupport)),ncol=1)
    
# HERE SET w0 to be the median levels of W for New York
W_pov <- as.matrix(YWXZ_obs$W[which(Xsupport_indices==5),])[,1]
se_pov <- sqrt( sum((W_pov - mean(W_pov)) ^2)/(n_obs-1)  )
if (include_NBH_MIN) {
  if (isMedianPov) {
    w0 <- c(median(W_pov),median(YWXZ_obs$W[which(Xsupport_indices==5),][,2]))
  } else {
    w0 <- c(median(W_pov) - se_pov,median(YWXZ_obs$W[which(Xsupport_indices==5),][,2]))
  }
} else { 
  # Only retain the first column of W
  YWXZ_obs$W <- as.matrix(YWXZ_obs$W[,1])
  if (isMedianPov) {
    w0 <- median(W_pov)
  } else {
    w0 <- median(W_pov) - se_pov
  }
}
cat("Setting counterfactual w0 value to (", w0, ").\n")
##################################################################################################    
    
####### Draw values for multiplier bootstrap #################################
rseed <- 4436766
set.seed(rseed)
xiBootDraws <- matrix(rnorm(n_obs*R_boot), nrow = R_boot, ncol = n_obs )
##############################################################################

######## Set nlopt options and random seed. ########################################
rseed <- 987698764
set.seed(rseed)
eq_tol <- 1.0e-04 # Tolerance for equality and inequality constraint violations

local_opts = list("algorithm" = "NLOPT_LN_COBYLA", "ftol_rel" = 1.0e-08, "ftol_abs" = 1.0e-08, "xtol_rel" = 1.0e-08,
                  "print_level" = print_val, tol_constraints_eq = eq_tol) # options for profile optimization routine.

global_opts <- list("algorithm" = "NLOPT_LN_COBYLA", "ftol_rel" = 1.0e-04, "ftol_abs" = 1.0e-04, "xtol_rel" = 1.0e-04,
                    "print_level" = print_val, "maxeval" = 1000, "stopval" = -Inf)

profile_opts = list("algorithm" = "NLOPT_LN_COBYLA", "local_opts" = local_opts, "xtol_rel" = 1.0e-08, tol_constraints_eq = eq_tol,
                    tol_constraints_ineq = c(eq_tol,eq_tol), "ftol_rel" = 1.0e-08, "ftol_abs" = 1.0e-08, "maxeval"=1000, "print_level" = print_val)
#####################################################################################

################ Create (t1,t2) pairs independently of simulations and theta #######################
tpoint_count <- 24
tpairs <- make_t_pairs(tpoint_count)
tpairs_by_n_obs <- array(tpairs,c(dim(tpairs),n_obs))
tpairs_by_n_obs <- aperm(tpairs_by_n_obs,c(3,2,1))
## n * 2 * J array
## tpairs_by_n_obs(i,k,j) is identical for all i giving the k-th element of the j-th tpair
## where J denotes the total number of tpairs, which is J = ((tpoint_count * tpoint_count+1)/2) - 1.
###################################################################################################

################ Create indicators for binning of conditioning variables X,Z. #############################
# First we deal with aspects of conditioning on xz, which do not depend on theta
# Generate unique numerical indices for each point of support on XZ
xz_indicators = matrix(0,nrow=n_obs,ncol=nSupportXZ)
for (i in 1:nSupportXZ) {
  xz_indicators[,i] = (XZsupport_indices == i)
}
# xz_indicators[i,j] indicates which of the j support point of XZ observation i falls into
# It has dimension n by the number of support points of XZ.
###########################################################################################################

################ Compute the RHS of the inequalities (3.2) - (3.4), which do not depend on theta
## The values are replicated across the number of tpairs to enable computing conditional
## probabilities without looping over XZ values, but instead through matrix comparisons.
tpairs_by_xz <- array(tpairs,c(dim(tpairs),nSupportXZ))
inequalityRHS_by_xz = t(pnorm(tpairs_by_xz[,2,]) - pnorm(tpairs_by_xz[,1,]))

###### COUNT NUMBER OF MOMENTS FOR SEIV GAUSSIAN AND SEIV GAUSSIAN MTO MODELS #########
J <- dim(tpairs_by_n_obs)[3] # number of tpairs
K <- dim(XZsupport)[1] # number of xz support points
momentcount <- J*K

################ Initialize values for profiled searches ##################################################
if (include_NBH_MIN) {
    pov_min_theta_TIV  <- c(-0.3704,0.9402,-0.2876,0.3417,0.3823,0.4967,0.2363,0.0688) # Triangular IV
    pov_min_theta_exog <- c(-0.4986,0.8914,-0.0791,-0.0213,0.2996,0.1589,0.2828,0.0676) # Ordered Probit
    v1 <- c(-0.02151177, 0.9022524, -0.812197, 0.5777164, 0.4645985, 0.4960496, 0.1024556, 0.2559215)
    v2 <- c(-0.3298451, 1.093974, -0.1597805, 0.2255501, 0.3784085, 0.7262483, 0.2446697, 0.3084423)
    v3 <- c(-0.8664188, 0.0850136, 0.3677914, -0.7703699, -0.4570376, -1.142988, -0.2163358, -0.2399566)
    v4 <- c(-0.06742408, 1.747111, 0.005923601, 0.9491897, 2.168834, 2.490019, 0.1951724, 0.3845195)
    v5 <- c(-0.0822951, 1.737675, 0.001447531, 0.9269304, 2.169322, 2.490015, 0.5314177, 0.394223)
    v6 <- c(-0.1394531, 1.611066, 0.1192675, 0.9171108, 2.073564, 2.490675, -0.02165676, 0.4310656)
    v7 <- c(-0.06742408, 1.747111, 0.005923601, 0.9491897, 2.168834, 2.490019, 0.1951724, 0.3845195)
    nDeterministicStartVecs <- 9
    deterministic_theta_start_vals = matrix(c(pov_min_theta_TIV,pov_min_theta_exog,v1,v2,v3,v4,v5,v6,v7),ncol=nDeterministicStartVecs)
} else {
    pov_theta_TIV <- c(-0.475, 0.9119, -0.1487, 0.2778, 0.1435, 0.2989, 0.0893) # Triangular IV
    pov_theta_exog <- c(-0.4935,0.8962,-0.0886,0.3068,0.1812,0.2800,0.0957) # Ordered Probit
    pov_theta_other <- c(-0.3648486, 1.041138, -0.3171777, -0.5532114, -0.2453283, 0.5207691, -0.4420153)
    # some minimizing vectors found in first draft using fewer moment inequalities
    v1 <- c(-1.035897, 2.446239, 0.9858357, 1.055011, 1.361132, 0.5089113, 0.2727134)
    v2 <- c(-0.4860253, 0.7054605, -0.2710354, 0.05271256, -0.190834, 0.1895751, -0.1207801)
    v3 <- c(-0.2601916, 0.6538879, -0.6571943, 0.0207307, -0.09569314, 0.3015204, 0.030403)
    v4 <- c(-0.3929452, 0.6495299, -0.6861087, -0.2063017, -0.1565981, 0.478618, -0.1424547)
    v5 <- c(-0.5559967, 0.7100101, -0.1282928, 0.03671564, -0.008380695, 0.1448245, 0.07808737)
    v6 <- c(-0.5882615, 1.943888, 0.6646634, 2.401796, 1.448784, 1.613387, 0.3040592)
    nDeterministicStartVecs <- 9
    deterministic_theta_start_vals <- matrix(c(pov_theta_TIV,pov_theta_exog,pov_theta_other,v1,v2,v3,v4,v5,v6),ncol=nDeterministicStartVecs)
}

nParams <- dim(deterministic_theta_start_vals)[1]
epsilon <- 0
if (echo) cat("Value of epsilon for discrepancy set to ", epsilon, ".\n", sep="")
universal_bound <- 2.5 # Set a universal upper and lower bound for all parameters
theta_lb <- rep(-universal_bound,nParams)
theta_ub <- rep(universal_bound,nParams)
if ((target_parameter == "CRP") | (target_parameter == "parameter_x")) {
  if (isPosRegion) { # restrict beta to be nonnegative
    theta_lb[3] <- 0
  } else { # restrict beta to be nonpositive
    theta_ub[3] <- 0
  }
}

BBC_SN_critical_value = sncv(momentcount, 0.05, n_obs)
BBC_SN_median_value = sncv(momentcount, 0.50, n_obs)

cat("BBC self-normalized critical value for ", momentcount, " moments is ", BBC_SN_critical_value, ".\n", sep="")
cat("BBC median self-normalized critical value for ", momentcount, " moments is ", BBC_SN_median_value, ".\n", sep="")

################# populate starting parameter values #################################
arglist <- list(xz_indicators=xz_indicators, tpairs_by_n=tpairs_by_n_obs,
                inequalityRHS_by_xz=inequalityRHS_by_xz, cut=0,
                theta_lb=theta_lb, theta_ub=theta_ub, alpha=alpha,
                opts=global_opts, profile_opts = profile_opts,
                Xsupport = Xsupport, Xsupport_indices = Xsupport_indices,
                Zsupport = Zsupport, XZsupport = XZsupport, XZsupport_indices = XZsupport_indices,
                nSupportXZ = nSupportXZ, nSupportX = nSupportX, nSupportZ = nSupportZ,
                tpairs=tpairs, last_profile_objective = NULL, last_profile_full_solution = NULL,
                studentize_moments = T, return_value_only = F, w0 = w0, x0 = x0, y0 = y0, # values needed for the parameter function, e.g. a marginal effect
                xiBootDraws = xiBootDraws, R = R_boot, # draws of xi for bootstrap
                discfun = GaussianDiscrepancyCpp, cv_type = cv_type, runDR = computeDRB,
                target_parameter = target_parameter, y = matrix(YWXZ_obs$Y,ncol=1), x=YWXZ_obs$X, w=YWXZ_obs$W, echo = echo,
                c_omega_bar = c_omega_bar, epsilon = epsilon, snCV = BBC_SN_critical_value, DRB_singleton_theta_set = DRB_singleton_theta_set)

if (target_parameter == "CRP") {
  component_range <- c(0.0,1.0)
} else {
  if (target_parameter == "CME")
    if (isPosRegion) {
      component_range <- c(0.0,0.8)
    } else {
      component_range <- c(-0.8,0.0)
    }
  else {
  if (target_parameter == "parameter_x") {
    if (isPosRegion) {
      component_range <- c(0.0,2.5)
    } else {
      component_range <- c(-2.5,0.0)
    }
  }
  }
}
arglist$target_lower_bound <- component_range[1]
arglist$target_upper_bound <- component_range[2]

#### SPECIFY THE FUNCTION TO PROFILE AS THE PARAM_FUN ELEMENT IN THE LIST ##################
# si is short for "search item" which specifies parameters by which to conduct a search.
if (target_parameter == "parameter_x") {
  # BETA
  si <- list( param_index = 3, component_range = component_range, solution = c(-Inf,Inf), discrepancy = c(Inf,Inf),
              full_param_vector = matrix(0,nrow=nParams,ncol=2), class = "component_search_item", target_parameter = target_parameter,
              param_fun = NULL, param_name = "parameter component 3", cv_type = cv_type)
} else {
if (target_parameter == "CME") {
  # CONDITIONAL MARGINAL EFFECT
  si <- list( param_index = 3, component_range = component_range, solution = c(-Inf,Inf), discrepancy = c(Inf,Inf),
              full_param_vector = matrix(0,nrow=nParams,ncol=2), class = "component_search_item", target_parameter = target_parameter,
              param_fun = ConMarginalEffect, param_name = "conditional marginal effect", cv_type = cv_type)
} else
if (target_parameter == "CRP") {
  # CONDITIONAL RESPONSE PROBABILITY
  si <- list( param_index = 3, component_range = c(0,1), solution = c(-Inf,Inf), discrepancy = c(Inf,Inf),
              full_param_vector = matrix(0,nrow=nParams,ncol=2), class = "component_search_item", target_parameter = target_parameter,
              param_fun = ConRP, param_name = "conditional choice probability", cv_type = cv_type)
}
}

if (is.null(si$param_fun)) {
  cat("Initiating search for parameter ", si$param_index, " over the range [", si$component_range[1], ",", si$component_range[2], "].\n", sep="")
  arglist$theta_lb[si$param_index] <- si$component_range[1]
  arglist$theta_ub[si$param_index] <- si$component_range[2]
} else {
  cat("Initiating search for ", si$param_name, " over the range [", si$component_range[1], ",", si$component_range[2], "].\n", sep="")
}
arglist$param_index <- si$param_index
arglist$param_fun <- si$param_fun
arglist$studentize_moments <- T
  
#-------------- CONDUCT START SEARCH OVER THE RANGE ------------------------
nRandomStartVecs <- 100 # Number of starting values for initial search.
theta_start_vals <- matrix(0,ncol = (nDeterministicStartVecs + nRandomStartVecs), nrow=nParams)
theta_start_vals[,1:nDeterministicStartVecs] <- deterministic_theta_start_vals
theta_param_vals <- rep(0,nDeterministicStartVecs + nRandomStartVecs) # Save the value of the function of theta or subvector for the given value of theta
if (echo) cat("Starting Values are as follows: ")
for (j in 1:nDeterministicStartVecs) { # replace any deterministic starting values in violation of constraints with values that fall within the required range
  if (is.null(si$param_fun)) {
    theta_param_vals[j] <- theta_start_vals[si$param_index,j]
    if ( (theta_param_vals[j] < si$component_range[1]) || (theta_param_vals[j] > si$component_range[2] )) {
      theta_start_vals[si$param_index,j] <- runif(1,si$component_range[1],si$component_range[2])
      theta_param_vals[j] <- theta_start_vals[si$param_index,j]
    }
  } else {
    theta_param_vals[j] <- si$param_fun(theta_start_vals[,j], moreargs = arglist)
    if ( (theta_param_vals[j] < si$component_range[1]) || 
         (theta_param_vals[j] > si$component_range[2])
         || any(theta_start_vals[,j] < arglist$theta_lb) || any(theta_start_vals[,j] > arglist$theta_ub)) {
      theta_start_vals[,j] <- draw_random_start_theta(theta_lb, theta_ub, si$param_fun, si$component_range[1], si$component_range[2], arglist)
      theta_param_vals[j] <- si$param_fun(theta_start_vals[,j], moreargs = arglist)
    }
  }
}

for (j in 1:nRandomStartVecs) {
  theta_start_vals[,nDeterministicStartVecs+j] <- draw_random_start_theta(theta_lb, theta_ub, si$param_fun, si$component_range[1], si$component_range[2], arglist)
  theta_param_vals[j] <- si$param_fun(theta_start_vals[,j], moreargs = arglist)
}
arglist$disc_vals <- rep(Inf,dim(theta_start_vals)[2])
arglist$disc_solns <- 0 * theta_start_vals

# -- Call MinDiscrepancyCpp to minimize the discrepancy function using nlopt with each different starting value --
search.time.nloptC <- system.time(retlist <- MinDiscrepancyCpp(theta_start_vals, theta_lb, theta_ub, arglist))
if (echo) cat("MinDiscrepancyCpp took ", search.time.nloptC[3], " seconds.\n", sep="")
disc_vals <- retlist$disc_vals
disc_solns <- retlist$disc_solns
param_vals <- retlist$param_vals

###################################### SEARCH FOR GLOBAL OPTIMUM WITHIN SEARCH RANGE CONCLUDED ###############################
min_value <- max( min(disc_vals), epsilon )

# FIND THE LOWEST AND HIGHEST VALUES OF THE *** SUBVECTOR or FUNCTION OF PARAMETERS ****

inset_indices <- which(disc_vals <= min_value )
if (echo) cat ("There were ", length(inset_indices), " searches that found a point in the analog set.\n")

min_gval <- min(param_vals[inset_indices])
min_gval_idx <- inset_indices[which.min(param_vals[inset_indices])]
min_gval_theta <- disc_solns[,min_gval_idx]

max_gval <- max(param_vals[inset_indices])
max_gval_idx <- inset_indices[which.max(param_vals[inset_indices])]
max_gval_theta <- disc_solns[,max_gval_idx]

# HERE PRINT HIGHEST AND LOWEST VALUE FOUND FROM THE INITIAL SEARCH

if (echo) {
  cat("Lowest value of function of interest found in the set was ", min_gval, " with parameter vector = (", sep="")
  cat(min_gval_theta, sep = ", ")
  cat(") and discrepancy ", disc_vals[min_gval_idx], "\n", sep ="")
  cat("Highest value of function of interest found in the set was ", max_gval, " with parameter vector = (", sep="")
  cat(max_gval_theta, sep = ", ")
  cat(") and discrepancy ", disc_vals[max_gval_idx], "\n", sep ="")
}

##### MINIMIZE CONSTRAINED DISCREPANCY FOR TARGET PARAMETER OVER A GRID OF POSSIBLE VALUES ####
  
arglist$profile_opts$stopval <- ifelse(min_value > epsilon,BBC_SN_median_value,epsilon)
arglist$profile_opts$local_opts$stopval <- arglist$profile_opts$stopval
  
search_grid_1 <- c(seq(max(si$component_range[1],min_gval),si$component_range[1],-search_grid_gap),si$component_range[1])
search_grid_2 <- c(seq(min(si$component_range[2],max_gval),si$component_range[2],search_grid_gap),si$component_range[2])  

g1 <- length(search_grid_1)
g2 <- length(search_grid_2)

if (g2 > 0) {
  search_grid <- c(rev(search_grid_1),search_grid_2[seq(2,g2)])
} else {
  search_grid <-rev(search_grid_1)
}

# Grid for search constructed.

# Matrices to save output at each point in the grid.
n_gpts <- length(search_grid)
si$grid_vals <- rep(Inf,n_gpts)
si$grid_vecs <- matrix(0, nrow = nParams, ncol = n_gpts)
si$grid_CVs <- rep(Inf,n_gpts)

#### Call ProfileDiscreapncyOnGridCpp to compute the profiled discrepancy function with target parameter fixed at each grid point in turn ####

search.time.nloptC <- system.time(ProfileDiscrepancyOnGridCpp(search_grid, g1, cbind(min_gval_theta, max_gval_theta), theta_lb, theta_ub, arglist$profile_opts$stopval, arglist, si))
if (echo) cat("ProfileDiscrepancyOnGridCpp took ", search.time.nloptC[3], " seconds.\n", sep="")

grid_vals <- si$grid_vals
grid_vecs <-si$grid_vecs
grid_CVs <- si$grid_CVs

theta_star <- grid_vecs[,which.min(grid_vals)] # full parameter vector at which minimum of profiled discrepancy achieved

cat("The discrepancy-minimizing vector of theta was found to be = (")
cat(theta_star, sep = ", ")
cat("),\n")
cat("which achieved a discrepancy of ", min(grid_vals), "\n", sep = "")
if (is.null(si$param_fun)) {
  cat("The corresponding value of the ", si$param_name, " was ", theta_star[si$param_index], "\n", sep = "")
} else {
  cat("The corresponding value of the ", si$param_name, " was ",
      si$param_fun(theta = theta_star, moreargs = arglist), "\n", sep = "")
}
cat("Grid search was conducted ranging from ", min(search_grid), " to ", max(search_grid), " in increments of ", search_grid_gap, "\n\n", sep="")

arglist$runDR <- F # Don't bootstrap when searching for the analog set estimate

cat("\n")
cat(" *************************** ANALOG SET *************************** ")
cat("\n") 

if(min(grid_vals) > epsilon) {
  is_empty_analog_set <- T
  cat("No points were found in the analog set.\n")
} else {
  is_empty_analog_set <- F
  min_analog_index <- min(which(grid_vals < epsilon))
  max_analog_index <- max(which(grid_vals < epsilon))
  cat("The minimum value of ", si$param_name, " on the grid in the analog set is ", search_grid[min_analog_index], "\n", sep="")
  cat("The maximum value of ", si$param_name, " on the grid in the analog set is ", search_grid[max_analog_index], "\n", sep="")
  # REFINE LOWER BOUND ESTIMATE IF NECESSARY
  if (component_range[1] < search_grid[min_analog_index]) {
    cat("---- Starting refined search for analog set lower bound ----- \n")
    LB_start <- max(component_range[1],search_grid[min_analog_index] - search_grid_gap)
    refined_analog_LB <- RefineBoundCpp(LB_start, search_grid[min_analog_index], T, 0.0005, epsilon, grid_vecs[,min_analog_index], arglist, si) 
  } else {
    cat("---- analog set lower bound coincides with lower bound of parameter component range ----- \n")
    refined_analog_LB <- list(bound = component_range[1], discrepancy = grid_vals[min_analog_index], theta = grid_vecs[,min_analog_index]) 
  }
  # REFINE UPPER BOUND ESTIMATE IF NECESSARY
  if (component_range[2] > search_grid[max_analog_index]) {
    cat("---- Starting refined search for analog set upper bound ----- \n")
    UB_start <- min(component_range[2],search_grid[max_analog_index] + search_grid_gap)
    refined_analog_UB <- RefineBoundCpp(UB_start, search_grid[max_analog_index], F, 0.0005, epsilon, grid_vecs[,max_analog_index], arglist, si)
  } else {
    cat("---- analog set upper bound coincides with upper bound of parameter component range ----- \n")
    refined_analog_UB <- list(bound = component_range[2], discrepancy = grid_vals[max_analog_index], theta = grid_vecs[,max_analog_index]) 
  }
  Analog_LB <- refined_analog_LB$bound
  Analog_UB <- refined_analog_UB$bound
  cat("Refined analog set interval: [", Analog_LB, ",", Analog_UB, "].\n", sep="")
  if (echo) cat("The value of theta achieving the LOWER bound discrepancy of ", refined_analog_LB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_analog_LB$theta, sep = ", ")
  if (echo) cat("),\n")
  if (echo) cat("The value of theta achieving the UPPER bound discrepancy of ", refined_analog_UB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_analog_UB$theta, sep = ", ")
  if (echo) cat("),\n")
}
if (echo) cat("\n")

cat("\n")
cat(" *************************** MEDIAN-CORRECTED SET ESTIMATE *************************** ")
cat("\n") 

if(min(grid_vals) > BBC_SN_median_value) {
  is_empty_SN_half <- T
  cat("No points were found in the median-corrected set estimate.\n")
} else {
  is_empty_SN_half <- F
  min_median_index <- min(which(grid_vals < BBC_SN_median_value))
  max_median_index <- max(which(grid_vals < BBC_SN_median_value))
  cat("The minimum value of ", si$param_name, " on the median-corrected interval estimate is ", search_grid[min_median_index], "\n", sep="")
  cat("The maximum value of ", si$param_name, " on the median-corrected interval estimate is ", search_grid[max_median_index], "\n", sep="")
  # REFINE LOWER BOUND ESTIMATE IF NECESSARY
  if (component_range[1] < search_grid[min_median_index]) {
    if (echo) cat("---- Starting refined search for median corrected lower bound ----- \n")
    LB_start <- max(component_range[1],search_grid[min_median_index] - search_grid_gap)
    refined_med_LB <- RefineBoundCpp(LB_start, search_grid[min_median_index], T, 0.0005, BBC_SN_median_value, grid_vecs[,min_median_index], arglist, si)
  }
  else {
    if (echo) cat("---- Median corrected lower bound coincides with lower bound of parameter component range ----- \n")
    refined_med_LB <- list(bound = component_range[1], discrepancy = grid_vals[min_median_index], theta = grid_vecs[,min_median_index]) 
  }
  # REFINE UPPER BOUND ESTIMATE IF NECESSARY
  if (component_range[2] > search_grid[max_median_index]) {
    if (echo) cat("---- Starting refined search for median corrected upper bound ----- \n")
    UB_start <- min(component_range[2],search_grid[max_median_index] + search_grid_gap)
    refined_med_UB <- RefineBoundCpp(UB_start, search_grid[max_median_index], F, 0.0005, BBC_SN_median_value, grid_vecs[,max_median_index], arglist, si)
  }
  else {
    if (echo) cat("---- Median corrected upper bound coincides with upper bound of parameter component range ----- \n")
    refined_med_UB <- list(bound = component_range[2], discrepancy = grid_vals[max_median_index], theta = grid_vecs[,max_median_index])
  }
  SN_50_LB <- refined_med_LB$bound
  SN_50_UB <- refined_med_UB$bound
  cat("Median corrected interval: [", SN_50_LB, ",", SN_50_UB, "].\n", sep="")
  if (echo) cat("The value of theta achieving the LOWER bound discrepancy of ", refined_med_LB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_med_LB$theta, sep = ", ")
  if (echo) cat("),\n")
  if (echo) cat("The value of theta achieving the UPPER bound discrepancy of ", refined_med_UB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_med_UB$theta, sep = ", ")
  if (echo) cat("),\n")
}
if (echo) cat("\n")

cat("\n")
cat(" *************************** SELF-NORMALIZED CONFIDENCE SET *************************** ")
cat("\n")  

if(min(grid_vals) > BBC_SN_critical_value) {
  is_empty_SN_95 <- T
  cat("No points were found in the SN confidence set.\n")
} else {
  is_empty_SN_95 <- F
  min95_index <- min(which(grid_vals < BBC_SN_critical_value))
  max95_index <- max(which(grid_vals < BBC_SN_critical_value))
  if (echo) cat("The minimum value of ", si$param_name, " on the grid in the 95% SN CI is ", search_grid[min95_index], "\n", sep="")
  if (echo) cat("The maximum value of ", si$param_name, " on the grid in the 95% SN CI is ", search_grid[max95_index], "\n", sep="")
  # REFINE LOWER BOUND ESTIMATE IF NECESSARY
  if (component_range[1] < search_grid[min95_index]) {
    if (echo) cat("---- Starting refined search for SN CI lower bound ----- \n")
    LB_start <- max(component_range[1],search_grid[min95_index] - search_grid_gap)
    refined_95_LB <- RefineBoundCpp(LB_start, search_grid[min95_index], T, 0.0005, BBC_SN_critical_value, grid_vecs[,min95_index], arglist, si)
  } else {
    if (echo) cat("---- SN CI lower bound coincides with lower bound of parameter component range ----- \n")
    refined_95_LB <- list(bound = component_range[1], discrepancy = grid_vals[min95_index], theta = grid_vecs[,min95_index])
  }
  # REFINE UPPER BOUND ESTIMATE IF NECESSARY
  if (component_range[2] > search_grid[max95_index]) {
    if (echo) cat("---- Starting refined search for SN CI upper bound ----- \n")
    UB_start <- min(component_range[2],search_grid[max95_index] + search_grid_gap)
    refined_95_UB <- RefineBoundCpp(UB_start, search_grid[max95_index], F, 0.0005, BBC_SN_critical_value, grid_vecs[,max95_index], arglist, si)
  } else {
    if (echo) cat("---- SN CI upper bound coincides with upper bound of parameter component range ----- \n")
    refined_95_UB <- list(bound = component_range[2], discrepancy = grid_vals[max95_index], theta = grid_vecs[,max95_index])
  }
  SN_95_LB <- refined_95_LB$bound
  SN_95_UB <- refined_95_UB$bound
  cat("Refined 95% SN interval: [", SN_95_LB, ",", SN_95_UB, "].\n", sep="")
  if (echo) cat("The value of theta achieving the LOWER bound discrepancy of ", refined_95_LB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_95_LB$theta, sep = ", ")
  if (echo) cat("),\n")
  if (echo) cat("The value of theta achieving the UPPER bound discrepancy of ", refined_95_UB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_95_UB$theta, sep = ", ")
  if (echo) cat("),\n")
}

##### COMPUTE DISCARD BOOTSTRAP RESAMPLING CONFIDENCE SET ###########################################
if (computeDRB) {
  cat("\n")
  cat(" *************************** DISCARD RESAMPLING BOOTSTRAP *************************** ")
  cat("\n")  
  arglist$runDR <- T
  if(min(grid_vals - grid_CVs) > 0) {
    is_empty_DRB_95 <- T
    cat("No points were found in the DR bootstrap confidence set.\n")
  } else {
    is_empty_DRB_95 <- F
    min95_index <- min(which(grid_vals < grid_CVs))
    max95_index <- max(which(grid_vals < grid_CVs))
    if (echo) cat("The minimum value of ", si$param_name, " on the grid in the 95% DRB CI is ", search_grid[min95_index], "\n", sep="")
    if (echo) cat("The maximum value of ", si$param_name, " on the grid in the 95% DRB CI is ", search_grid[max95_index], "\n", sep="")
    # REFINE LOWER BOUND ESTIMATE IF NECESSARY
    if (component_range[1] < search_grid[min95_index]) {
      if (echo) cat("---- Starting refined search for CI lower bound ----- \n")
      LB_start <- max(component_range[1],search_grid[min95_index] - search_grid_gap)
      refined_95_LB <- RefineBoundCpp(LB_start, search_grid[min95_index], T, 0.0005, BBC_SN_critical_value, grid_vecs[,min95_index], arglist, si)
  } else {
    if (echo) cat("---- 95% DRB CI lower bound coincides with lower bound of parameter component range ----- \n")
    refined_95_LB <- list(bound = component_range[1], discrepancy = grid_vals[min95_index], theta = grid_vecs[,min95_index])
  }
  # REFINE UPPER BOUND ESTIMATE IF NECESSARY
  if (component_range[2] > search_grid[max95_index]) {
    if (echo) cat("---- Starting refined search for CI upper bound ----- \n")
    UB_start <- min(component_range[2],search_grid[max95_index] + search_grid_gap)
    refined_95_UB <- RefineBoundCpp(UB_start, search_grid[max95_index], F, 0.0005, BBC_SN_critical_value, grid_vecs[,max95_index], arglist, si)
  } else {
    if (echo) cat("---- 95% DRB CI upper bound coincides with upper bound of parameter component range ----- \n")
    refined_95_UB <- list(bound = component_range[2], discrepancy = grid_vals[max95_index], theta = grid_vecs[,max95_index])
  }
  DRB_95_LB <- refined_95_LB$bound
  DRB_95_UB <- refined_95_UB$bound
  cat("Refined 95% DRB interval: [", DRB_95_LB, ",", DRB_95_UB, "].\n", sep="")
  if (echo) cat("The value of theta achieving the LOWER bound discrepancy of ", refined_95_LB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_95_LB$theta, sep = ", ")
  if (echo) cat("),\n")
  if (echo) cat("The value of theta achieving the UPPER bound discrepancy of ", refined_95_UB$discrepancy, " was:\n (", sep = "")
  if (echo) cat(refined_95_UB$theta, sep = ", ")
  if (echo) cat("),\n")
  }
}

cat("\n\n")
cat(" *************************** End Results *************************** ")
cat("\n\n")
cat(" --- Search over ", ifelse(isPosRegion, "nonnegative", "nonpositive"), " region of parameter space --- \n", sep="")

if (is_empty_analog_set) {
  cat("Analog set is empty.\n")
} else {
  cat("Analog estimate interval: [", Analog_LB, ",", Analog_UB, "].\n", sep="")
}

if (is_empty_SN_half) {
  cat("SN half-median-unbiased interval is empty.\n")
} else {
  cat("SN half-median-unbiased interval: [", SN_50_LB, ",", SN_50_UB, "].\n", sep="")
}

if (is_empty_SN_95) {
  cat("SN 95% confidence interval is empty.\n")
} else {
  cat("SN 95% confidence interval: [", SN_95_LB, ",", SN_95_UB, "].\n", sep="")
}

if (computeDRB) {
  if (is_empty_DRB_95) {
    cat("DRB 95% confidence interval is empty.\n")
  } else {
    cat("DRB 95% confidence interval: [", DRB_95_LB, ",", DRB_95_UB, "].\n", sep="")
  }
}

DisplayFunctionEvaluations()

cat("\n ***** Execution took ", (proc.time() - start_time)["elapsed"], " seconds ***** ", sep="")

sink()






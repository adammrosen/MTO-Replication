#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include <string> 
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]

#include <time.h>

using namespace Rcpp;
using namespace arma;

static int fcount_DRBoot = 0;
static int fcount_evaluations = 0;
static int fcount_fulldisc_evaluations = 0;
static int fcount_profiled_evaluations = 0;
static int fcount_refining_evaluations = 0;
static int fcount_OmegaBarObjective = 0;

/*##################################################################################################################################
# FUNCTION capcon_thresholds_Cpp computes the threshold function define in Section 3.1 of the MTO paper.
# INPUT ARGUMENTS
#   y: n*1 vector of ordered outcomes
#   w: n*k_w vector of realizations of endogenous variables
#   x: n*k_x vector of realizations of included exogenous variables
#   beta_vec: parameter vector beta to use when computing thresholds
#   gamma_vec: parameter vector gamma to use when computing thresholds.
#   cutoffs: parameter vector (c_1,...,c_J)
# RETURNS
#   n x 2 matrix of upper and lower endpoints of interval in which each U_i must fall if
#   population parameters correspond to beta_vec, gamma_vec, and cutoffs.
#######################################################################################################################################*/
mat capcon_thresholds_Cpp(uvec y, mat x, mat w, vec beta_vec, vec gamma_vec, vec cutoffs) {
  vec extended_cutoffs(cutoffs.n_elem + 2);
  extended_cutoffs(0) = -datum::inf;
  extended_cutoffs(extended_cutoffs.n_elem-1) = datum::inf;
  for (uword i=1; i<(extended_cutoffs.n_elem-1); i++) {
    extended_cutoffs(i) = cutoffs(i-1);
  }
  mat cc_thresholds(y.n_elem,2);
  cc_thresholds.col(0) = extended_cutoffs.elem(y) - x * gamma_vec - w * beta_vec;
  cc_thresholds.col(1) = extended_cutoffs.elem(y+1) - x * gamma_vec - w * beta_vec;
  return cc_thresholds;
}

/*##################################################################################################################################
# FUNCTION capcon_indicators_Cpp forms indicator functions for the events in equations (2.7) - (2.8)
# INPUT ARGUMENTS
#   cc_thresholds: Thresholds for each U_i computed from capcon_thresholds_Cpp 
#   tpairs: K x 2 matrix of (s,t) pairs in the inequalities described in Section 3.2.
# RETURNS
#   cc_indicators: an n x K matrix where cc_indicators(i,j) indicates whether observation i satisfies t1 <= c(y,x,w) <= c(y+1,x,w) <= t2
#   for the j-th t1-t2 pair.
#######################################################################################################################################*/
umat capcon_indicators_Cpp(mat cc_thresholds, mat tpairs) {
  umat cc_indicators(cc_thresholds.n_rows, tpairs.n_rows, fill::zeros);
  for (uword i=0; i<cc_thresholds.n_rows; i++) {
    for (uword j=0; j<tpairs.n_rows; j++) {
      cc_indicators(i,j) = ( (tpairs(j,0) <= cc_thresholds(i,0)) && (cc_thresholds(i,1) <= tpairs(j,1)) );
    }
  }
  return cc_indicators;
}

/*##################################################################################################################################
# FUNCTION GaussianMFunValsCpp computes values of each moment function at each observation i. The sample mean of each of these functions
#   across observation i is the sample analog of the population mean of that function, which is required to be less than
#   or equal to zero.  That is, the mean of these moment functions comprise the moment inequalities that restrict theta.
# INPUT ARGUMENTS
#   theta: The parameter vector at which to evaluate the moment functions.
#   moreargs: A list containing additional variables, data, and other quantities pre-computed in R.
# RETURNS
#   moment_indicators: an n x J matrix where moment_indicators(i,j) indicates whether observation i satisfies t1 <= c(y,x,w) <= c(y+1,x,w) <= t2
#   for the j-th t1-t2 pair and value of (x,z).  This is the value of omega(x_i,z_i,t,theta)
#######################################################################################################################################*/
mat GaussianMFunValsCpp(vec theta, List moreargs) {
  
  /* **** UNPACK data *********/
  uvec y = as<uvec>(moreargs["y"]);
  mat w = as<mat>(moreargs["w"]);
  mat x = as<mat>(moreargs["x"]);
  /* ***************************/
  
  /* **** UNPACK moreargs **********************************************/
  mat xz_indicators = as<mat>(moreargs["xz_indicators"]);
  mat tpairs = as<mat>(moreargs["tpairs"]);
  cube tpairs_by_n = as<cube>(moreargs["tpairs_by_n"]);
  mat inequalityRHS_by_xz = as<mat>(moreargs["inequalityRHS_by_xz"]);
  vec theta_lb = as<vec>(moreargs["theta_lb"]);
  vec theta_ub = as<vec>(moreargs["theta_ub"]);
  /* ******************************************************************/
  
  /* *********** INITIALIZE variables **********************************************/
  int dimbeta = w.n_cols; // Number of columns of w
  int n = w.n_rows; // Sample size
  int J = tpairs_by_n.n_slices; // Number of (t1,t2) pairs
  int Kx = as<int>(moreargs["nSupportX"]);
  int Kz = as<int>(moreargs["nSupportZ"]);
  int K = Kx*Kz;
  
  vec cutoffs(2); 
  for (int i = 0; i < 2; i++) cutoffs(i) = theta(i);
  vec beta_vec(dimbeta); 
  for (int i = 0; i < dimbeta; i++) beta_vec(i) = theta(i+2);
  vec gamma_vec(theta.n_elem-(dimbeta+2)); 
  for (uword i = 0; i < theta.n_elem-(dimbeta+2); i++) gamma_vec(i) = theta(i+dimbeta+2);
  /* ******************************************************************/
  
  /* *********************** Compute containment indicators ***************/
  umat cc_indicators_long = capcon_indicators_Cpp(capcon_thresholds_Cpp(y, x, w, beta_vec, gamma_vec, cutoffs), tpairs);
  umat cc_indicators = umat(n,J);
  vec con_upper_probs = vec(J); // upper probabilities on containments
  for (int i=0; i<J; i++) {
    cc_indicators.col(i) = cc_indicators_long.col(i);
    con_upper_probs(i) = inequalityRHS_by_xz(1,i);
  }
  // Element i,j of cc_indicators indicates whether observation i satisfies t1 <= c(y,x,w) <= c(y+1,x,w) <= t2
  // for the j-th t1-t2 pair.  This is an indicator for the event omega(x_i,z_i,t,theta)
  /* ******************************************************************/
  
  /* *############## COMPUTE INFERENCE PARAMETER ESTIMATES ZETA_HAT ##################################### */
  /* Compute indicators for estimation of inference parameters zeta. */
  mat moment_indicators = mat(n,K*J,fill::zeros);
  int column_index = 0;
  for (int k = 0; k < K; k++) { // k indexes a support point of (X,Z)
    for (int j = 0; j < J; j++) {
      moment_indicators.col(column_index++) = ((cc_indicators.col(j) - (con_upper_probs(j)*ones(n,1))) % xz_indicators.col(k)) ;
    }
  }
  return(moment_indicators);
}

/*##################################################################################################################################
# FUNCTION GaussianDiscrepancyCpp the discrepancy function \hat{Q} defined in (3.7).
#   Uses the function GaussianMFunValsCpp defined above.
# INPUT ARGUMENTS
#   theta: The parameter vector at which to evaluate the moment functions.
#   moreargs: A list containing additional variables, data, and other quantities pre-computed in R.
# RETURNS
#   sqrt(n) times the maximium of the studentized moments, corresponding to \hat{Q} defined in (3.7).
#######################################################################################################################################*/
// [[Rcpp::export]]
double GaussianDiscrepancyCpp(vec theta, List moreargs) {
  mat moment_indicators = GaussianMFunValsCpp(theta,moreargs);
  int n = moment_indicators.n_rows;
  rowvec mu_hat = mean(moment_indicators,0);
  rowvec sigma_hat = stddev(moment_indicators);
  rowvec studentized_moments = mu_hat/sigma_hat;
  return(std::sqrt(n) * max(studentized_moments));
}

/*##############################################################################################################################
# Function GaussianDiscrepancy wraps GaussianDiscrepancyCpp with an argument format amenable to use NLOPT solvers
# INPUT ARGUMENTS
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   value of the discrepancy function computed by GaussianDiscrepancyCpp.
#######################################################################################################################################*/
double GaussianDiscrepancy(unsigned k, const double *theta, double *grad, void *fdata) {
  fcount_evaluations++;
  List moreargs = *((List*)fdata);
  vec theta_vec = vec(k);
  for (unsigned j=0; j<k; j++) theta_vec(j) = theta[j];
  return(GaussianDiscrepancyCpp(theta_vec, moreargs));
}

/*##############################################################################################################################
# FUNCTION BootstrapVstar
# INPUT ARGUMENTS
#   theta: The parameter vector at which to evaluate the moment functions.
#   moreargs: A list containing additional variables, data, and other quantities pre-computed in R.
#   boot_index: The index of the bootstrap repetition being computed.
# RETURNS
#   value of the bootstrap process computed at theta in the bootstrap repetition indicated by boot_index.
#######################################################################################################################################*/
vec BootstrapVstar(vec theta, List moreargs, int boot_index) {
  mat moment_indicators = GaussianMFunValsCpp(theta, moreargs);
  int n = moment_indicators.n_rows;
  int J = moment_indicators.n_cols;
  rowvec mu_hat = mean(moment_indicators);
  rowvec sigma_hat = stddev(moment_indicators);
  mat rep_mu_hat = mat(n,J); // repeat mu_hat in a matrix accordingly
  rep_mu_hat.each_row() = mu_hat; // n x J matrix with each row equal to mu_hat
  mat rep_sigma_hat = mat(n,J); // repeat sigma_hat in a matrix accordingly
  rep_sigma_hat.each_row() = sigma_hat; // n x J matrix with each col equal to sigma_hat
  rowvec xiVec = (as<mat>(moreargs["xiBootDraws"])).row(boot_index);
  return((xiVec * ((moment_indicators - rep_mu_hat)/rep_sigma_hat) / sqrt(n)).t());
}

/*##############################################################################################################################
# FUNCTION DRBootstrapDiscrepancy computes the DR Bootstrap statistic given in Appendix F, equation (F.6) of
#   the online supplement for Chesher, Rosen, and Siddique (2022) following BBC18 equation (4.5).
# INPUT ARGUMENTS
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   value of the DR Bootstrap test statistic #######################################################################################################################################*/
double DRBootstrapDiscrepancy(unsigned k, const double *theta, double *grad, void *fdata) {
  fcount_DRBoot++;
  List moreargs = *((List*)fdata);
  int r = as<int>(moreargs["r"]);
  vec theta_vec = vec(k);
  for (unsigned j=0; j<k; j++) theta_vec(j) = theta[j];
  vec vstar_theta = BootstrapVstar(theta_vec, moreargs, r);
  mat moment_indicators = GaussianMFunValsCpp(theta_vec,moreargs);
  int n = moment_indicators.n_rows;
  int J = moment_indicators.n_cols;
  rowvec mu_hat = mean(moment_indicators);
  rowvec sigma_hat = stddev(moment_indicators);
  colvec tvec = sqrt(n) * (mu_hat/sigma_hat).t();
  double Mn = as<double>(moreargs["Mn"]);
  colvec rep_max_tvec = colvec(J);
  rep_max_tvec.fill(max(tvec) - Mn);
  ucolvec Psi = find(tvec > rep_max_tvec); // This is the \hat\Psi_{\theta} set in BBC
  return(max(vstar_theta.elem(Psi))); 
}

/*##############################################################################################################################
# FUNCTION OmegaBarObjective computes the function of theta the statistic used to compute \bar{w}_n maximized.
#   See the def in the online supplelement to Chesher, Rosen, and Siddique (2022) as defined in (F.4) following the prescription from BBC18 page 14.
# INPUT ARGUMENTS
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   value of the function that is to be maximized with respect to theta.
########################################################################################################################################*/
double OmegaBarObjective(unsigned k, const double *theta, double *grad, void *fdata) {
  fcount_OmegaBarObjective++;
  List moreargs = *((List*)fdata);
  int s = as<int>(moreargs["OmegaBarIteration"]);
  vec theta_vec(k);
  for (unsigned j = 0; j < k; j++) theta_vec(j) = theta[j];
  vec vstar_theta = BootstrapVstar(theta_vec, moreargs, s);
  return(max(abs(vstar_theta)));
}

/*#####################################################################################################################################
# FUNCTION DisplayFunctionEvaluations writes the number of function evaluations from various steps in the computations to the output console
# INPUT ARGUMENTS -- none --
# RETURNS -- none --
#######################################################################################################################################*/
// [[Rcpp::export]]
void DisplayFunctionEvaluations() {
  Rcpp::Rcout << "Full Parameter Discrepancy Evaluations: " << fcount_fulldisc_evaluations << std::endl;
  Rcpp::Rcout << "Profiled Discrepancy Evaluations: " << fcount_profiled_evaluations << std::endl;
  Rcpp::Rcout << "Interval Boundary Profiled Discrepancy Evaluations: " << fcount_refining_evaluations << std::endl;
  Rcpp::Rcout << "DR Bootstrap Evaluations: " << fcount_DRBoot << std::endl;
  Rcpp::Rcout << "Omega Bar Objective Function Evaluations: " << fcount_OmegaBarObjective << std::endl;
}

/*#####################################################################################################################################
# FUNCTION DisplaySearchResult writes the results of minimization of a criterion function with respect to theta to the output console
# INPUT ARGUMENTS 
#   theta_arr: pointer to a double from which the minimizing parameter vector found is reference.
#   min_value: minimum value of the objective function found, achieved by theta_arr
#   k: number of components of theta_arr
# RETURNS -- none --
#######################################################################################################################################*/
void DisplaySearchResult(double *theta_arr, double min_value, int k) {
  Rcpp::Rcout << "Found minimum at theta = (";
  for (int i=0; i<k; i++) {
    if (i==(k-1)) {
      Rcpp::Rcout << theta_arr[i] << ").\n";
    } else {
      Rcpp::Rcout << theta_arr[i] << ", ";
    }
  }
  Rcpp::Rcout << "Minimum value was " << std::setprecision(8) << min_value << "." << std::endl;
}

/*#####################################################################################################################################
# FUNCTION DisplaySearchResult same as above (function overloaded) but writing additional output.
# INPUT ARGUMENTS 
#   p_index: parameter index if the criterion is minimized subject to the p_index component of theta being constrained to set value.
#   g_index: the constraint value for subvector inference.
#   theta_arr: pointer to a double from which the minimizing parameter vector found is reference.
#   evaluation: number of funtion evaluatios conducted in the search for a minimum.
#   min_value: minimum value of the objective function found, achieved by theta_arr
#   k: number of components of theta_arr
# RETURNS -- none --
#######################################################################################################################################*/
void DisplaySearchResult(int p_index, double p_value, double *theta_arr, int evaluations, double min_value, int k) {
  Rcpp::Rcout << "Found minimum at theta = (";
  for (int i=0; i<k; i++) {
    if (i==(k-1)) {
      Rcpp::Rcout << theta_arr[i] << ").\n";
    } else {
      Rcpp::Rcout << theta_arr[i] << ", ";
    }
  }
  Rcpp::Rcout << "Minimum value was " << std::setprecision(8) << min_value << " after " << evaluations << " function evaluations." << std::endl;
}

/*#####################################################################################################################################
# FUNCTION gpdf computes the mean zero unit variance gaussian pdf 
# INPUT ARGUMENTS 
#   x: value at which to compute the pdf
# RETURNS
#   value of the pdf evaluated at x.
#######################################################################################################################################*/
double gpdf(double x) {
  static const double inv_sqrt_2pi = 0.3989422804014327;
  return inv_sqrt_2pi * exp(-0.5 * (x * x) );
}

/*#####################################################################################################################################
# FUNCTION gcdf computes the mean zero unit variance gaussian cdf 
# INPUT ARGUMENTS 
#   x: value at which to compute the cdf
# RETURNS
#   value of the cdf evaluated at x.
#######################################################################################################################################*/
double gcdf(double value) {
   static const double inv_sqrt_2 = 0.707106781187;
   return 0.5 * erfc(-value * inv_sqrt_2);
}

/*#########################################################################################################################################################
# FUNCTION CME_constraint computes the difference between the conditional marginal effect and the value imposed by an equality or inequality constraint.
#   The arguments are specified so as to adhere to NLOPT's requirements for imposing equality or inequality constraints.
#   Values at which the CME are evaluated and the value of the constraint are passed through fdata, which points to a list.
# INPUT ARGUMENTS 
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   value of the conditional marginal effect minus the value imposed by the constraint.
#######################################################################################################################################*/
double CME_constraint(unsigned k, const double *theta, double *grad, void *fdata) {
  List argl = *((List*)fdata);
  double gvalue = argl["gvalue"];
  vec w0 = as<vec>(argl["w0"]);
  vec x0 = as<vec>(argl["x0"]);
  int y0 = as<int>(argl["y0"]);
  vec wx = join_vert(w0,x0);
  double wbeta_plus_xgamma = 0;
  for (unsigned j=2; j < k; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
  vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
  return( (gpdf(cutoffs[y0] - wbeta_plus_xgamma) - gpdf(cutoffs[y0+1] - wbeta_plus_xgamma)) * theta[2] - gvalue );
}

/*#####################################################################################################################################
# FUNCTION DRB_CME_constraint is used by the DR Bootstrap to impose the restriction that parameter vector theta is among those that
#   achieve the sample minimum, when the constraint set in the bootstrap is not simply taken to be the singleton set containing only
#   the sample minimum.
# INPUT ARGUMENTS 
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   The difference between the discrepancy function and the sample minimum.
#######################################################################################################################################*/
double DRB_CME_constraint(unsigned k, const double *theta, double *grad, void *fdata) {
  List argl = *((List*)fdata);
  double Tnvalue = argl["last_profile_objective"];
  return ( (GaussianDiscrepancy(k, theta, NULL, &argl)) - Tnvalue);
}

/*#####################################################################################################################################
# FUNCTION CME_bound_constraints is used to impose upper and lower bound inequality constraints on the parameter space for conditional marginal effects. 
#   The arguments are specified so as to adhere to NLOPT's requirements for imposing multiple inequality constraints.
# INPUT ARGUMENTS 
#   m: number of constraints being imposed, length of result.
#   result: output of the constraint should be stored here. This argument is passed by reference.Tthe function itself has no return value.
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS -- none --
#######################################################################################################################################*/
void CME_bound_constraints(unsigned m, double* result, unsigned k, const double *theta, double *grad, void *fdata) {
  List argl = *((List*)fdata);
  double target_lower_bound = as<double>(argl["target_lower_bound"]);
  double target_upper_bound = as<double>(argl["target_upper_bound"]);
  vec w0 = as<vec>(argl["w0"]);
  vec x0 = as<vec>(argl["x0"]);
  int y0 = as<int>(argl["y0"]);
  vec wx = join_vert(w0,x0);
  double wbeta_plus_xgamma = 0;
  for (unsigned j=2; j < k; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
  vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
  double cme = (gpdf(cutoffs[y0] - wbeta_plus_xgamma) - gpdf(cutoffs[y0+1] - wbeta_plus_xgamma)) * theta[2];
  result[0] = target_lower_bound - cme;
  result[1] = cme - target_upper_bound; 
}

/*#####################################################################################################################################
# FUNCTION CME computes the value of the conditional marginal effect at the values of w0, x0, y0 specified in the list argl.
# INPUT ARGUMENTS 
#   theta: parameter vector at which to compute the conditional marginal effect.
#   argl: A list that includes values of w0, x0, y0 at which to compute the conditional marginal effect.
# RETURNS
#   value of the conditional marginal effect.
#######################################################################################################################################*/
// [[Rcpp::export]]
double CME(vec theta, List argl) {
  vec w0 = as<vec>(argl["w0"]);
  vec x0 = as<vec>(argl["x0"]);
  int y0 = as<int>(argl["y0"]);
  vec wx = join_vert(w0,x0);
  double wbeta_plus_xgamma = 0;
  for (unsigned j=2; j < theta.n_elem; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
  vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
  return((gpdf(cutoffs[y0] - wbeta_plus_xgamma) - gpdf(cutoffs[y0+1] - wbeta_plus_xgamma)) * theta[2]);
}

/*#########################################################################################################################################################
# FUNCTION ConRProb_constraint computes the difference between a counterfactual response probability and the value imposed by an equality or inequality constraint.
#   The arguments are specified so as to adhere to NLOPT's requirements for imposing equality or inequality constraints.
#   Values at which the CRP are evaluated and the value of the constraint are passed through fdata, which points to a list.
# INPUT ARGUMENTS 
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   value of the counterfactual response probability minus the value imposed by the constraint.
#######################################################################################################################################*/
 double ConRProb_constraint(unsigned k, const double *theta, double *grad, void *fdata) {
    List argl = *((List*)fdata);
    double gvalue = argl["gvalue"];
    vec w0 = as<vec>(argl["w0"]);
    vec x0 = as<vec>(argl["x0"]);
    int y0 = as<int>(argl["y0"]);
    vec wx = join_vert(w0,x0);
    double wbeta_plus_xgamma = 0;
    for (unsigned j=2; j < k; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
    vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
    return((gcdf(cutoffs[y0+1] - wbeta_plus_xgamma) - gcdf(cutoffs[y0] - wbeta_plus_xgamma))  - gvalue);
}

/*#####################################################################################################################################
# FUNCTION DRB_ConRProb_constraint is used by the DR Bootstrap to impose the restriction that parameter vector theta is among those that
#   achieve the sample minimum, when the constraint set in the bootstrap is not simply taken to be the singleton set containing only
#   the sample minimum.
# INPUT ARGUMENTS 
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS
#   The difference between the discrepancy function and the sample minimum.
#######################################################################################################################################*/
double DRB_ConRProb_constraint(unsigned k, const double *theta, double *grad, void *fdata) {
  List argl = *((List*)fdata);
  double Tnvalue = argl["last_profile_objective"];
  return ( (GaussianDiscrepancy(k, theta, NULL, &argl)) - Tnvalue);
}

/*#####################################################################################################################################
# FUNCTION CRP_bound_constraints is used to impose upper and lower bound inequality constraints on the parameter space for counterfactual 
#   response probabilities. The arguments are specified so as to adhere to NLOPT's requirements for imposing multiple inequality constraints.
# INPUT ARGUMENTS 
#   m: number of constraints being imposed, length of result.
#   result: output of the constraint should be stored here. This argument is passed by reference.Tthe function itself has no return value.
#   k: length of parameter vector
#   theta: parameter vector (c1,c2,beta,gamma) as an array of doubles (i.e. point to a double at which the array begins).
#   gradient: -- not used -- Required as an argument for use with NLOPT, but set to NULL.
#   fdata: Pointer to void for format required by NLOPT.  Here used to pass moreargs which is then typecast to a list.
# RETURNS -- none --
#######################################################################################################################################*/
void CRP_bound_constraints(unsigned m, double* result, unsigned k, const double *theta, double *grad, void *fdata) {
  List argl = *((List*)fdata);
  double target_lower_bound = as<double>(argl["target_lower_bound"]);
  double target_upper_bound = as<double>(argl["target_upper_bound"]);
  vec w0 = as<vec>(argl["w0"]);
  vec x0 = as<vec>(argl["x0"]);
  int y0 = as<int>(argl["y0"]);
  vec wx = join_vert(w0,x0);
  double wbeta_plus_xgamma = 0;
  for (unsigned j=2; j < k; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
  vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
  double crp = gcdf(cutoffs[y0+1] - wbeta_plus_xgamma) - gcdf(cutoffs[y0] - wbeta_plus_xgamma);
  result[0] = target_lower_bound - crp;
  result[1] = crp - target_upper_bound; 
}

/*#####################################################################################################################################
# FUNCTION ConRProb computes the value of the counterfactual response probability at the values of w0, x0, y0 specified in the list argl.
# INPUT ARGUMENTS 
#   theta: parameter vector at which to compute the conditional marginal effect.
#   argl: A list that includes values of w0, x0, y0 at which to compute the counterfactual response probability.
# RETURNS
#   value of the counterfactual response probability.
#######################################################################################################################################*/
// [[Rcpp::export]]
double ConRProb(vec theta, List argl) {
  vec w0 = as<vec>(argl["w0"]);
  vec x0 = as<vec>(argl["x0"]);
  int y0 = as<int>(argl["y0"]);
  vec wx = join_vert(w0,x0);
  double wbeta_plus_xgamma = 0;
  for (unsigned j=2; j < theta.n_elem; j++) wbeta_plus_xgamma += (wx(j-2) * theta[j]);
  vec cutoffs = {-datum::inf, theta[0], theta[1], datum::inf};
  return(gcdf(cutoffs[y0+1] - wbeta_plus_xgamma) - gcdf(cutoffs[y0] - wbeta_plus_xgamma));
}

/*#####################################################################################################################################
# FUNCTION DisplayArray writes to the output console the values stored in an array of doubles.
# INPUT ARGUMENTS 
#   theta_arr: a pointer to the beginning of the array representing the full parameter vector.
#   k: the number of elements in the array.
# RETURNS
#   -- none --
#######################################################################################################################################*/
void DisplayArray(double* theta_arr, int k) {
  Rcpp::Rcout << "(";
  for (int i=0; i<k; i++) {
    if (i==(k-1)) Rcpp::Rcout << theta_arr[i] << ").\n";
    else Rcpp::Rcout << theta_arr[i] << ", ";
  }
}

/*#####################################################################################################################################
# FUNCTION EstimateOmegaBarN 
# INPUT ARGUMENTS 
#   theta_arr: a pointer to the beginning of the array representing the full parameter vector
#   lb: lower bounds on parameter vector
#   ub: upper bounds on parameter vector
#   k: length of parameter vector
#   moreargs: A list containing additional variables, data, and other quantities pre-computed in R.
# RETURNS
#   A bootstrap-based estimator for \bar{\omega}_n as described on BBC18 page 14.
#######################################################################################################################################*/
double EstimateOmegaBarN(double* theta_arr, double* lb, double* ub, int k, List moreargs) {
  int R = as<mat>(moreargs["xiBootDraws"]).n_rows;
  int n = as<mat>(moreargs["xiBootDraws"]).n_cols;
  int c = as<int>(moreargs["c_omega_bar"]);
  double* theta_boot = new double [k];
  for (int j=0; j < k; j++) theta_boot[j] = theta_arr[j]; 
  List profile_opts = as<List>(moreargs["profile_opts"]);
  std::string target_parameter = as<std::string>(moreargs["target_parameter"]);
  vec OmegaBarStat = vec(R);
  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, k);
  if (target_parameter == "CME") {
    nlopt_add_equality_constraint(opt, CME_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
  } else if (target_parameter == "CRP") {
    nlopt_add_equality_constraint(opt, ConRProb_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
  }
  double maxf;
  ivec exit_codes = ivec(R);
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
  nlopt_set_max_objective(opt, OmegaBarObjective, &moreargs);
  nlopt_set_ftol_rel(opt, as<double>(profile_opts["ftol_rel"]));
  nlopt_set_ftol_abs(opt, as<double>(profile_opts["ftol_abs"]));
  nlopt_set_xtol_rel(opt, as<double>(profile_opts["xtol_rel"]));
  nlopt_set_maxeval(opt, as<int>(profile_opts["maxeval"]));
  nlopt_set_stopval(opt, as<double>(profile_opts["stopval"]));
  for (int s = 0; s < R; s++) {
    moreargs["OmegaBarIteration"] = s;
    int nlopt_return_code = nlopt_optimize(opt, &(theta_boot[0]), &maxf);
    if (nlopt_return_code < 0) {
      Rcpp::Rcout << "nlopt failed! Exit code " << nlopt_return_code << endl;
    }
    OmegaBarStat(s) = maxf; // Compute the statistic whose high level quantile will estimate \bar{\omega}_n following BBC18 page 14
  }
  vec p = vec(1);
  p(0) = 1-pow(n,-c);
  vec q = quantile(OmegaBarStat,p);
  nlopt_destroy(opt);
  delete[] theta_boot;
  return(q[0]);
}

/*#####################################################################################################################################
# FUNCTION BootstrapCV computes the discard resampling critical value follow BBC18.
# INPUT ARGUMENTS 
#   val: Value of the target parameter that is being tested. Which target parameter is passed through moreargs.
#   lb_arr: lower bounds on parameter vector passed as an array (pointer to a double).
#   ub_arr: upper bounds on parameter vector passed as an array (pointer to a double).
#   k: length of parameter vector
#   moreargs: A list containing additional variables, data, and other quantities pre-computed in R and C++.
# RETURNS
#   The discard resampling bootstrap critical value for a level alpha test.
#######################################################################################################################################*/
double BootstrapCV(double val, double* theta_hat, double* lb_arr, double* ub_arr, int k, List moreargs) {
  int R = as<int>(moreargs["R"]);
  bool echo = as<bool>(moreargs["echo"]);
  bool DRB_singleton_theta_set = as<bool>(moreargs["DRB_singleton_theta_set"]);
  double alpha = as<double>(moreargs["alpha"]);
  std::string target_parameter = as<std::string>(moreargs["target_parameter"]);
  int p_index = (as<int>(moreargs["param_index"]) - 1);
  List profile_opts = as<List>(moreargs["profile_opts"]);
  double Mn = log((as<mat>(moreargs["w"])).n_rows) * EstimateOmegaBarN(theta_hat, lb_arr, ub_arr, k, moreargs);
  vec DRstatvals = vec(R);
  DRstatvals.fill(datum::inf);
  double* theta_arr = new double [k];

  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, k); 
  
  if ( !DRB_singleton_theta_set ) { /* -- If only using theta_hat from sample, no need to optimize criterion -- */
    /* -- Set up optimization problem -- Impose constraints -- */
    if (target_parameter == "parameter_x") { // fix parameter component
      lb_arr[p_index] = theta_hat[p_index];
      ub_arr[p_index] = theta_hat[p_index];
    } 
    else if (target_parameter == "CME") {
      if (echo) Rcpp::Rcout << "DR Bootstrap is imposing equality constraints for Tn(theta) = " << as<double>(moreargs["last_profile_objective"]) << " and CME = " << val << endl;
      nlopt_add_equality_constraint(opt, DRB_CME_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
      nlopt_add_equality_constraint(opt, CME_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
    }
    else if (target_parameter == "CRP") {
      if (echo) Rcpp::Rcout << "DR Bootstrap is imposing equality constraints for Tn(theta) = " << as<double>(moreargs["last_profile_objective"]) << " and CRP = " << val << endl;
      nlopt_add_equality_constraint(opt, DRB_ConRProb_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
      nlopt_add_equality_constraint(opt, ConRProb_constraint, &moreargs, as<double>(profile_opts["tol_constraints_eq"]));
    }
    nlopt_set_lower_bounds(opt, lb_arr);
    nlopt_set_upper_bounds(opt, ub_arr);
    nlopt_set_ftol_rel(opt, as<double>(profile_opts["ftol_rel"]));
    nlopt_set_ftol_abs(opt, as<double>(profile_opts["ftol_abs"]));
    nlopt_set_xtol_rel(opt, as<double>(profile_opts["xtol_rel"]));
    nlopt_set_maxeval(opt, as<int>(profile_opts["maxeval"]));
    nlopt_set_stopval(opt, as<double>(profile_opts["stopval"]));
  }
  
  double minf = 0.0;
  int cv_rank = floor(R * alpha);
  vec top_bootstrap_draws = vec(cv_rank + 1);
  ivec exit_codes = ivec(R);
  double p = 1-alpha;
  double DRBcv;
  moreargs["Mn"] = Mn;
  
  for (int r = 0; r < R; r++) {
    moreargs["r"] = r;
      for (int j = 0; j < k; j++) theta_arr[j] = theta_hat[j];
      if (DRB_singleton_theta_set) {
        minf = DRBootstrapDiscrepancy(k, theta_hat, NULL, &moreargs); // Compute Boostrap Discrepancy at theta_hat
      } 
      else {
        nlopt_set_min_objective(opt, DRBootstrapDiscrepancy, &moreargs);
        nlopt_optimize(opt, &(theta_arr[0]), &minf);
      }
      DRstatvals(r) = minf;
      if (r > cv_rank) {
        if (minf > top_bootstrap_draws(cv_rank)) {
          top_bootstrap_draws(cv_rank) = minf;
          top_bootstrap_draws = sort(top_bootstrap_draws,"descend");
          nlopt_set_stopval(opt, top_bootstrap_draws(cv_rank));
          //for (int m=0; m <= cv_rank; m++) Rcpp::Rcout << top_bootstrap_draws(m) << std::endl;
        }
      } else {
          top_bootstrap_draws(r) = minf; 
          if  (r == cv_rank) {
            top_bootstrap_draws = sort(top_bootstrap_draws,"descend");
            //for (int m=0; m <= cv_rank; m++) Rcpp::Rcout << top_bootstrap_draws(m) << std::endl;
            nlopt_set_stopval(opt, top_bootstrap_draws(cv_rank));
          }
      }
  }
  nlopt_destroy(opt);
  DRBcv = top_bootstrap_draws(cv_rank);
  if (echo) Rcpp::Rcout << "DR Bootstrap critical value for level " << p << " is " << DRBcv << std::endl;
  return(DRBcv);
}

/*#####################################################################################################################################
# FUNCTION MinDiscrepancyCpp minimizes the discrepancy function with respect to the full parameter vector using
#   each column of theta_start_vals as a starting value. Note that minimization automatically stops when a parameter vector
#   is found such that the discrepancy function is zero or negative, such that all inequalities are satisfied. The parameter
#   found is then known to be in the analog set estimate, and all other confidence sets and set estimates.
# INPUT ARGUMENTS 
#   theta_start_vals: A K x J matrix in which each column 1,...,J represents a different starting value from which to
#     start search for a minimizer of the discrepancy wrt to the K-variate parameter theta.
#   theta_lb: lower bounds on parameter vector passed as an Armadillo vec.
#   theta_ub: upper bounds on parameter vector passed as an Armadillo vec.
#   argl: A list containing additional variables, data, and other quantities pre-computed in R.
# RETURNS
#   A list comprising the following elements:
#       $disc_solns: A K x J matrix in which each column represents the value at which minimization stopped when search
#         started from the j-th column of theta_start_vals
#       $disc_vals: A J-vector whose j-th entry is minimum value of the discrepancy function obtained in the j-th search.
#       $param_vals: The value of the target parameter at the value of theta obtained in the j-th search.
#######################################################################################################################################*/
//[[Rcpp::export]]
List MinDiscrepancyCpp(mat theta_start_vals, vec theta_lb, vec theta_ub, List argl) {
  bool echo = as<bool>(argl["echo"]);
  double tol_array[2] =  {1e-8, 1e-8}; // For imposing boundary constraints
  if (echo) Rcpp::Rcout << "Running MinDiscrepancyCpp" << std::endl;
  std::string target_parameter = as<std::string>(argl["target_parameter"]);
  clock_t start, end;
  List global_opts = as<List>(argl["opts"]);
  List retlist;
  int k = theta_start_vals.n_rows;
  double* theta_arr = new double [k];
  double* lb_arr = new double [k];
  double* ub_arr = new double [k];
  vec min_fvals(theta_start_vals.n_cols);
  min_fvals.fill(datum::inf);
  vec param_vals(theta_start_vals.n_cols);
  param_vals.fill(datum::inf);
  mat min_fvecs(k,theta_start_vals.n_cols, fill::zeros);
  for (int j=0; j < k; j++) {
    lb_arr[j] = theta_lb(j);
    ub_arr[j] = theta_ub(j);
  }
  double minf = 0.0;
  unsigned j = 0;
  while (j < theta_start_vals.n_cols) {
    int exit_code;
    if (echo) Rcpp::Rcout << "Beginning search using starting parameter vector, " << j+1 << " of " << theta_start_vals.n_cols << "." << std::endl;
    for (int h=0; h < k; h++) theta_arr[h] = theta_start_vals(h,j);
    if (echo) Rcpp::Rcout << "Starting value is ";
    if (echo) DisplayArray(theta_arr, k);
    double startfval = GaussianDiscrepancy(k, theta_arr, NULL, &argl);
    if (echo) Rcpp::Rcout << "With function value " << startfval << std::endl;
    fcount_evaluations = 0;
    start = clock();
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, k);
    nlopt_set_min_objective(opt, GaussianDiscrepancy, &argl);
    nlopt_set_lower_bounds(opt, lb_arr);
    nlopt_set_upper_bounds(opt, ub_arr);
    nlopt_set_ftol_rel(opt, as<double>(global_opts["ftol_rel"]));
    nlopt_set_ftol_abs(opt, as<double>(global_opts["ftol_abs"]));
    nlopt_set_xtol_rel(opt, as<double>(global_opts["xtol_rel"]));
    nlopt_set_maxeval(opt, as<int>(global_opts["maxeval"]));
    nlopt_set_stopval(opt, as<double>(global_opts["stopval"]));
    if (target_parameter == "CME") {
      nlopt_add_inequality_mconstraint(opt, 2, CME_bound_constraints, &argl, tol_array);
    } else if (target_parameter == "CRP") {
      nlopt_add_inequality_mconstraint(opt, 2, CRP_bound_constraints, &argl, tol_array);
    }
    exit_code = nlopt_optimize(opt, &(theta_arr[0]), &minf);
    end = clock();
    min_fvals(j) = minf;
    for (int h=0; h < k; h++) min_fvecs(h,j) = theta_arr[h]; 
    if (echo) DisplaySearchResult(theta_arr, minf, theta_start_vals.n_rows);
    if (echo) {
      if (target_parameter == "CME") {
        param_vals(j) = CME(min_fvecs.col(j), argl);
        Rcpp::Rcout << "The corresponding value of the CME is " << param_vals(j) << endl;
      } else if (target_parameter == "CRP") {
        param_vals(j) = ConRProb(min_fvecs.col(j), argl);
        Rcpp::Rcout << "The corresponding value of the CRP is " << param_vals(j) << endl;
      }
    }
    if (echo) Rcpp::Rcout << "NLOPT took " << double(end-start)/CLOCKS_PER_SEC << " seconds to run " << fcount_evaluations << " function evaluations. NLopt exit code: " << exit_code << "." << std::endl;
    fcount_fulldisc_evaluations += fcount_evaluations;
    nlopt_destroy(opt);
    j++;
  }
  if (echo) Rcpp::Rcout << "MinDiscrepancyCpp completed a total of " << fcount_fulldisc_evaluations << " function evaluations." << std::endl;
  retlist["disc_solns"] = min_fvecs;
  retlist["disc_vals"] = min_fvals;
  retlist["param_vals"] = param_vals;
  delete[] theta_arr;
  delete[] lb_arr;
  delete[] ub_arr;
  return(retlist);
}
/* ------------------------------- END OF FUNCTION MinDiscrepancyCpp ---------------------------------   */

/*#####################################################################################################################################
# FUNCTION ProfileDiscrepancyCpp evaluates the profile discrepancy function \hat{T}_n(g) by minimizing the discrepancy function
#   \hat{Q}(theta) subject to g(theta) = gvalue.
# INPUT ARGUMENTS 
#   gvalue: The value at which g(theta) is fixed.
#   argl: A list containing additional variables, data, and other quantities pre-computed in R.
#   si: Stands for "search item", a list containing various information from R about the search being conducted.
# RETURNS
#   A list comprising the following elements:
#       $objective: The value of the discrepancy function obtained at the constrained minimum.
#       $solution: The parameter vector theta at which g(theta) = gvalue that achievesthe constrained minimum.
#       $DR_boot_cv: The value of the discard resampling critical value for the given gvalue.
#######################################################################################################################################*/
List ProfileDiscrepancyCpp(double gvalue, List argl, List si) {
  bool echo = as<bool>(argl["echo"]);
  if (echo) Rcpp::Rcout << "Running ProfileDiscrepancyCpp" << std::endl;
  std::string target_parameter = as<std::string>(si["target_parameter"]);
  std::string param_name = as<std::string>(si["param_name"]);
  List profile_opts = as<List>(argl["profile_opts"]);
  List retlist;
  int p_index = as<int>(si["param_index"]);
  bool runDR = as<bool>(argl["runDR"]);
  clock_t start, end;
  vec theta_lb = as<vec>(argl["theta_lb"]);
  vec theta_ub = as<vec>(argl["theta_ub"]);
  vec theta_start = as<vec>(argl["theta_start"]);
  int k = theta_lb.n_elem;
  double* theta_arr = new double [k];
  double* lb_arr = new double [k];
  double* ub_arr = new double [k];
  for (int j=0; j < k; j++) {
    theta_arr[j] = theta_start(j);
    lb_arr[j] = theta_lb(j);
    ub_arr[j] = theta_ub(j);
  }
  double minfval = 0.0;
  
  /* -----------Compute Profile Discrepancy for g(theta) value **************/
  start = clock();
  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, k);
  if (target_parameter == "parameter_x") {
    lb_arr[p_index-1] = gvalue;
    ub_arr[p_index-1] = gvalue;
    if (echo) Rcpp::Rcout << "Searching for profiled discrepancy at theta[" << p_index << "] = " << gvalue << std::endl;
  } 
  else {
    if (echo) Rcpp::Rcout << "Searching for profiled discrepancy at " << param_name << " = " << gvalue << endl;
    argl["gvalue"] = gvalue;
    if (target_parameter == "CME") {
      nlopt_add_equality_constraint(opt, CME_constraint, &argl, as<double>(profile_opts["tol_constraints_eq"]));
    }
    else if (target_parameter == "CRP") {
      nlopt_add_equality_constraint(opt, ConRProb_constraint, &argl, as<double>(profile_opts["tol_constraints_eq"]));
    }
  }
  nlopt_set_min_objective(opt, GaussianDiscrepancy, &argl);
  nlopt_set_lower_bounds(opt, lb_arr);
  nlopt_set_upper_bounds(opt, ub_arr);
  nlopt_set_ftol_rel(opt, as<double>(profile_opts["ftol_rel"]));
  nlopt_set_ftol_abs(opt, as<double>(profile_opts["ftol_abs"]));
  nlopt_set_xtol_rel(opt, as<double>(profile_opts["xtol_rel"]));
  nlopt_set_maxeval(opt, as<int>(profile_opts["maxeval"]));
  nlopt_set_stopval(opt, as<double>(profile_opts["stopval"]));
  fcount_evaluations = 0;
  int exit_code = nlopt_optimize(opt, &(theta_arr[0]), &minfval);
  nlopt_destroy(opt);
  end = clock();
  if (echo) {
    Rcpp::Rcout << "Search concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds with exit code " << exit_code << "." << std::endl;
    DisplaySearchResult(p_index, gvalue, theta_arr, fcount_evaluations, minfval, k);
  }
  fcount_refining_evaluations += fcount_evaluations;
  argl["last_profile_objective"] = minfval;
  retlist["objective"] = minfval;
  rowvec argmin = rowvec(k);
  for (int j=0; j < k; j++) argmin(j) = theta_arr[j];
  retlist["solution"] = argmin;
  retlist["DR_boot_cv"] = 0;

  if ((minfval > 0) && runDR) { // If the value of the objective function is less than 0, no need to bootstrap
      if (echo) {
        if (runDR) Rcpp::Rcout << "Computing Discard Resampling Bootstrap Critical Value." << std::endl;
      }
      start = clock();
      double bootCV = BootstrapCV(gvalue, theta_arr, lb_arr, ub_arr, k, argl);
      end = clock();
      if (echo) Rcpp::Rcout << "Bootstrapping " << as<int>(argl["R"]) << " iterations concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds with CV of " << bootCV << std::endl << std::endl;
      retlist["DR_boot_cv"] = bootCV;
    }
  delete[] theta_arr;
  delete[] lb_arr;
  delete[] ub_arr;
  return(retlist);
}
/* ------------------------------- END OF FUNCTION ProfileDiscrepancyCpp ---------------------------------   */

/*#####################################################################################################################################
# FUNCTION ProfileDiscrepancyOnGridCpp evaluates the profile discrepancy function \hat{T}_n(g) by minimizing the discrepancy function
#   \hat{Q}(theta) subject to g(theta) = gvalue for each gvalue on the grid of values "search_grid".
# INPUT ARGUMENTS 
#   search_grid: Grid of values at which to evaluated the profiled discrepancy function \hat{T}(.).
#   g_optindex_r: Grid index of the value at which the smallest value of the objective function was attained in the full parameter search.
#   theta_start: K x 2 matrix indicating the parameter vector at which the smallest and largest values of g such that g(theta) <= 0 have been found.
#   theta_lb: lower bounds on parameter vector passed as an Armadillo vec.
#   theta_ub: upper bounds on parameter vector passed as an Armadillo vec.
#   stopval: objective function value at which to stop searching for a lower value, e.g. zero.
#   arglist: A list containing additional variables, data, and other quantities pre-computed in R.
#   si: Stands for "search item", a list containing various information from R about the search being conducted.
# RETURNS
#   The number of seconds needed for execution.
#   In addition, the following quantities are passed by reference in the si list:
#     $grid_vals: solution found for the constrained optimization problem at each point in search_grid.
#     $grid_CVs: DRB critical value computed for each point in search_grid (simply zero if \hat{T}(g) = 0 at a given point g.)
#     $grid_vecs: A matrix in which each column j represents the parameter vector solving the associated constrained minimization
#                 for the j-th components of search_grid.
#######################################################################################################################################*/
//[[Rcpp::export]]
double ProfileDiscrepancyOnGridCpp(vec search_grid, int g_optindex_r, mat theta_start, vec theta_lb, vec theta_ub, double stopval, List arglist, List si) {
  bool echo = as<bool>(arglist["echo"]);
  Rcpp::Rcout << "Running ProfileDiscrepancyOnGridCpp" << std::endl;
  std::string target_parameter = as<std::string>(si["target_parameter"]);
  std::string param_name = as<std::string>(si["param_name"]);
  List profile_opts = as<List>(arglist["profile_opts"]);
  int p_index = as<int>(si["param_index"]);
  bool runDR = as<bool>(arglist["runDR"]);
  clock_t start, end;
  clock_t total_start, total_end;
  total_start = clock();
  int n_gpts = search_grid.n_elem;
  int k = theta_start.n_rows;
  double* theta_arr = new double [k];
  double* lb_arr = new double [k];
  double* ub_arr = new double [k];
  for (int j=0; j < k; j++) {
    theta_arr[j] = theta_start(j,0);
    lb_arr[j] = theta_lb(j);
    ub_arr[j] = theta_ub(j);
  }
  double minfval = 0.0;
  ivec exit_codes = ivec(n_gpts); 
  int i = 0;
  int grid_index = g_optindex_r;
  bool below_g_start = true;
  vec grid_vals = as<vec>(si["grid_vals"]);
  vec grid_CVs = as<vec>(si["grid_CVs"]);
  mat grid_vecs = as<mat>(si["grid_vecs"]);
  
  /* -------- while loop over grid points --------- */
  while (i < n_gpts) { 
    if (below_g_start) {
      grid_index--;
    } else {
      grid_index++;
    }
    if (grid_index == -1) {
      grid_index = g_optindex_r;
      below_g_start = false;
      for (int j=0; j < k; j++) theta_arr[j] = theta_start(j,1);
    }
    double gvalue = search_grid(grid_index);
    start = clock();
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, k);
    nlopt_set_min_objective(opt, GaussianDiscrepancy, &arglist);
    /*********Compute Profile Discrepancy for g(theta) value **************/
    if (target_parameter == "parameter_x") {
      lb_arr[p_index-1] = gvalue;
      ub_arr[p_index-1] = gvalue;
      Rcpp::Rcout << "Searching for profiled discrepancy at theta[" << p_index << "] = " << gvalue << endl;
    } 
    else if (target_parameter == "CME") {
      Rcpp::Rcout << "Searching for profiled discrepancy at " << param_name << " = " << gvalue << endl;
      arglist["gvalue"] = gvalue;
      nlopt_add_equality_constraint(opt, CME_constraint, &arglist, as<double>(profile_opts["tol_constraints_eq"]));
    } 
    else if (target_parameter == "CRP") {
      Rcpp::Rcout << "Searching for profiled discrepancy at " << param_name << " = " << gvalue << endl;
      arglist["gvalue"] = gvalue;
      nlopt_add_equality_constraint(opt, ConRProb_constraint, &arglist, as<double>(profile_opts["tol_constraints_eq"]));
    }
    nlopt_set_lower_bounds(opt, lb_arr);
    nlopt_set_upper_bounds(opt, ub_arr);
    nlopt_set_ftol_rel(opt, as<double>(profile_opts["ftol_rel"]));
    nlopt_set_ftol_abs(opt, as<double>(profile_opts["ftol_abs"]));
    nlopt_set_xtol_rel(opt, as<double>(profile_opts["xtol_rel"]));
    nlopt_set_maxeval(opt, as<int>(profile_opts["maxeval"]));
    nlopt_set_stopval(opt, stopval);
    fcount_evaluations = 0;
    exit_codes[grid_index] = nlopt_optimize(opt, &(theta_arr[0]), &minfval);
    end = clock();
    if (echo) Rcpp::Rcout << "Search concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds with exit code " << exit_codes[grid_index] << "." << std::endl;
    DisplaySearchResult(p_index, gvalue, theta_arr, fcount_evaluations, minfval, k);
    fcount_profiled_evaluations += fcount_evaluations;
    grid_vals(grid_index) = minfval;
    arglist["last_profile_objective"] = minfval;
    for (int j=0; j < k; j++) {
      grid_vecs(j,grid_index) = theta_arr[j]; 
    }
    nlopt_destroy(opt);
    
    if (minfval <= 0) { // If the value of the objective function is less than 0, no need to bootstrap the CV
      grid_CVs(grid_index) = 0;
    } else if (runDR) {
      if (echo) {
        Rcpp::Rcout << "Computing bootstrap critical value." << std::endl;
        if (runDR) Rcpp::Rcout << "Discard Resampling Bootstrap (DR)." << std::endl;
      }
      start = clock();
      double bootCV = BootstrapCV(gvalue, theta_arr, lb_arr, ub_arr, k, arglist);
      end = clock();
      if (echo) Rcpp::Rcout << "Bootstrapping " << as<int>(arglist["R"]) << " iterations concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds with CV of " << bootCV << std::endl;
      grid_CVs(grid_index) = bootCV;
    }
    i++;
  }
  /* -------- end of while loop over grid points --------- */
  si["grid_vals"] = grid_vals;
  si["grid_CVs"] = grid_CVs;
  si["grid_vecs"] = grid_vecs;
  delete[] theta_arr;
  delete[] lb_arr;
  delete[] ub_arr;
  total_end = clock();
  return(double(total_end-total_start)/CLOCKS_PER_SEC);
}

/** ------------------- FUNCTION RefineBoundCpp -------------------------
 
 
######## Function RefineBoundCpp ####################################################################
#   This function is used to refine a previously performed grid search to find the boundary of a non-rejection interval
#   for a univariate functional of a partially identified parameter vector. It takes as an argument a value on a previously
#   explored grid found to be rejected by the inclusion criteria, and the closest point on the grid not rejected by the
#   inclusion criteria.  The function iteratively tests values from the rejected point closer and closer to the non-rejected
#   point until a point is found that is not rejected.  It stops at the first point not rejected and returns that value.
#   Input arguments:
#   outerbound: The value at which to start the search.
#   innerbound: The value at which to stop the search, an extreme point of those already checked for inclusion.
#   fromLower: T for a lower bound search, F for an upper bound search.
#   increment: The difference between points to incrementally test until there is non-rejection.
#   critical_value: The point is rejected if the profiled discrepancy is found to exceed this value.
#   param_start: The parameter vector at which to begin the profiled search at the initial outerbound.
#   echo: Boolean indicating whether to display output regarding sequential tests.
#   argl: arglist variable list used for the profile discrepancy search.
################################################################################################################*/
//[[Rcpp::export]]
List RefineBoundCpp(double outerbound, double innerbound, bool fromLower, double increment, double default_cv, vec param_start, List argl, List si) {
  bool echo = as<bool>(argl["echo"]);
  bool runDR = as<bool>(argl["runDR"]); // Note: default_cv is ignored if runDR == true.
  bool first_iteration = true;
  clock_t start, end;
  bool blnReject = true;
  int sgnIsLower = 2*fromLower - 1;
  double val = outerbound;
  argl["theta_start"] = param_start;
  List profile_nloptr_retobj;
  List returnlist;
  while (blnReject && ( (sgnIsLower * (innerbound - val)) > 0) ) {
    val += (sgnIsLower * increment);
    if (echo) Rcpp::Rcout << "Testing value " << val << " with starting value from last iteration's optimum." << endl;
    start = clock();
    profile_nloptr_retobj = ProfileDiscrepancyCpp(val, argl, si);
    end = clock();
    if (echo) Rcpp::Rcout << "Test concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds. Discrepancy of " << as<double>(profile_nloptr_retobj["objective"]) << " obtained." << endl;
    blnReject = (as<double>(profile_nloptr_retobj["objective"]) > (runDR ? (as<double>(profile_nloptr_retobj["DR_boot_cv"])) : default_cv) );
    if (blnReject && (!first_iteration)) { // Check if using the inner point optimum as a starting value could produce a lower discrepancy
      if (echo) Rcpp::Rcout << "----- EXECUTING INSIDE IF CONDITION ------ " << endl;
      argl["theta_start"] = param_start;
      if (echo) Rcpp::Rcout << "Testing value " << val << " with inner point starting value." << endl;
      start = clock();
      List profile_nloptr_retobj2 = ProfileDiscrepancyCpp(val, argl, si);
      end = clock();
      if (echo) Rcpp::Rcout << "Test concluded in " << double(end-start)/CLOCKS_PER_SEC << " seconds. Discrepancy of " << as<double>(profile_nloptr_retobj2["objective"]) << " obtained." << endl;
        blnReject = (as<double>(profile_nloptr_retobj2["objective"]) > (runDR ? (as<double>(profile_nloptr_retobj2["DR_boot_cv"])) : default_cv));
        if (!blnReject) profile_nloptr_retobj = profile_nloptr_retobj2;
      }
    argl["theta_start"] = profile_nloptr_retobj["solution"];
    first_iteration = false;
    }
  returnlist["bound"] = val;
  returnlist["theta"] = profile_nloptr_retobj["solution"];
  returnlist["discrepancy"] = profile_nloptr_retobj["objective"];
  return(returnlist);
}

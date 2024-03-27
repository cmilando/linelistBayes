// -------------------------------------------------------------------------------
// Author: Li Tenglong
// Annotations by CWM
// Date started: 12.20.2023
// Notes:
// * anytime you see something in double square brackets, 
//   its important for the compiler
//
// Functions in this file: (all in this file until you make a package)
//  [x] findmiss - finds missing values in a vector
//  [x] dnb      - get the PMF value for a given x in a standard negative binomial distribution
//  [x] pnb      - get the CDF value for a given x in a standard negative binomial distribution
//  [x] rnb      - get a random draw from a negative binomial distribution
//  [x] get_mu_vec - function to get a vector of mu, line-specific estimates of mean of right-trunc neg.binom
//  [x] dummy    - create a matrix of indicator values for is_week_i and is_weekend
//  [ ] logLikNB - calculating the log-likelihood of the right-truncated negative binomial distr.
//  [ ] prop     - function to make proportion of counts from 0 to maxdelay
//  [ ] lambda   - create a lambda function to compute mean of poisson dist.
//  [ ] getr     - curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
//  [ ] backnow  - main function to do Gibbs MCMC, impute missing delays, back-calculation, then now-casting
//
// -------------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

// -------------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector findmiss(NumericVector x){
  // Description:
  //   function that returns the indices of NA values of a vector
  //
  // Arguments:
  //   x: a numeric vector, -inf to inf
  //
  // Returns: 
  //   the indices of x for which is_na(x) returns True
  IntegerVector y = seq(0, x.size() - 1);
  LogicalVector z = is_na(x);
  return y[z];
}

// -------------------------------------------------------------------------------
double dnb(double x, double r, double m){
  // Description:
  //   dnb function used by mapply for likelihood.
  //   These functions provide information about the 
  //   Negative binomial distribution with dispersion r and mean m
  //   R:: implementation chosen for speed, returns a scalar
  //
  // Arguments:
  //   x: a double random variable, 0 to inf
  //   r: the dispersion parameter, higher r means tighter around mu 
  //   m: the mean, mu
  //
  // Returns: 
  //   the Probability density function value at x 
  //   of a distribution given by r and m
  
  //   implied in the function:
  //   log: true --> probabilities p are given as log(p)
  double y = R::dnbinom_mu(x, r, m, true);
  return y;
}

// -------------------------------------------------------------------------------
double pnb(int x, double r, double m){
  // Description:
  //   pnb function used by mapply for likelihood.
  //   These functions provide information about the 
  //   Negative binomial distribution with dispersion r and mean m
  //   R:: implementation chosen for speed, returns a scalar
  //
  // Arguments:
  //   x: a double random variable, 0 to inf
  //   r: the dispersion parameter, size. higher r means tighter around mu 
  //   m: the mean, mu > 0
  //
  // Returns: 
  //   the Cumulative density function value at x 
  //   of a distribution given by r and m
  
  // implied in the function:
  //   lower: true --> Calculate the probability of the region where 
  //                 the random variable is less than or equal to x
  //   log: true --> probabilities p are given as log(p)
  // 
  double y = R::pnbinom_mu(x, r, m, true, true);
  return y;
}

// -------------------------------------------------------------------------------
double rnb(double r, double p){
  // Description:
  //   rnb function used by mapply for likelihood.
  //   These functions provide information about the 
  //   Negative binomial distribution with dispersion r and mean m
  //   R:: implementation chosen for speed, returns a scalar
  //
  // Arguments:
  //   r: the dispersion parameter, size
  //   p: a probability, ...?
  //
  // Returns: 
  //   a single random value in this distribution, given p?
  
  double y = R::rnbinom(r, p);
  return y;
}

// -------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector get_mu_vec(NumericMatrix x12, NumericVector beta) {
  // Description:
  //   Function to compute mu=exp(X1*beta + X2*gamma), the parameter representing the mean value 
  //
  // Arguments:
  //   x12:  matrix of indicator variables for is_week_i and is_weekend, n x p
  //   beta: existing week-wise beta coefficients and a single weekend parameter, p x 1
  //         this includes what is called beta in the text
  //         as well as gamma, the single parameter for weekend
  //         but shortened to just beta for this script
  //
  // Returns: 
  //   a vector, exp(x12 %*% beta)
  
  // Conversions into armadillo syntax for ease of processing in 
  // Armadillo linear algebra
  // the a_ represents the same variable as above but in armadillo
  arma::mat a_x12  = as<arma::mat>(x12);  // dim(x1) is n observations x n weeks
  arma::vec a_beta = as<arma::vec>(beta); // beta are corresponding model coefficients

  // matrix multiplication, X1*B + X2*gamma
  arma::mat a_mu_mat = a_x12 * a_beta;

  // wrap() is conversion back into an R object
  NumericMatrix mu_mat = wrap(a_mu_mat);

  // getting all the rows (_) of the first column (0) of lc
  // it MUST have one 1 column
  if (mu_mat.ncol() != 1) {
        stop("ncol of mu_vector must = 1");
  }
  NumericVector mu_vec = mu_mat(_, 0); 

  // taking the exponent and returning
  NumericVector exp_mu_vec  = exp(mu_vec);

  return exp_mu_vec;
}

// -------------------------------------------------------------------------------
// Create Matrix of dummy variables
// [[Rcpp::export]]
NumericMatrix dummy(IntegerVector week, IntegerVector weekend) {
  // Description:
  //   
  //
  // Arguments
  //   week:     d$week of the report, allows you to assume that things don't change within this window
  //   weekend:  d$weekend, 0 or 1, indicator of weekend (1)
  //
  // Returns:
  //   cov: a matrix of dummy variables where the first n columns are is_week_i
  //        and the last column is is_weekend. the rows are data for individuals
  //
  
  int nw = max(week);   // maximum number of weeks
  int nc = nw + 1;      // adds an extra column for weekends
  int n  = week.size(); // this is the rows of line-list data (so 1 per person) 
  NumericMatrix cov(n, nc); // the final matrix has a row for each person and a column for each week + 
                            // one additional column for is_weekend

  for (int i = 0; i < nc; ++i) {
    if (i < nw) {
      cov(_,i) = ifelse(week == (i+1), 1, 0); // see above
    } else {
      cov(_,i) = weekend; // see above
    }
  }
  
  return cov;
}


// -------------------------------------------------------------------------------
// [[Rcpp::export]]
double logLikNB(NumericVector delay_vec, NumericMatrix x12, 
            NumericVector disp, NumericVector betaplus, int maxdelay) {
  // Description:
  //   TL: faster version of like for large sample
  //   CM: this is the log-likelihood function for the right-truncated negative binomial distribution
  //
  // Arguments
  //   delay_vec:  vector, reporting delay vector
  //   x12:        matrix, cov, the indicator values for is_week_i & is_weekend
  //   disp:       vector, the flag for each outcome, =1 if cd is null, otherwise =1 and =2 if some other condition is met
  //                       I think this is the dispersion associated with changing cd? unclear
  //   betaplus:   vector, current row of parameter estimates, betas(gammas) and .... other params???
  //   maxdelay:   integer, maximum reporting delay
  //
  // Returns:
  //   the logLiklihood of the right-truncated neg. binom distribution given the current param set
  
  // ----------------
  // get the number of columns of the indicator matrix
  int nbeta = x12.ncol();        // number of parameters estimated for is_week_i and is_weekend
  int n     = delay_vec.size();  // number of people in the line-list dataset
  NumericVector s_vec;       // vector of the size parameters
  double r1;                 // either the first index of the two sizes for the two dispersion or not
  
  // ----------------
  // splits betaplus into two vectors, from 0:(nbeta-1), and then nbeta:extra
  // beta = the beta coefficients for each week and is_weekend
  NumericVector beta = betaplus[seq(0, nbeta - 1)];   

  // the remainder are the other parameters being estimated which are ...?
  // TODO: The additional parameters that are being estimated, which are ...?
  // if cd is defined then this is just one extra
  // if cd is not defined then there are two extra
  NumericVector r    = betaplus[seq(nbeta, betaplus.size() - 1)]; // these are the additional params, either 1 or 2
  if (! ((r.size() == 1) || (r.size() == 2))) {
    stop("size of r is not in c(1, 2)");
  }

  // ----------------
  // and then mu_vec is the expected value for the mean for each person's reporting delay distribution
  // get mu, = exp(X1*beta + X2*gamma), one value for each person, E(delay)
  NumericVector mu_vec   = get_mu_vec(x12, beta); 

  // ----------------
  // TODO: ???
  // (basically, if this is type = 0, if cd is not defined)
  if (max(disp) == 1){ // originally this was max but we probably want all() right?
    // ### ERROR CHECKING ###
    // if(! (all(disp == 1))) {
    //   stop("all disp != 1");
    // }
    if(! (r.size() == 1)) {
      stop("r.size() must = 1");
    }
    // ####
    r1 = as<double>(r); // type conversion to a scalar
    s_vec = disp * r1;  // you know that disp = 1 everywhere, so this is just repeated r right?
  } else {
    // ### ERROR CHECKING ###
    // these is the case where cd is defined and where disp == 2
    if(! (r.size() == 2)) {
      stop("r.size() must = 2");
    }
    //if(! ((max(disp) == 2) && (min(disp) == 1))) {
    //  stop("all disp not in 1 or 2");
    //}
    // ####
    r1 = r[0];            // set r1 equal to the first dispersion parameter
    s_vec = rep(r1, n);   // right, why not do this above?
    LogicalVector select = (disp==2);
    double r2 = r[1];     // the two dispersions
    s_vec[select] = r2;
  }
  
  // ************
  // Now get the log-liklihood of the right-truncated negative binomial distribution
  // This is manifested as the ratio between P(x[i]) and sum{P(x[i] <= K)} for each combination of variables
  // (i.e., size = s[i] and mean = mu[i])

  // (1) First, the numerator
  //     for each reporting delay value (x[i] in delay_vec)
  //     calculate the P(x[i]), the probability of that delay
  //     given a negative binomial distribution defined by size = s[i] and mean = mu[i]
  NumericVector nbinom_pmf_vec = mapply(delay_vec, s_vec, mu_vec, dnb);

  // (2) Then, calculate the denominator(?)
  //     create a vector that looks like delay_vec but is instead the max reporting delay
  IntegerVector maxdelay_vec  = rep(maxdelay, n);

  // and then get the CDF value for max reporting delay for a distribution given size? = s, and mean = m1
  // this varies for every value because the size parameter and mu are different for each distribution
  NumericVector max_nbinom_cdf = mapply(maxdelay_vec, s_vec, mu_vec, pnb);

  // then, the result (log-likelihood of all the right-truncated negative binomial distribution^s^) 
  // is the difference of all summed numerators and denominators
  double result = sum(nbinom_pmf_vec) - sum(max_nbinom_cdf);
  // ************
  
  //
  return result;
}



// -------------------------------------------------------------------------------
NumericVector prop(NumericVector x, NumericVector onset, int maxdelay, int cd) {
  // Description:
  //    function to make proportion of counts from 0 to maxdelay
  //
  // Arguments
  //
  // Returns:
  //
  
  //
  LogicalVector v = (x <= maxdelay) & (onset >= cd);
  NumericVector x1 = x[v];
  int dem = x1.size();
  NumericVector p1 (maxdelay);

  //
  for (int i = 0; i < maxdelay; ++i){
    p1[i] = sum(x1 == (maxdelay - i));
  }

  NumericVector p = p1 / dem;
  NumericVector p2 = cumsum(p);
  NumericVector result = 1 - p2;

  return result;
}


// -------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector lambda(NumericVector curve, NumericVector si) {
  // Description:
  //    create a lambda function to compute mean of poisson dist.
  //    function used in getr and updater
  //    curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
  //
  // Arguments
  //
  // Returns:
  //
  int k = si.size();
  int n = curve.size();
  NumericVector c;
  NumericVector sic;
  NumericVector result (n-1);
  for (int i = 0; i < (n-1); ++i){
    if (i < k) {
      c = curve[seq(0, i)];
      sic = rev(si[seq(0, i)]);
      result[i] = sum(c * sic);
    } else {
      c = curve[seq(i - k + 1, i)];
      sic = rev(si);
      result[i] = sum(c * sic);
    }
  }
  return result;
}

// -------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector getr(NumericVector curve, NumericVector si, int size){
  // Description:
  //    curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
  //    Compute the mean only
  //
  // Arguments
  //    curve: i.e., the current backcaulation, dim of (nd + maxdelay)   
  //    si: serial interval of 14?
  //    size: NB_size (default is 6)
  // Returns:
  //
  
  int n = curve.size();
  int nr = n - size - 1;
  NumericVector result (nr);
  double shape;
  double scale;
  NumericVector incid;
  NumericVector dem = lambda(curve, si);
  NumericVector dem1;
  for (int i = 0; i < nr; ++i){
    incid     = curve[seq(i + 1, i + size + 1)];  //daily counts for days i to i + size
    shape     = sum(incid) + 1;
    dem1      = dem[seq(i, i + size)];
    scale     = 1 / (sum(dem1) + 0.2);
    result[i] = shape * scale ; 
  }
  return result;
}

// -------------------------------------------------------------------------------
// #' @export
// [[Rcpp::export]]
List backnow_cm(NumericVector outcome, NumericVector days, 
             IntegerVector week, IntegerVector weekend, 
             int iter, double sigma, int maxdelay, 
             NumericVector si, int size, 
             int workerID,
             int printProgress,
             Nullable<int> cd = R_NilValue){
  // Description:
  //
  //
  // Arguments:
  //   outcome:  d$delay, which = report - onset
  //   days:     d$report, the day that the report is given on, center on min(report), so 1, ...., max()
  //   week:     d$week of the report, allows you to assume that things don't change within this window
  //   weekend:  d$weekend, 0 or 1, indicator of weekend (1)
  //   iter:     the number of iterations to do in the big loop
  //   sigma:    set to 0.2, the standard deviation of the normal distribution used for parameter draws
  //   maxdelay: the maxmimum reporting delay, 20 days
  //   si:       the PDF of the serial interval by day, stands in for the generation time between infection and onset
  //   size:     ?? 6, seems like tau, the sliding window size in days
  //   workerID: workerID number
  //   printProgress: prints progress using worker ID
  //   cd:       ?? I believe its a placeholder for incorrect reporting delays?
  //              Note: cd must be earlier than the first day of nowcasting period
  //
  // Returns:
  //   Back:
  //   R:

  // --------------------------------
  // VARIABLE DEFINITIONS
  int cday;                                 // ??? unclear, a placeholder for the incorrect reporting delay?
  int type;                                 // ??? a flag for something, either 0 or 1
  int nc;                                   // ??? n parameters to define, 1 for every week + either 2 or 3
  int mi;                                   // ???
  int dv;                                   // ???
  int nm;                                   // ???

  double ratio;                             // ??? what is ratio
  double decision;                          // ??? what is decision
  double r1;                                // ???
  double m;                                 // ???

  NumericMatrix cov;                        // ??? 
  NumericVector par0;                       // create a placeholder for initial parameter estimates
  NumericVector oldpar;                     // initialize the loop oldpar vector
  NumericVector newpar;                     // initialize the loop newpar vector
  NumericVector mm;                         // ???
  NumericVector par1;                       // ???
  NumericVector mean;                       // ???
  NumericVector prob;                       // ???
  NumericVector p;                          // ???
  NumericVector s;                          // ???

  // Variables that depend on inputs
  NumericVector outcome1 = clone(outcome);           // ??? creates a copy of the outcome: why?
  int n  = outcome1.size();                          // the length of the outcome vector
  NumericVector dpind = rep(1.0, n);                 // ??? dispersion param changes? specific flag for each outcome
  int nw = max(week);                                // total number of weeks reporting data exist over
  int nd = max(days);                                // total number of days reporting data exist over
  IntegerVector x = seq(0, maxdelay);                // ???
  NumericVector dayseq = as<NumericVector>(x);       // ???
  NumericMatrix rt (iter, nd + maxdelay - size - 1); // the r(t) variable at the end, minus size??
  
  // -------------------------------- 
  // INITIALIZE VARIABLE SIZES & FLAGS
  // Ok so `nc` is the number of parameters that are being estimated
  // there are one for every week, and 3 extra if c is specified, 2 otherwise
  // perhaps this is the weekly value of the reporting delay, which doesn't change much within week
  if (cd.isNotNull()){      // if cd (which is ...???) is specified as input
    cday = as <int> (cd);   // make sure that its an integer
    nc   = nw + 3;          // set nc to nw + 3, which is + 1 for is_weekend then + 2 for TODO:... reasons???
    LogicalVector kk = (days > cday);  // for any report that has a day > reporting day max, set kk = T, otherwise kk = F
    dpind[kk] = 2;          // for any delay where kk=2, set dpind flag to 2    
    type      = 1;          // this type of model is called type 1
  } else {
    cday = -maxdelay;       // otherwise if cd is not specified, set it equal to -1 * maxdelay
    nc   = nw + 2;          // and nc now equals nw + 2, so +1 for is_weekend then +1 for ....?
    type = 0;               // this type of model is called type 0
  }

  // So now you can define parameter and back
  // matrix of Bayesian parameter estimates, and how they evolve over iterations
  NumericMatrix parameter (iter, nc);       
  // matrix of back-calculated counts; and how they evolve over iterations
  NumericMatrix back (iter, nd + maxdelay); 

  // --------------------------------
  // FILL IN MISSING
  // find all the missing delays (aka when report - onset = NA, because onset was missing)
  IntegerVector missind = findmiss(outcome1);  
  // the number of missing indices
  int nmiss = missind.size();                  
  // randomly sample from 1 to maxdelay, with replacement
  NumericVector miss0 = as<NumericVector>(sample(maxdelay, nmiss, true)); 
  // fill these values to delay vector
  outcome1[missind] = miss0; 
  
  // --------------------------------
  // INITIALIZE BAYESIAN PARAMETER ESTIMATES
  par0 = runif(nc);        // create a vector of `nc` uniform random variables between 0 and 1
  parameter(0,_) = par0;   // set the first row of the parameter estimation to be this runif() variable

  // --------------------------------
  // CREATE INDICATOR MATRIX AND MISSING MATRIX
  // creates a matrix of indicator values, first n columns are is_week_i, last column is is_weekend
  cov = dummy(week, weekend); 

  // Creates a separate indicator value matrix for the subset of days that are/were missing data
  NumericMatrix misscov (nmiss, cov.ncol());
  for (int i=0;i < nmiss; ++i){
    mi = missind[i];
    misscov(i,_) = cov(mi,_);
  }
  
  // --------------------------------
  std::string name_base = "./tmp/w" + std::to_string(workerID);
  std::string old_name = name_base + "-" + "0" + ".txt";
  std::string new_name = old_name;
  std::ofstream file(old_name.c_str()); // Open the file
  file.close();

  // --------------------------------  
  ////////////////////////////////////////////////////////////////////
  /// The MCMC loop //////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //Rcpp::Rcout << "ITER# ";
  for (int i = 1; i < iter; ++i) { 
    
    // Print for every 500
    if(printProgress == 1) {
      if(i % 500 == 0) {
        new_name = name_base + "-" + std::to_string(i) + ".txt";
        //Rcpp::Rcout << i << "\t";
        if (std::rename(old_name.c_str(), new_name.c_str()) != 0) {
          // If std::rename returns a non-zero value, the rename operation failed
          Rcpp::stop("Failed to rename the file.");
        }
        old_name = new_name;
      }
    }

    // initialize the starting old parameters
    oldpar = parameter(i - 1, _);

    // --------------------------------
    // LOOP PART 1: The "Gibbs" Sampler, although this is not Gibbs, its metropolis-within-gibbs
    // This runs once for each parameter. Future iterations could look at randomizing this or making
    // some similarity beween 
    for (int k = 0; k < nc; ++k) {

      // for the variable that is i, ...., nw, 1 for each week and then 1 for is_weekend
      // THIS IS THE WHOLE REASON WE HAVE TO DO GIBBS?! otherwise we could just sample independently for each week
      if (k < nw + 1) {
        
        // make a copy of the old parameter values
        newpar    = clone(oldpar);

        // make a draw from a normal distribution with mean=oldpar[k] and sd=sigma
        newpar[k] = R::rnorm(oldpar[k], sigma);

        // compute the exp of the difference of the log-likelihoods? aka the ratio?
        ratio     = exp(logLikNB(outcome1, cov, dpind, newpar, maxdelay) - 
                        logLikNB(outcome1, cov, dpind, oldpar, maxdelay));

        // same as decision = ifelse(ratio > 1, 1, ratio)
        decision  = (ratio>1)?1:ratio;

        // ok so then this first gets a random uniform variable R::runif(0,1)
        // and if this is less than decision, make the parameter value = newpark[k],
        // but if its >= decision, keep the oldparameter value
        parameter(i,k) = (R::runif(0,1) < decision)?newpar[k]:oldpar[k];

        // updated oldpar[k]
        oldpar[k] = parameter(i, k);
      
      // for the extra 1 or 2 variables, which are ......
      } else {
        newpar    = clone(oldpar);
        newpar[k] = exp(R::rnorm(log(oldpar[k]), sigma));
        
        if (newpar[k] < 100) {
          ratio = (newpar[k]/oldpar[k])* 
            exp(logLikNB(outcome1, cov, dpind, newpar, maxdelay) - 
                logLikNB(outcome1, cov, dpind, oldpar, maxdelay));
        } else {
          ratio = 0;
        }
        
        decision = (ratio>1)?1:ratio;
        parameter(i,k) = (R::runif(0,1)<decision)?newpar[k]:oldpar[k];
        oldpar[k] = parameter(i,k);
      }
    }
    
    // --------------------------------
    // LOOP Part 2: Impute the missing delays
    // Rprintf(">> Part 2: Imputation of missing delays\n");
    // TODO: Why would you not just do this on the last iteration?

    // TYPE == 1 means there are two types of reporting delay sizes, aka r1 and r2
    if (type == 1) {
      par1      = oldpar[seq(0, nc - 3)];
      r1        = oldpar[nc - 2];
      double r2 = oldpar[nc - 1];
      mm        = get_mu_vec(misscov, par1);
      mean      = unique(mm);
      dv        = mean.size();
      //
      for (int k = 0 ; k < dv; ++k) {
        m = mean[k];
        LogicalVector select1 = (mm==m)&(dpind==1);
        LogicalVector select2 = (mm==m)&(dpind==2);
        int nm1 = sum(select1);
        int nm2 = sum(select2);

        if (nm1 > 0){
          prob = Rcpp::dnbinom_mu(dayseq, r1, m);
          p = prob / sum(prob);
          s = sample(dayseq, nm1, true, p); 
          IntegerVector ind = missind[select1];
          outcome1[ind] = s;
        }

        if (nm2 > 0){
          prob = Rcpp::dnbinom_mu(dayseq, r2, m);
          p = prob / sum(prob);
          s = sample(dayseq, nm2, true, p); 
          IntegerVector ind = missind[select2];
          outcome1[ind] = s;
        }
      }
    } else { // type != 1, is simpler, means there is just one size
      par1 = oldpar[seq(0, nc - 2)];
      r1   = oldpar[nc - 1];
      mm   = get_mu_vec(misscov, par1);
      mean = unique(mm);
      dv   = mean.size();
      for (int k = 0; k < dv; ++k){
        m = mean[k];
        LogicalVector select = (mm == m);
        nm = sum(select);
        prob = Rcpp::dnbinom_mu(dayseq, r1, m);
        p = prob / sum(prob);
        s = sample(dayseq, nm, true, p); 
        IntegerVector ind = missind[select];
        outcome1[ind] = s;
      }
    }
    
    // --------------------------------
    // LOOP Part 3: Back-calculation
    // seems like the upper limit of this is nd + maxdelay
    NumericVector backc (nd + maxdelay);     // initialize an empty vector for backc
    NumericVector backd = days - outcome1;   // initialize a vector for backd
    for (int j = 0; j <  (nd + maxdelay); ++j){
      backc[j] = sum(backd == (j - maxdelay + 1));
    }
    
    // --------------------------------
    // LOOP Part 4: Nowcasting
    NumericVector weights = prop(outcome1, backd, maxdelay, cday);
    NumericVector back1 = backc[seq(nd, nd + maxdelay - 1)];
    NumericVector check0 (back1.size());
    LogicalVector l = (back1 == 0);
    check0[l] = 1;
    NumericVector back2 = back1 + check0;  
    NumericVector trunc = mapply(back2, weights, rnb);
    NumericVector now   = back1 + trunc;
    NumericVector now1  = now - check0;
    LogicalVector check = (now1 < 0);
    now1[check] = 0;
    backc[seq(nd, nd + maxdelay - 1)] = now1; // set the rest of the backc vector = to now1
    back(i,_) = backc;
    rt(i,_) = getr(backc, si, size);
  } ///
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // Final print
  if(printProgress == 1) {
    //if(i % 500 == 0) {
      new_name = name_base + "-" + std::to_string(iter) + ".txt";
      //Rcpp::Rcout << i << "\t";
      if (std::rename(old_name.c_str(), new_name.c_str()) != 0) {
        // If std::rename returns a non-zero value, the rename operation failed
        Rcpp::stop("Failed to rename the file.");
      }
      old_name = new_name;
    //}
  }
  
  // --------------------------------
  // OUTPUT
  List output = List::create(Named("Back") = back, Named("R") = rt);
  return output;
}













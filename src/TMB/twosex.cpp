// CurrentAssessment is a TMB model that mimics the current Alaskan ADMB Sablefish stock assessment model
/*
 * Statistical, separable age-structured population model for sablefish
 * Alaska Fisheries Science Center, October 2008
 * D. Hanselman:dana.hanselman@noaa.gov
 * UPDATED (and commented)  by D. Goethel: daniel.goethel@noaa.gov   (10/15/20)
 * Translated from ADMB to TMB by C.Marsh: craig.marsh10@gmail.com  (11/10/22)
 * Tips
 * - TMB indicies start at 0 (i.e., like C++) where as ADMB starts at 1 (i.e., like R)
 * - Object labels should start with the transformation that is assumed for example natural log of F for longline should follow ln_F_ll
 *   and for logistic proportion something like logis_prop. This is to aid readability and keep syntax consistent
 *
 * - Key when reading object names
 *      - dom = domestic
 *      - ll = Longline
 *      - trwl = Trawl
 *      - srv = survey
 *      - lgth = Length
 *      - nmfs = Trawl survey
 */


/* General model description
 Single Area (combined across GOA, BS, AI)
 Split sexes, seperate weight at ages, age-length keys and unsexed
 All fishery catch include  catch and discards
 LL catch is adjusted (increased) for whale depredation (CPUE is not ajusted)
 DOM LL survey RPN/RPW is adjusted (increased) for whale depredation
 All fisheries and surveys are across entire model domain (GOA, BS, AI)
 */

#include <TMB.hpp>

// Zerofun functions see - for a discussion on this https://github.com/kaskr/adcomp/issues/7
template<class Type>
Type posfun(Type x, Type eps, Type &pen) {
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  Type xp = -(x/eps-1);
  return CppAD::CondExpGe(x,eps,x,
                          eps*(1/(1+xp+pow(xp,2)+pow(xp,3)+pow(xp,4)+pow(xp,5))));
}


// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// dinverse_gaussian
template<class Type>
Type dinverse_gaussian(Type x, Type mean, Type cv, int give_log=0){
  //return sqrt(lambda/(2*M_PI*pow(x,3))) * exp( -1.0 * lambda*pow(x-mean,2) / (2*pow(mean,2)*x) );
  Type sd = cv * mean;
  Type lambda = pow(mean,3.0) / pow(sd,2.0);
  Type logres = 0.5*(log(lambda) - 3.0*log(x) - log(2.0*M_PI)) - ( lambda*pow(x-mean,2.0) / (2.0*pow(mean,2.0)*x) );
  if(give_log) return logres; else return exp(logres);
}

// Zerofun functions
template <class Type>
Type ZeroFun(Type x, Type delta) {
  if (x >= delta)
    return x;

  return delta / (2.0 - (x / delta));
}
/*
 *  an index folding method to help get the folded index given dim1 and dim2
 *  @param ndx_1 (starts at 0 goes to (n_regions - 1))
 *  @param ndx2
 *  @param n_elements_dim1
 *  @return an index to look up the tagged partition which covers both these indicies
 */
int get_folded_ndx(int ndx_1, int ndx2, int n_elements_dim1) {
  return ndx2 * n_elements_dim1 + ndx_1;
}
/*
 * Updated the multinom fun because TMB's version cannot deal with p = 0
 * where as R's can
 */
template <class Type>
Type dmultinom_upd(vector<Type> x, vector<Type> p, int give_log=0)
{
  p /= p.sum(); // double check its nomalised
  Type logres = lgamma(x.sum() + Type(1));
  Type p_delta;
  Type delta = 0.00001;
  for(int i = 0; i < p.size(); ++i) {
    p_delta = ZeroFun(p(i), delta); // delta for the zerofun function to stop p = 0
    logres += x(i)*log(p_delta) - lgamma(x(i)+Type(1));
  }
  if(give_log)
    return logres;
  else
    return exp(logres);
}
/*
 *  Evaluate the density or log-likelihood for the dirichlet-multinomial distribution using linear-parameterisation from Thorson et.al 2017
 *  @param n input samples size
 *  @param theta over-dispersion parameter
 *  @param p_hat vector of fitted proportions
 *  @param p vector of observed proportions
 */
template <class Type>
Type ddirichletmulti(vector<Type> p_obs, vector<Type> p_hat, Type n, Type theta, int give_log=0) {
  Type s1 = 0.0;
  Type s2 = 0.0;
  Type loglik = 0.0;
  for(int ndx = 0; ndx < p_obs.size(); ndx++){
    s1 += lgamma(n * p_obs(ndx) + 1);
    s2 += lgamma(n * p_obs(ndx) + theta * n * p_hat(ndx)) - lgamma(theta * n * p_hat(ndx));
  }
  loglik = lgamma(n + 1) - s1 + lgamma(theta * n) - lgamma(n + theta * n) + s2;
  if(give_log)
    return loglik;
  else
    return exp(loglik);
}
/*
 *  an index folding method to help get the region and release year for tagged fish in the tagged partition
 *  @param region_ndx (starts at 0 goes to (n_regions - 1))
 *  @param release_year_ndx (0:n_years_to_retain_tagged_cohorts_for - 1). 0 is current year release, 1 is tag fished release last year, 2 is tag release 2 years before current
 *  @return an index to look up the tagged partition which covers both these indicies
 */
int get_tag_release_event_ndx(int region_ndx, int release_event_year_ndx, int n_regions) {
  return release_event_year_ndx * n_regions + region_ndx;
}
/*
 *  square of a variable
 */
template <class Type>
Type square(Type x){return x*x;}
/*
 *  square of a variable overloaded for vector input
 */
template <class Type>
vector<Type> square(vector<Type> x) {
  return x*x;
}

// Return the natural logarithm of one plus the specified value.
template <class Type>
Type log1p(Type x) {
  return log(1 + x);
}
//Calculates the log of 1 plus the exponential of the specified value without overflow.
template <class Type>
Type log1p_exp(Type a) {
  if (a > 0.0) {
    return a + log1p(exp(-a));
  }
  return log1p(exp(a));
}

// transform Y -Inf-Inf -> X bound lb - ub
template <class Type>
Type invlogit_general(Type& Y, Type& lb, Type& ub) {
  return(lb + (ub - lb) * (1 / (1 + exp(-Y))));
}
// transform Y -Inf-Inf -> X bound lb - ub
template <class Type>
vector<Type> invlogit_general(vector<Type>& Y, Type& lb, Type& ub) {
  vector<Type> X(Y.size());
  for(int i = 0; i < X.size(); ++i) {
    X(i) = lb + (ub - lb) * (1 / (1 + exp(-Y(i))));
  }
  return(X);
}

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}


/*
 * Geometric mean
 */
template <class Type>
Type geo_mean(vector<Type>& x){
  return exp((log(x).sum())/x.size());
}

template <class Type>
void get_covar(vector<Type> theta, matrix<Type>& covar, int n, int type) {
  // Source
  // https://github.com/glmmTMB/glmmTMB/blob/9cb23d77afc042be5bced3713c40529023d9ef7e/glmmTMB/src/glmmTMB.cpp#L361
  if (type == 1) {
    // case: diag_covstruct
    //std::cerr << "diag_covstruct" << std::endl;
    Type sd = exp(theta(0));
    for(int i = 0; i < n; ++i)
      covar(i,i) = sd * sd;
  } else if (type == 2) {
    // case: us_covstruct
    //std::cerr << "us_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    // for covariance
    //density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        covar(i,j) = nldens.cov()(i,j) * sd(i) * sd(j);
      }
    }
  } else if (type == 3){
    // case: cs_covstruct Compound Symmetry: Heterogenous.
    //std::cerr << "cs_covstruct" << std::endl;

    vector<Type> sd = exp(theta.head(n));
    Type corr_transf = theta(n);
    Type a = Type(1) / (n - Type(1)); // Correcting for violations of sphericity
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : rho * sd(i) * sd(j)) ;

  } else if (type == 4){
    // case: toep_covstruct Toeplitz: Heterogenous.
    //std::cerr << "toep_covstruct" << std::endl;
    vector<Type> sd = exp(theta.head(n));
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        covar(i,j) = (i==j ? Type(1) * sd(i) * sd(i) : parms( (i > j ? i-j : j-i) - 1 ) * sd(i) * sd(j));
  } else if (type == 5){
    //std::cerr << "ar1_covstruct" << std::endl;
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    vector<Type> sd(n);
    sd.fill(exp(theta(0)));
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    for(int i=0; i < n; i++){
      for(int j=0; j < n; j++){
        covar(i,j) = pow(phi, abs(i-j)) * sd(i) * sd(i);
      }
    }
  }

  return;
}



/*
 * rmultinomm - for simulate call
 */
template <class Type>
vector<Type> rmultinom(vector<Type> prob, Type N) {
  vector<Type> sim_X(prob.size());
  sim_X.setZero();
  // Now simulate using the uniform random variable
  Type rng_uniform;
  Type cumulative_expect;
  while(N > 0) {
    rng_uniform = runif(Type(0),Type(1));
    //std::cout << rng_uniform << " ";
    cumulative_expect = 0.0;
    for (unsigned i = 0; i < prob.size(); ++i) {
      cumulative_expect += prob[i];
      if (cumulative_expect >= rng_uniform) {
        sim_X[i] += 1.0;
        break;
      }
      //sim_X[prob.size() - 1] += 1.0;
    }
    N -= 1;
  }
  //std::cout << "\n";
  return(sim_X);
}

/*
 *  Simulate a single draw from the dirichlet-multinomial distribution using linear-parameterisation from Thorson et.al 2017
 *  @param n input samples size
 *  @param theta over-dispersion parameter
 *  @param fitted_props vector of fitted proportions
 */
template <class Type>
vector<Type> rdirichletmulti(vector<Type> fitted_props, Type& n, Type& theta) {
  vector<Type> dirichlet_draw(fitted_props.size());
  //Type n_eff = (1 + theta * n) / (1 + theta); // calculate the effective sample size
  for(int ndx = 0; ndx < fitted_props.size(); ndx++)
    dirichlet_draw(ndx) = rgamma(fitted_props(ndx) * theta * n, (Type)1.0);// shape, rate = 1.0

  Type dirich_total = dirichlet_draw.sum();
  dirichlet_draw /= dirich_total;
  return(rmultinom(dirichlet_draw, n));
  // Should we simulate the final multinomial sample using n or n_eff?
  // I did some simulations contained in my Thesis/Demonstration/Distribution code which
  // shows based on Thorson 2017 paper from his Pearson residual calculations that we should simualte using n not n_eff
}
/*
 * inverse centred log transform
 * Up to a constant so standardise
 */
template <class Type>
vector<Type> inv_crl(vector<Type>& y){
  return exp(y) / (exp(y)).sum();
}

// Beverton-Holt SR relationship function
template<class Type>
Type BevertonHolt(Type SSB, Type B0, Type h) {
  Type ssb_ratio = SSB / B0;
  Type part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
  return (ssb_ratio / part_2);
}

// Beverton-Holt SR relationship function without equilibrium assumptions
template<class Type>
Type BevertonHoltNoEquil(Type a, Type b, Type SSB) {
  Type Rt = (a + SSB) / (SSB * b);
  return Rt;
}

/*
 * Simplex functionality
 */
// simplex_jacobian - calculate the jacobian for estiamting the logistic simplex rather than the simplex
template <class Type>
Type simplex_jacobian(vector<Type> zk, vector<Type> X_orig) {
  std::cout << "Simplex Jacobian\n";
  Type Jacobian = 1.0;
  for (int k = 0; k < (zk.size() - 1); ++k) { // because I store zk in a matrix with xk in simplex_restore() the zk vector has one more null/empty element than it should have so we don't iterate over this.
    if (k == 0) {
      Jacobian *= (zk[k] * (1 - zk[k]));
    } else {
      Jacobian *= (zk[k] * (1 - zk[k])) * (1 - sum(vector<Type>(X_orig.segment(0,k))));
    }
  }
  return Jacobian;
}

// This class takes the transformed (logistic simplex) values, and their original unit vector and gives the
// recent update value.
// returns a matrix row 1 (C++ call: row(0)) is the zk - needed for the jacobian
// teh second row (C++ call: row(1)) is xk the untransformed values
template <class Type>
matrix<Type> Simplex_restore(vector<Type> logit_unit_vec, vector<Type> last_unit_vec) {
  std::cout << "Simplex_restore\n";
  int n_pars = last_unit_vec.size();
  matrix<Type> matrix_return(2,n_pars);
  Type total = 0.0;
  for (int k = 0; k < (logit_unit_vec.size()); ++k) {
    matrix_return(0,k) = invlogit(logit_unit_vec[k] + log(Type(1) / ((Type(logit_unit_vec.size()) + Type(1)) - Type(k + 1))));

    if (k == 0)
      matrix_return(1,k) = matrix_return(0,k);
    else
      matrix_return(1,k) = (Type(1) - sum(vector<Type>(last_unit_vec.segment(0,k)))) * matrix_return(0,k);
    total += matrix_return(1,k);
    //std::cerr << "k = " << k + 1 << " zk = " << matrix_return(0,k) << " xk = " << matrix_return(1,k) << std::endl; // Output warning
  }
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;
  matrix_return(1,n_pars - 1) = 1 - total;
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;

  return matrix_return;
}

/*
 * centred log transform
 */
template <class Type>
vector<Type> crl(vector<Type>& x) {
  return log(x / geo_mean(x));
}

// Calculate spawner per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param paa vector of proportions at age
 * @param propZ_ssb scalar interpolating SSB
 * @param ages vector of ages
 * @return SPR
 */

template<class Type>
Type get_SPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& paa, Type propZ_ssb, vector<Type>& ages) {
  Type Na = 1;
  vector<Type> Za = natural_mortality + sel * F;
  vector<Type> Sa = exp(-Za);
  Type SPR = Na * exp(-Za(0) * propZ_ssb) * waa(0) * paa(0);
  for(unsigned age_iter = 1; age_iter < (ages.size() * 4); ++age_iter) {
    if(age_iter >= ages.size() ) {
      SPR += Na * exp(-Za(ages.size() - 1) * propZ_ssb) * waa(ages.size() - 1) * paa(ages.size() - 1);
      Na *= Sa(ages.size() - 1);
    } else {
      SPR += Na * exp(-Za(age_iter) * propZ_ssb) * waa(age_iter) * paa(age_iter);
      Na *= Sa(age_iter);
    }
  }
  return SPR;
}

// Calculate Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_YPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type Na = 1;
  Type YPR = (sel(0) * F) / (natural_mortality + sel(0) * F) * (Type(1) - exp(-(natural_mortality + sel(0) * F))) * Na * waa(0);
  for(unsigned age_iter = 1; age_iter < ages.size() * 4; ++age_iter) {
    if(age_iter >= ages.size() ) {
      Na *= exp(-(natural_mortality + sel(ages.size() - 1) * F));
      YPR += (sel(ages.size() - 1) * F) / (natural_mortality + sel(ages.size() - 1) * F) * (Type(1) - exp(-(natural_mortality + sel(ages.size() - 1) * F))) * Na * waa(ages.size() - 1);
    } else {
      Na *= exp(-(natural_mortality + sel(age_iter) * F));
      YPR += (sel(age_iter) * F) / (natural_mortality + sel(age_iter) * F) * (Type(1) - exp(-(natural_mortality + sel(age_iter) * F))) * Na * waa(age_iter);

    }
  }
  return YPR;
}
// Calculate derivative for Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return derivative YPR
 */
template<class Type>
Type get_dYPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type ln_F = log(F);
  Type h = 0.001;
  Type v = -get_YPR(exp(ln_F + 2.0 * h), sel, natural_mortality, waa, ages) + 8.0 * get_YPR(exp(ln_F + h), sel, natural_mortality, waa, ages) - 8.0 * get_YPR(exp(ln_F - h), sel, natural_mortality, waa, ages) + get_YPR(exp(ln_F - 2.0 * h), sel, natural_mortality, waa, ages);
  Type g = v / (12.0 * h);
  return(g / F);
}
// Calculate Equilbrium SSB for a given F
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_SSBe(Type& F, vector<Type>& N_init, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages, vector<Type>& paa, Type propZ_ssb, int n_runs, int maxAgePlusGroup) {

  vector<Type> Fa = F * sel;
  vector<Type> Za = Fa + natural_mortality;
  vector<Type> Sa = exp(-Za);
  array<Type> N_equilibrium(ages.size(), n_runs);
  N_equilibrium.col(0) = N_init;

  // Run annual cycle
  for(int year_ndx = 1; year_ndx < n_runs; ++year_ndx) {
    // Initial recruitment
    N_equilibrium(0, year_ndx) = N_init(0);
    // Ageing + Z
    for(int age_ndx = 1; age_ndx < ages.size() ; ++age_ndx)
      N_equilibrium(age_ndx, year_ndx) = N_equilibrium(age_ndx - 1, year_ndx - 1) * Sa(age_ndx - 1);
    if(maxAgePlusGroup == 1) {
      N_equilibrium(ages.size() - 1, year_ndx) = N_equilibrium(ages.size()  - 2, year_ndx - 1) * Sa(ages.size()  - 2) +
        N_equilibrium(ages.size()  - 1, year_ndx - 1) * Sa(ages.size()  - 1);
    }
  }

  // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
  Type ssbe = 0;
  for(int age_ndx = 0; age_ndx < ages.size(); ++age_ndx)
    ssbe += N_equilibrium(age_ndx, n_runs - 1) * exp(-Za(age_ndx) * propZ_ssb) * paa(age_ndx) * waa(age_ndx);

  return(ssbe);
}

/*
*  SetupSelectivities.h
*
*  This script has functions that are used to map selectivity parameters and control variables
*  to develop an array (age x year) for each fishery and survey.
*
*/


// logistic ogive function parametersised with a_50 and a_to95
template <class Type>
vector<Type> logistic_ogive(vector<Type>& ages, Type& sel_50, Type& sel_95) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + pow(Type(19.0), (sel_50 - ages[age]) / sel_95));
  }
  return logis;
}
// Alternative logistic ogive function parametersised with a_50 and delta
template <class Type>
vector<Type> logistic_alt_ogive(vector<Type>& ages, Type& sel_50, Type& delta) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + exp(-delta * (ages[age] - sel_50)));
  }
  return logis;
}
// Double normal Punt et. al 1996 gamma parameterization
template <class Type>
vector<Type> double_normal_ogive(vector<Type>& ages, Type& sel_50, Type& delta) {
  //std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> d_norm(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    d_norm[age] = (pow(ages[age] / sel_50, sel_50 / (0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50)))*exp((sel_50-ages[age])/(0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50))));
  }
  return d_norm;
}

/*
 *  BuildSelectivity is the function that will build a selectivity array that can have different time-blocks and a different ogive per block given parameters and control variables
 *  @param sel_params dim: parameters (untransformed) where n_time_blocks is the number of time-blocks n_time_blocks x n_params
 *  @param sel_type vector of length = n_time_blocks. THis will define how many columns are in sel_params
 *  0 = Specifies whether it is a logistic (= 0)
 *  1 = Specifies the gamma formuation which mimics the double normal with two parameters
 *  2 = power function not scaled to have a max = 1 like in the assessment
 *  3 = Alternative logistic parameterisation
 *  4 = exponential decay selectivity
 *  5 = double normal with 3 parameters
 *  @param ages vector of ages.
 *  @param sel_array is a prepopulated array that will be fulled with the correct . Specifies which time-block selectivity to which year.
 *  @param normalise if true then it will scale the selectivity to have max value = 1
 */
template <class Type>
void BuildSelectivity(array<Type> sel_params, vector<int> sel_type, vector<Type> ages, array<Type>& sel_array, bool normalise) {
  int n_time_blocks = sel_params.dim(0);
  int n_ages = ages.size();

  for(int i = 0; i < n_time_blocks; ++i) {
    if(sel_type(i) == 0) {
      // logistic, expect two parameters
      sel_array.col(i) = logistic_alt_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 1) {
      // Double normal
      sel_array.col(i) = double_normal_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 2) {
      // Power function used for the GOA trawl selectivity
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) = (1.0 / pow(ages[age_ndx],sel_params(i, 0)));
      // scale by max TODO: DELETE THIS after validated it shouldn't be here
      // A max call can cause AD issues. Or cause a split in the AD graph which I don't like at all.... consider using the expoential decay i.e., sel_type == 4 instead
    }  else if (sel_type(i) == 3) {
      // Alternative logistic parameterisation
      sel_array.col(i) = logistic_ogive(ages, sel_params(i, 0), sel_params(i, 1));
    } else if (sel_type(i) == 4) {
      // Exponential decay function
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) = exp(-1.0 * ages[age_ndx] * sel_params(i, 0));
    } else if (sel_type(i) == 5) {
      // Double normal with 3 parameters
      // param order mu, sigma_r, sigma_l
      Type delta = 5.0;
      Type tmp = log(0.5);
      Type stmp = 0.0;
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
      stmp = 1.0 / (1.0 + exp(-delta * (ages[age_ndx] - sel_params(i, 0))));
      sel_array(age_ndx, i) = stmp * exp(tmp * square((ages[age_ndx] - sel_params(i, 0)) / sel_params(i, 1))) + (1.0 - stmp) *
        exp(tmp * square((ages[age_ndx] - sel_params(i, 0)) / sel_params(i, 2)));
      }
    }
    // scale by max TODO: DELETE THIS after validated it shouldn't be here
    // A max call can cause AD issues.
    if(normalise) {
      Type max_sel = max(vector<Type>(sel_array.col(i)));
      for(int age_ndx = 0; age_ndx < n_ages; ++age_ndx)
        sel_array(age_ndx, i) /= max_sel;
    }
  }
}




template<class Type>
Type objective_function<Type>::operator() () {
   using namespace density;

  // Input parameters
  // model dimensions
  DATA_VECTOR(ages);                                // assumes min(ages) >= 1, also assumes the last age is a plus group
  DATA_VECTOR(years);                               // annual years
  DATA_VECTOR(length_bins);                         // Length bins, the last length bin value is the minimum for a length plus group
  DATA_INTEGER(n_projections_years);                // number of years to project the model beyond max(years)

  int n_years = years.size();
  int n_projyears = n_years + n_projections_years;
  int n_ages = ages.size();
  int n_length_bins = length_bins.size();


  // Biology parameters
  DATA_INTEGER(n_init_rec_devs);                    // Number of initial recruitment devs parameters "init_ln_rec_dev" Note: should cannot be greater than n_ages - 2 (we don't apply it to first age or plus group)

  // this will effect the expected size of the parameter 'ln_rec_dev', if global_rec_devs = 1. then ln_rec_dev.size() = n_years + n_ages + 1 else (n_years + n_ages + 1). with the first (n_years + n_ages + 1) corresponding to region 1 and so in block
  DATA_ARRAY(M);                                    // Natural Mortality: dim = n_ages x n_projyears
  DATA_ARRAY(maturity);                             // Proportion ages mature: dim = n_ages x n_projyears
  DATA_ARRAY(male_mean_weight_by_age);              // male_mean_weight_by_age (tonnes): dim = n_ages x n_projyears
  DATA_ARRAY(female_mean_weight_by_age);            // female_mean_weight_by_age (tonnes): dim = n_ages x n_projyears

  DATA_ARRAY(male_age_length_transition);           // Proportion at among length bins for each age for male: dim = n_ages x n_lengths x n_years
  DATA_ARRAY(female_age_length_transition);         // Proportion at among length bins for each age for female: dim = n_ages x n_lengths x n_years

  // This needs to change in future assessments, this should be dealt with by selectivities
  DATA_VECTOR(proportion_male);                     // YUCK!!!! proportion of males in RPN. A value for each year.
  DATA_VECTOR(proportion_male2);                    // YUCK!!!! proportion of males in RPN. A value for each year.

  DATA_SCALAR(sigma_R);                             // standard deviation for recruitment;
  DATA_VECTOR(spawning_time_proportion);            // proportion of time within a year that spawning occurs needed for each year length = n_projyears, bound between 0 and 1

  // Fishing stuff
  DATA_INTEGER(catch_likelihood);                   // 0 is the ADMB formulation, 1 == normal with catch_sd
  DATA_SCALAR(catch_sd);                            // only used if catch_likelihood == 1
  DATA_SCALAR(prop_F_hist);                         // Proportion of ll_F_avg that is applied during initialization
  DATA_VECTOR(ll_fishery_catch);                    // Observed catch for Longline fishery. length = n_years
  DATA_VECTOR(trwl_fishery_catch);                  // Observed catch for Trawl fishery. length = n_years
  DATA_SCALAR(loglik_wgt_ll_catch);                 // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_trwl_catch);               // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_Fs);                       // Log-likelihood multiplier (Craig's not a fan of these)

  // Selectivity indicator switches
  DATA_IVECTOR(ll_sel_type);                        // Selectivity type for each row of ln_ll_sel_m_pars and ln_ll_sel_f_pars
  DATA_IVECTOR(ll_sel_by_year_indicator);           // Selectivity time-block to apply in each model year
  DATA_IVECTOR(trwl_sel_type);                      // Selectivity type for each row of ln_trwl_sel_m_pars and ln_ll_sel_f_pars
  DATA_IVECTOR(trwl_sel_by_year_indicator);         // Selectivity time-block to apply in each model year

  // Survey stuff
  DATA_IVECTOR(srv_dom_ll_sel_type);                // Selectivity type for each row of ln_ll_sel_m_pars and ln_ll_sel_f_pars
  DATA_IVECTOR(srv_dom_ll_sel_by_year_indicator);   // Selectivity time-block to apply in each model year
  DATA_IVECTOR(srv_dom_ll_q_by_year_indicator);     // Catchability time-block to apply when deriving model predictions each year

  // NMFS bottom trawl survey
  DATA_IVECTOR(srv_nmfs_trwl_sel_type);                // Selectivity type for each row of ln_ll_sel_m_pars and ln_ll_sel_f_pars
  DATA_IVECTOR(srv_nmfs_trwl_sel_by_year_indicator);   // Selectivity time-block to apply in each model year
  DATA_IVECTOR(srv_nmfs_trwl_q_by_year_indicator);     // Catchability time-block to apply when deriving model predictions each year

  // Fishery longline CPUE (fishery-dependent)
  DATA_IVECTOR(ll_cpue_q_by_year_indicator);       // Catchability time-block to apply when deriving model predictions each year

  // Observational stuff
  DATA_MATRIX(ageing_error_matrix);                 // Ageing error/missclassification matrix n_ages x n_ages
  DATA_VECTOR(log_likelihood_weights);              // YUCK!!! as soon as you can delete this and 'weight' likelihoods using the respective dispersion parameters
  // 0 = ll_age, 1 = trwl lgth male, 2 = trwl lgth female, 3 = srv_dom_ll_bio

  // Longline catch at age (sex Aggregated)
  DATA_IVECTOR(ll_catchatage_indicator);            // length(ll_catchatage_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_ll_catchatage = sum(ll_catchatage_indicator);
  DATA_ARRAY(obs_ll_catchatage);                   // Longline fishery composition observations dim = n_ages x ll_catchatage_comp
  DATA_ARRAY_INDICATOR(keep_ll_catchatage_comp, obs_ll_catchatage); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(ll_catchatage_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(ll_catchatage_comp_likelihood);             // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  array<Type> pred_ll_catchatage(n_ages, n_ll_catchatage); // Sex disaggregated predicted catch at age
  DATA_SCALAR(loglik_wgt_ll_catchatage);               // Log-likelihood multiplier (Craig's not a fan of these)

  // Longline catch at length (sex disaggregated)
  DATA_IVECTOR(ll_catchatlgth_indicator);            // length(ll_catchatlgth_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_ll_catchatlgth = sum(ll_catchatlgth_indicator);
  DATA_ARRAY(obs_ll_catchatlgth_m);                   // Trawl fishery composition observationsn_years dim = n_length_bins x ll_catchatlgth_comp
  DATA_ARRAY_INDICATOR(keep_ll_catchatlgth_m_comp, obs_ll_catchatlgth_m);
  DATA_ARRAY(obs_ll_catchatlgth_f);                   // Trawl fishery composition observations dim = n_length_bins x ll_catchatlgth_comp
  DATA_ARRAY_INDICATOR(keep_ll_catchatlgth_f_comp, obs_ll_catchatlgth_f);
  DATA_INTEGER(ll_catchatlgth_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(ll_catchatlgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  DATA_SCALAR(loglik_wgt_ll_catchatlgth_m);                   // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_ll_catchatlgth_f);                   // Log-likelihood multiplier (Craig's not a fan of these)

  array<Type> pred_ll_catchatlgth_m(n_length_bins, n_ll_catchatlgth); // Sex disaggregated predicted catch at age
  array<Type> pred_ll_catchatlgth_f(n_length_bins, n_ll_catchatlgth); // Sex disaggregated predicted catch at age

  // Trawl catch at length (sex disaggregated)
  DATA_IVECTOR(trwl_catchatlgth_indicator);            // length(trwl_catchatlgth_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_trwl_catchatlgth = sum(trwl_catchatlgth_indicator);
  DATA_ARRAY(obs_trwl_catchatlgth_m);                   // Trawl fishery composition observations dim = n_length_bins x ll_catchatlgth_comp
  DATA_ARRAY_INDICATOR(keep_trwl_catchatlgth_m_comp, obs_trwl_catchatlgth_m);
  DATA_ARRAY(obs_trwl_catchatlgth_f);                   // Trawl fishery composition observations dim = n_length_bins x ll_catchatlgth_comp
  DATA_ARRAY_INDICATOR(keep_trwl_catchatlgth_f_comp, obs_trwl_catchatlgth_f);
  DATA_INTEGER(trwl_catchatlgth_covar_structure);             // placeholder - 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(trwl_catchatlgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  DATA_SCALAR(loglik_wgt_trwl_catchatlgth_m);                   // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_trwl_catchatlgth_f);                   // Log-likelihood multiplier (Craig's not a fan of these)

  array<Type> pred_trwl_catchatlgth_m(n_length_bins, n_trwl_catchatlgth); // Sex disaggregated predicted catch at age
  array<Type> pred_trwl_catchatlgth_f(n_length_bins, n_trwl_catchatlgth); // Sex disaggregated predicted catch at age

  // Survey biomass for the Domestic Longline survey
  DATA_IVECTOR(srv_dom_ll_bio_indicator);               // length(srv_dom_ll_bio_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_dom_ll_bio = sum(srv_dom_ll_bio_indicator);
  DATA_VECTOR(obs_dom_ll_bio);                          // Survey Domestic longline biomass obs. length =  n_dom_ll_bio
  DATA_VECTOR(se_dom_ll_bio);                           // Survey Domestic longline biomass standard errors. length =  n_dom_ll_bio
  DATA_VECTOR_INDICATOR(keep_obs_dom_ll_bio, obs_dom_ll_bio);
  DATA_INTEGER(srv_dom_ll_bio_likelihood);                  // 0 = ADMB, 1 = lnorm
  vector<Type> pred_dom_ll_bio(n_dom_ll_bio);           // Sex aggregated predicted
  pred_dom_ll_bio.setZero();                            // initialise vector to be filled with zeros

  // Survey NMFS GOA trawl survey
  DATA_IVECTOR(srv_nmfs_trwl_bio_indicator);               // length(srv_nmfs_trwl_bio_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_nmfs_trwl_bio = sum(srv_nmfs_trwl_bio_indicator);
  DATA_VECTOR(obs_nmfs_trwl_bio);                          // Survey NMFS GOA trawl biomass obs. length =  n_nmfs_trwl_bio
  DATA_VECTOR(se_nmfs_trwl_bio);                           // Survey NMFS GOA trawl biomass standard errors. length =  n_nmfs_trwl_bio
  DATA_VECTOR_INDICATOR(keep_obs_nmfs_trwl_bio, obs_nmfs_trwl_bio);
  DATA_INTEGER(srv_nmfs_trwl_bio_likelihood);                  // 0 = ADMB, 1 = lnorm
  vector<Type> pred_nmfs_trwl_bio(n_nmfs_trwl_bio);        // Sex aggregated predicted
  pred_nmfs_trwl_bio.setZero();                           // initialise vector to be filled with zeros

  // Longline CPUE
  DATA_IVECTOR(ll_cpue_indicator);                      // length(ll_cpue_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_ll_cpue = sum(ll_cpue_indicator);
  DATA_VECTOR(obs_ll_cpue);                              // Fishery longline CPUE obs. length =  n_ll_cpue
  DATA_VECTOR(se_ll_cpue);                               // Fishery longline CPUE standard errors. length =  n_ll_cpue
  DATA_VECTOR_INDICATOR(keep_obs_ll_cpue, obs_ll_cpue);
  DATA_INTEGER(ll_cpue_likelihood);                      // 0 = ADMB, 1 = lnorm
  vector<Type> pred_ll_cpue(n_ll_cpue);                  // Sex aggregated predicted
  pred_ll_cpue.setZero();                               // initialise vector to be filled with zeros

  // Survey age for the Domestic Longline
  DATA_IVECTOR(srv_dom_ll_age_indicator);               // length(srv_dom_ll_age_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_srv_dom_ll_age = sum(srv_dom_ll_age_indicator);
  DATA_ARRAY(obs_srv_dom_ll_age);                          // Survey Domestic longline age obs. length =  n_srv_dom_ll_age
  DATA_ARRAY_INDICATOR(keep_obs_srv_dom_ll_age, obs_srv_dom_ll_age);
  DATA_INTEGER(srv_dom_ll_age_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(srv_dom_ll_age_comp_likelihood);            // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  array<Type> pred_srv_dom_ll_age(n_ages, n_srv_dom_ll_age); // Sex aggregated predicted catch at age
  DATA_SCALAR(loglik_wgt_srv_dom_ll_age);                   // Log-likelihood multiplier (Craig's not a fan of these)

  // Survey length Comp from the Domestic Longline survey
  DATA_IVECTOR(srv_dom_ll_lgth_indicator);               // length(srv_dom_ll_lgth_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_srv_dom_ll_lgth = sum(srv_dom_ll_lgth_indicator);
  DATA_ARRAY(obs_srv_dom_ll_lgth_m);                          // Survey Domestic longline length obs. length =  n_srv_dom_ll_lgth
  DATA_ARRAY(obs_srv_dom_ll_lgth_f);                          // Survey Domestic longline length obs. length =  n_srv_dom_ll_lgth
  DATA_ARRAY_INDICATOR(keep_obs_srv_dom_ll_lgth_m, obs_srv_dom_ll_lgth_m);
  DATA_ARRAY_INDICATOR(keep_obs_srv_dom_ll_lgth_f, obs_srv_dom_ll_lgth_f);
  DATA_INTEGER(srv_dom_ll_lgth_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(srv_dom_ll_lgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's multinomial. //0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  array<Type> pred_srv_dom_ll_lgth_m(n_length_bins, n_srv_dom_ll_lgth); // Sex disaggregated predicted catch at age
  array<Type> pred_srv_dom_ll_lgth_f(n_length_bins, n_srv_dom_ll_lgth); // Sex disaggregated predicted catch at age
  DATA_SCALAR(loglik_wgt_srv_dom_ll_lgth_m);                   // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_srv_dom_ll_lgth_f);                   // Log-likelihood multiplier (Craig's not a fan of these)

  // NMFS bottom trawl age frequency
  DATA_IVECTOR(srv_nmfs_trwl_age_indicator);               // length(srv_nmfs_trwl_age_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_srv_nmfs_trwl_age = sum(srv_nmfs_trwl_age_indicator);
  DATA_ARRAY(obs_srv_nmfs_trwl_age);                          // Survey Domestic longline age obs. length =  n_srv_nmfs_trwl_age
  DATA_ARRAY_INDICATOR(keep_obs_srv_nmfs_trwl_age, obs_srv_nmfs_trwl_age);
  DATA_INTEGER(srv_nmfs_trwl_age_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(srv_nmfs_trwl_age_comp_likelihood);            // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  array<Type> pred_srv_nmfs_trwl_age(n_ages, n_srv_nmfs_trwl_age); // Sex aggregated predicted catch at age
  DATA_SCALAR(loglik_wgt_srv_nmfs_trwl_age);                   // Log-likelihood multiplier (Craig's not a fan of these)

  // NMFS bottom trawl length Comp
  DATA_IVECTOR(srv_nmfs_trwl_lgth_indicator);               // length(srv_nmfs_trwl_lgth_indicator) = n_years.  1 = calculate catch at age in this year, 0 = don't calculate catch at age
  int n_srv_nmfs_trwl_lgth = sum(srv_nmfs_trwl_lgth_indicator);
  DATA_ARRAY(obs_srv_nmfs_trwl_lgth_m);                          // Survey Domestic longline length obs. length =  n_srv_nmfs_trwl_lgth
  DATA_ARRAY(obs_srv_nmfs_trwl_lgth_f);                          // Survey Domestic longline length obs. length =  n_srv_nmfs_trwl_lgth
  DATA_ARRAY_INDICATOR(keep_obs_srv_nmfs_trwl_lgth_m, obs_srv_nmfs_trwl_lgth_m);
  DATA_ARRAY_INDICATOR(keep_obs_srv_nmfs_trwl_lgth_f, obs_srv_nmfs_trwl_lgth_f);
  DATA_INTEGER(srv_nmfs_trwl_lgth_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(srv_nmfs_trwl_lgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's Multinomial, 2 = dirichlet-multinomial (not implemented)
  array<Type> pred_srv_nmfs_trwl_lgth_m(n_length_bins, n_srv_nmfs_trwl_lgth); // Sex disaggregated predicted catch at age
  array<Type> pred_srv_nmfs_trwl_lgth_f(n_length_bins, n_srv_nmfs_trwl_lgth); // Sex disaggregated predicted catch at age
  DATA_SCALAR(loglik_wgt_srv_nmfs_trwl_lgth_m);                   // Log-likelihood multiplier (Craig's not a fan of these)
  DATA_SCALAR(loglik_wgt_srv_nmfs_trwl_lgth_f);                   // Log-likelihood multiplier (Craig's not a fan of these)

  // Estimable parameters
  PARAMETER(ln_mean_rec);                           // Unfish equil recruitment (logged) (estimated)
  PARAMETER_VECTOR(ln_rec_dev);                     // Recruitment deviations they include years before the asssessment starts to final year: length = n_years
  PARAMETER_VECTOR(ln_init_rec_dev);                // Recruitment deviations to apply during initialization they include years before the assessment starts: length = n_init_rec_devs

  // will be moved to the parameter section later
  PARAMETER_ARRAY(ln_ll_sel_pars);                       // log selectivity parameters for longline male, dim = time-blocks,  max(sel parameters), n_sex
  PARAMETER_ARRAY(ln_trwl_sel_pars);                     // log selectivity parameters for Trawl male, dim = time-blocks,  max(sel parameters), n_sex

  PARAMETER(ln_ll_F_avg);                                // log average longline Fishing mortality
  PARAMETER_VECTOR(ln_ll_F_devs);                        // Annual fishing mortality deviation
  PARAMETER(ln_trwl_F_avg);                              // log average trawl Fishing mortality
  PARAMETER_VECTOR(ln_trwl_F_devs);                      // Annual fishing mortality deviations

  PARAMETER_VECTOR(ln_ll_cpue_q);                        // log catchabilities parameters for srv_dom_ll observation

  PARAMETER_VECTOR(ln_srv_nmfs_trwl_q);                   // log catchabilities parameters for srv_dom_ll observation
  PARAMETER_ARRAY(ln_srv_nmfs_trwl_sel_pars);             // log selectivity parameters for NMFS bottom trawl survey male, dim = time-blocks,  max(sel parameters), n_sex
  PARAMETER_VECTOR(ln_srv_dom_ll_q);                      // log catchabilities parameters for srv_dom_ll observation
  PARAMETER_ARRAY(ln_srv_dom_ll_sel_pars);                // log selectivity parameters for domestic longline survey male, dim = time-blocks,  max(sel parameters), n_sex


  // Note: about ln_init_rec_dev - it is the opposite order to how ADMB model is formulated. I am sorry but it was easier to code.
  //       the first init_rec_dev corresponds to recruitment for age class 2 which would have arrived in styr - 1, and so on.
  // Untransform parameters
  Type mean_rec = exp(ln_mean_rec);
  vector<Type> rec_dev = exp(ln_rec_dev);
  Type F_hist = exp(ln_ll_F_avg);
  Type init_F_hist = F_hist * prop_F_hist;
  array<Type> ll_sel_pars(ln_ll_sel_pars.dim);
  ll_sel_pars = exp(ln_ll_sel_pars);
  array<Type> trwl_sel_pars(ln_trwl_sel_pars.dim);
  trwl_sel_pars = exp(ln_trwl_sel_pars);
  array<Type> srv_dom_ll_sel_pars(ln_srv_dom_ll_sel_pars.dim);
  srv_dom_ll_sel_pars = exp(ln_srv_dom_ll_sel_pars);
  vector<Type> srv_dom_ll_q = exp(ln_srv_dom_ll_q);
  vector<Type> ll_cpue_q = exp(ln_ll_cpue_q);
  array<Type> srv_nmfs_trwl_sel_pars(ln_srv_nmfs_trwl_sel_pars.dim);
  srv_nmfs_trwl_sel_pars = exp(ln_srv_nmfs_trwl_sel_pars);
  vector<Type> srv_nmfs_trwl_q = exp(ln_srv_nmfs_trwl_q);


  // Initialise consistently used variables throughout the code
  int model_type = 0; //"Assessment";
  int year_ndx;
  int age_ndx;
  int len_ndx;
  int min_age = 0;
  int n_regions = 1;
  while(min_age < ages(0)){
    min_age++;
  }
  Type m_plus_group = 0.0;
  Type f_plus_group = 0.0;

  // Declare Derived quantities
  vector<Type> SSB(n_projyears);
  vector<Type> annual_recruitment(n_projyears);
  vector<Type> init_natage_m(n_ages);                     // Initial numbers at age Males
  vector<Type> init_natage_f(n_ages);                     // Initial numbers at age Females
  array<Type> weight_maturity_prod_f(n_ages, n_projyears);// Female weight and proportion mature, used to calcualte SSBs etc
  array<Type> natage_m(n_ages, n_projyears + 1);          // Male numbers at age at the beginning of the year from start year to end of last projection year
  array<Type> natage_f(n_ages, n_projyears + 1);          // Female numbers at age at the beginning of the year from start year to end of last projection year
  array<Type> natlength_m(n_length_bins, n_years);        // Male numbers at length at the beginning of the year from start year to end year
  array<Type> natlength_f(n_length_bins, n_years);        // Female numbers at length at the beginning of the year from start year to end year

  array<Type> Z_m(n_ages, n_projyears);                   // Male total mortality at age from start year to end year
  array<Type> Z_f(n_ages, n_projyears);                   // Female total mortality at age from start year to end year
  array<Type> S_m(n_ages, n_projyears);                   // Male Survival at age from start year to end year
  array<Type> S_f(n_ages, n_projyears);                   // Female Survival at age from start year to end year
  array<Type> S_m_mid(n_ages, n_projyears);               // Male Survival at age from start year to end year
  array<Type> S_f_mid(n_ages, n_projyears);               // Female Survival at age from start year to end year
  array<Type> F_ll_m(n_ages, n_projyears);                // Male Fishing mortality Longline at age from start year to end year
  array<Type> F_ll_f(n_ages, n_projyears);                // Female Fishing mortality Longline at age from start year to end year
  array<Type> F_trwl_m(n_ages, n_projyears);              // Male Fishing mortality Trawl at age from start year to end year
  array<Type> F_trwl_f(n_ages, n_projyears);              // Female Fishing mortality Trawl at age from start year to end year
  array<Type> catchatage_ll_m(n_ages, n_years);           // Male Catch at age Longline at age from start year to end year
  array<Type> catchatage_ll_f(n_ages, n_years);           // Female Catch at age Longline at age from start year to end year
  array<Type> catchatage_trwl_m(n_ages, n_years);         // Male Catch at age Trawl at age from start year to end year
  array<Type> catchatage_trwl_f(n_ages, n_years);         // Female Catch at age Trawl at age from start year to end year
  vector<Type> annual_F_ll(n_years);                      // Fishing mortality for longline gear for each model year
  vector<Type> annual_ll_catch_pred(n_years);             // Fishing mortality for longline gear for each model year
  vector<Type> annual_F_trwl(n_years);                    // Fishing mortality for trawl gear for each model year
  vector<Type> annual_trwl_catch_pred(n_years);           // Fishing mortality for longline gear for each model year
  array<Type> sel_ll_f(n_ages, ln_ll_sel_pars.dim(0));              // Longline selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_ll_m(n_ages, ln_ll_sel_pars.dim(0));              // Longline selectivity Male. dim: n_ages x n_projyears
  array<Type> sel_trwl_f(n_ages, ln_trwl_sel_pars.dim(0));            // Trawl selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_trwl_m(n_ages, ln_trwl_sel_pars.dim(0));            // Trawl selectivity Male. dim: n_ages x n_projyears
  array<Type> sel_srv_dom_ll_f(n_ages, ln_srv_dom_ll_sel_pars.dim(0));              // Longline selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_srv_dom_ll_m(n_ages, ln_srv_dom_ll_sel_pars.dim(0));              // Longline selectivity Male. dim: n_ages x n_projyears
  array<Type> sel_srv_nmfs_trwl_f(n_ages, ln_srv_nmfs_trwl_sel_pars.dim(0));              // Longline selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_srv_nmfs_trwl_m(n_ages, ln_srv_nmfs_trwl_sel_pars.dim(0));              // Longline selectivity Male. dim: n_ages x n_projyears
  



  Type alpha = 0.0;                                       // alpha for the stock recruit relationship
  Type beta = 0.0;                                        // beta for the stock recruit relationship
  Type Bzero = 0.0;
  Type sigma_R_sq = sigma_R * sigma_R;
  // Initialise Derived quantities
  SSB.setZero();
  pred_dom_ll_bio.setZero();
  pred_ll_cpue.setZero();
  init_natage_m.setZero();
  init_natage_f.setZero();
  annual_F_ll.fill(0.0);
  annual_F_trwl.fill(0.0);
  Z_f.fill(0.0);
  Z_m.fill(0.0);
  S_f.fill(0.0);
  S_m.fill(0.0);
  /*
   * Calculate some initial Derived quantities
   */
  weight_maturity_prod_f = maturity * female_mean_weight_by_age;
  /*
   * Build selectivity objects - The BuildSelectivity() method can be found in SetupSelectivities.hpp
   */
  BuildSelectivity(ll_sel_pars.col(0), ll_sel_type, ages, sel_ll_m, false);
  BuildSelectivity(ll_sel_pars.col(1), ll_sel_type, ages, sel_ll_f, false);
  BuildSelectivity(trwl_sel_pars.col(0), trwl_sel_type, ages, sel_trwl_m, true);
  BuildSelectivity(trwl_sel_pars.col(1), trwl_sel_type, ages, sel_trwl_f, true);

  BuildSelectivity(srv_dom_ll_sel_pars.col(0), srv_dom_ll_sel_type, ages, sel_srv_dom_ll_m, false);
  BuildSelectivity(srv_dom_ll_sel_pars.col(1), srv_dom_ll_sel_type, ages, sel_srv_dom_ll_f, false);

  BuildSelectivity(srv_nmfs_trwl_sel_pars.col(0), srv_nmfs_trwl_sel_type, ages, sel_srv_nmfs_trwl_m, true);
  BuildSelectivity(srv_nmfs_trwl_sel_pars.col(1), srv_nmfs_trwl_sel_type, ages, sel_srv_nmfs_trwl_f, true);



  // F, Z and survivorship
  annual_F_ll = exp(ln_ll_F_avg + ln_ll_F_devs);
  annual_F_trwl = exp(ln_trwl_F_avg + ln_trwl_F_devs);
  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    F_ll_m.col(year_ndx) = annual_F_ll(year_ndx) * sel_ll_m.col(ll_sel_by_year_indicator(year_ndx));
    F_ll_f.col(year_ndx) = annual_F_ll(year_ndx) * sel_ll_f.col(ll_sel_by_year_indicator(year_ndx));
    F_trwl_m.col(year_ndx) = annual_F_trwl(year_ndx) * sel_trwl_m.col(trwl_sel_by_year_indicator(year_ndx));
    F_trwl_f.col(year_ndx) = annual_F_trwl(year_ndx) * sel_trwl_f.col(trwl_sel_by_year_indicator(year_ndx));
    Z_f.col(year_ndx) = F_ll_f.col(year_ndx) + F_trwl_f.col(year_ndx) + M.col(year_ndx);
    Z_m.col(year_ndx) = F_ll_m.col(year_ndx) + F_trwl_m.col(year_ndx) + M.col(year_ndx);
    S_f.col(year_ndx) = exp(-1.0 * Z_f.col(year_ndx));
    S_m.col(year_ndx) = exp(-1.0 * Z_m.col(year_ndx));
    S_f_mid.col(year_ndx) = exp(-0.5 * Z_f.col(year_ndx));
    S_m_mid.col(year_ndx) = exp(-0.5 * Z_m.col(year_ndx));
  }

  vector<Type> nll(25); // slots
  vector<Type> nll_weighted(25); // slots
  nll.setZero();

  /* nll components
   * 0 - ll- fishery age comp
   * 1 - trwl-fishery length male comp
   * 2 - trwl-fishery length female comp
   * 3 - srv Domestic ll Biomass
   * 4 - srv Japanese ll Biomass
   * 5 - LL Fishery CPUE
   * 6 - srv Domestic ll Age
   * 7 - srv Domestic ll Length male
   * 8 - srv Domestic ll Length female
   * 9 - srv Japanese ll Age
   * 10 - srv Japanese ll Length male
   * 11 - srv Japanese ll Length female
   * 12 - srv GOA trwl Age
   * 13 - srv GOA trwl  Length male
   * 14 - srv GOA trwl  Length female
   * 15 - ll-fishery length male comp
   * 16 - ll-fishery length female comp
   * 17 - srv GOA trwl biomass index
   * 18 - srv Japanese Fishery longline biomass index
   * 19 - srv Japanese Fishery longline Length Frequency
   * 20 - Longline fishery catch Sum of squares
   * 21 - Trawl fishery catch Sum of squares
   * 22 - Recruitment penalty/hyper prior if model is hierachical
   * 23 - F_penalty_ll - to make F_devs identifiable and a positive definite hessian
   * 24 - F_penalty_trwl - to make F_devs identifiable and a positive definite hessian
   */

  /*
   * Initialise the partition (age structure)
   * at equilibrium only M and R0
   */
  init_natage_m(0) = exp(ln_mean_rec)/2.0; //mean_rec * exp((sigma_R*sigma_R)/2)/2;
  init_natage_f(0) = init_natage_m(0);

  for(age_ndx = 1; age_ndx < n_ages; age_ndx++) {
    // include recruitment variation in intial age-comp
    init_natage_f(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0)) * age_ndx) / 2.0;
    init_natage_m(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0)) * age_ndx) / 2.0;
  }
  // plus group
  init_natage_f(n_ages - 1) = (exp(ln_mean_rec - (M(n_ages - 1, 0)) * (n_ages - 1)) / (1.0 - exp(-(M(n_ages - 1, 0))))) / 2.0;
  init_natage_m(n_ages - 1) = (exp(ln_mean_rec - (M(n_ages - 1, 0)) * (n_ages - 1)) / (1.0 - exp(-(M(n_ages - 1, 0))))) / 2.0;
  // Calculate B0
  for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
    Bzero += init_natage_f(age_ndx) * pow(exp(-M(age_ndx, 0)), spawning_time_proportion(0)) * weight_maturity_prod_f(age_ndx, 0);

  /*
   * Initialise the partition (age structure) with initial F and initial age devs
   */
  init_natage_m(0) = exp(ln_mean_rec + ln_rec_dev(0) + sigma_R_sq/2.0)/2.0; //mean_rec * exp((sigma_R*sigma_R)/2)/2;
  init_natage_f(0) = init_natage_m(0);

  for(age_ndx = 1; age_ndx < n_ages; age_ndx++) {
    if(age_ndx <= n_init_rec_devs) {
      // include recruitment variation in intial age-comp
      init_natage_f(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0) + init_F_hist * sel_ll_f(age_ndx, 0)) * age_ndx + ln_init_rec_dev(age_ndx - 1) + sigma_R_sq/2) / 2;
      init_natage_m(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0) + init_F_hist * sel_ll_m(age_ndx, 0)) * age_ndx + ln_init_rec_dev(age_ndx - 1) + sigma_R_sq/2) / 2;
    } else {
      // assume the last initial age deviation for these
      init_natage_f(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0) + init_F_hist * sel_ll_f(age_ndx, 0)) * age_ndx + ln_init_rec_dev(ln_init_rec_dev.size() - 1) + sigma_R_sq/2) / 2;
      init_natage_m(age_ndx) = exp(ln_mean_rec - (M(age_ndx - 1, 0) + init_F_hist * sel_ll_m(age_ndx, 0)) * age_ndx + ln_init_rec_dev(ln_init_rec_dev.size() - 1) + sigma_R_sq/2) / 2;
    }
  }
  // plus group
  init_natage_f(n_ages - 1) = (exp(ln_mean_rec - (M(n_ages - 1, 0) + init_F_hist * sel_ll_f(n_ages - 1, 0)) * (n_ages - 1)) / (1.0 - exp(-(M(n_ages - 1, 0) + init_F_hist * sel_ll_f(n_ages - 1, 0))))) / 2.0;
  init_natage_m(n_ages - 1) = (exp(ln_mean_rec - (M(n_ages - 1, 0) + init_F_hist * sel_ll_m(n_ages - 1, 0)) * (n_ages - 1)) / (1.0 - exp(-(M(n_ages - 1, 0) + init_F_hist * sel_ll_m(n_ages - 1, 0))))) / 2.0;
  // set numbers at age
  natage_m.col(0) = init_natage_m;
  natage_f.col(0) = init_natage_f;
  /*
   * Run the annual cycle
   */
  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    //fill in recruitment in current year - a slight ineffieciency as we have already done this for year year_ndx = 0 above but meh!
    natage_m(0, year_ndx) = exp(ln_mean_rec + ln_rec_dev(year_ndx) + sigma_R_sq/2)/2; //
    natage_f(0, year_ndx) = natage_m(0, year_ndx);
    annual_recruitment(year_ndx) = natage_m(0, year_ndx) + natage_f(0, year_ndx);
    // Fill in N in following year, Survival from prev since uses i and N used i+1
    m_plus_group = natage_m(n_ages - 1, year_ndx);
    f_plus_group = natage_f(n_ages - 1, year_ndx);
    //std::cout << "year ndx = " << year_ndx << "\n";
    for(age_ndx = 0; age_ndx < (n_ages - 1); age_ndx++) {
      //std::cout << "age ndx = " << age_ndx << "\n";
      natage_m(age_ndx + 1, year_ndx + 1) =  natage_m(age_ndx, year_ndx) * S_m(age_ndx, year_ndx);
      natage_f(age_ndx + 1, year_ndx + 1) =  natage_f(age_ndx, year_ndx) * S_f(age_ndx, year_ndx);
      SSB(year_ndx) += natage_f(age_ndx, year_ndx) * pow(S_f(age_ndx, year_ndx) ,spawning_time_proportion(year_ndx)) * weight_maturity_prod_f(age_ndx, year_ndx);
    }

    // plus group accumulation
    natage_m(n_ages - 1, year_ndx + 1) +=  m_plus_group * S_m(n_ages - 1, year_ndx);
    natage_f(n_ages - 1, year_ndx + 1) +=  f_plus_group * S_f(n_ages - 1, year_ndx);
    SSB(year_ndx) += f_plus_group * pow(S_f(n_ages - 1, year_ndx), spawning_time_proportion(year_ndx)) * weight_maturity_prod_f(n_ages - 1, year_ndx);

    // Calculate numbers at length
    natlength_m.col(year_ndx) = (male_age_length_transition.col(year_ndx).matrix().transpose() * natage_m.col(year_ndx).matrix()).col(0);
    natlength_f.col(year_ndx) = (female_age_length_transition.col(year_ndx).matrix().transpose() * natage_f.col(year_ndx).matrix()).col(0);

    // Calculate Catch at age
    catchatage_ll_m.col(year_ndx) = F_ll_m.col(year_ndx) / Z_m.col(year_ndx) * natage_m.col(year_ndx) * (1.0 - S_m.col(year_ndx));
    catchatage_ll_f.col(year_ndx) = F_ll_f.col(year_ndx) / Z_f.col(year_ndx) * natage_f.col(year_ndx) * (1.0 - S_f.col(year_ndx));
    annual_ll_catch_pred(year_ndx) = (catchatage_ll_f.col(year_ndx) * female_mean_weight_by_age.col(year_ndx) + catchatage_ll_m.col(year_ndx) * male_mean_weight_by_age.col(year_ndx)).sum();
    catchatage_trwl_m.col(year_ndx) = F_trwl_m.col(year_ndx) / Z_m.col(year_ndx) * natage_m.col(year_ndx) * (1.0 - S_m.col(year_ndx));
    catchatage_trwl_f.col(year_ndx) = F_trwl_f.col(year_ndx) / Z_f.col(year_ndx) * natage_f.col(year_ndx) * (1.0 - S_f.col(year_ndx));
    annual_trwl_catch_pred(year_ndx) = (catchatage_trwl_f.col(year_ndx) * female_mean_weight_by_age.col(year_ndx) + catchatage_trwl_m.col(year_ndx) * male_mean_weight_by_age.col(year_ndx)).sum();
  }

  /*
   * Calculate expected/Predicted observation values
   * and also calculate objective score contribution to nll
   */
  // Calculate Catch at age
  int ll_catchatage_ndx = 0;
  int ll_catchatlgth_ndx = 0;
  int trwl_catchatlgth_ndx = 0;
  int srv_dom_ll_bio_ndx = 0;

  int ll_cpue_ndx = 0;
  int srv_dom_ll_age_ndx = 0;
  int srv_dom_ll_lgth_ndx = 0;

  int srv_nmfs_trwl_age_ndx = 0;
  int srv_nmfs_trwl_lgth_ndx = 0;
  int srv_nmfs_trwl_bio_ndx = 0;


  vector<Type> temp_numbers_at_age(n_ages);
  vector<Type> temp_numbers_at_lgth(n_length_bins);
  vector<Type> temp_numbers_at_age_after_ageing_error(n_ages);


  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    // Longline Fishery age comp
    if(ll_catchatage_indicator(year_ndx) == 1) {
      temp_numbers_at_age = ((catchatage_ll_m.col(year_ndx) / sum(vector<Type>(catchatage_ll_m.col(year_ndx)))) + (catchatage_ll_f.col(year_ndx) / sum(vector<Type>(catchatage_ll_f.col(year_ndx))))) / 2.0;
      temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
      pred_ll_catchatage.col(ll_catchatage_ndx) = temp_numbers_at_age_after_ageing_error / sum(temp_numbers_at_age_after_ageing_error);
      // Log-likelihood contribution
      if(ll_catchatage_comp_likelihood == 0) {
        // ADMB's multinomial method
        Type sample_size = sum(obs_ll_catchatage.col(ll_catchatage_ndx));
        temp_numbers_at_age = obs_ll_catchatage.col(ll_catchatage_ndx) / sample_size;
        nll(0) -= sample_size * sum(vector<Type>((temp_numbers_at_age + 0.001) * log(pred_ll_catchatage.col(ll_catchatage_ndx) + 0.001)));
      } else if(ll_catchatage_comp_likelihood == 1) {
        // TMB's multinomial method
        nll(0) -= dmultinom_upd(obs_ll_catchatage.col(ll_catchatage_ndx).vec(), pred_ll_catchatage.col(ll_catchatage_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        Type sample_size = sum(obs_ll_catchatage.col(ll_catchatage_ndx));
        vector<Type> temp_observed_age = rmultinom(pred_ll_catchatage.col(ll_catchatage_ndx).vec(), sample_size);
        obs_ll_catchatage.col(ll_catchatage_ndx) = temp_observed_age;
      }
      ++ll_catchatage_ndx;
    }

    // Longline Fishery length frequency
    if(ll_catchatlgth_indicator(year_ndx) == 1) {
      pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx) = (male_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_ll_m.col(year_ndx).matrix()).col(0);
      pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx) = (female_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_ll_f.col(year_ndx).matrix()).col(0);
      // normalise to sum = 1
      pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx) /= sum(pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx));
      pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx) /= sum(pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx));
      // Log-likelihood contribution
      if(ll_catchatlgth_comp_likelihood == 0) {
        // ADMB's multinomial method
        // males
        Type sample_size = sum(obs_ll_catchatlgth_m.col(ll_catchatlgth_ndx));
        temp_numbers_at_lgth = obs_ll_catchatlgth_m.col(ll_catchatlgth_ndx) / sample_size;
        nll(15) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx) + 0.001)));
        // females
        sample_size = sum(obs_ll_catchatlgth_f.col(ll_catchatlgth_ndx));
        temp_numbers_at_lgth = obs_ll_catchatlgth_f.col(ll_catchatlgth_ndx) / sample_size;
        nll(16) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx) + 0.001)));
      } else if(ll_catchatage_comp_likelihood == 1) {
        // TMB's multinomial method
        nll(15) -= dmultinom_upd(obs_ll_catchatlgth_m.col(ll_catchatlgth_ndx).vec(), pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx).vec(), true);
        nll(16) -= dmultinom_upd(obs_ll_catchatlgth_f.col(ll_catchatlgth_ndx).vec(), pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        // Males first
        Type sample_size = sum(obs_ll_catchatlgth_m.col(ll_catchatlgth_ndx));
        vector<Type> temp_observed_lgth = rmultinom(pred_ll_catchatlgth_m.col(ll_catchatlgth_ndx).vec(), sample_size);
        obs_ll_catchatlgth_m.col(ll_catchatlgth_ndx) = temp_observed_lgth;
        // Females second
        sample_size = sum(obs_ll_catchatlgth_f.col(ll_catchatlgth_ndx));
        temp_observed_lgth = rmultinom(pred_ll_catchatlgth_f.col(ll_catchatlgth_ndx).vec(), sample_size);
        obs_ll_catchatlgth_f.col(ll_catchatlgth_ndx) = temp_observed_lgth;
      }
      ++ll_catchatlgth_ndx;
    }

    // Trawl Fishery length frequency
    if(trwl_catchatlgth_indicator(year_ndx) == 1) {
      pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx) = (male_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_trwl_m.col(year_ndx).matrix()).col(0);
      pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx) = (female_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_trwl_f.col(year_ndx).matrix()).col(0);
      // normalise to sum = 1
      pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx) /= sum(pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx));
      pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx) /= sum(pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx));
      // Log-likelihood contribution
      if(trwl_catchatlgth_comp_likelihood == 0) {
        // ADMB's multinomial method
        // males
        Type sample_size = sum(obs_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx));
        temp_numbers_at_lgth = obs_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx) / sample_size;
        nll(1) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx) + 0.001)));
        // females
        sample_size = sum(obs_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx));
        temp_numbers_at_lgth = obs_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx) / sample_size;
        nll(2) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx) + 0.001)));

      } else if(ll_catchatage_comp_likelihood == 1) {
        // TMB's multinomial method
        nll(1) -= dmultinom_upd(obs_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx).vec(), pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx).vec(), true);
        nll(2) -= dmultinom_upd(obs_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx).vec(), pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        // Males first
        Type sample_size = sum(obs_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx));
        vector<Type> temp_observed_lgth = rmultinom(pred_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx).vec(), sample_size);
        obs_trwl_catchatlgth_m.col(trwl_catchatlgth_ndx) = temp_observed_lgth;
        // Females second
        sample_size = sum(obs_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx));
        temp_observed_lgth = rmultinom(pred_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx).vec(), sample_size);
        obs_trwl_catchatlgth_f.col(trwl_catchatlgth_ndx) = temp_observed_lgth;
      }
      ++trwl_catchatlgth_ndx;
    }

    // Domestic longline survey biomass
    if(srv_dom_ll_bio_indicator(year_ndx) == 1) {
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
        pred_dom_ll_bio(srv_dom_ll_bio_ndx) += proportion_male(year_ndx) * natage_m(age_ndx, year_ndx) * S_m_mid(age_ndx, year_ndx) * sel_srv_dom_ll_m(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx) + (1.0 - proportion_male(year_ndx)) * natage_f(age_ndx, year_ndx) * S_f_mid(age_ndx, year_ndx) * sel_srv_dom_ll_f(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx);
      // account for catchability and times 2 ???
      pred_dom_ll_bio(srv_dom_ll_bio_ndx)  *= 2 * srv_dom_ll_q(srv_dom_ll_q_by_year_indicator(year_ndx));
      if(srv_dom_ll_bio_likelihood == 0) {
        nll(3) += square((log(obs_dom_ll_bio(srv_dom_ll_bio_ndx) + 0.0001) - log(pred_dom_ll_bio(srv_dom_ll_bio_ndx) + 0.0001) ))/ (2.0 * square(se_dom_ll_bio(srv_dom_ll_bio_ndx) / obs_dom_ll_bio(srv_dom_ll_bio_ndx)));
        // Simulate lognormal observations
        // you may want include a bias correction for the mean of -0.5 sigma^2 so the distribution has expectation = predicted value.
        SIMULATE {
          obs_dom_ll_bio(srv_dom_ll_bio_ndx) = exp(rnorm(log(pred_dom_ll_bio(srv_dom_ll_bio_ndx) + 0.0001), (se_dom_ll_bio(srv_dom_ll_bio_ndx) / obs_dom_ll_bio(srv_dom_ll_bio_ndx))));
        }
      } else {
        nll(3) -= dlnorm(obs_dom_ll_bio(srv_dom_ll_bio_ndx), log(pred_dom_ll_bio(srv_dom_ll_bio_ndx)) - 0.5 * se_dom_ll_bio(srv_dom_ll_bio_ndx) * se_dom_ll_bio(srv_dom_ll_bio_ndx), se_dom_ll_bio(srv_dom_ll_bio_ndx), true);
        // Simulate lognormal observations
        SIMULATE {
          obs_dom_ll_bio(srv_dom_ll_bio_ndx) = exp(rnorm(log(pred_dom_ll_bio(srv_dom_ll_bio_ndx)) - 0.5 * se_dom_ll_bio(srv_dom_ll_bio_ndx) * se_dom_ll_bio(srv_dom_ll_bio_ndx), se_dom_ll_bio(srv_dom_ll_bio_ndx)));
        }
      }

      ++srv_dom_ll_bio_ndx;
    }

    // longline Fishery CPUE
    if(ll_cpue_indicator(year_ndx) == 1) {
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
        pred_ll_cpue(ll_cpue_ndx) += natage_m(age_ndx, year_ndx) * S_m_mid(age_ndx, year_ndx) * sel_ll_m(age_ndx, ll_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx) + natage_f(age_ndx, year_ndx) * S_f_mid(age_ndx, year_ndx) * sel_ll_f(age_ndx, ll_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx);
      // account for catchability and times 2 ???
      pred_ll_cpue(ll_cpue_ndx)  *= ll_cpue_q(ll_cpue_q_by_year_indicator(year_ndx));
      //nll(5) += square((log(obs_ll_cpue(ll_cpue_ndx) + 0.0001) - log(pred_ll_cpue(ll_cpue_ndx) + 0.0001) ))/ (2.0 * square(se_ll_cpue(ll_cpue_ndx) / obs_ll_cpue(ll_cpue_ndx)));
      if(ll_cpue_likelihood == 0) {
        nll(5) += square((log(obs_ll_cpue(ll_cpue_ndx) + 0.0001) - log(pred_ll_cpue(ll_cpue_ndx) + 0.0001) ))/ (2.0 * square(se_ll_cpue(ll_cpue_ndx) / obs_ll_cpue(ll_cpue_ndx)));
        // Simulate lognormal observations
        // you may want include a bias correction for the mean of -0.5 sigma^2 so the distribution has expectation = predicted value.
        SIMULATE {
          obs_ll_cpue(ll_cpue_ndx) = exp(rnorm(log(pred_ll_cpue(ll_cpue_ndx) + 0.0001), (se_ll_cpue(ll_cpue_ndx) / obs_ll_cpue(ll_cpue_ndx))));
        }
      } else {
        nll(5) -= dlnorm(obs_ll_cpue(ll_cpue_ndx), log(pred_ll_cpue(ll_cpue_ndx)) - 0.5 * se_ll_cpue(ll_cpue_ndx) * se_ll_cpue(ll_cpue_ndx), se_ll_cpue(ll_cpue_ndx), true);
        // Simulate lognormal observations
        SIMULATE {
          obs_ll_cpue(ll_cpue_ndx) = exp(rnorm(log(pred_ll_cpue(ll_cpue_ndx)) - 0.5 * se_ll_cpue(ll_cpue_ndx) * se_ll_cpue(ll_cpue_ndx), se_ll_cpue(ll_cpue_ndx)));
        }
      }
      ++ll_cpue_ndx;
    }

    // Survey Domestic Longline age frequency
    if(srv_dom_ll_age_indicator(year_ndx) == 1) {
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++) {
        temp_numbers_at_age(age_ndx) = proportion_male(year_ndx) * natage_m(age_ndx, year_ndx) * sel_srv_dom_ll_m(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) + (1.0 - proportion_male(year_ndx)) * natage_f(age_ndx, year_ndx) * sel_srv_dom_ll_f(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx));
      }
      temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
      pred_srv_dom_ll_age.col(srv_dom_ll_age_ndx) = temp_numbers_at_age_after_ageing_error / sum(temp_numbers_at_age_after_ageing_error);
      // Log-likelihood contribution
      if(srv_dom_ll_age_comp_likelihood == 0) {
        Type sample_size = sum(obs_srv_dom_ll_age.col(srv_dom_ll_age_ndx));
        temp_numbers_at_age = obs_srv_dom_ll_age.col(srv_dom_ll_age_ndx) / sample_size;
        nll(6) -= sample_size * sum(vector<Type>((temp_numbers_at_age + 0.001) * log(pred_srv_dom_ll_age.col(srv_dom_ll_age_ndx) + 0.001)));
      } else if(ll_catchatage_comp_likelihood == 1) {
        nll(6) -= dmultinom_upd(obs_srv_dom_ll_age.col(srv_dom_ll_age_ndx).vec(), pred_srv_dom_ll_age.col(srv_dom_ll_age_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        Type sample_size = sum(obs_srv_dom_ll_age.col(srv_dom_ll_age_ndx));
        vector<Type> temp_observed_age = rmultinom(pred_srv_dom_ll_age.col(srv_dom_ll_age_ndx).vec(), sample_size);
        obs_srv_dom_ll_age.col(srv_dom_ll_age_ndx) = temp_observed_age;
      }
      ++srv_dom_ll_age_ndx;
    }

    // Survey Domestic Longline Length frequency
    if(srv_dom_ll_lgth_indicator(year_ndx) == 1) {
      temp_numbers_at_lgth = natage_m.col(year_ndx) * sel_srv_dom_ll_m.col(srv_dom_ll_sel_by_year_indicator(year_ndx));
      pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx) = male_age_length_transition.col(year_ndx).matrix().transpose() * temp_numbers_at_lgth.matrix();
      temp_numbers_at_lgth = natage_f.col(year_ndx) * sel_srv_dom_ll_f.col(srv_dom_ll_sel_by_year_indicator(year_ndx));
      pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) = female_age_length_transition.col(year_ndx).matrix().transpose() * temp_numbers_at_lgth.matrix();

      //pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) = (female_age_length_transition.col(year_ndx).matrix().transpose() * (natage_f.col(year_ndx) * sel_srv_dom_ll_f.col(srv_dom_ll_sel_by_year_indicator(year_ndx))).matrix()).col(0);
      // normalise to sum = 1
      pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx) /= sum(pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx));
      pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) /= sum(pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx));
      // Log-likelihood contribution
      if(srv_dom_ll_lgth_comp_likelihood == 0) {
        // ADMB's multinomial method
        // males
        Type sample_size = sum(obs_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx));
        temp_numbers_at_lgth = obs_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx) / sample_size;
        nll(7) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx) + 0.001)));
        // females
        sample_size = sum(obs_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx));
        temp_numbers_at_lgth = obs_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) / sample_size;
        nll(8) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) + 0.001)));

      } else if(srv_dom_ll_lgth_comp_likelihood == 1) {
        // TMB's multinomial method
        nll(7) -= dmultinom_upd(obs_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx).vec(), pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx).vec(), true);
        nll(8) -= dmultinom_upd(obs_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx).vec(), pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        // Males first
        Type sample_size = sum(obs_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx));
        vector<Type> temp_observed_lgth = rmultinom(pred_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx).vec(), sample_size);
        obs_srv_dom_ll_lgth_m.col(srv_dom_ll_lgth_ndx) = temp_observed_lgth;
        // Females second
        sample_size = sum(obs_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx));
        temp_observed_lgth = rmultinom(pred_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx).vec(), sample_size);
        obs_srv_dom_ll_lgth_f.col(srv_dom_ll_lgth_ndx) = temp_observed_lgth;
      }
      ++srv_dom_ll_lgth_ndx;
    }
    

    // Survey GOA Trawl Age frequency
    if(srv_nmfs_trwl_age_indicator(year_ndx) == 1) {
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++) {
        temp_numbers_at_age(age_ndx) = proportion_male(year_ndx) * natage_m(age_ndx, year_ndx) * sel_srv_nmfs_trwl_m(age_ndx, srv_nmfs_trwl_sel_by_year_indicator(year_ndx)) + (1.0 - proportion_male(year_ndx)) * natage_f(age_ndx, year_ndx) * sel_srv_nmfs_trwl_f(age_ndx, srv_nmfs_trwl_sel_by_year_indicator(year_ndx));
      }
      temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
      pred_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx) = temp_numbers_at_age_after_ageing_error / sum(temp_numbers_at_age_after_ageing_error);
      // Log-likelihood contribution
      if(srv_nmfs_trwl_age_comp_likelihood == 0) {
        // ADMB's multinomial method
        Type sample_size = sum(obs_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx));
        temp_numbers_at_age = obs_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx) / sample_size;
        nll(12) -= sample_size * sum(vector<Type>((temp_numbers_at_age + 0.001) * log(pred_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx) + 0.001)));
      } else if(ll_catchatage_comp_likelihood == 1) {
        nll(12) -= dmultinom_upd(obs_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx).vec(), pred_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        Type sample_size = sum(obs_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx));
        vector<Type> temp_observed_age = rmultinom(pred_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx).vec(), sample_size);
        obs_srv_nmfs_trwl_age.col(srv_nmfs_trwl_age_ndx) = temp_observed_age;
      }
      ++srv_nmfs_trwl_age_ndx;
    }
    // Survey GOA Trawl Length frequency
    if(srv_nmfs_trwl_lgth_indicator(year_ndx) == 1) {
      temp_numbers_at_lgth = natage_m.col(year_ndx) * sel_srv_nmfs_trwl_m.col(srv_nmfs_trwl_sel_by_year_indicator(year_ndx));
      pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx) = male_age_length_transition.col(year_ndx).matrix().transpose() * temp_numbers_at_lgth.matrix();
      temp_numbers_at_lgth = natage_f.col(year_ndx) * sel_srv_nmfs_trwl_f.col(srv_nmfs_trwl_sel_by_year_indicator(year_ndx));
      pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) = female_age_length_transition.col(year_ndx).matrix().transpose() * temp_numbers_at_lgth.matrix();

      //pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) = (female_age_length_transition.col(year_ndx).matrix().transpose() * (natage_f.col(year_ndx) * sel_srv_nmfs_trwl_f.col(srv_nmfs_trwl_sel_by_year_indicator(year_ndx))).matrix()).col(0);
      // normalise to sum = 1
      pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx) /= sum(pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx));
      pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) /= sum(pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx));
      // Log-likelihood contribution
      if(srv_nmfs_trwl_lgth_comp_likelihood == 0) {
        // ADMB's multinomial method
        // males
        Type sample_size = sum(obs_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx));
        temp_numbers_at_lgth = obs_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx) / sample_size;
        nll(13) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx) + 0.001)));
        // females
        sample_size = sum(obs_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx));
        temp_numbers_at_lgth = obs_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) / sample_size;
        nll(14) -= sample_size * sum(vector<Type>((temp_numbers_at_lgth + 0.001) * log(pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) + 0.001)));

      } else if(srv_nmfs_trwl_lgth_comp_likelihood == 1) {
        nll(13) -= dmultinom_upd(obs_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx).vec(), pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx).vec(), true);
        nll(14) -= dmultinom_upd(obs_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx).vec(), pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx).vec(), true);
      }
      // Simulate Multinomial observations. These will be integers i.e., numbers at age
      // so will contain both proportions and N_eff
      SIMULATE {
        // Males first
        Type sample_size = sum(obs_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx));
        vector<Type> temp_observed_lgth = rmultinom(pred_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx).vec(), sample_size);
        obs_srv_nmfs_trwl_lgth_m.col(srv_nmfs_trwl_lgth_ndx) = temp_observed_lgth;
        // Females second
        sample_size = sum(obs_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx));
        temp_observed_lgth = rmultinom(pred_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx).vec(), sample_size);
        obs_srv_nmfs_trwl_lgth_f.col(srv_nmfs_trwl_lgth_ndx) = temp_observed_lgth;
      }
      ++srv_nmfs_trwl_lgth_ndx;
    }

    // NMFS GOA trawl survey biomass
    if(srv_nmfs_trwl_bio_indicator(year_ndx) == 1) {
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
        pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) += proportion_male(year_ndx) * natage_m(age_ndx, year_ndx) * S_m_mid(age_ndx, year_ndx) * sel_srv_nmfs_trwl_m(age_ndx, srv_nmfs_trwl_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx) + (1.0 - proportion_male(year_ndx)) * natage_f(age_ndx, year_ndx) * S_f_mid(age_ndx, year_ndx) * sel_srv_nmfs_trwl_f(age_ndx, srv_nmfs_trwl_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx);
      // account for catchability and times 2 ???
      pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx)  *= 2 * srv_nmfs_trwl_q(srv_nmfs_trwl_q_by_year_indicator(year_ndx));
      if(srv_nmfs_trwl_bio_likelihood == 0) {
        nll(17) += square((log(obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) + 0.0001) - log(pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) + 0.0001) ))/ (2.0 * square(se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) / obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx)));
        // Simulate lognormal observations
        // you may want include a bias correction for the mean of -0.5 sigma^2 so the distribution has expectation = predicted value.
        SIMULATE {
          obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) = exp(rnorm(log(pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) + 0.0001), (se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) / obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx))));
        }
      } else {
        nll(17) -= dlnorm(obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx), log(pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx)) - 0.5 * se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) * se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx), se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx), true);
        // Simulate lognormal observations
        SIMULATE {
          obs_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) = exp(rnorm(log(pred_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx)) - 0.5 * se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx) * se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx), se_nmfs_trwl_bio(srv_nmfs_trwl_bio_ndx)));
        }
      }
      ++srv_nmfs_trwl_bio_ndx;
    }

  }

  /*
   * Additional objective function components that are not observations
   */
  // Sum of squares for catch using the norm2 function same which is (sum(x^2))
  if(catch_likelihood == 0) {
    nll(20) = (square(vector<Type>(log(ll_fishery_catch+0.01)-log(annual_ll_catch_pred+0.01)))).sum();
    nll(21) = (square(vector<Type>(log(trwl_fishery_catch+0.8)-log(annual_trwl_catch_pred+0.8)))).sum();
  } else {
    for(year_ndx = 0; year_ndx < ll_fishery_catch.size(); ++year_ndx) {
      nll(20) -= dnorm(log(ll_fishery_catch(year_ndx)), log(annual_ll_catch_pred(year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd, true);
      nll(21) -= dnorm(log(trwl_fishery_catch(year_ndx)), log(annual_trwl_catch_pred(year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd, true);
    }
  }

  SIMULATE {
    // Set simulated catch to predicted
    if(catch_likelihood == 0) {
      trwl_fishery_catch = annual_trwl_catch_pred;
      ll_fishery_catch = annual_ll_catch_pred;
    } else {
      for(year_ndx = 0; year_ndx < ll_fishery_catch.size(); ++year_ndx) {
        ll_fishery_catch(year_ndx) =  exp(rnorm(log(annual_ll_catch_pred(year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd));
        trwl_fishery_catch(year_ndx) =  exp(rnorm(log(annual_trwl_catch_pred(year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd));
      }
    }
  }
  // Recruitment penalty
  for(year_ndx = 0; year_ndx < n_init_rec_devs; ++year_ndx)
    nll(22) += square(ln_init_rec_dev(year_ndx) + sigma_R_sq / 2.)/(2.* sigma_R_sq);
  for(year_ndx = 0; year_ndx < ln_rec_dev.size(); ++year_ndx)
    nll(22) += square(ln_rec_dev(year_ndx) + sigma_R_sq / 2.)/(2.* sigma_R_sq);
  nll(22) += (ln_rec_dev.size() + n_init_rec_devs) * log(sigma_R);
  // F-dev penalty
  nll(23) = square(ln_ll_F_devs).sum();
  nll(24) = square(ln_trwl_F_devs).sum();

  // Apply Log likelihood weights Yuck!!
  nll_weighted = nll;
  nll_weighted(0) *= loglik_wgt_ll_catchatage;            // 0 - ll- fishery age comp
  nll_weighted(1) *= loglik_wgt_trwl_catchatlgth_m;       // 1 - trwl-fishery length male comp
  nll_weighted(2) *= loglik_wgt_trwl_catchatlgth_f;       // 2 - trwl-fishery length female comp
  nll_weighted(3) *= 1.0;                                 // 3 - srv Domestic ll Biomass
  nll_weighted(4) *= 1.0;                                 // 4 - srv Japanese ll Biomass
  nll_weighted(5) *= 1.0;                                 // 5 - LL Fishery CPUE
  nll_weighted(6) *= loglik_wgt_srv_dom_ll_age;           // 6 - srv Domestic ll Age
  nll_weighted(7) *= loglik_wgt_srv_dom_ll_lgth_m;        // 7 - srv Domestic ll Length male
  nll_weighted(8) *= loglik_wgt_srv_dom_ll_lgth_f;        // 8 - srv Domestic ll Length female

  nll_weighted(12) *= loglik_wgt_srv_nmfs_trwl_age;       // 12 - srv GOA trwl Age
  nll_weighted(13) *= loglik_wgt_srv_nmfs_trwl_lgth_m;    // 13 - srv GOA trwl  Length male
  nll_weighted(14) *= loglik_wgt_srv_nmfs_trwl_lgth_f;    // 14 - srv GOA trwl  Length female
  nll_weighted(15) *= loglik_wgt_ll_catchatlgth_m;        // 15 - ll-fishery length male comp
  nll_weighted(16) *= loglik_wgt_ll_catchatlgth_f;        // 16 - ll-fishery length female comp
  nll_weighted(17) *= 1.0;                                // 17 - srv GOA trwl biomass index
  nll_weighted(18) *= 1.0;                                // 18 - srv Japanese Fishery longline biomass index

  nll_weighted(20) *= loglik_wgt_ll_catch;                // 20 - Longline fishery catch Sum of squares
  nll_weighted(21) *= loglik_wgt_trwl_catch;              // 21 - Trawl fishery catch Sum of squares
  nll_weighted(22) *= 1.0;                                // 22 - Recruitment penalty/hyper prior if model is hierachical
  nll_weighted(23) *= loglik_wgt_Fs;                      // 23 - Longline F penalty
  nll_weighted(24) *= loglik_wgt_Fs;                      // 24 - Trawl F penalty

  vector<Type> depletion = SSB / Bzero * 100;
  /*
   * Report section
   */
  REPORT(nll);
  REPORT(nll_weighted);
  REPORT(mean_rec);
  REPORT(rec_dev);
  REPORT(Bzero);
  REPORT(init_natage_f);
  REPORT(init_natage_m);
  REPORT(ln_init_rec_dev);
  REPORT(SSB);
  REPORT(depletion);
  REPORT(natage_f);
  REPORT(natage_m);
  REPORT(natlength_m);
  REPORT(natlength_f);
  REPORT(init_F_hist);
  REPORT(annual_F_trwl);
  REPORT(annual_F_ll);
  REPORT(annual_recruitment);
  REPORT( ln_rec_dev );
  REPORT( sigma_R );
  REPORT(S_f);
  REPORT(S_m);
  REPORT(Z_f);
  REPORT(Z_m);
  REPORT(weight_maturity_prod_f);
  REPORT(catchatage_ll_m);
  REPORT(catchatage_trwl_m);
  REPORT(catchatage_ll_f);
  REPORT(catchatage_trwl_f);
  REPORT(annual_trwl_catch_pred);
  REPORT(annual_ll_catch_pred);

  // Selectivity outputs
  REPORT(sel_ll_m);
  REPORT(sel_ll_f);
  REPORT(sel_trwl_f);
  REPORT(sel_trwl_m);
  REPORT(sel_srv_dom_ll_f);
  REPORT(sel_srv_dom_ll_m);

  REPORT(sel_srv_nmfs_trwl_f);
  REPORT(sel_srv_nmfs_trwl_m);


  REPORT(srv_dom_ll_sel_pars);

  REPORT(srv_nmfs_trwl_sel_pars);
  REPORT(ll_sel_pars);
  REPORT(trwl_sel_pars);

  // Catchability coeffecients

  REPORT(ll_cpue_q);

  REPORT( srv_nmfs_trwl_q );
  REPORT(srv_dom_ll_q);

  REPORT(F_ll_m);
  REPORT(F_ll_f);
  REPORT(F_trwl_m);
  REPORT(F_trwl_f);

  //
  REPORT(ages);
  REPORT( min_age );
  REPORT(years);
  REPORT(length_bins);
  REPORT(n_projections_years);
  REPORT( n_regions );
  // Report model expected/predicted values
  REPORT(pred_ll_catchatage);
  REPORT(pred_trwl_catchatlgth_m);
  REPORT(pred_trwl_catchatlgth_f);
  REPORT(pred_dom_ll_bio);

  REPORT(pred_ll_cpue);
  REPORT(pred_srv_dom_ll_age);
  REPORT(pred_srv_dom_ll_lgth_m);
  REPORT(pred_srv_dom_ll_lgth_f);

  REPORT(pred_srv_nmfs_trwl_age);
  REPORT(pred_srv_nmfs_trwl_lgth_m);
  REPORT(pred_srv_nmfs_trwl_lgth_f);
  REPORT(pred_ll_catchatlgth_m);
  REPORT(pred_ll_catchatlgth_f);
  REPORT(pred_nmfs_trwl_bio);

  // Report model observed values
  REPORT(obs_ll_catchatage);
  REPORT(obs_trwl_catchatlgth_m);
  REPORT(obs_trwl_catchatlgth_f);

  REPORT(obs_srv_dom_ll_age);
  REPORT(obs_srv_dom_ll_lgth_m);
  REPORT(obs_srv_dom_ll_lgth_f);

  REPORT(obs_srv_nmfs_trwl_age);
  REPORT(obs_srv_nmfs_trwl_lgth_m);
  REPORT(obs_srv_nmfs_trwl_lgth_f);
  REPORT(obs_ll_catchatlgth_m);
  REPORT(obs_ll_catchatlgth_f);
  REPORT(obs_nmfs_trwl_bio);
  REPORT(obs_dom_ll_bio);
  REPORT(obs_ll_cpue);

  REPORT(ll_fishery_catch);
  REPORT(trwl_fishery_catch);

  // REport standard errors
  REPORT(se_dom_ll_bio);

  REPORT(se_ll_cpue);

  REPORT(se_nmfs_trwl_bio);
  // Report model type to help R functions
  REPORT(model_type);

  // Likelihood indicators
  REPORT(ll_catchatage_comp_likelihood);
  REPORT( ll_catchatlgth_comp_likelihood);
  REPORT(trwl_catchatlgth_comp_likelihood);
  REPORT(srv_dom_ll_bio_likelihood);

  REPORT(srv_nmfs_trwl_bio_likelihood);
  REPORT(ll_cpue_likelihood);

  REPORT(srv_dom_ll_age_comp_likelihood);
  REPORT(srv_dom_ll_lgth_comp_likelihood);

  REPORT(srv_nmfs_trwl_age_comp_likelihood);
  REPORT(srv_nmfs_trwl_lgth_comp_likelihood);

  // observation time indicators
  REPORT(ll_catchatage_indicator);
  REPORT(ll_catchatlgth_indicator);
  REPORT(ll_cpue_indicator);
  REPORT(trwl_catchatlgth_indicator);
  REPORT(srv_dom_ll_age_indicator);
  REPORT(srv_dom_ll_lgth_indicator);
  REPORT(srv_dom_ll_bio_indicator);
  REPORT(srv_nmfs_trwl_age_indicator);
  REPORT( srv_nmfs_trwl_lgth_indicator);
  REPORT(srv_nmfs_trwl_bio_indicator);


  // Likelihood weights
  REPORT(loglik_wgt_ll_catch);
  REPORT(loglik_wgt_trwl_catch);
  REPORT(loglik_wgt_ll_catchatage);
  REPORT(loglik_wgt_trwl_catchatlgth_m);
  REPORT(loglik_wgt_trwl_catchatlgth_f);
  REPORT(loglik_wgt_trwl_catch);
  REPORT(loglik_wgt_srv_dom_ll_age);
  REPORT(loglik_wgt_srv_dom_ll_lgth_m);
  REPORT(loglik_wgt_srv_dom_ll_lgth_f);

  REPORT(loglik_wgt_srv_nmfs_trwl_age);
  REPORT(loglik_wgt_srv_nmfs_trwl_lgth_m);
  REPORT(loglik_wgt_srv_nmfs_trwl_lgth_f);

  REPORT(loglik_wgt_ll_catchatlgth_m);
  REPORT(loglik_wgt_ll_catchatlgth_f);

  // REMOVE these objects once we have validated
  // I created them for reporting interim calculations

  // standard error reports
  ADREPORT(Bzero);
  ADREPORT(mean_rec);
  ADREPORT(depletion);
  ADREPORT(SSB);
  ADREPORT(annual_F_trwl);
  ADREPORT(annual_F_ll);
  ADREPORT(annual_recruitment);
  ADREPORT(sel_ll_m);
  ADREPORT(sel_ll_f);
  ADREPORT(sel_trwl_f);
  ADREPORT(sel_trwl_m);
  ADREPORT(sel_srv_dom_ll_f);
  ADREPORT(sel_srv_dom_ll_m);

  ADREPORT(sel_srv_nmfs_trwl_f);
  ADREPORT(sel_srv_nmfs_trwl_m);


  return nll_weighted.sum();
}



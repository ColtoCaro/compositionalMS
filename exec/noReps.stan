


data {
int<lower=0> npep ;
int<lower=0> N ;
int<lower=0> nprotch ;

real SN90 ;
real SN95 ;
real y[N] ;
real IS[npep] ;
real SN[npep] ;

//real lowIS[npep] ;
//real lowSN[npep] ;

int<lower=1,upper=npep> opep[N] ;
int<lower=1,upper=nprotch> oprotch[N] ;

vector[2] betamean ;
}

parameters{

//real beta_pep[npep] ;  // pep effect
//vector[2] beta[nprotch] ;  // channel/prot averages independent
real beta[nprotch] ;
real deviate1[nprotch] ;
//real deviate2[nprotch] ;

real<lower=0> tau ; // peptide effect var
real<lower=0> sigma ; // experimental error
//vector<lower=0>[2] bsds ; // random int and slope SD
//real<lower=0> bsds ;
//real<lower=-1,upper=1> rho1 ; //cor between slope and int
//real<lower=-1,upper=1> rho2 ;
//real<lower=-1,upper=1> rho3 ;

//vector[2] mu_beta ;
//cholesky_factor_corr[2] L ;
//vector<lower=0>[2] priorsd ;

}

transformed parameters {
//vector[2] beta[nprotch] ;
//for (i in 1:nprotch){
//  beta[i] = mu_beta + bsds .* (L * alpha[i]) ;  //int_SSN*beta gives the protch mean
//}

//matrix[2,2] Gcorr ;
//Gcorr[1,1] = 1 ;
//Gcorr[1,2] = rho1 ;
//Gcorr[1,3] = rho3 ;
//Gcorr[2,1] = rho1 ;
//Gcorr[2,2] = 1 ;
//Gcorr[2,3] = rho2 ;
//Gcorr[3,1] = rho3 ;
//Gcorr[3,2] = rho2 ;
//Gcorr[3,3] = 1 ;

real beta_slope1[nprotch] ;
//real beta_slope2[nprotch] ;
for (i in 1:nprotch){
  beta_slope1[i] = beta[i] + deviate1[i] ;
  //beta_slope2[i] = beta[i] + deviate2[i] ;
}

}

model {
// First do parameters that effect the protein level
//mu_beta ~ normal(0,10) ;

//bsds ~ normal(0, 5) ; //This is a half-normal
//L ~ lkj_corr_cholesky(1) ;

//tau ~ normal(0, 2) ;
for (i in 1:nprotch){
  //beta[i] ~ multi_normal(betamean, quad_form_diag(Gcorr, bsds)) ;  //
beta[i] ~ normal(0, 5) ;
deviate1[i] ~ normal(0, 2) ;
//deviate2[i] ~ normal(0, 2) ;
}

//Now compute peptide deviations
  //peptide prior
//beta_pep ~ normal(0, 5) ;


//finally add it together
sigma ~ normal(0, 5) ;
for (i in 1:N)
  y[i] ~ normal( beta[oprotch[i]] + SN[opep[i]] * beta_slope1[oprotch[i]] , sigma) ;  //

}  // end model statement


generated quantities {

  vector[nprotch] ppred;
  //vector[nprotch] ppred95;
  for (i in 1:nprotch){
  ppred[i] = beta[i] + SN95 * beta_slope1[i];  //
  //ppred95[i] =  beta[i] + SN95 * beta_slope1[i] ;
}

}




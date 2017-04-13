// Stan code for all heierarchical global protein models and ptms

data{
  int<lower=0> N_ ;  //number of data points
  int<lower=0> n_b ; //number of biologicaly unique proteins
  int<lower=0> n_c ;  //number of protein*conditions (columns that aren't replicates)
  int<lower=0> n_gc ;  //number of group conditions per protein
  int<lower=0> n_t ; // number of column/tag effects
  int<lower=0> n_p ; // number of ptm peptides (0 if not a ptm experiment)
  int<lower=0> n_ptm ; // number of ptm's

  int<lower=0, upper=n_c> techID[N_] ; //ID for technical replicates
  int<lower=0, upper=n_t> tag[N_] ; //ID for tag*tenplex combinations
  int<lower=0, upper=n_b> bioID[N_] ; //ID for biological replicates
  int<lower=0, upper=n_gc> condID[N_] ; //ID for group conditions
  int<lower=0, upper=n_ptm> ptm[N_] ; //ID for ptms (determines variance parameter)
  int<lower=0, upper=n_p> ptmPep[N_] ; //ID for ptm peptides

  int<lower=0> useCov ;  //Indicator for use of covariate (either 0 or 1)
  real<lower=0, upper=1> covariate[N_] ; // covariate
  real pp ; // prediction percentile for the covariate

  real lr[N_] ; // log-ratio observations

}

transformed data {
  int n_bRaw ;
  int n_gcRaw ;
  n_bRaw = (n_b == 0) ? 0 : n_b - 1 ;
  n_gcRaw = (n_gc == 0) ? 0 : n_gc - 1 ;
  }

parameters{
  real beta[n_c] ;
  vector[n_bRaw] beta_bRaw ;
  vector[n_gcRaw] beta_gcRaw ;

  real deviate[useCov * n_c] ; // The correlated slope

  real alpha[n_p] ; // ptm means

  vector[n_t - 1] beta_tRaw ; // tag*tenplex effect

  real<lower = 0> sigma[n_ptm + 1] ; //experimental error
  real<lower = 0> tau[useCov] ; // variance that determines the amount of
                                // correlation between means and slopes

}

transformed parameters{
  real beta_t[n_t] ;  // always have tag effects
  real slope[useCov*n_c] ;
  real beta_b[n_b] ;  //sometimes have biological replicates
  real beta_gc[n_gc] ;  //sometimes have grouped conditions
  real betaP[useCov*n_c] ;  //predicted protein level at pp*covariate

  if(useCov > 0){
  for (i in 1:n_c){
    slope[i] = beta[i] + deviate[i] ;
    betaP[i] = beta[i] + pp*slope[i] ;
  }
  }

  if(n_b > 0){
  for(i in 1:(n_b-1)){
    beta_b[i] = beta_bRaw[i] ;
  }
  beta_b[n_b] = -sum(beta_bRaw) ;
  }

  if(n_gc > 0){
  for(i in 1:(n_gc-1)){
    beta_gc[i] = beta_gcRaw[i] ;
  }
  beta_gc[n_gc] = -sum(beta_gcRaw) ;
  }

  for(i in 1:(n_t-1)){
    beta_t[i] = beta_tRaw[i] ;
  }
  beta_t[n_t] = -sum(beta_tRaw) ;

} // end transform parameters

model{
  //first set parameters that apply to all models
  for(i in 1:n_c){
    beta[i] ~ normal(0, 10) ;
  }

  for(i in 1:(n_t-1)){
    beta_tRaw[i] ~ normal(0, 5) ;
  }

  for(i in 1:(n_ptm + 1)){
  sigma[i] ~ normal(0, 5) ;
  }

  //set ptm distributions
  if(n_ptm >0){
    for(i in 1:n_p){
      alpha[i] ~ normal(0, 10) ;
    }
  }

  //Now work on the different mean protein models
  if(useCov == 0){
  // base model
  if(n_b == 0 && n_gc == 0){

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]], sigma[1]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }
    }
  } // end base model

  // bioRep model
  if(n_b > 0 && n_gc == 0){
    for(i in 1:(n_b-1)){
      beta_bRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_b[bioID[i]] + beta_t[tag[i]], sigma[1]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }
    }
  } // end bioRep model

    // Condition model
  if(n_b == 0 && n_gc > 0){
    for(i in 1:(n_gc - 1)){
      beta_gcRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_gc[condID[i]] + beta_t[tag[i]], sigma[1]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }
    }
  } // end Condition model

    // Full model
  if(n_b > 0 && n_gc > 0){

    for(i in 1:(n_b - 1)){
      beta_bRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:(n_gc - 1)){
      beta_gcRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_b[bioID[i]] + beta_gc[condID[i]] +
      beta_t[tag[i]] , sigma[1]) ;
      }
    if(ptm[i] > 0){
        lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }

    }
  } // end full model

  } //end no covariate

  //Now repeat with covariate use
  if(useCov == 1){
    tau ~ normal(0, 5) ;
    for(i in 1:n_c){
      deviate[i] ~ normal(0, tau) ;
    }
  // base model
  if(n_b == 0 && n_gc == 0){

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_t[tag[i]] + covariate[i]*slope[techID[i]], sigma[1]) ;
    }
      if(ptm[i] > 0){
        lr[i] ~ normal(betaP[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]],
        sigma[1 + ptm[i]]) ;
      }
  }

  } // end base model

  // bioRep model
  if(n_b > 0 && n_gc == 0){
    for(i in 1:(n_b-1)){
      beta_bRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_b[bioID[i]] + beta_t[tag[i]]
      + covariate[i]*slope[techID[i]], sigma[1]) ;
    }
  if(ptm[i] > 0){
        lr[i] ~ normal(betaP[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }
    }

  } // end bioRep model

    // Condition model
  if(n_b == 0 && n_gc > 0){
    for(i in 1:(n_gc - 1)){
      beta_gcRaw[i] ~ normal(0, 10) ;
    }


    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[techID[i]] + beta_gc[condID[i]] + beta_t[tag[i]]
      + covariate[i]*slope[techID[i]], sigma[1]) ;
    }

    if(ptm[i] > 0){
        lr[i] ~ normal(betaP[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
    }
      }
  } // end Condition model

    // Full model
  if(n_b > 0 && n_gc > 0){

    for(i in 1:(n_b - 1)){
      beta_bRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:(n_gc - 1)){
      beta_gcRaw[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
        lr[i] ~ normal(beta[techID[i]] + beta_b[bioID[i]] + beta_gc[condID[i]] +
        beta_t[tag[i]] + covariate[i]*slope[techID[i]], sigma[1]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(betaP[techID[i]] + beta_t[tag[i]] + alpha[ptmPep[i]], sigma[1 + ptm[i]]) ;
      }

    }
  } // end full model

  } //end with covariate




} //end stan program





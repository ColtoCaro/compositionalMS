// Stan code for all heierarchical global protein models and ptms

data{
  int<lower=0> N_ ;  //number of data points
  int<lower=0> n_b ; //number of biologicaly unique proteins
  int<lower=0> n_c ;  //number of protein*condition combinations
  int<lower=0> n_t ;  //number of tag/plex combinations
  int<lower=0> n_p ; // number of ptm peptides (0 if not a ptm experiment)
  int<lower=0> n_ptm ; // number of ptm's
  int<lower=0> n_nc[n_c] ; //number of bio id's within each condition
  int<lower=0> max_nc ; //maximum number of bioReps

  int<lower=0, upper=n_c> condID[N_] ; //ID for condition*prot
  int<lower=0, upper=n_b> bioID[N_] ; //ID for biological replicates (nested
                                      //within condID)
  //int<lower=0, upper=n_t> tagID[N_] ;
  int<lower=0, upper=n_ptm> ptm[N_] ; //ID for ptms (determines variance parameter)
  int<lower=0, upper=n_p> ptmPep[N_] ; //ID for ptm peptides
  int<lower=0, upper=n_b> condToBio[n_c, max(max_nc, 1)] ; //Mapping from bio to cond

  int<lower=0> useCov ;  //Indicator for use of covariate (either 0 or 1)
  real<lower=0> covariate[N_] ; // covariate
  real pp ; // prediction percentile for the covariate

  real lr[N_] ; // log-ratio observations

}

transformed data{
  int<lower=0> bioInd ;
  bioInd = (n_b == 0) ? 0 : 1 ;
}

parameters{
  real beta[n_c * (1-bioInd)] ;
  real beta_b[n_b] ;

//  real deviate_c[useCov * n_c * (1-bioInd)] ; // The correlated slope
  real deviate_b[useCov * n_b] ; // The correlated slope per bioRep

  real alpha[n_p] ; // ptm means

  real<lower = 0> sigma_raw[n_c * (1-bioInd)] ; // experimental error
  real<lower = 0> sigma_rawb[n_b] ; // experimental error
  real<lower = 0> scale ; //heierarchical variance scale
  real<lower = 0> xi[n_ptm] ; // vc's for ptms
  real<lower = 0> delta[useCov] ; //slope multiplier
}

transformed parameters{
//  real slope_c[useCov*n_c* (1-bioInd)] ;
  real betaP_c[useCov*n_c* (1-bioInd)] ;  //predicted protein level
  real slope_b[useCov*n_b] ;
  real betaP_b[useCov*n_b] ;  //predicted protein level at pp*covariate
  real<lower = 0> sigma[n_c * (1-bioInd)] ;
  real<lower = 0> sigmab[n_b] ;

//create the real variance parameters
if(n_b == 0){
  for(i in 1:n_c){
    sigma[i] = scale*sigma_raw[i] ;
  }
}else{
  for(i in 1:n_b){
    sigmab[i] = scale*sigma_rawb[i] ;
  }
}


// Now take care of covariate parameters
  if(n_b > 0 && useCov > 0){
      for (i in 1:n_b){
        //slope_b[i] = beta_b[i] + deviate_b[i] ;
        betaP_b[i] = beta_b[i] * (1 + pp*delta[useCov]);
      }

  }

  if(n_b == 0 && useCov > 0){
      for (i in 1:n_c){
        //slope_c[i] = rho*beta[i]*10 + sqrt(1-rho^2)*deviate_c[i]*10 ;
        betaP_c[i] = beta[i] * (1 + pp*delta[useCov]) ;
      }
  }


} // end transform parameters

model{
  //first set parameters that apply to all models
  scale ~ normal(0, 5) ;

  //set ptm distributions
  if(n_ptm > 0){
    for(i in 1:n_ptm){
      xi[i] ~ cauchy(0, 5) ;
    }
    for(i in 1:n_p){
      alpha[i] ~ normal(0, 10) ;
    }
  }

  //Now work on the different mean protein models
  if(useCov == 0){
  // base model
  if(n_b == 0){

   for(i in 1:n_c){
      beta[i] ~ normal(0, 10) ;
      sigma_raw[i] ~ inv_gamma(1, 1) ;
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[condID[i]] , sigma[condID[i]]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[condID[i]]  + alpha[ptmPep[i]],
          xi[ptm[i]]) ;
      }
    }
  } // end base model

  // bioRep model
  if(n_b > 0){
    for(i in 1:n_b){
      beta_b[i] ~ normal(0, 10) ;
      sigma_rawb[i] ~ inv_gamma(1, 1) ;
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta_b[bioID[i]] , sigmab[bioID[i]]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta_b[bioID[i]] + alpha[ptmPep[i]],
        xi[ptm[i]]) ;
      }
    }
  } // end bioRep model

  } //end no covariate

  //Now repeat with covariate use
  if(useCov == 1){
    //tau ~ normal(0, 5) ;
    //slope ~ normal(0,5) ;
  // base model
  if(n_b == 0){

   for(i in 1:n_c){
      //deviate_c[i] ~ normal(0, 1) ;

      beta[i] ~ normal(0, 10) ;
      sigma_raw[i] ~ inv_gamma(1, 1) ;
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[condID[i]]*covariate[i], sigma[condID[i]]) ;
    }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[condID[i]] + alpha[ptmPep[i]],
        xi[ptm[i]]) ;
      }
  }

  } // end base model

  // bioRep model
  if(n_b > 0){
    for(i in 1:(n_b)){
      deviate_b[i] ~ normal(0, 5) ;
      beta_b[i] ~ normal(0, 10) ;
      sigma_rawb[i] ~ inv_gamma(1, 1) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta_b[bioID[i]]
      + covariate[i]*slope_b[bioID[i]], sigmab[bioID[i]]) ;
    }
  if(ptm[i] > 0){
        lr[i] ~ normal(betaP_b[bioID[i]] + alpha[ptmPep[i]],
        xi[ptm[i]]) ;
      }
    }

  } // end bioRep model


  } //end with covariate




} //end model statement

generated quantities{
  real avgCond[n_c*bioInd] ;

  if(useCov == 0 && bioInd ==1){
    for(i in 1:n_c){
      avgCond[i] = 0;
      for(j in 1:n_nc[i]){
        avgCond[i] = avgCond[i] + beta_b[condToBio[i, j]]/n_nc[i] ;
      }
    }
  }

  if(useCov == 1 && bioInd ==1){
    for(i in 1:n_c){
      avgCond[i] = 0;
      for(j in 1:n_nc[i]){
        avgCond[i] = avgCond[i] + betaP_b[condToBio[i, j]]/n_nc[i] ;
      }
    }
  }


} // end stan program





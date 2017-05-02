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
  int<lower=0, upper=n_t> tagID[N_] ;
  int<lower=0, upper=n_ptm> ptm[N_] ; //ID for ptms (determines variance parameter)
  int<lower=0, upper=n_p> ptmPep[N_] ; //ID for ptm peptides
  int<lower=0, upper=n_b> condToBio[n_c, max(max_nc, 1)] ; //Mapping from bio to cond

  int<lower=0> useCov ;  //Indicator for use of covariate (either 0 or 1)
  real<lower=0, upper=1> covariate[N_] ; // covariate
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

  real deviate_c[useCov * n_c * (1-bioInd)] ; // The correlated slope
  real deviate_b[useCov * n_b] ; // The correlated slope per bioRep

  real alpha[n_p] ; // ptm means

  real<lower = 0> sigmaK ; //constant of proportionality
  real<lower = 0> sigma[n_t] ; //experimental error
  real<lower = 0> tau[useCov] ; // variance that determines the amount of
                                // correlation between means and slopes
}

transformed parameters{
  real slope_c[useCov*n_c* (1-bioInd)] ;
  real betaP_c[useCov*n_c* (1-bioInd)] ;  //predicted protein level
  real slope_b[useCov*n_b] ;
  real betaP_b[useCov*n_b] ;  //predicted protein level at pp*covariate


  if(n_b > 0 && useCov > 0){
      for (i in 1:n_b){
        slope_b[i] = beta_b[i] + deviate_b[i] ;
        betaP_b[i] = beta_b[i] + pp*slope_b[i] ;
      }

  }
  if(n_b == 0 && useCov > 0){
      for (i in 1:n_c){
        slope_c[i] = beta[i] + deviate_c[i] ;
        betaP_c[i] = beta[i] + pp*slope_c[i] ;
      }
  }


} // end transform parameters

model{
  //first set parameters that apply to all models

  for(i in 1:(n_t)){
  sigma[i] ~ normal(0, 5) ;
  }

  //set ptm distributions
  if(n_ptm > 0){
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
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[condID[i]] ,
      sigmaK*fabs(beta[condID[i]])+sigma[tagID[i]]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta[condID[i]]  + alpha[ptmPep[i]],
          sigma[tagID[i]]) ;
      }
    }
  } // end base model

  // bioRep model
  if(n_b > 0){
    for(i in 1:n_b){
      beta_b[i] ~ normal(0, 10) ;
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta_b[bioID[i]] ,
                     sigmaK*fabs(beta_b[bioID[i]]) + sigma[tagID[i]]) ;
      }
      if(ptm[i] > 0){
        lr[i] ~ normal(beta_b[bioID[i]] + alpha[ptmPep[i]], sigma[tagID[i]]) ;
      }
    }
  } // end bioRep model

  } //end no covariate

  //Now repeat with covariate use
  if(useCov == 1){
    tau ~ normal(0, 5) ;

  // base model
  if(n_b == 0){

   for(i in 1:n_c){
      deviate_c[i] ~ normal(0, 5) ;
      beta[i] ~ normal(0, 10) ;
    }
    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta[condID[i]]  +
      covariate[i]*slope_c[condID[i]], sigma[tagID[i]]) ;
    }
      if(ptm[i] > 0){
        lr[i] ~ normal(betaP_c[condID[i]] + alpha[ptmPep[i]],
        sigma[tagID[i]]) ;
      }
  }

  } // end base model

  // bioRep model
  if(n_b > 0){
    for(i in 1:(n_b)){
      deviate_b[i] ~ normal(0, 5) ;
      beta_b[i] ~ normal(0, 10) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
      lr[i] ~ normal(beta_b[bioID[i]]
      + covariate[i]*slope_b[bioID[i]], sigma[tagID[i]]) ;
    }
  if(ptm[i] > 0){
        lr[i] ~ normal(betaP_b[bioID[i]] + alpha[ptmPep[i]], sigma[tagID[i]]) ;
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





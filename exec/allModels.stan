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
  int<lower=0, upper=n_b> condToBio[n_c, max(max_nc, 1)] ; //Mapping from cond to bio
  int<lower=0, upper=n_c> bioToCond[n_b] ; //Mapping from bio to cond

  int<lower=0> useCov ;  //Indicator for use of covariate (either 0 or 1)
  real<lower=0> covariate[N_] ; // covariate
  //real pp ; // prediction percentile for the covariate

  real lr[N_] ; // log-ratio observations

  real<lower=0> pop_sd ; //user input population prior standard deviation
  int<lower=0, upper=1> simpleMod ;  //Indicator to use a simple model

}

transformed data{
  int<lower=0> bioInd ;
  bioInd = (simpleMod == 1) ? 0 : 1 ;
}

parameters{
  real beta[n_c] ;
  real beta_b[n_b * (bioInd)] ;

  real alpha[n_p] ; // ptm means

  real<lower = 0> sigma_raw[n_c * (1 - bioInd)] ; // experimental error
  real<lower = 0> sigma_rawb[n_b * bioInd] ; // experimental error
  real<lower = 0> scale ; //heierarchical variance scale
//  real<lower = 0> tau[bioInd] ; //population level variance
  real<lower = 0> xi[n_ptm] ; // vc's for ptms
  real<lower = 0> delta[useCov] ; //slope multiplier
}

transformed parameters{
  real betaP_c[useCov*n_c] ;  //predicted protein level
  real betaP_b[useCov * n_b * bioInd] ;  //predicted protein level at pp
  real<lower = 0> sigma[n_c * (1 - bioInd)] ;
  real<lower = 0> sigmab[n_b * bioInd] ;

//create the real variance parameters
if(bioInd == 0){
  for(i in 1:n_c){
    sigma[i] = scale*sigma_raw[i] ;
  }
}

if(bioInd == 1){
  for(i in 1:n_b){
    sigmab[i] = scale*sigma_rawb[i] ;
    //beta_b[i] =  beta[bioToCond[i]] + beta_rawb[i] ;
  }
}


// //Now take care of covariate parameters
if(bioInd == 1 && useCov > 0){
    for (i in 1:n_b){
      betaP_b[i] = beta_b[i] * (1 + delta[useCov]);
    }

}

if(useCov > 0){
    for (i in 1:n_c){
      betaP_c[i] = beta[i] * (1 + delta[useCov]) ;
    }
}


} // end transform parameters

 model{
  //first set parameters that apply to all models
  scale ~ normal(0, 5) ;
  beta ~ normal(0, 10) ;

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
  if(simpleMod == 1){

   for(i in 1:n_c){
      sigma_raw[i] ~ inv_gamma(2, 1) ;
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
  if(bioInd > 0){
//    tau ~ normal(0, 5) ;
    for(i in 1:n_b){
      //beta_rawb[i] ~ normal(0, 10) ;
      beta_b[i] ~ normal(beta[bioToCond[i]], pop_sd) ;
      sigma_rawb[i] ~ inv_gamma(2, 1) ;
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

  //base model
 if(bioInd == 0){

  for(i in 1:n_c){
     sigma_raw[i] ~ inv_gamma(2, 1) ;
   }
 for(i in 1:N_){
   if(ptm[i] == 0){
   lr[i] ~ normal(beta[condID[i]] * (1 + covariate[i]*delta[useCov]),
     sigma[condID[i]]) ;
 }
 if(ptm[i] > 0){
   lr[i] ~ normal(beta[condID[i]] + alpha[ptmPep[i]],
   xi[ptm[i]]) ;
 }
 }

  } // end base model

 //  // bioRep model
  if(bioInd > 0){
    for(i in 1:(n_b)){
      beta_b[i] ~ normal(beta[bioToCond[i]], pop_sd) ;
      sigma_rawb[i] ~ inv_gamma(2, 1) ;
    }

    for(i in 1:N_){
      if(ptm[i] == 0){
        lr[i] ~ normal(beta_b[bioID[i]] * (1 + covariate[i]*delta[useCov]),
        sigmab[bioID[i]]) ;
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
  real avgCond[n_c * bioInd] ;

if(bioInd == 1){
  if(useCov == 0){
    for(i in 1:n_c){
      avgCond[i] = 0;
      for(j in 1:n_nc[i]){
        avgCond[i] = avgCond[i] + beta_b[condToBio[i, j]]/n_nc[i] ;
      }
    }
  }

  if(useCov == 1){
    for(i in 1:n_c){
      avgCond[i] = 0;
      for(j in 1:n_nc[i]){
        avgCond[i] = avgCond[i] + betaP_b[condToBio[i, j]]/n_nc[i] ;
      }
    }
  }
} // end if bioInd == 1

} // end generated quantities and stan program





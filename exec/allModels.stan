// Stan code for all heierarchical global protein models and ptms

data{
  int<lower=0> N_ ;  //number of data points
  int<lower=0> n_b ; //number of biologicaly unique proteins
  int<lower=0> n_c ;  //number of protein*condition combinations
  int<lower=0> n_t ;  //number of tag/plex combinations
  int<lower=0> n_p ; // number of ptm peptides (0 if not a ptm experiment)
  int<lower=0> n_ptm ; // number of ptm's
  int<lower=0> n_nc[n_c] ; //number of bio id's within each condition
  int<lower=0> n_pep[n_b] ; //number of peptides per each biorep
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
  real sn[N_] ; // channel signal-noise observations
  real rsn[N_] ; // reference signal-noise observations

  real<lower=0> pop_sd ; //user input population prior standard deviation
  int<lower=0, upper=1> simpleMod ;  //Indicator to use a simple model
  int<lower=0, upper=2> varPool;  //Integer denoting variance pooling.  0 = none, 1 = partial, 2 = complete

}

transformed data{
  int<lower=0> bioInd ;
  int<lower=0> multVar ;
  int<lower=0> scaleInd ;
  bioInd = (simpleMod == 1) ? 0 : 1 ;
  multVar = (varPool == 2) ? 0 : 1 ;  //indicator for multiple variance components
  scaleInd = (varPool == 1) ? 1 : 0 ;  //indicator for multiple variance components
}

parameters{
  real beta[n_c * (1 - bioInd)] ;
  real beta_b[n_b * (bioInd)] ;

  real alpha[n_p] ; // ptm means

  real<lower = 0> sigma_raw[n_c * (1 - bioInd) * (multVar)] ; // experimental error
  real<lower = 0> sigma1[(1 - multVar)] ; // experimental error

  real<lower = 0> sigma_rawb[n_b * bioInd * multVar] ; // experimental error: deprecated
  real<lower = 0> scale[scaleInd] ; //heierarchical variance scale

  real<lower = 0> mu ; // coefficient for channel s/n to modulate variance
  real<lower = 0> nu ; // coefficient for reference s/n to modulate variance

//  real<lower = 0> tau[bioInd] ; //population level variance
  real<lower = 0> xi[n_ptm] ; // vc's for ptms
  real<lower = 0> delta[useCov] ; //slope multiplier
}

transformed parameters{
  real betaP_c[useCov*n_c] ;  //predicted protein level
  real betaP_b[useCov * n_b * bioInd] ;  //predicted protein level at pp
  real<lower = 0> sigma[n_c * (1 - bioInd) * multVar] ;
  real<lower = 0> sigmab[n_b * bioInd] ;
  //real<lower = 0> sigma1[(1 - multVar)] ;

//create the real variance parameters
if(bioInd == 0){
  if(multVar == 1){
    for(i in 1:n_c){
      if(varPool == 0){ //Fit variance with no pooling
          sigma[i] = 5 * sigma_raw[i] ;
        }else{ //Do partial pooling
          if(n_pep[n_c] == 1){
            sigma[i] = scale[1] ;
          }else{
            sigma[i] = scale[1]*sigma_raw[i] ;
          }
        }
     }
  }//End if(multVar) == 1

}//End if bioInd ==0

if(bioInd == 1){
  for(i in 1:n_b){
    if(n_pep[n_b] == 1){
      sigmab[i] = scale[1] ;
    }else{
      sigmab[i] = scale[1]*sigma_rawb[i] ;
      //beta_b[i] =  beta[bioToCond[i]] + beta_rawb[i] ;
    }
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
  if(varPool == 1){
    scale ~ normal(0, 5) ;
    mu ~ normal(0, 5) ;
    nu ~ normal(0, 5) ;
  }

  //beta ~ normal(0, 10) ;  Simple model has been deprecated

  //set ptm distributions
  if(n_ptm > 0){
    for(i in 1:n_ptm){
      xi[i] ~ normal(0, 5) ;
    }
    for(i in 1:n_p){
      alpha[i] ~ normal(0, 10) ;
    }
  }

  //Now work on the different mean protein models
  if(useCov == 0){
  // base model
  if(simpleMod == 1){
   beta ~ normal(0, 10) ;
   if(multVar){
               for(i in 1:n_c){
                 sigma_raw[i] ~ inv_gamma(2, 1) ;
              }
              for(i in 1:N_){
                if(ptm[i] == 0){
                lr[i] ~ normal(beta[condID[i]] , sigma[condID[i]] + mu * 1/pow(sn[i], 2) + nu * 1/pow(rsn[i], 2)) ;
                }
                if(ptm[i] > 0){
                  lr[i] ~ normal(beta[condID[i]]  + alpha[ptmPep[i]],
                    xi[ptm[i]]) ;
                }
              }
   }else{
          sigma1 ~ inv_gamma(2, 5) ;

          for(i in 1:N_){
            if(ptm[i] == 0){
            lr[i] ~ normal(beta[condID[i]] , sigma1) ;
            }
            if(ptm[i] > 0){
              lr[i] ~ normal(beta[condID[i]]  + alpha[ptmPep[i]],
                xi[ptm[i]]) ;
            }
          }
   }

  } // end base model

 // bioRep model
  if(bioInd > 0){
//    tau ~ normal(0, 5) ;
    for(i in 1:n_b){
      beta_b[i] ~ normal(0, 10) ;
      //beta_b[i] ~ normal(beta[bioToCond[i]], pop_sd) ;
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





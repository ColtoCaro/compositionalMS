// Stan code for all heierarchical global protein models and ptms

data{
  int<lower=0> N_ ;  //number of data points
  int<lower=0> n_b ; //number of biologicaly unique proteins
  int<lower=0> n_t ;  //number of tag/plex combinations
  
  int<lower=0, upper=n_b> bioID[N_] ; //ID for biological replicates (nested
                                      //within condID)
  int<lower=0, upper=n_t> tagID[N_] ;
  
  real lr[N_] ; // observations

}

parameters{
  real beta_b[n_b] ;
  real<lower=0> scale;
  //real<lower = 0> sigma[n_b] ; //[n_t] experimental error
  real<lower = 0> sigma_raw[n_b] ;
  }

transformed parameters{
  real<lower = 0> sigma[n_b];
  for(i in 1:n_b){
   sigma[i] = sigma_raw[i]*scale ;
  }
}

model{
  //first set parameters that apply to all models
  scale ~ normal(0,4) ;
  for(i in 1:n_b){
  sigma_raw[i] ~ cauchy(0, 1) ;
  }

  
  // bioRep model
    for(i in 1:n_b){
      beta_b[i] ~ normal(0, 10) ;
    }
    for(i in 1:N_){
      lr[i] ~ normal(beta_b[bioID[i]] , sigma[bioID[i]]) ; //
    }

} //end model statement







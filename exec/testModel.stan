
data{
  int<lower=0> N ;
  real y_[N] ;
}

parameters{
real beta ;
}

model{
beta ~ normal(0, 10) ;
for (i in 1:N){
y_[i] ~ normal(beta, 1) ;
}
}

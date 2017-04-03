//Stan code to figure out the += step

data{
  int<lower = 0> N;
  vector<lower=0, upper=1>[N] theta;
}

parameters {
real<lower=0,upper=.5> phi;
real<lower=0.1> lambda;
}

model {
lambda ~ pareto(0.1, 1.5);

theta ~ beta(lambda * phi, lambda * (1 - phi));
}



data{
  int<lower=0> m;
  int<lower=0> p;
  vector[m]  I;
  matrix[m,p]  X;
  vector[m]    y;
  vector<lower=0>[m]  D;
}

parameters{
  vector[p]  beta;
  real<lower=0> sig2v;
  vector[m] theta;
}

transformed parameters{
  vector[m] sig2v_vec = (1/sig2v) * I;
  vector[m] D_vec = 1 ./ D;
}

model{
  #sig2v ~ uniform(0,1);
  #beta ~ normal(0,5);
  theta ~ multi_normal_prec(X*beta,diag_matrix(sig2v_vec));
  y ~ multi_normal_prec(theta, diag_matrix(D_vec));
}


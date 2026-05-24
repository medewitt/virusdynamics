functions {
  vector ode_base(real t,
                  vector y,
                  real lambda,
                  real d,
                  real beta,
                  real a,
                  real k,
                  real u) {
    vector[3] dydt;

    dydt[1] = lambda - d * y[1] - beta * y[1] * y[3];
    dydt[2] = beta * y[1] * y[3] - a * y[2];
    dydt[3] = k * y[2] - u * y[3];

    return dydt;
  }
}

data {
  int<lower=1> N;                 // number of observation times
  real t0;                        // initial time
  array[N] real ts;               // observation times (ts[n] > t0)
  vector[3] y0;                   // initial state [x, y, v]
  array[N] vector[3] measures;    // observed states
}

parameters {
  real<lower=0> lambda;
  real<lower=0> d;
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> k;
  real<lower=0> u;
  real<lower=0> sigma;
}

transformed parameters {
  array[N] vector[3] y_hat = ode_rk45(ode_base, y0, t0, ts,
                                      lambda, d, beta, a, k, u);
}

model {
  // Weakly informative priors on the rate constants
  lambda ~ normal(0, 1e5);
  d ~ normal(0, 1);
  beta ~ normal(0, 1e-6);
  a ~ normal(0, 1);
  k ~ normal(0, 100);
  u ~ normal(0, 10);
  sigma ~ normal(0, 1);

  for (n in 1:N) {
    measures[n] ~ normal(y_hat[n], sigma);
  }
}

generated quantities {
  // Basic reproductive ratio of the virus
  real R0 = beta * lambda * k / (a * d * u);
}

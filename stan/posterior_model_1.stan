data {
    int <lower=0> N; // number of observations
    vector <lower=0,upper=180> [N] y;       // survival times
    array[N] int <lower=0,upper=1> censor;  //censoring indicator

    real <lower=0> mu_k;          // mean for k
    real <lower=0> sigma_k;       // std for k
    real <lower=0> mu_lambda;     // mean for lambda
    real <lower=0> sigma_lambda;  // std for lambda
}

parameters {
    real <lower=0.1,upper=10> k;       // shape parameter
    real <lower=0.1,upper=100> lambda;  // scale parameter
}

model {
    k ~ lognormal(mu_k, sigma_k);
    lambda ~ lognormal(mu_lambda, sigma_lambda);
    for (n in 1:N) {
        if (censor[n] == 0) {
            target += weibull_lpdf(y[n] | k, lambda);
        } else {
            target += weibull_lccdf(y[n] | k, lambda);
        }
    }
}

generated quantities {
    array [N] real y_sim;   // survival times
    vector[N] log_lik;      // likelihood
    for (n in 1:N) {
        if (censor[n] == 0) {
            log_lik[n] = weibull_lpdf(y[n] | k, lambda);
        } else {
            log_lik[n] = weibull_lccdf(y[n] | k, lambda);
        }
        y_sim[n] = weibull_rng(k, lambda);
    }
}
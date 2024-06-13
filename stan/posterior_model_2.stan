data {
    int <lower=0> N; // number of observations
    int <lower=0> M; // number of features
    matrix <lower=0,upper=1> [N, M] X;      // covariates matrix
    vector <lower=0,upper=180> [N] y;       // survival times
    array[N] int <lower=0,upper=1> censor;  //censoring indicator

    real <lower=0> mu_k;          // mean for k
    real <lower=0> sigma_k;       // std for k
    real <lower=0> mu_theta;      // mean for theta
    real <lower=0> sigma_theta;   // std for theta
}

parameters {
    real<lower=0> k;              // shape parameter
    vector <lower=0> [M] theta;   // regression coefficient
}

model {
    k ~ lognormal(mu_k, sigma_k);
    theta ~ lognormal(mu_theta, sigma_theta);
    for (n in 1:N) {
        if (censor[n] == 0) {
            target += weibull_lpdf(y[n] | k, exp(X[n] * theta));
        } else {
            target += weibull_lccdf(y[n] | k, exp(X[n] * theta));
        }
    }
}

generated quantities {
    array [N] real y_sim;   // survival times
    vector[N] log_lik;      // likelihood
    for (n in 1:N) {
        if (censor[n] == 0) {
            log_lik[n] = weibull_lpdf(y[n] | k, exp(X[n] * theta));
        } else {
            log_lik[n] = weibull_lccdf(y[n] | k, exp(X[n] * theta));
        }
        y_sim[n] = weibull_rng(k, exp(X[n] * theta));
    }
}
data {
    int <lower=0> N; // number of observations
    int <lower=0> M; // number of features
    matrix [N, M] X; // covariates matrix

    real <lower=0> mu_k;          // mean for k
    real <lower=0> sigma_k;       // std for k
    real <lower=0> mu_theta;      // mean for theta
    real <lower=0> sigma_theta;   // std for theta
}

generated quantities {
    real <lower=0> k = normal_rng(mu_k, sigma_k);   // shape parameter
    vector <lower=0> [M] theta;                     // regression coefficient
    for (m in 1:M) {
        theta[m] = normal_rng(mu_theta, sigma_theta);
    }
    array [N] real y_sim; // simulated survival times
    for (n in 1:N) {
        y_sim[n] = weibull_rng(k, exp(X[n] * theta));
    }
}
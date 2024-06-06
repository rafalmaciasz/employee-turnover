data {
    int <lower=0> N; // number of observations
    int <lower=0> K; // number of features
    matrix [N, K] X; // covariates matrix
}

generated quantities {
    real <lower=0> alpha = normal_rng(1, 0.2);
    vector <lower=0> [K] beta;
    for (k in 1:K) {
        beta[k] = normal_rng(2, 0.2);
    }
    array [N] real y_sim; // survival times
    for (n in 1:N) {
        y_sim[n] = weibull_rng(alpha, X[n] * beta);
    }
}
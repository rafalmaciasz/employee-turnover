data {
    int <lower=0> N; // number of observations
    int <lower=0> M; // number of features
    matrix [N, M] X; // covariates matrix

    real <lower=0> mu_k;          // mean for k
    real <lower=0> sigma_k;       // std for k
    real <lower=0> mu_lambda;     // mean for lambda
    real <lower=0> sigma_lambda;  // std for lambda
}

generated quantities {
    real <lower=0> lambda = normal_rng(mu_lambda, sigma_lambda);
    vector <lower=0> [M] k;
    for (m in 1:M) {
        k[m] = normal_rng(mu_k, sigma_k);
    }
    array [N] real y_sim; // simulated survival times
    for (n in 1:N) {
        y_sim[n] = weibull_rng(X[n] * k, lambda);
    }
}
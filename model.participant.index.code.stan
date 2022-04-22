data {
    // global data (values shared across at least two datasets)
    int<lower=1> N_verb;
    int<lower=1> N_scenario;
    int<lower=1> N_polarity_tense;
    int<lower=1> N_participant;

    // norming data
    int<lower=1> N_norming;
    vector[N_norming] y_n;
    int<lower=1, upper=N_scenario> scenario_n[N_norming];
    int<lower=1, upper=N_participant> participant_n[N_norming];

    // templatic validation data
    int<lower=1> N_templatic;
    vector[N_templatic] y_t;
    int<lower=1, upper=N_verb> verb_t[N_templatic];
    int<lower=1, upper=N_polarity_tense> polarity_tense_t[N_templatic];
    int<lower=1, upper=N_participant> participant_t[N_templatic];
    
    // contentful validation data
    int<lower=1> N_contentful;
    vector[N_contentful] y_c;
    int<lower=1, upper=N_scenario> scenario_c[N_contentful];
    int<lower=1, upper=N_verb> verb_c[N_contentful];
    int<lower=1, upper=N_polarity_tense> polarity_tense_c[N_contentful];
    int<lower=1, upper=N_participant> participant_c[N_contentful];
}

parameters {

    //
    // FIXED EFFECTS
    //

    // single fixed effect for Beta precision (should we do by polarity*tense?)
    real a_prec;

    // polarity*tense fixed effects for Beta mean
    vector[N_polarity_tense] a_mean;

    //
    // RANDOM EFFECTS AND HYPERPRIORS
    //

    // non-centered params for by-participant random intercepts
    // for the precision
    real<lower=0> tau_a_prec_p;
    vector[N_participant] z_a_prec_p;

    // non-centered params for random intercepts for the mean
    real<lower=0> tau_a_p;
    real<lower=0> tau_a_v;
    real<lower=0> tau_a_s;
    vector[N_participant] z_a_p;
    vector[N_verb] z_a_v;
    vector[N_scenario] z_a_s;

    matrix[N_polarity_tense, N_verb] z_B_v;
    cholesky_factor_corr[N_polarity_tense] L_Omega_B_v;

    matrix[N_polarity_tense, N_scenario] z_B_s;
    cholesky_factor_corr[N_polarity_tense] L_Omega_B_s;

    matrix[N_polarity_tense, N_participant] z_B_p;
    cholesky_factor_corr[N_polarity_tense] L_Omega_B_p;

    vector<lower = 0, upper = pi()/2>[N_polarity_tense] tau_B_v_unif;
    vector<lower = 0, upper = pi()/2>[N_polarity_tense] tau_B_s_unif;
    vector<lower = 0, upper = pi()/2>[N_polarity_tense] tau_B_p_unif;
}

transformed parameters {
    // beta parameters
    vector<lower=0>[N_norming] alpha_n;
    vector<lower=0>[N_norming] beta_n;
    vector<lower=0>[N_templatic] alpha_t;
    vector<lower=0>[N_templatic] beta_t;
    vector<lower=0>[N_contentful] alpha_c;
    vector<lower=0>[N_contentful] beta_c;

    // mean, precision parameters for betas
    vector<lower=0, upper=1>[N_norming] mu_n;
    vector<lower=0>[N_norming] prec_n;
    vector<lower=0, upper=1>[N_templatic] mu_t;
    vector<lower=0>[N_templatic] prec_t;
    vector<lower=0, upper=1>[N_contentful] mu_c;
    vector<lower=0>[N_contentful] prec_c;

    // by-verb random intercepts + slopes for polarity*tense
    // for the mean
    matrix[N_polarity_tense, N_verb] B_v;
    vector[N_verb] a_v;

    // by-scenario random intercepts + slopes for polarity*tense
    // for the mean
    matrix[N_polarity_tense, N_scenario] B_s;
    vector[N_scenario] a_s;

    // by-participant random intercepts + slopes for polarity*tense
    // for the mean
    matrix[N_polarity_tense, N_participant] B_p;
    vector[N_participant] a_p;

    // by-participant random intercepts for the precision
    vector[N_participant] a_prec_p;

    // scale vectors for by-verb and by-scenario random effects
    vector[N_polarity_tense] tau_B_v;
    vector[N_polarity_tense] tau_B_s;
    vector[N_polarity_tense] tau_B_p;

    tau_B_v = 2.5 * tan(tau_B_v_unif); // tau_B_v ~ cauchy(0, 2.5)
    tau_B_s = 2.5 * tan(tau_B_s_unif); // tau_B_s ~ cauchy(0, 2.5)
    tau_B_p = 2.5 * tan(tau_B_p_unif); // tau_B_p ~ cauchy(0, 2.5)

    // B_V ~ multi_normal(0, Sigma)
    // Sigma = diag(tau) * Omega * diag(tau)
    // Omega = L_Omega * L_Omega^T
    B_v = diag_pre_multiply(tau_B_v, L_Omega_B_v) * z_B_v; 
    B_s = diag_pre_multiply(tau_B_s, L_Omega_B_s) * z_B_s;
    B_p = diag_pre_multiply(tau_B_p, L_Omega_B_p) * z_B_p;

    a_prec_p = tau_a_prec_p * z_a_prec_p;

    a_v = tau_a_v * z_a_v;
    a_s = tau_a_s * z_a_s;
    a_p = tau_a_p * z_a_p;

    //
    // NORMING
    //
    mu_n = inv_logit(a_s[scenario_n] + a_p[participant_n]);
    prec_n = exp(a_prec + a_prec_p[participant_n]);
    alpha_n = mu_n .* prec_n;
    beta_n = (1 - mu_n) .* prec_n;

    //
    // TEMPLATIC VALIDATION
    //
    for (i in 1:N_templatic) {
        // - fixed effects for polarity*tense
        // - by-verb random intercepts + slopes for polarity*tense
        // - by-participant random intercepts + slopes for polarity*tense
        mu_t[i] = inv_logit(a_mean[polarity_tense_t[i]] +
                            a_v[verb_t[i]] +
                            B_v[polarity_tense_t[i], verb_t[i]] +
                            a_p[participant_t[i]] +
                            B_p[polarity_tense_t[i], participant_t[i]]);
        prec_t[i] = exp(a_prec + a_prec_p[participant_t[i]]);
    }
    alpha_t = mu_t .* prec_t;
    beta_t = (1 - mu_t) .* prec_t;

    //
    // CONTENTFUL VALIDATION
    //
    for (i in 1:N_contentful) {
        // - fixed effects for polarity*tense
        // - by-verb random intercepts + slopes for polarity*tense
        // - by-scenario random intercepts, slopes for polarity*tense
        // - by-participant random intercepts, slopes for polarity*tense
        mu_c[i] = inv_logit(a_mean[polarity_tense_c[i]] +
                            a_v[verb_c[i]] +
                            B_v[polarity_tense_c[i], verb_c[i]] + 
                            a_s[scenario_c[i]] +
                            B_s[polarity_tense_c[i], scenario_c[i]] + 
                            a_p[participant_c[i]] +
                            B_p[polarity_tense_c[i], participant_c[i]]);
        prec_c[i] = exp(a_prec + a_prec_p[participant_c[i]]);
    }
    alpha_c = mu_c .* prec_c;
    beta_c = (1 - mu_c) .* prec_c;
}


model {

    //
    // LIKELIHOODS
    //

    y_n ~ beta(alpha_n, beta_n);
    y_t ~ beta(alpha_t, beta_t);
    y_c ~ beta(alpha_c, beta_c);

    //
    // FIXED EFFECTS PRIORS
    //
    a_prec ~ normal(0,1);
    a_mean ~ normal(0,1);

    //
    // RANDOM EFFECTS PRIORS
    //
    // random intercepts for the mean
    z_a_p ~ normal(0,1);
    z_a_v ~ normal(0,1);
    z_a_s ~ normal(0,1);
    tau_a_p ~ exponential(1);
    tau_a_v ~ exponential(1);
    tau_a_s ~ exponential(1);

    // by-participant random intercepts for precision
    z_a_prec_p ~ normal(0,1);
    tau_a_prec_p ~ exponential(1);

    // priors for the Cholesky factorization of
    // covariance matrices for by-verb and by-scenario
    // random effects
    to_vector(z_B_s) ~ std_normal();
    L_Omega_B_s ~ lkj_corr_cholesky(2);

    to_vector(z_B_v) ~ std_normal();
    L_Omega_B_v ~ lkj_corr_cholesky(2);

    to_vector(z_B_p) ~ std_normal();
    L_Omega_B_p ~ lkj_corr_cholesky(2);
}

generated quantities {
    // log likelihoods (needed for WAIC/PSIS calculations)
    vector[N_norming] ll_n;
    vector[N_contentful] ll_c;
    vector[N_templatic] ll_t;
    vector[N_norming + N_contentful + N_templatic] ll_all;

    // covariance matrices for by-scenario and by-verb random effects
    cov_matrix[N_polarity_tense] Sigma_B_v;
    cov_matrix[N_polarity_tense] Sigma_B_s;
    cov_matrix[N_polarity_tense] Sigma_B_p;

    // means, sds of beta parameters
    // (both alphas + betas and means, precisions)
    real alpha_n_mean;
    real beta_n_mean;
    real alpha_c_mean;
    real beta_c_mean;
    real alpha_t_mean;
    real beta_t_mean;

    real mu_n_mean;
    real prec_n_mean;
    real mu_t_mean;
    real prec_t_mean;
    real mu_c_mean;
    real prec_c_mean;

    real alpha_n_sd;
    real beta_n_sd;
    real alpha_c_sd;
    real beta_c_sd;
    real alpha_t_sd;
    real beta_t_sd;

    real mu_n_sd;
    real prec_n_sd;
    real mu_c_sd;
    real prec_c_sd;
    real mu_t_sd;
    real prec_t_sd;

    Sigma_B_v = diag_pre_multiply(tau_B_v, L_Omega_B_v) * diag_pre_multiply(tau_B_v, L_Omega_B_v)';
    Sigma_B_s = diag_pre_multiply(tau_B_s, L_Omega_B_s) * diag_pre_multiply(tau_B_s, L_Omega_B_s)';
    Sigma_B_p = diag_pre_multiply(tau_B_p, L_Omega_B_p) * diag_pre_multiply(tau_B_p, L_Omega_B_p)';

    alpha_n_mean = mean(alpha_n);
    beta_n_mean = mean(beta_n);
    alpha_c_mean = mean(alpha_c);
    beta_c_mean = mean(beta_c);
    alpha_t_mean = mean(alpha_t);
    beta_t_mean = mean(beta_t);

    alpha_n_sd = sd(alpha_n);
    beta_n_sd = sd(beta_n);
    alpha_c_sd = sd(alpha_c);
    beta_c_sd = sd(beta_c);
    alpha_t_sd = sd(alpha_t);
    beta_t_sd = sd(beta_t);

    mu_n_mean = mean(mu_n);
    prec_n_mean = mean(prec_n);
    mu_c_mean = mean(mu_c);
    prec_c_mean = mean(prec_c);
    mu_t_mean = mean(mu_t);
    prec_t_mean = mean(prec_t);

    mu_n_sd = sd(mu_n);
    prec_n_sd = sd(prec_n);
    mu_c_sd = sd(mu_c);
    prec_c_sd = sd(prec_c);
    mu_t_sd = sd(mu_t);
    prec_t_sd = sd(prec_t);

    for (i in 1:N_norming) {
        ll_n[i] = beta_lpdf(y_n[i] | alpha_n[i], beta_n[i]);
    }
    for (i in 1:N_contentful) {
        ll_c[i] = beta_lpdf(y_c[i] | alpha_c[i], beta_c[i]);
    }
    for (i in 1:N_templatic) {
        ll_t[i] = beta_lpdf(y_t[i] | alpha_t[i], beta_t[i]);
    }
    ll_all = append_row(ll_n, append_row(ll_c, ll_t));
}

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
    matrix<lower=0, upper=1>[N_templatic, N_polarity_tense] polarity_tense_mat_t;
    int<lower=1, upper=N_participant> participant_t[N_templatic];
    
    // contentful validation data
    int<lower=1> N_contentful;
    vector[N_contentful] y_c;
    int<lower=1, upper=N_scenario> scenario_c[N_contentful];
    int<lower=1, upper=N_verb> verb_c[N_contentful];
    int<lower=1, upper=N_polarity_tense> polarity_tense_c[N_contentful];
    matrix<lower=0, upper=1>[N_contentful, N_polarity_tense] polarity_tense_mat_c;
    int<lower=1, upper=N_participant> participant_c[N_contentful];
}

parameters {

    //
    // FIXED EFFECTS
    //

    // polarity*tense fixed effects for Beta mean
    vector[N_polarity_tense] B_pt_mu;

    // single fixed effect for Beta precision
    real<lower=0> B0_prec;

    //
    // RANDOM EFFECTS PARAMETERS
    //

    // corr_matrix[N_polarity_tense] Omega_B_v;
    matrix[N_polarity_tense, N_verb] z_B_v;
    cholesky_factor_corr[N_polarity_tense] L_Omega_B_v;

    // corr_matrix[N_polarity_tense] Omega_B_s;
    matrix[N_polarity_tense, N_scenario] z_B_s;
    cholesky_factor_corr[N_polarity_tense] L_Omega_B_s;

    // vector[N_polarity_tense] tau_B_v;
    vector<lower = 0, upper = pi()/2>[N_polarity_tense] tau_B_v_unif;

    // vector[N_polarity_tense] tau_B_s;
    vector<lower = 0, upper = pi()/2>[N_polarity_tense] tau_B_s_unif;
}

transformed parameters {
    // beta parameters
    vector[N_norming] alpha_n;
    vector[N_norming] beta_n;
    vector[N_templatic] alpha_t;
    vector[N_templatic] beta_t;
    vector[N_contentful] alpha_c;
    vector[N_contentful] beta_c;

    // mean, precision parameters for betas
    vector[N_norming] mu_n;
    vector[N_templatic] mu_t;
    vector[N_contentful] mu_c;

    // by-verb random intercepts + slopes for polarity*tense
    // (intercept first)
    matrix[N_polarity_tense, N_verb] B_v;

    // by-scenario random intercepts + slopes for polarity*tense
    // (intercept first)
    matrix[N_polarity_tense, N_scenario] B_s;

    // scale vectors for by-verb and by-scenario random effects
    vector[N_polarity_tense] tau_B_v;
    vector[N_polarity_tense] tau_B_s;

    tau_B_v = 2.5 * tan(tau_B_v_unif); // tau_B_v ~ cauchy(0, 2.5)
    tau_B_s = 2.5 * tan(tau_B_s_unif); // tau_B_s ~ cauchy(0, 2.5)

    // B_V ~ multi_normal(0, Sigma)
    // Sigma = diag(tau) * Omega * diag(tau)
    // Omega = L_Omega * L_Omega^T
    B_v = diag_pre_multiply(tau_B_v, L_Omega_B_v) * z_B_v; 
    B_s = diag_pre_multiply(tau_B_s, L_Omega_B_s) * z_B_s;

    //
    // NORMING
    //
    for (i in 1:N_norming) {
        // - by-scenario random intercepts
        mu_n[i] = inv_logit(B_s[1, scenario_n[i]]);
    };
    alpha_n = mu_n * exp(B0_prec);
    beta_n = (1 - mu_n) * exp(B0_prec);

    //
    // TEMPLATIC VALIDATION
    //
    for (i in 1:N_templatic) {
        // - polarity*tense fixed effects
        // - by-verb random intercepts + slopes for polarity*tense
        mu_t[i] = inv_logit(B_pt_mu[polarity_tense_t[i]] + 
                            polarity_tense_mat_t[i, ] * B_v[, verb_t[i]]);
    }
    alpha_t = mu_t * exp(B0_prec);
    beta_t = (1 - mu_t) * exp(B0_prec);

    //
    // CONTENTFUL VALIDATION
    //
    for (i in 1:N_contentful) {
        // - polarity*tense fixed effects
        // - by-verb random intercepts + slopes for polarity*tense
        // - by-scenario random intercepts, slopes for polarity*tense
        mu_c[i] = inv_logit(B_pt_mu[polarity_tense_c[i]] + 
                            polarity_tense_mat_c[i, ] * B_v[, verb_c[i]] + 
                            polarity_tense_mat_c[i, ] * B_s[, scenario_c[i]]);
    }
    alpha_c = mu_c * exp(B0_prec);
    beta_c = (1 - mu_c) * exp(B0_prec);
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

    B_pt_mu ~ normal(0,1);
    B0_prec ~ normal(0,1);

    //
    // RANDOM EFFECTS PRIORS
    //

    // priors for the Cholesky factorization of
    // covariance matrices for by-verb and by-scenario
    // random effects
    to_vector(z_B_s) ~ std_normal();
    L_Omega_B_s ~ lkj_corr_cholesky(1);

    to_vector(z_B_v) ~ std_normal();
    L_Omega_B_v ~ lkj_corr_cholesky(1);
}

generated quantities {
    // covariance matrices for by-scenario and by-verb random effects
    cov_matrix[N_polarity_tense, N_polarity_tense] Sigma_B_v;
    cov_matrix[N_polarity_tense, N_polarity_tense] Sigma_B_s;

    // means, sds of beta parameters
    // (both alphas + betas and means)
    real alpha_n_mean;
    real beta_n_mean;
    real alpha_c_mean;
    real beta_c_mean;
    real alpha_t_mean;
    real beta_t_mean;

    real mu_n_mean;
    real mu_t_mean;
    real mu_c_mean;

    real alpha_n_sd;
    real beta_n_sd;
    real alpha_c_sd;
    real beta_c_sd;
    real alpha_t_sd;
    real beta_t_sd;

    real mu_n_sd;
    real mu_c_sd;
    real mu_t_sd;

    Sigma_B_v = diag_pre_multiply(tau_B_v, L_Omega_B_v) * diag_pre_multiply(tau_B_v, L_Omega_B_v)';
    Sigma_B_s = diag_pre_multiply(tau_B_s, L_Omega_B_s) * diag_pre_multiply(tau_B_s, L_Omega_B_s)';

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
    mu_c_mean = mean(mu_c);
    mu_t_mean = mean(mu_t);

    mu_n_sd = sd(mu_n);
    mu_c_sd = sd(mu_c);
    mu_t_sd = sd(mu_t);
}
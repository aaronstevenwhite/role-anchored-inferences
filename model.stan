data {
    // NOTE: role and participant have not yet been incorporated

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
    int<lower=0, upper=N_polarity_tense> polarity_tense_t[N_templatic];
    matrix<lower=0, upper=1>[N_templatic, N_polarity_tense] polarity_tense_mat_t;
    int<lower=1, upper=N_participant> participant_t[N_templatic];
    
    // contentful validation data
    int<lower=1> N_contentful;
    vector[N_contentful] y_c;
    int<lower=1, upper=N_scenario> scenario_c[N_contentful];
    int<lower=1, upper=N_verb> verb_c[N_contentful];
    int<lower=0, upper=N_polarity_tense> polarity_tense_c[N_templatic];
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
    real B0_prec;


    //
    // RANDOM EFFECTS
    //

    // by-verb random intercepts + slopes for polarity*tense
    // (intercept first)
    vector[N_polarity_tense] B_v[N_verb];

    // by-scenario random intercepts + slopes for polarity*tense
    // (intercept first)
    vector[N_polarity_tense] B_s[N_scenario];

    // by-participant random intercepts
    // (one each for mean and precision terms)
    // vector[N_participant] B0_p[2];

    // by-participant random slopes for polarity*tense
    // vector[N_polarity_tense] B_pt_p[N_verb];

    //
    // HYPERPRIOR PARAMETERS
    //

    corr_matrix[N_polarity_tense] Omega_B_v;
    vector[N_polarity_tense] tau_B_v;

    corr_matrix[N_polarity_tense] Omega_B_s;
    vector[N_polarity_tense] tau_B_s;
}

transformed parameters {

    // mean, precision parameters for betas
    vector[N_norming] mu_n;
    vector[N_norming] prec_n;
    vector[N_templatic] mu_t;
    vector[N_templatic] prec_t;
    vector[N_contentful] mu_c;
    vector[N_contentful] prec_c;

    //
    // NORMING
    //

    for (i in 1:N_norming) {
        // - by-scenario random intercept
        mu_n[i] = inv_logit(B_s[scenario_n[i]][1]);
    };
    prec_n = exp(rep_vector(B0_prec, N_norming));

    //
    // TEMPLATIC VALIDATION
    //
    for (i in 1:N_templatic) {
        // - polarity*tense fixed effects
        // - by-verb random intercepts, slopes for polarity*tense
        mu_t[i] = inv_logit(B_pt_mu[polarity_tense_t[i]] + 
                            polarity_tense_mat_t[i] * B_v[verb_t[i]]);
    }
    prec_t = exp(rep_vector(B0_prec, N_templatic));

    //
    // CONTENTFUL VALIDATION
    //
    for (i in 1:N_contentful) {
        // - polarity*tense fixed effects
        // - by-verb random intercepts, slopes for polarity*tense
        // - by-scenario random intercepts, slopes for polarity*tense
        mu_c[i] = inv_logit(B_pt_mu[polarity_tense_c[i]] + 
                            polarity_tense_mat_c[i] * B_v[verb_c[i]] + 
                            polarity_tense_mat_c[i] * B_s[scenario_c[i]]);
    }
    prec_c = exp(rep_vector(B0_prec, N_contentful));
}


model {

    //
    // LIKELIHOODS
    //

    y_n ~ beta(mu_n .* prec_n, (1 - mu_n) .* prec_n);
    y_t ~ beta(mu_t .* prec_t, (1 - mu_t) .* prec_t);
    y_c ~ beta(mu_c .* prec_c, (1 - mu_c) .* prec_c);

    //
    // FIXED EFFECTS PRIORS
    //

    B_pt_mu ~ normal(0,10);
    B0_prec ~ normal(0,10);

    //
    // RANDOM EFFECTS PRIORS
    //

    // We enforce that by-scenario and by-verb random ints+slopes be zero-centered
    B_s ~ multi_normal(rep_row_vector(0, N_polarity_tense), quad_form_diag(Omega_B_v, tau_B_v));
    B_v ~ multi_normal(rep_row_vector(0, N_polarity_tense), quad_form_diag(Omega_B_s, tau_B_s));

    //
    // RANDOM EFFECTS HYPERPRIORS
    // (for covariance)
    //

    tau_B_v ~ cauchy(0, 2.5);
    Omega_B_v ~ lkj_corr(1);

    tau_B_s ~ cauchy(0, 2.5);
    Omega_B_s ~ lkj_corr(1);
}


# Author: Will Gantt
# 

# Load RStan
library(rstan);
library(loo);
options(mc.cores = parallel::detectCores());
rstan_options(auto_write = TRUE);

# NOTE: you must update these before running
REPO_PATH <- "/home/wgantt/role-anchored-inferences";
OUTPUT_PATH <- "/data/wgantt/mi/stan_models/desire/random-slopes";

PREPROCESS_SCRIPT <- file.path(REPO_PATH, "preprocessing/preprocess.r");
MODEL_NAME <- file.path(REPO_PATH, "model.participant.stan");
MODEL_SAVE_PATH <- file.path(OUTPUT_PATH, "desire-slopes.rds");
DIAGNOSTIC_FILE_PATH <- paste(MODEL_NAME, "desire.diagnostic", sep=".");

# Preprocess the data
source(PREPROCESS_SCRIPT);

# Fit the model
data <- data.w
fit <- stan(file = MODEL_NAME,
            data = data,
            open_progress = TRUE,
            refresh = 20,
            verbose = TRUE,
            diagnostic_file = DIAGNOSTIC_FILE_PATH,
            seed = 1337,
            iter = 12000
);

# Save the model
saveRDS(fit, MODEL_SAVE_PATH);

# Extract the pointwise log-likelihood and compute WAIC
ll_all <- extract_log_lik(fit, parameter_name = "ll_all");
fit.waic <- waic(ll_all);

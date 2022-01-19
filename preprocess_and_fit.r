# Author: Will Gantt
# 

# MegaIntensionality CSVs after preprocessing
w_n <- read.csv("want_norming.csv");
w_t <- read.csv("want_templatic.csv");
w_c <- read.csv("want_contentful.csv");

b_n <- read.csv("believe_norming.csv");
b_t <- read.csv("believe_templatic.csv");
b_c <- read.csv("believe_contentful.csv");

# Global data
N_verb <- length(unique(w_n$verb));
N_scenario <- length(unique(w_n$scenario));
N_polarity_tense <- length(unique(w_c$polarity)) * length(unique(w_c$tense));
N_participant <- length(unique(w_n$participant)) + length(unique(w_c$participant)) + length(unique(w_t$participant));
N_norming <- nrow(w_n);
N_templatic <- nrow(w_t);
N_contentful <- nrow(w_c);

# Norming data
y_n.w <- w_n$response;
scenario_n.w <- as.integer(as.factor(w_n$scenario));
participant_n.w <- as.integer(as.factor(w_n$participant));


# Templatic validation data
y_t.w <- w_t$response;
verb_t.w <- as.integer(as.factor(w_t$verb));
w_t$polarity.tense <- relevel(interaction(w_t$polarity, w_t$tense), ref="positive.past");
polarity_tense_t.w <- as.integer(w_t$polarity.tense);
polarity_tense_mat_t.w <- model.matrix(~ 1 + polarity.tense, data=w_t);
participant_t.w <- as.integer(as.factor(w_t$participant));


# Contenful validation data
y_c.w <- w_c$response;
scenario_c.w <- as.integer(as.factor(w_c$scenario));
verb_c.w <- as.integer(as.factor(w_c$verb));
w_c$polarity.tense <- relevel(interaction(w_c$polarity, w_c$tense), ref="positive.past");
polarity_tense_c.w <- as.integer(w_c$polarity.tense);
polarity_tense_mat_c.w <- model.matrix(~ 1 + polarity.tense, data=w_c);
participant_c.w <- as.integer(as.factor(w_c$participant));

data.w <- list(
  N_verb=N_verb,
  N_scenario=N_scenario,
  N_polarity_tense=N_polarity_tense,
  N_participant=N_participant,
  N_norming=N_norming,
  N_templatic=N_templatic,
  N_contentful=N_contentful,
  y_n=y_n.w,
  scenario_n=scenario_n.w,
  participant_n=participant_n.w,
  y_t=y_t.w,
  verb_t=verb_t.w,
  polarity_tense_t=polarity_tense_t.w,
  polarity_tense_mat_t=polarity_tense_mat_t.w,
  participant_t=participant_t.w,
  y_c=y_c.w,
  scenario_c=scenario_c.w,
  verb_c=verb_c.w,
  polarity_tense_c=polarity_tense_c.w,
  polarity_tense_mat_c=polarity_tense_mat_c.w,
  participant_c=participant_c.w
);

# Load RStan
library(rstan);
options(mc.cores = parallel::detectCores());
rstan_options(auto_write = TRUE);

# Fit the model
model_name <- "model.stan";
fit <- stan(file = model_name,
            data = data.w,
            open_progress = TRUE,
            refresh = 20,
            verbose = TRUE,
            diagnostic_file = paste(model_name, ".diagnostic", sep=""),
            seed = 1337,
            iter = 5000
);

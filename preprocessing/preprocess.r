# Author: Will Gantt
# 

# MegaIntensionality CSVs after preprocessing
# NOTE: you will need to change these paths
w_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_norming.csv");
w_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_templatic.csv");
w_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_contentful.csv");

b_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_norming.csv");
b_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_templatic.csv");
b_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_contentful.csv");

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

y_n.b <- b_n$response;
scenario_n.b <- as.integer(as.factor(b_n$scenario));
participant_n.b <- as.integer(as.factor(b_n$participant));



# Templatic validation data
y_t.w <- w_t$response;
verb_t.w <- as.integer(as.factor(w_t$verb));
w_t$polarity.tense <- relevel(interaction(w_t$polarity, w_t$tense), ref="positive.past");
polarity_tense_t.w <- as.integer(w_t$polarity.tense);
polarity_tense_mat_t.w <- model.matrix(~ 1 + polarity*tense, data=w_t);
participant_t.w <- as.integer(as.factor(w_t$participant));

y_t.b <- b_t$response;
verb_t.b <- as.integer(as.factor(b_t$verb));
b_t$polarity.tense <- relevel(interaction(b_t$polarity, b_t$tense), ref="positive.past");
polarity_tense_t.b <- as.integer(b_t$polarity.tense);
polarity_tense_mat_t.b <- model.matrix(~ 1 + polarity*tense, data=b_t);
participant_t.b <- as.integer(as.factor(b_t$participant));


# Contenful validation data
y_c.w <- w_c$response;
scenario_c.w <- as.integer(as.factor(w_c$scenario));
verb_c.w <- as.integer(as.factor(w_c$verb));
w_c$polarity.tense <- relevel(interaction(w_c$polarity, w_c$tense), ref="positive.past");
polarity_tense_c.w <- as.integer(w_c$polarity.tense);
polarity_tense_mat_c.w <- model.matrix(~ 1 + polarity*tense, data=w_c);
participant_c.w <- as.integer(as.factor(w_c$participant));

y_c.b <- b_c$response;
scenario_c.b <- as.integer(as.factor(b_c$scenario));
verb_c.b <- as.integer(as.factor(b_c$verb));
b_c$polarity.tense <- relevel(interaction(b_c$polarity, b_c$tense), ref="positive.past");
polarity_tense_c.b <- as.integer(b_c$polarity.tense);
polarity_tense_mat_c.b <- model.matrix(~ 1 + polarity*tense, data=b_c);
participant_c.b <- as.integer(as.factor(b_c$participant));

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

data.b <- list(
  N_verb=N_verb,
  N_scenario=N_scenario,
  N_polarity_tense=N_polarity_tense,
  N_participant=N_participant,
  N_norming=N_norming,
  N_templatic=N_templatic,
  N_contentful=N_contentful,
  y_n=y_n.b,
  scenario_n=scenario_n.b,
  participant_n=participant_n.b,
  y_t=y_t.b,
  verb_t=verb_t.b,
  polarity_tense_t=polarity_tense_t.b,
  polarity_tense_mat_t=polarity_tense_mat_t.b,
  participant_t=participant_t.b,
  y_c=y_c.b,
  scenario_c=scenario_c.b,
  verb_c=verb_c.b,
  polarity_tense_c=polarity_tense_c.b,
  polarity_tense_mat_c=polarity_tense_mat_c.b,
  participant_c=participant_c.b
);

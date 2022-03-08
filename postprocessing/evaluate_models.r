library(loo)
library(assertthat)

# Path will have to be changed to match local repo
source('/home/wgantt/role-anchored-inferences/preprocessing/preprocess.r');

w_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_norming.csv");
w_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_templatic.csv");
w_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_contentful.csv");

b_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_norming.csv");
b_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_templatic.csv");
b_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_contentful.csv");

# Computes WAIC for each level of a categorical predictor (verb or scenario)
# for a fitted Stan model from this repository
by_predictor_waic <- function(fit, inference_type, predictor, output_csv_root=NULL, output_path=NULL) {
	# Extract pointwise LL matrices for norming data and for
	# contentful and templatic validation data
	ll_n <- extract_log_lik(fit, parameter_name = "ll_n");
	ll_c <- extract_log_lik(fit, parameter_name = "ll_c");
	ll_t <- extract_log_lik(fit, parameter_name = "ll_t");

	# Are we working with belief or desire inferences?
	if (inference_type == "w") {
		data <- data.w
	} else if (inference_type == "b") {
		data <- data.b
	} else {
		stop("Unrecognized inference type!");
	}

	# Do we want to compute WAIC for scenarios or for verbs?
	if (predictor == "scenario") {
		predictor_n <- data$scenario_n;
		predictor_c <- data$scenario_c;
		predictor_t <- NULL;
		N_levels <- data$N_scenario;
		# arbitrarily take level names from contentful desire data
		# (levels are the same across datasets)
		level <- levels(as.factor(w_c$scenario));
	} else if (predictor == "verb") {
		predictor_n <- NULL;
		predictor_c <- data$verb_c;
		predictor_t <- data$verb_t;
		N_levels <- data$N_verb;
		level <- levels(as.factor(w_c$verb));
	} else {
		stop("Unrecognized predictor name!");
	}

	# Compute WAIC for each level of the predictor
	# (scenario or verb)
	waic_n <- NULL;
	waic_c <- NULL;
	waic_t <- NULL;
	for (l in 1:N_levels) {
		# The contentful validation data includes both verb and
		# scenario, so whichever predictor we're concerned with,
		# we can compute WAIC for it on the contentful data
		waic_c <- rbind(waic_c, waic(ll_c[ , predictor_c == l]));

		# The norming data has only scenario and lacks verb
		if (!is.null(predictor_n)) {
			waic_n <- rbind(waic_n, waic(ll_n[ , predictor_n == l]));
		}

		# The templatic data has only verb and lacks scenario
		if (!is.null(predictor_t)) {
			waic_t <- rbind(waic_t, waic(ll_t[ , predictor_t == l]));
		}
	}

	# If a root name for the ouptut files was provided,
	# create them
	if (!is.null(output_csv_root)) {
		if (is.null(output_path)) {
			stop("No output path specified!");
		} else if (!is.dir(output_path)) {
			stop("Output directory does not exist!");
		}

		# strip off "estimates" and "pointwise" columns
		# before generating output dataframes
		waic_c_df <- waic_c[ , -c(1,2)];
		waic_c_df <- cbind(level, waic_c_df);
		output_file_name <- paste(output_csv_root, inference_type, "contentful", predictor, "waic", "csv", sep=".");
		write.csv(waic_c_df, file.path(output_path, output_file_name), row.names=FALSE);
		if (!is.null(waic_n)) {
			waic_n_df <- waic_n[ , -c(1,2)];
			waic_n_df <- cbind(level, waic_n_df);
			output_file_name <- paste(output_csv_root, inference_type, "norming", predictor, "waic", "csv", sep=".");
			write.csv(waic_n_df, file.path(output_path, output_file_name), row.names=FALSE);
		}
		if (!is.null(waic_t)) {
			waic_t_df <- waic_t[ , -c(1,2)];
			waic_t_df <- cbind(level, waic_t_df);
			output_file_name <- paste(output_csv_root, inference_type, "templatic", predictor, "waic", "csv", sep=".");
			write.csv(waic_t_df, file.path(output_path, output_file_name), row.names=FALSE);
		}
	}

	# Return vectors of WAICs per predictor level, per dataset
	return(list(norming=waic_n, contentful=waic_c, templatic=waic_t));
}

# Computes either WAIC or LOO-CV-PSIS for (verb, polarity, tense) triples
# on the templatic and contentful validation data (norming data is homogenous
# with respect to tense and polarity, and so is excluded)
by_verb_polarity_tense <- function(fit, inference_type, output_csv_root=NULL, output_path=NULL, func="waic") {
	# Extract pointwise LL matrices for
	# contentful and templatic validation data
	ll_c <- extract_log_lik(fit, parameter_name = "ll_c");
	ll_t <- extract_log_lik(fit, parameter_name = "ll_t");

	# What criterion should we use for model evaluation: WAIC or LOO-CV-PSIS?
	if (func == "waic") {
		print("Using WAIC to evaluate the fit.")
		criterion <- waic;
	} else if (func == "psis") {
		print("Using PSIS to evaluate the fit. This may take a while.")
		criterion <- loo;
		# For PSIS, we need to specify the relative effective
		# sample size of each observation; otherwise, Stan will
		# issue warnings about overly optimistic error estimates
		ll_c_all_chains <- extract_log_lik(fit, parameter_name = "ll_c", merge_chains = FALSE);
		ll_t_all_chains <- extract_log_lik(fit, parameter_name = "ll_t", merge_chains = FALSE);
		r_eff_c <- relative_eff(ll_c_all_chains);
		r_eff_t <- relative_eff(ll_t_all_chains);

	} else {
		stop("Unrecognized evaluation method! Valid options are 'waic' and 'psis'.");
	}

	# Are we working with belief or desire inferences?
	if (inference_type == "w") {
		data_c <- w_c;
		data_t <- w_t;
	} else if (inference_type == "b") {
		data_c <- b_c;
		data_t <- b_t;
	} else {
		stop("Unrecognized inference type! Valid options are 'w' and 'b'");
	}

	verb.levels <- levels(as.factor(data_c$verb));
	polarity.levels <- levels(as.factor(data_c$polarity));
	tense.levels <- levels(as.factor(data_c$tense));

	# Compute criterion (WAIC or LOO-CV-PSIS) for each
	# verb-polarity-tense tuple for both templatic and
	# contentful validations.
	criterion_c <- NULL;
	criterion_t <- NULL;
	for (verb in verb.levels) {
		verb.mask.c <- (data_c$verb == verb);
		verb.mask.t <- (data_t$verb == verb);
		for (polarity in polarity.levels) {
			polarity.mask.c <- (data_c$polarity == polarity);
			polarity.mask.t <- (data_t$polarity == polarity);
			for (tense in tense.levels) {
				# Compute the criterion for the current split for
				# templatic and contentful validation data
				tense.mask.c <- (data_c$tense == tense);
				tense.mask.t <- (data_t$tense == tense);
				mask.c <- (verb.mask.c & polarity.mask.c & tense.mask.c);
				mask.t <- (verb.mask.t & polarity.mask.t & tense.mask.t);
				if (func == "waic") {
					split.criterion.c <- criterion(ll_c[, mask.c]);
					split.criterion.t <- criterion(ll_t[, mask.t]);
				} else {
					split.criterion.c <- criterion(ll_c[, mask.c], r_eff=r_eff_c[mask.c]);
					split.criterion.t <- criterion(ll_t[, mask.t], r_eff=r_eff_t[mask.t]);
				}
				split.criterion.c <- append(c(verb=verb, polarity=polarity, tense=tense), split.criterion.c);
				split.criterion.t <- append(c(verb=verb, polarity=polarity, tense=tense), split.criterion.t);
				criterion_c <- rbind(criterion_c, split.criterion.c);
				criterion_t <- rbind(criterion_t, split.criterion.t);
			}
		}
	}

	# If a root name for the ouptut files was provided, create them
	if (!is.null(output_csv_root)) {
		if (is.null(output_path)) {
			stop("No output path specified!");
		} else if (!is.dir(output_path)) {
			stop("Output directory does not exist!");
		}

		if (func == "waic") {
			# strip off "estimates" and "pointwise" columns
			# from WAIC output before generating final dataframes
			# (this information is overly fine-grained and not what we need)
			criterion_c_df <- criterion_c[ , -c(4,5)];
			criterion_t_df <- criterion_t[ , -c(4,5)];
		} else {
			# additionally remove "psis object" and "diagnostics" columns
			criterion_c_df <- criterion_c[ , -c(4,5,6,7)];
			criterion_t_df <- criterion_t[ , -c(4,5,6,7)];
		}
		output_file_name_c <- paste(output_csv_root, inference_type, "contentful", func, "csv", sep=".");
		output_file_name_t <- paste(output_csv_root, inference_type, "templatic", func, "csv", sep=".");
		write.csv(criterion_c_df, file.path(output_path, output_file_name_c), row.names=FALSE);
		write.csv(criterion_t_df, file.path(output_path, output_file_name_t), row.names=FALSE);
	}

	# Return vectors of WAICs/PSIS estimates per predictor level, per dataset
	return(list(contentful=criterion_c, templatic=criterion_t));

}

# Runs by_predictor_waic for all predictors.
# Still have to specify inference type, since models
# are fit to *either* belief inferences *or* desire inferences.
# Currently doesn't actually return anything; simply assumes
# you will be writing to CSVs.
run_all <- function(fit, inference_type, output_csv_root=NULL, output_path=NULL) {
	predictors <- c("verb", "scenario");
	for (p in predictors) {
		by_predictor_waic(fit, inference_type, p, output_csv_root, output_path);
	}
}

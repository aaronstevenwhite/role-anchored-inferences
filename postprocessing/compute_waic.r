library(loo)
library(assertthat)

# Path will have to be changed to match local repo
source('/home/wgantt/role-anchored-inferences/preprocessing/preprocess.r');

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

library(loo)

w_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_norming.csv");
w_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_templatic.csv");
w_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/want_contentful.csv");

b_n <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_norming.csv");
b_t <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_templatic.csv");
b_c <- read.csv("/home/wgantt/role-anchored-inferences/data/preprocessed/believe_contentful.csv");

dump_likelihoods <- function(fit, inference_type, output_path, output_csv_root) {
	"
	Dumps likelihoods from a fitted Stan model to CSVs.

	parameters:
	-----------
	fit: the model fit for which the likelihoods are to be dumped
	inference_type: the type of inference for which this model was fit
	  (use 'w' for want/desire; 'b' for belief)
	output_path: the directory to which the output CSVs will be written
	output_csv_root: a base filename to use for the output CSVs

	returns:
	--------
	Nothing. Writes one CSV each for norming, templatic, and contentful data 
          to disk at the location specified.
	"
	# Extract num_examples x num_samples likelihood matrices
	# for each data type (norming, contentful, templatic)
	ll_n_df <- as.data.frame(t(extract_log_lik(fit, parameter_name = "ll_n")));
	ll_c_df <- as.data.frame(t(extract_log_lik(fit, parameter_name = "ll_c")));
	ll_t_df <- as.data.frame(t(extract_log_lik(fit, parameter_name = "ll_t")));

	# Column names for the likelihood samples
	format_col_name <- function(x) paste("sample", x);
	ll_n_col_names <- sapply(1:ncol(ll_n_df), format_col_name);
	ll_c_col_names <- sapply(1:ncol(ll_c_df), format_col_name);
	ll_t_col_names <- sapply(1:ncol(ll_t_df), format_col_name);
	names(ll_n_df) <- ll_n_col_names;
	names(ll_c_df) <- ll_c_col_names;
	names(ll_t_df) <- ll_t_col_names;

	# Load the relevant data
	if (inference_type == "w") {
		data_n <- w_n;
		data_t <- w_t;
		data_c <- w_c;
	} else if (inference_type == "b") {
		data_n <- b_n;
		data_t <- b_t;
		data_c <- b_c;
	} else {
		stop("Unrecognized inference type!");
	}

	# Join likelihoods for each sample to corresponding dataframe
	ll_n_df <- cbind(data_n, ll_n_df);
	ll_c_df <- cbind(data_c, ll_c_df);
	ll_t_df <- cbind(data_t, ll_t_df);

	# Output the augmented dataframes
	output_ll_n_file_name <- paste(output_csv_root, "norming.csv", sep="-");
	output_ll_c_file_name <- paste(output_csv_root, "contentful.csv", sep="-");
	output_ll_t_file_name <- paste(output_csv_root, "templatic.csv", sep="-");

	write.csv(ll_n_df, file.path(output_path, output_ll_n_file_name), row.names=FALSE);
	write.csv(ll_c_df, file.path(output_path, output_ll_c_file_name), row.names=FALSE);
	write.csv(ll_t_df, file.path(output_path, output_ll_t_file_name), row.names=FALSE);
}

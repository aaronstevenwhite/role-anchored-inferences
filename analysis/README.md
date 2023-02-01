# Model Analysis

This directory contains things related to model analysis. Right now, this just includes:

1. CSVs containing per-(verb,polarity,tense) WAIC estimates for the contentful and templatic validation data for both belief and desire inferences (`{belief,desire}/*.waic.csv`)
2. CSVs containing per-(verb,polarity,tense) estimates of WAIC *differences* between random slopes and random intercepts models on this same data (`{belief,desire}/*.comparison.csv`). This difference is expressed by three columns:
  - `best_model`: for a given tuple, the model that performed best on the examples for that tuple, as judged by WAIC evaluated on those examples only (either "intercepts" or "slopes").
  - `elpd_diff`: point estimate of the *difference* in expected log
     pointwise predictive density ([ELPD](https://mc-stan.org/loo/reference/loo-glossary.html)) between the two mdoels ("slopes" and "intercepts"), as
     computed by the [loo_compare](https://mc-stan.org/loo/reference/loo_compare) function from the `loo` package.
  - `se_elpd_diff`: The standard error of `elpd_diff`, also output by `loo_compare`.

The functions for producing both types of CSVs can be found in `postprocessing/evaluate_models.r`.

Data from further analyses, such as those focused on the posteriors for specific parameters, may be added later.

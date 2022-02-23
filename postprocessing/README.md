# Postprocessing and Analysis

This directory contains:
- CSVs of samples for the by-verb and by-scenario random slopes for polarity and tense, for both the belief and desire models. These slopes are dummy coded, so the posterior samples have to be post-processed to collect statistics on the true value associated with each slope. That's the primary purpose of the 'Clean Stan Model Outputs' Jupyter Notebook
- CSVs of the raw summary statistics for the by-verb and by-scenario random slopes. The 'Clean Stan Model Outputs' notebook does some postprocessing on these as well, though these are less interesting, as they are statistics about the raw dummy coded parameter values, not the true slopes themselves.
- A silly png of the contrasts, which I consulted while doing the post-processing

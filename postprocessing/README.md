# Postprocessing and Analysis

This directory contains:
- `evaluate_models.r`: functions for computing information criteria at different granularities, and for comparing
   models on the basis of those criteria. This is currently the only relevant file in this directory. The rest is
   obsolete.
- `Clean Stan Model Outputs.ipynb`: Previously, predictors were dummy coded rather than index coded, so the posterior samples had to be post-processed to collect statistics on the true value associated with each slope. That was the primary purpose of the 'Clean Stan Model Outputs' Jupyter Notebook, but the latest models have been index-coded, so this can be ignored.

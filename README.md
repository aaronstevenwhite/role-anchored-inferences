# Role-Anchored Doxastic and Bouletic Inference Models

This repository contains code for generalized linear models of the doxastic ("belief") and bouletic ("desire") inference judgment data in the [MegaIntensionality dataset](http://megaattitude.io/projects/mega-intensionality/), part of the [MegaAttitude project](http://megaattitude.io/). The models are implemented in Stan and invoked via the RStan R package. Preprocessing code is in a mix of Python (in a Jupyter Notebook) and R. See [below](#workflow) for instructions on how to fit the models.

## Workflow

The general preprocessing and model fitting workflow is as follows:

1. Generate pre-processed CSVs from the raw MegaIntensionality data by running the MegaIntensionality Preprocessing Notebook.
2. Select the model you want to fit. There are four, which differ from each other depending on whether: 1) they are fit to doxastic inference judgments or to bouletic inference judgments, and 2) they use by-scenario random slopes (for all combinations of polarity and tense of the antecedent frame) or just by-scenario random intercepts. Scripts for running each of these models are provided in the root of this repository (the `fit-*.r` files).
3. Start a new session in the R terminal (preferably inside a terminal multiplexer, such as [tmux](https://github.com/tmux/tmux) or [screen](https://www.gnu.org/software/screen/), as model fitting will take a long time).
4. Assuming you are running the R terminal from the root of the repository, run `source(<model_script>)`, where `<model_script>` is the R script for fitting the model you selected in step 3. Note that this will handle all further preprocessing, assuming you have run (1).

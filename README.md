# Role-Anchored Doxastic and Bouletic Inference Models

This repository contains code for generalized linear models of the doxastic and bouletic inference judgment data in the [MegaIntensionality dataset](http://megaattitude.io/projects/mega-intensionality/), part of the [MegaAttitude project](http://megaattitude.io/). The models are implemented in Stan and invoked via the RStan R package. Preprocessing code is in a mix of Python (in a Jupyter Notebook) and R. See [below](#workflow) for instructions on how to fit the models.

## Workflow

The general preprocessing and model fitting workflow is as follows:

1. Generate pre-processed CSVs from the raw MegaIntensionality data by running the MegaIntensionality Preprocessing Notebook.
2. Start a new session in the R terminal (preferably inside a terminal multiplexer, such as [tmux](https://github.com/tmux/tmux) or [screen](https://www.gnu.org/software/screen/), as model fitting will take a long time).
3. From the R terminal, source `preprocess_and_fit.r`, which performs further transformations on the preprocessed data and fits the model.

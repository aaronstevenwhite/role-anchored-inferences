# Role-Anchored Doxastic and Bouletic Inference Models

This repository contains code for mixed effects models of the doxastic ("belief") and bouletic ("desire") inference judgment data in the [MegaIntensionality dataset](http://megaattitude.io/projects/mega-intensionality/), part of the [MegaAttitude project](http://megaattitude.io/). The models are implemented in Stan and invoked via the RStan R package. Preprocessing code is in a mix of Python (in a Jupyter Notebook) and R. See [below](#workflow) for instructions on how to fit the models.

## Workflow

The general preprocessing and model fitting workflow is as follows:

1. Generate pre-processed CSVs from the raw MegaIntensionality norming and validation data by running the `MegaIntensionality Preprocessing.ipynb` Notebook in the `preprocessing` directory. (NOTE: you should *not* have reason to run `preprocess.r`; this is a separate script used by the model fitting scripts below.) The MegaItensionality data is stored as zip files in `data/raw` in this repo, and these must of course be unzipped before use. Alternatively, you can use the already-preprocessed versions of the files in `data/preprocessed`.
2. Select the model you want to fit. There are eight, which differ from each other depending on whether: 1) they are fit to doxastic inference judgments or to bouletic inference judgments, 2) they use by-scenario random slopes (for all combinations of polarity and tense of the antecedent frame) or just by-scenario random intercepts; and 3) the categorical variable coding used: a dummy coding (k-1 values for k levels) or an index coding (k values for k levels). Scripts for running each of these models are provided in the root of this repository (the `fit-*.r` files). The latest models use index coding, and I have only been testing with these models lately, so you are strongly encouraged to use them (the `fit-*-index-code.r` scripts). Below, I indicate for each script the raw Stan file (`.stan`) containing the model code that was used to fit it (you can also find this information by inspecting the script itself). Apologies for the cumbersome naming of the Stan files; it's a holdover from earlier stages of development that I haven't bothered to change.
    - `fit-{belief,desire}-slopes.r`: `model.participant.stan`
    - `fit-{belief,desire}-intercepts.r`: `model.participant.scenario.int.stan`
	- `fit-{belief,desire}-slopes-index-code.r`: `model.participant.index.code.stan`
	- `fit-{belief,desire}-intercepts-index-code.r`: `model.participant.scenario.int.index.code.stan`
3. Start a new session in the R terminal. Unless debugging, you should almost certainly do this inside a terminal multiplexer, such as [tmux](https://github.com/tmux/tmux) or [screen](https://www.gnu.org/software/screen/), as model fitting will take a long time.
4. Inside the script you have chosen, you will need to modify the variables for the paths to (1) your local copy of this repo (`REPO_PATH`), and (2) the path to the directory where you want the model to be saved (`OUTPUT_PATH`). Ideally, these would be command line arguments, but they are not right now.
5. Assuming you are running the R terminal from the root of the repository, run `source(<model_script>)`, where `<model_script>` is the R script for fitting the model you selected in step 2. Note that this will handle all further preprocessing, assuming you have run (1). Each of these scripts will save pointwise log likelihood estimates for all observations from the norming data, as well as from the contentful and templatic validation data. They also compute WAIC for the fit across all three datasets on the basis of these log likelihoods. Any further analysis you must do yourself, after the model has been fit.

## Saved models and samples

If you are working on the FACTSlab machines, you can find the latest model fits in `/data/wgantt/mi/stan_models/{belief,desire}/random-{slopes,intercepts}-index-code/`. These are *very large* files and take a non-trivial amount of time to load. This can be done with the [`readRDS`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS) function (e.g. `saved_fit <- readRDS('/path/to/saved/fit.rds'`)>).

For each model, I have also saved the likelihood estimates for each example in the norming, contentful, and templatic data across all (post-warmup) posterior samples (merged across chains) in gzipped archives in these same directories. Each archive contains one CSV per dataset (norming, contentful, templatic). These too are large files (the contentful and templatic CSVs are about 1.6G each when extracted), since these models all had 20,000 post-warmup draws, and there is a column for each draw, for each example. These CSVs were created using `postprocessing/dump_likelihoods.r`.

## Other Notes

Currently, all the analysis code is focused on computing information criteria (ICs) at various granularities and for comparing models based on ICs. Presently, there is no code for analysis of model parameters or for plotting, which are shortcomings that need to be remedied.

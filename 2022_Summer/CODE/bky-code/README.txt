Source Code for Blevins, Khwaja, and Yang (2016)
================================================

Estimation
----------

List of main files for estimation:

* canadafastfood.mat - Matlab data file
* estimation.sh - Shell script for controlling estimation
* specX_first.m - First stage estimation for specification X (1, 2, or 3)
* specX_simulate.m - Forward simulation for specification X (1, 2, or 3)
* specX_second.m - Second stage estimation for specification X (1, 2, or 3)
* specX_bootstrap.m - Bootstrap code for specification X (1, 2, or 3)
* results - directory where output files will be stored

To estimate specification X, run the following scripts in order:

1. specX_first
2. specX_simulate
3. specX_second

Then, to estimate bootstrap standard errors, run `specX_bootstrap`.
Alternatively, run the shell script `estimation.sh` to carry out all
four steps (including bootstrap standard errors) for all specifications.


Model Fit and Counterfactuals
-----------------------------

List of main files for model fit and counterfactual simulations:

* results/specX_param.mat - Estimated parameters for specification X.
* model_fit.sh - Shell script for carrying out model fit simulations.
* counterfactual.sh - Shell script for carrying out counterfactual simulations.

Run the appropriate shell scripts to carry out the model fit and
counterfactual simulations, produce graphs, etc.  This relies on being
able to run Stata and Matlab from the command line.  The
`model_fit.sh` script generates all of the subfigures that constitute
Figure 3 in the paper.  The `counterfactutal.sh` script produces the
counterfactual results in Tables 8, 9, and 10 in the paper.

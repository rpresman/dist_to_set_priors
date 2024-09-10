# Distance-to-Set Priors and Bayesian Constraint Relaxation

This is the repo with supporting code for the paper ["Distance-to-Set Priors and Bayesian Constraint Relaxation"](https://proceedings.mlr.press/v206/presman23a/presman23a.pdf) (Presman and Xu, 2023), presented at AISTATS 2023. This README will provide an overview of what files are located in the repo for each experiment.

## Regression over the $\ell_2$-ball

The following files contain supporting `stan` code:

* `hmc1_squared.stan` contains the squared distance-to-set prior
* `hmc2_linear.stan` contains the unsquared distance-to-set prior

The results in this example can be obtained by running the code found in the file `optim_ridge.R`.

## Sampling a Lower-Dimensional Surface

The following files contain support `stan` code:

* `relaxed_sqd_vmf_sampling.stan` contains the distance-to-set prior
* `relaxed_vmf_sampling.stan` contains the constraint relaxatio method of Duan et al. (2022)

The results in this example can be obtained by running the code found in the file `robust_vmf.R`.

## Real Data Case Study

The results in this example can be obtained by running the code found in the file `stoch_dominance.R`. This file contains all supporting functions, including the HMC sampler, along with the data needed to run the this example.



# bayesian-inference

Jupyter Notebooks for the publication Bayesian Inference for Integrating Yarrowia lipolytica Multi-omics Datasets with Metabolic Modeling

## What is in this Repository?

This repository contains all the necessary Jupyter Notebooks and code to replicate the results obtained in our paper. This also serves as a jumping off point for expanding the work done here to your own models if necessary.

## What you will need

A conda environment is provided in this repository for the sake of easier reproducability. There are not many major packages required to run the code, but the code requires Ensemble Modeling with Linear-Logarithmic Kinetics (emll) [source](https://github.com/pstjohn/emll) which you can either install independently or download through the given environment.

### Conda Environment

To install the working environment type the following commands

```
conda env create -f environment.yml
```
followed by
```
source activate BMI
```

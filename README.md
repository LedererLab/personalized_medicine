# Euclidean Distance Ridge Estimator

This repository provides the data and implementations of the methods described in [Tuning parameter calibration for prediction in personalized medicine](To Do).

# Usage 

The file `PAV.md` contains a function `PAVedr` to select a tuning parameter within a set of tuning parameters for ridge regression that minimizes the prediction error for an specific individual covariate. The aforementioned paper contains detailed descriptions of these methods.


# Simulations

We provide an example code in `SimulationStudy.md` for a comparison of averaged individual prediction errors with 10-fold cross-validation. This program requires `R 3.6.1` or earlier.

# Real Data Analyses

**Kidney Transplant Data Analysis :**

The processed data is in `personalized_medicine/data/E-GEOD-33070_clinical.RData`, 
`personalized_medicine/data/E-GEOD-33070_gene.RData`, and 
`personalized_medicine/data/index_for_parameters.RData`, which are preprocessed by Kristoffer Herland Hellton and Shih-Ting Huang. The original raw kidney transplant data was collected and provided by [Expression Levels of Obesity-Related Genes Are Associated with Weight Change in Kidney Transplant Recipients](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059962). The programs for our statistical analysis is in RealDataAnalysis.md.

# Repository Authors 

* Shih-Ting Huang, Ph.D. student in Mathematical Statistics, Ruhr-University Bochum

* Yannick Düren, Ph.D. student in Mathematical Statistics, Ruhr-University Bochum

* Kristoffer Herland Hellton, Researcher, Norwegian Computing Center

* Johannes Lederer, Professor in Mathematical Statistics, Ruhr-University Bochum

# Other files

**RidgeCv.md** : K-fold cross-validation for ridge regression.

**RidgeEstimator.md** : Computing ridge estimator.

## Supported Languages and platforms

All of the codes in this repository are written in R and supports all plarforms which are
 supported by R itself.

## Dependencies

This repository does not depend on any R libraries or external sources.

## Licensing

The HDIM package is licensed under the MIT license. To
view the MIT license please consult `LICENSE.txt`.

## References
[Tuning parameter calibration for prediction in personalized medicine](To Do)


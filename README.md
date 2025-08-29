# Regression Discontinuity Tree

Welcome to the GitHub page of Regression Discontinuity Tree. 
This page provides the codes for this method proposed by Reguly (2025): Discovering Heterogeneous Treatment Effects in Regression Discontinuity Designs


## Use of the codes

This repository contains codes and a description of how to estimate heterogeneous treatment effects in RDD. Codes are built under MatLab, version 2024a.

1. Add the downloaded/forked/cloned **algorithm** folder to your path in MatLab.
2. Open **tryout.m**, which gives a toy example of how to use the package.
  - The main function is `runTree_uni()` function.
  - However, you need to carefully specify your `sample` object as well as `optCART` object. The comments should guide you on how to set them.
3. `tostring()` function will draw you the resulting tree.
 
**Note:** This is the second version, containing nonparametric estimator.

## Replication Reguly (2025)

1. **algorithm** folder contains necessary functions to run regression discontinuity tree. One needs to add this folder to the path when using any other script.
2. **simulations** folder contain codes to replicate simulation results reported in Section 4. **sharp** folder contains `runRDD_sharp.m` file which runs the monte carlo simulations. Caution -- it takes considerable amount of time  and uses multiple cores. (1day+ with Apple M1 Max chip with 64GB memory and 9 cores) Note if one does not run with 50,000 observations or reduces the MC iteration number it reduces the time significantly. After running the simulations, one can evaluate the results with `evaluate_results/createTables.m`.
3. **empirical_application** folder contains codes to replicate Section 5. It also contains a further readme file to navigate through the replication.

*All feedbacks are welcomed!* Please write to: agoston.reguly@uni-corvinus.hu

## Python implementation

A python implementation of the **parametric** version of the paper is done by Firat Yaman. Check [here](https://www.firat-yaman.de/Codes.html).

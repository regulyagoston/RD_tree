# Regression Discontinuity Tree

Welcome to the GitHub page of Regression Discontinuity Tree. 
This page provides the codes for this method proposed by [Reguly (2021): Heterogeneous Treatment Effects in Regression Discontinuity Designs](https://arxiv.org/abs/2106.11640).

Regression Discontinuity Tree extends Athey-Imbens' (2016) causal trees to regression discontinuity designs.


## Use of the codes

This repository contains codes and a description of how to estimate heterogeneous treatment effects in RDD. Codes are built under MatLab, version 2016b.

1. Add the downloaded/forked/cloned **codes** folder to your path in MatLab.
2. Open **tryout.m**, which gives a toy example of how to use the package.
  - The main function is `runTree()` function.
  - However, you need to carefully specify your `sample` object as well as `optCART` object. The comments should guide you on how to set them.
3. `tostring()` function will draw you the resulting tree.
4. After line 86 you can find parts of the `runTree()` function if you want to have a more detailed picture of what the algorithm does.

 
**Note:** This is the first version, many refinements are needed. Description of the codes is work-in-progress.

*All feedbacks are welcomed!* Please write to: reguly_agoston-at-uni-corvinus.hu

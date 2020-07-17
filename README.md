# ABTCK

Augmented Bayesian treed co-kriging


## MIT License

Copyright (C) 2019 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This file is part of the Github repository, ABTCK, Augmented Bayesian treed co-kriging

## Description 

MATLAB scripts implementing our proposed Bayesian procedure,'Augmented Bayesian treed co-kriging'. This procedure can be used for the statistical analysis of computer models when the available simulations are in multiple fidelity, they are not necessarily based on hierarchically nested experimental designs, and they may present local features, non-stationarities, or discontinuities.

This is a MATLAB code has been used to produce the numerical results and the statistical analysis reported in the paper  

* "Bayesian analysis of multifidelity computer models with local features and non-nested experimental designs: Application to the WRF model", by Alexandros Konomi and Georgios Karagianniswhich is submitted in Technometrics.

## Requirements 

MATLAB compiler (R2017a or later)

## Files 

* ABCK_M: has part of the imputation techniques used in the paper

* cov_functions: has all the covaraince functions neccessary for the GP 

* function_fold: has all the functions used in the simulation study

* Multi_GP_operations: has all the MCMC operations for the co-kriging Gaussian process hyperparamter updates

* Multi_prediction: has the functions used for different prediction strategies

* Multi_tree_nested: has all the MCMC treed update of the model when the design is nested

* Multi_tree_NON_nestedB: has all the MCMC treed update of the model when the design is nested

* Random_var: has different algorithms for generating data from distributions necessary in the BTC

* help_tree_operations: some operations that help to define the treed form

* Results: This is an empty folder where results are saved

* FirstSimulationStudy: Contain the example 1 of the paper (section 3.1). It also contain different variations of that example. The user can choose to change the function.  

* SecondSinulationStudy: Contain the example 2 of the paper (section 3.1). Here the user can also change the function and the setting

* Online_ExampleP: Contain the example in the supplementary material S.3 of the paper.

* Application_WRF2: Contain the application to WRF dataset in the paper (section 4). The user can change the settings for testing/training datasets. 







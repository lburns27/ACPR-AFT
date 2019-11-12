# ACPR-AFT
(A)CPR-AFT contains methods for feature selection and survival prediction utlizing continuum power regression (CPR) in conjuction with a parametric and semi-parametric accelerated failure time (AFT) model.

This repository contains R code for the various algorithms and functions involved in (A)CPR-AFT. The code is organized as outlined below.

* **CPR & LOOCV based on PRESS**: Functions for running CPR and LOOCV based on PRESS
  * doCPR: function for running CPR given desired alpha
  * press1: function for running LOOCV based on press
  
* **Algorithm 1: CPR-AFT**: Runs Steps 1-4 of CPR-AFT algorithm 
  * Required: Run **CPR & LOOCV based on PRESS** first
  * Input: dataset named 'dat' (column 1 = survival time, column 2 = censoring indicator, columns 3+ = gene expression)
  * Output options: CPR components, CPR coefficients, or VIP (see R code for details)

* **Algorithm 2: ACPR-AFT**: Runs Steps 1-6 of ACPR-AFT algorithm 
  * Required: Run **CPR & LOOCV based on PRESS** first; R packages (flexsurv, lss)
  * sAFT: Runs semi-parametric AFT version of Algorithm 2 (censored times adjusted based on sAFT)
  * genF: Runs generalized-F version of Algorithm 2 (censored times adjusted based on genF)
  * Input: dataset named 'dat' (column 1 = survival time, column 2 = censoring indicator, columns 3+ = gene expression)
  * Output options: CPR components, CPR coefficients, or VIP (see R code for details)

* **Algorithm 3: Survival Prediction Algorithm**: Runs Steps 1-6 of the Survival Prediction Algorithm.
  * Required: Run **CPR & LOOCV based on PRESS** first; R packages (lss, flexsurv, survival, survivalROC, concreg, Hmisc)
  * Input: filtered dataset based on marginal screening named 'dat' (column 1 = survival time, column 2 = censoring indicator, columns 3+ = gene expression)
  * Runs (A)CPR-AFT across random training/test set splits, calculates PI for both training and test tests in each case, and computes various measures of prediction accuracy.
  * Output options: Final output are the summarized measures of prediction accuracy (intermediate results also saved and stored throughout)

* **Youden & AUC**: Computes Youden & AUC values based on gene ranking by CPR coefficients and VIP resulting from (A)CPR-AFT algorithm
  * Required: R packages (MESS)
  * Inputs: 
    * data frame with n rows (number of genes) and 2 columns ($cprCoeff and $vip).  Note, these are obtained by running the (A)CPR-AFT algorithms (algorithms 1 & 2)
    * effectGenes: number of significant genes
  * Output option: Specificity, Sensitivity, Youden, AUC & AUC standard deviation for the CPR coefficients and VIP 

* **Simulations**: Contains R code for simulating data
  * Scheme 1: Univariate aproach; genes linked to survival one at a time
  * Scheme 2: Multivariate approach; incorporates correlation between features
  * For both schemes, there are options to simulate from the following models: LN, LL1, LL2, W1, W2

# Copyright & Citations
Copyright © Lauren Spirko-Burns and Karthik Devarajan 

Spirko-Burns, L., Devarajan, K. Supervised dimension reduction for large-scale “omics" data with censored survival outcomes under possible non-proportional hazards. Under review. bioRxiv 586529; doi: https://doi.org/10.1101/586529. COBRA pre-print series, Article 119 (March 2019). http://biostats.bepress.com/cobra/art119.

Spirko, L. (2017). Variable Selection and Supervised Dimension Reduction for Large-Scale Genomic Data with Censored Survival Outcomes. Ph.D. Dissertation. Department of Statistical Science, Temple University, Philadelphia.

# License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">ACPR-AFT</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Lauren Spirko Burns and Karthik Devarajan</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

# 

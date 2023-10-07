# MVBCF Paper
Materials for replicating the results in "Bayesian Causal Forests for Multivariate Outcomes: Application to Irish Data From an International Large Scale Education Assessment."

The preprint of this article can be found at https://browse.arxiv.org/pdf/2303.04874.pdf

The TIMSS data used in the study is available for download at https://timss2019.org/international-database/

### Code Scripts and Description:

MVBCF_Code.cpp - RCPP Implementation of Multivariate BCF.

MVBCF_RI_Code.cpp - RCPP Implementation of Multivariate BCF (With Random Intercepts).

<br/>

DGP1_Github.R - Code for evaluating model performance on Data Generating Process 1.

DGP2_Github.R - Code for evaluating model performance on Data Generating Process 2.

DGP3_Github.R - Code for evaluating model performance on Data Generating Process 3.

<br/>

Desk_Github.R - Code for estimating treatment effects of "Has Study Desk" treatment.

Hungry_Github.R - Code for estimating treatment effects of "Often Hungry" treatment.

Absent_Github.R - Code for estimating treatment effects of "Often Absent" treatment.

<br/>

Desk_Github_OOS.R - Code for 10-fold cross validation results with "Has Study Desk" treatment.

Hungry_Github_OOS.R - Code for 10-fold cross validation results with "Often Hungry" treatment.

bsent_Github_OOS.R - Code for 10-fold cross validation results with "Often Absent" treatment.

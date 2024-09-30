# MVBCF Paper
Materials for replicating the results in "Bayesian Causal Forests for Multivariate Outcomes: Application to Irish Data From an International Large Scale Education Assessment."

The preprint of this article can be found at https://browse.arxiv.org/pdf/2303.04874.pdf

The TIMSS data used in the study is available for download at https://timss2019.org/international-database/

An R package implementing the model from the paper is available at https://nathan-mcjames.github.io/mvbcf/ 

![alt text](https://github.com/Nathan-McJames/MVBCF_Paper/blob/main/Pictures/paper_plot.svg?raw=true)

### Code Scripts and Description:

MVBCF_Code.cpp - RCPP Implementation of Multivariate BCF.

MVBCF_RI_Code.cpp - RCPP Implementation of Multivariate BCF (With Random Intercepts).

<br/>

GitHub_DGP1.R - Code for evaluating model performance on Data Generating Process 1.

GitHub_DGP2.R - Code for evaluating model performance on Data Generating Process 2.

GitHub_DGP3.R - Code for evaluating model performance on Data Generating Process 3.

<br/>

GitHub_Desk.R - Code for estimating treatment effects of "Has Study Desk" treatment.

GitHub_Hungry.R - Code for estimating treatment effects of "Often Hungry" treatment.

GitHub_Absent.R - Code for estimating treatment effects of "Often Absent" treatment.

<br/>

GitHub_Desk_OOS.R - Code for 10-fold cross validation results with "Has Study Desk" treatment.

GitHub_Hungry_OOS.R - Code for 10-fold cross validation results with "Often Hungry" treatment.

GitHub_Absent_OOS.R - Code for 10-fold cross validation results with "Often Absent" treatment.

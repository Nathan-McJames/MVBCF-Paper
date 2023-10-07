# MVBCF Paper
Materials for replicating the results in "Bayesian Causal Forests for Multivariate Outcomes: Application to Irish Data From an International Large Scale Education Assessment."

The preprint of this article can be found at https://browse.arxiv.org/pdf/2303.04874.pdf

The TIMSS data used in the study is available for download at https://timss2019.org/international-database/

## Code Scripts and Description:

MVBCF_Code.R - RCPP Implementation of Multivariate BCF.\n
MVBCF_RI_Code.R - RCPP Implementation of Multivariate BCF (With Random Intercepts).

DGP1_Github.R - Code for evaluating model performance on Data Generating Process 1.\n
DGP2_Github.R - Code for evaluating model performance on Data Generating Process 2.\n
DGP3_Github.R - Code for evaluating model performance on Data Generating Process 3.

Desk_Github.R - Code for estimating treatment effects of "Has Study Desk" treatment.\n
Hungry_Github.R - Code for estimating treatment effects of "Often Hungry" treatment.\n
Absent_Github.R - Code for estimating treatment effects of "Often Absent" treatment.

Desk_Github_OOS.R - Code for 10-fold cross validation results with "Has Study Desk" treatment.\n
Hungry_Github_OOS.R - Code for 10-fold cross validation results with "Often Hungry" treatment.\n
bsent_Github_OOS.R - Code for 10-fold cross validation results with "Often Absent" treatment.

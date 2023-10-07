#import libraries
library(mvtnorm)
library(bartCause)
library(skewBART)
library(tidyverse)
library(dplyr)
library(magrittr)
library(zeallot)
library(dbarts)
library(bcf)
library(scoringRules)
library(matrixStats)

#Define Helper Functions
in_cred<-function(samples, value, interval)
{
  upper_quantile<-1-(1-interval)/2
  lower_quantile<-0+(1-interval)/2
  
  q1<-quantile(samples, lower_quantile)
  q2<-quantile(samples, upper_quantile)
  
  in_cred<-ifelse(value>=q1 & value<=q2, T, F)
}


#Set random seed
seed_val<-round(runif(1, min=1, max=9999999))
set.seed(seed_val)


#Train Data
n<-sample(c(500, 1000), 1)

X1<-runif(n)
X2<-runif(n)
X3<-runif(n)
X4<-runif(n)
X5<-runif(n)
X6<-rbinom(n, 1, 0.5)
X7<-rbinom(n, 1, 0.5)
X8<-rbinom(n, 1, 0.5)
X9<-sample(c(0, 1, 2, 3, 4), n, replace=T)
X10<-sample(c(0, 1, 2, 3, 4), n, replace=T)

X<-cbind(X1, X2, X3, X5, X6, X7, X8, X9, X10)

Mu1<-(11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
Mu2<-(9*sin(pi*X1*X2)+22*(X3-0.5)^2+0*X4+8*X6+X9)*10+300

Tau1<-(2*X4+2*X5)*10
Tau2<-(1*X4+3*X5)*10

true_propensity<-X4

Z<-rbinom(n, 1, true_propensity)

Y<-cbind(Mu1+Z*Tau1, Mu2+Z*Tau2) + mvtnorm::rmvnorm(n, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow=2, byrow=T))


#Test Data
n_test<-1000

X1_test<-runif(n_test)
X2_test<-runif(n_test)
X3_test<-runif(n_test)
X4_test<-runif(n_test)
X5_test<-runif(n_test)
X6_test<-rbinom(n_test, 1, 0.5)
X7_test<-rbinom(n_test, 1, 0.5)
X8_test<-rbinom(n_test, 1, 0.5)
X9_test<-sample(c(0, 1, 2, 3, 4), n_test, replace=T)
X10_test<-sample(c(0, 1, 2, 3, 4), n_test, replace=T)

X_test<-cbind(X1_test, X2_test, X3_test, X5_test, X6_test, X7_test, X8_test, X9_test, X10_test)

Mu1_test<-(11*sin(pi*X1_test*X2_test)+18*(X3_test-0.5)^2+10*X4_test+12*X6_test+X9_test)*10+300
Mu2_test<-(9*sin(pi*X1_test*X2_test)+22*(X3_test-0.5)^2+0*X4_test+8*X6_test+X9_test)*10+300

Tau1_test<-(2*X4_test+2*X5_test)*10
Tau2_test<-(1*X4_test+3*X5_test)*10

true_propensity_test<-X4_test

Z_test<-rbinom(n_test, 1, true_propensity_test)

Y_test<-cbind(Mu1_test+Z_test*Tau1_test, Mu2_test+Z_test*Tau2_test) + mvtnorm::rmvnorm(n_test, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow=2, byrow=T))

#estimate of propensity score
p_mod<-bart(x.train = X, y.train = Z, x.test = X_test, k=3)
p<-colMeans(pnorm(p_mod$yhat.train))
p_test<-colMeans(pnorm(p_mod$yhat.test))

#adding propensity score to matrix
X2<-X
X2_test<-X_test
X<-cbind(X, p)
X_test<-cbind(X_test, p_test)
Z2<-cbind(Z,Z)

#set some parameters
n_tree_mu<-50
n_tree_tau<-20
n_iter<-1000
n_burn<-500

#Run the mvbcf model
sourceCpp(file = "~/your_directory/MVBCF_Code.cpp")

mu_val<-1
tau_val<-0.375
v_val<-1
wish_val<-1
min_val<-1

mvbcf_mod <- fast_bart(X,
                       Y,
                       Z2,
                       X2,
                       X_test, #here is the test data for the mu part of the model
                       X2_test, #here is the test data for the tau part of the model
                       0.95,
                       2,
                       0.25,
                       3,
                       diag((mu_val)^2/n_tree_mu, 2),
                       diag((tau_val)^2/n_tree_tau, 2),
                       v_val,
                       diag(wish_val, 2),
                       n_iter,
                       n_tree_mu,
                       n_tree_tau,
                       min_val)


mvbcf_tau_preds1<-rowMeans(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)])
mvbcf_ate1<-mean(mvbcf_tau_preds1)
mvbcf_tau_preds2<-rowMeans(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)])
mvbcf_ate2<-mean(mvbcf_tau_preds2)





#run the bart model
colnames(X)<-rep("", ncol(X))
colnames(X_test)<-rep("", ncol(X_test))
colnames(X)<-paste0("V", 1:ncol(X))
colnames(X_test)<-paste0("V", 1:ncol(X_test))
bart_mod1<-bartc(Y[,1], Z, X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau)
bart_tau_preds1<-colMeans(predict(bart_mod1, X_test, type="icate"))
bart_ate1<-mean(bart_tau_preds1)


bart_mod2<-bartc(Y[,2], Z, X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau)
bart_tau_preds2<-colMeans(predict(bart_mod2, X_test, type="icate"))
bart_ate2<-mean(bart_tau_preds2)


#Fit Univariate BCF Model to Outcome 1
bcfmod1<-bcf(Y[,1], Z, X2, 
             pihat=p, nburn=n_burn, nsim=n_iter-n_burn, 
             ntree_control = n_tree_mu,
             ntree_moderate = n_tree_tau,
             n_chains=1,
             n_cores=1,
             n_threads=1,
             use_muscale = F,
             use_tauscale = F,
             sd_control = sd(Y[,1]),
             sd_moderate = 0.375*sd(Y[,1]),
             save_tree_directory = "/your_directory/",
             include_pi = "control",
             random_seed = seed_val)


pred_out1<-predict(object=bcfmod1,
                   x_predict_control=X_test,
                   x_predict_moderate=X2_test,
                   pi_pred=p_test,
                   z_pred=Z_test,
                   n_chains=1,
                   n_cores=1,
                   n_threads=1,
                   save_tree_directory = "/your_directory/")


bcf_tau_preds1<-colMeans(pred_out1$tau)
bcf_ate1<-mean(bcf_tau_preds1)

#Fit Univariate BCF Model to Outcome 2
bcfmod2<-bcf(Y[,2], Z, X2, 
             pihat=p, nburn=n_burn, nsim=n_iter-n_burn, 
             ntree_control = n_tree_mu,
             ntree_moderate = n_tree_tau,
             n_chains=1,
             n_cores=1,
             n_threads=1,
             use_muscale = F,
             use_tauscale = F,
             sd_control = sd(Y[,2]),
             sd_moderate = 0.375*sd(Y[,2]),
             save_tree_directory = "/your_directory/",
             include_pi = "control",
             random_seed = seed_val)


pred_out2<-predict(object=bcfmod2,
                   x_predict_control=X_test,
                   x_predict_moderate=X2_test,
                   pi_pred=p_test,
                   z_pred=Z_test,
                   n_chains=1,
                   n_cores=1,
                   n_threads=1,
                   save_tree_directory = "/your_directory/")

bcf_tau_preds2<-colMeans(pred_out2$tau)
bcf_ate2<-mean(bcf_tau_preds2)

#fit multivariate skew bart
hypers <- Hypers(X = cbind(X, Z), Y = Y, num_tree = n_tree_mu+n_tree_tau)
opts <- Opts(num_burn = n_burn, num_save = n_iter-n_burn)
X_test_z<-rbind(cbind(X_test, 0*Z_test), cbind(X_test, 0*Z_test+1))

fitted_Multiskewbart <- MultiskewBART(X = cbind(X, Z), Y = Y, test_X = X_test_z, hypers=hypers, opts=opts, do_skew = F)

z_0_preds<-fitted_Multiskewbart$y_hat_test[1:n_test,,]
z_1_preds<-fitted_Multiskewbart$y_hat_test[-c(1:n_test),,]
z_1_0_preds<-z_1_preds-z_0_preds

mvbart_icates1<-rowMeans(z_1_0_preds[,1,])
mvbart_ate1<-mean(mvbart_icates1)
mvbart_icates2<-rowMeans(z_1_0_preds[,2,])
mvbart_ate2<-mean(mvbart_icates2)



#get metrics
mvbcf_pehe1<-sqrt(mean((Tau1_test-mvbcf_tau_preds1)^2))
bart_pehe1<-sqrt(mean((Tau1_test-bart_tau_preds1)^2))
bcf_pehe1<-sqrt(mean((Tau1_test-bcf_tau_preds1)^2))
mvbart_pehe1<-sqrt(mean((Tau1_test-mvbart_icates1)^2))
mvbcf_bias1<-mvbcf_ate1-mean(Tau1_test)
bart_bias1<-bart_ate1-mean(Tau1_test)
bcf_bias1<-bcf_ate1-mean(Tau1_test)
mvbart_bias1<-mvbart_ate1-mean(Tau1_test)



mvbcf_tau_501<-mean(diag(apply(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)], 1, in_cred, Tau1_test, 0.5)))
bcf_tau_501<-mean(diag(apply(pred_out1$tau, 2, in_cred, Tau1_test, 0.5)))
bart_tau_501<-mean(diag(apply(predict(bart_mod1, X_test, type="icate"), 2, in_cred, Tau1_test, 0.5)))
mvbart_tau_501<-mean(diag(apply(z_1_0_preds[,1,], 1, in_cred, Tau1_test, 0.5)))

mvbcf_tau_951<-mean(diag(apply(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)], 1, in_cred, Tau1_test, 0.95)))
bcf_tau_951<-mean(diag(apply(pred_out1$tau, 2, in_cred, Tau1_test, 0.95)))
bart_tau_951<-mean(diag(apply(predict(bart_mod1, X_test, type="icate"), 2, in_cred, Tau1_test, 0.95)))
mvbart_tau_951<-mean(diag(apply(z_1_0_preds[,1,], 1, in_cred, Tau1_test, 0.95)))





mvbcf_pehe2<-sqrt(mean((Tau2_test-mvbcf_tau_preds2)^2))
bart_pehe2<-sqrt(mean((Tau2_test-bart_tau_preds2)^2))
bcf_pehe2<-sqrt(mean((Tau2_test-bcf_tau_preds2)^2))
mvbart_pehe2<-sqrt(mean((Tau2_test-mvbart_icates2)^2))
mvbcf_bias2<-mvbcf_ate2-mean(Tau2_test)
bart_bias2<-bart_ate2-mean(Tau2_test)
bcf_bias2<-bcf_ate2-mean(Tau2_test)
mvbart_bias2<-mvbart_ate2-mean(Tau2_test)


mvbcf_tau_502<-mean(diag(apply(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)], 1, in_cred, Tau2_test, 0.5)))
bcf_tau_502<-mean(diag(apply(pred_out2$tau, 2, in_cred, Tau2_test, 0.5)))
bart_tau_502<-mean(diag(apply(predict(bart_mod2, X_test, type="icate"), 2, in_cred, Tau2_test, 0.5)))
mvbart_tau_502<-mean(diag(apply(z_1_0_preds[,2,], 1, in_cred, Tau2_test, 0.5)))

mvbcf_tau_952<-mean(diag(apply(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)], 1, in_cred, Tau2_test, 0.95)))
bcf_tau_952<-mean(diag(apply(pred_out2$tau, 2, in_cred, Tau2_test, 0.95)))
bart_tau_952<-mean(diag(apply(predict(bart_mod2, X_test, type="icate"), 2, in_cred, Tau2_test, 0.95)))
mvbart_tau_952<-mean(diag(apply(z_1_0_preds[,2,], 1, in_cred, Tau2_test, 0.95)))

df<-data.frame(mvbcf_pehe1,
               bart_pehe1,
               bcf_pehe1,
               mvbart_pehe1,
               mvbcf_bias1,
               bart_bias1,
               bcf_bias1,
               mvbart_bias1,
               
               
               mvbcf_tau_501,
               bcf_tau_501,
               bart_tau_501,
               mvbart_tau_501,
               
               mvbcf_tau_951,
               bcf_tau_951,
               bart_tau_951,
               mvbart_tau_951,
               
               
               mvbcf_pehe2,
               bart_pehe2,
               bcf_pehe2,
               mvbart_pehe2,
               mvbcf_bias2,
               bart_bias2,
               bcf_bias2,
               mvbart_bias2,
               
               
               mvbcf_tau_502,
               bcf_tau_502,
               bart_tau_502,
               mvbart_tau_502,
               
               mvbcf_tau_952,
               bcf_tau_952,
               bart_tau_952,
               mvbart_tau_952,
               
               
               seed_val,
               n
               
)



write.table(df,"~/your_directory/DGP2_Github.csv", row.names = F, append=T, sep = ",", col.names=F)





